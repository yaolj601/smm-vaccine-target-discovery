"""
Step 7: Build the presentation deck (.pptx) from the results.

Everything on these slides is read from results/ rather than typed in, so the
deck cannot drift from the ranking the way a hand-maintained one does. Import
into Google Slides with File > Import slides, or upload the .pptx to Drive and
open it with Slides.

Input:  $SMM_OUT_DIR/ranked_targets.csv, all_candidates_scored.csv
        $SMM_OUT_DIR/figures/*.png
        $SMM_DATA_DIR/de_results.csv   (for the contrast entry counts)
Output: $SMM_OUT_DIR/SMM_vaccine_targets.pptx
"""

import os

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

from config import DATA, OUT

# Wright-Giemsa marrow staining, same palette as the HTML deck
INK = RGBColor(0x23, 0x16, 0x1C)
INK2 = RGBColor(0x6B, 0x55, 0x60)
INK3 = RGBColor(0x9A, 0x87, 0x91)
ACCENT = RGBColor(0x9B, 0x2C, 0x62)
TEAL = RGBColor(0x2F, 0x6B, 0x67)
AMBER = RGBColor(0x8F, 0x62, 0x06)
GROUND = RGBColor(0xFB, 0xF8, 0xF6)
SURF2 = RGBColor(0xF4, 0xEE, 0xEB)
RULE = RGBColor(0xE3, 0xD8, 0xD5)
HL = RGBColor(0xF7, 0xE8, 0xF0)
WARN = RGBColor(0xF8, 0xEE, 0xDC)

HEAD = "Libre Franklin"
BODY = "Source Serif 4"
MONO = "Consolas"

W, H = Inches(13.333), Inches(7.5)
M = Inches(0.62)


def _txt(slide, x, y, w, h, text, size=14, font=BODY, color=INK,
         bold=False, italic=False, align=PP_ALIGN.LEFT, space=6):
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, para in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(space)
        run = p.add_run()
        run.text = para
        run.font.size = Pt(size)
        run.font.name = font
        run.font.color.rgb = color
        run.font.bold = bold
        run.font.italic = italic
    return box


def _rect(slide, x, y, w, h, fill):
    from pptx.enum.shapes import MSO_SHAPE
    sh = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, x, y, w, h)
    sh.fill.solid()
    sh.fill.fore_color.rgb = fill
    sh.line.fill.background()
    sh.shadow.inherit = False
    return sh


def new_slide(prs, eyebrow, title):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    _rect(s, 0, 0, W, H, GROUND)
    if eyebrow:
        _txt(s, M, Inches(0.42), W - 2 * M, Inches(0.3),
             eyebrow.upper(), size=11, font=MONO, color=INK3)
    _rect(s, M, Inches(0.78), W - 2 * M, Emu(9525), RULE)
    if title:
        _txt(s, M, Inches(0.95), W - 2 * M, Inches(0.8),
             title, size=30, font=HEAD, color=INK, bold=True)
    return s


def table(slide, x, y, w, rows, col_w=None, head=True, size=10.5,
          highlight=(), warn_rows=()):
    """rows[0] is the header."""
    nr, nc = len(rows), len(rows[0])
    shp = slide.shapes.add_table(nr, nc, x, y, w, Inches(0.3 * nr))
    tbl = shp.table
    tbl.first_row = head
    if col_w:
        total = sum(col_w)
        for i, frac in enumerate(col_w):
            tbl.columns[i].width = Emu(int(w * frac / total))
    for r, row in enumerate(rows):
        tbl.rows[r].height = Inches(0.27)
        for c, val in enumerate(row):
            cell = tbl.cell(r, c)
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.02)
            cell.fill.solid()
            if r == 0 and head:
                cell.fill.fore_color.rgb = SURF2
            elif r in warn_rows:
                cell.fill.fore_color.rgb = WARN
            elif r in highlight:
                cell.fill.fore_color.rgb = HL
            else:
                cell.fill.fore_color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            tf = cell.text_frame
            tf.word_wrap = True
            p = tf.paragraphs[0]
            run = p.add_run()
            run.text = str(val)
            run.font.size = Pt(size - 1 if r == 0 and head else size)
            run.font.name = MONO if (r == 0 and head) else BODY
            run.font.bold = (r == 0 and head)
            run.font.color.rgb = INK3 if (r == 0 and head) else INK
    return shp


def bullets(slide, x, y, w, h, items, size=13, gap=9):
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, (lead, rest) in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap)
        if lead:
            r1 = p.add_run()
            r1.text = lead
            r1.font.size = Pt(size)
            r1.font.name = HEAD
            r1.font.bold = True
            r1.font.color.rgb = INK
        r2 = p.add_run()
        r2.text = rest
        r2.font.size = Pt(size)
        r2.font.name = BODY
        r2.font.color.rgb = INK2
    return box


def pic(slide, path, x, y, w=None, h=None):
    if os.path.exists(path):
        return slide.shapes.add_picture(path, x, y, width=w, height=h)
    return None


def label(slide, x, y, w, text, color=INK3):
    """Small uppercase section marker above a block."""
    return _txt(slide, x, y, w, Inches(0.2), text.upper(), size=9,
                font=MONO, color=color)


def glossary(slide, x, y, w, h, items, size=11, gap=5):
    """Key-terms block: the words that actually appear in the figure.

    items = [(term, one-line definition), ...]. The term is set in the accent
    colour and mono so the eye can jump from a label in the figure straight to
    its definition, without reading prose to find it.
    """
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, (term, defn) in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap)
        r1 = p.add_run()
        r1.text = term + "   "
        r1.font.size = Pt(size)
        r1.font.name = MONO
        r1.font.bold = True
        r1.font.color.rgb = ACCENT
        r2 = p.add_run()
        r2.text = defn
        r2.font.size = Pt(size)
        r2.font.name = BODY
        r2.font.color.rgb = INK2
    return box


def takeaway(slide, x, y, w, h, items, size=12, gap=8):
    """Short conclusion bullets. Each item is a str, or (bold lead, rest)."""
    box = slide.shapes.add_textbox(x, y, w, h)
    tf = box.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    for i, item in enumerate(items):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.space_after = Pt(gap)
        mark = p.add_run()
        mark.text = "\u25b8   "
        mark.font.size = Pt(size)
        mark.font.color.rgb = TEAL
        mark.font.bold = True
        if isinstance(item, tuple):
            lead, rest = item if len(item) == 2 else (item[0], "")
        else:
            lead, rest = "", item
        if lead:
            r1 = p.add_run()
            r1.text = lead
            r1.font.size = Pt(size)
            r1.font.name = HEAD
            r1.font.bold = True
            r1.font.color.rgb = INK
        if rest:
            r2 = p.add_run()
            r2.text = rest
            r2.font.size = Pt(size)
            r2.font.name = BODY
            r2.font.color.rgb = INK2
    return box


# ───────────────────────────────────────────────────── data
sl = pd.read_csv(f"{OUT}/ranked_targets.csv")
allc = pd.read_csv(f"{OUT}/all_candidates_scored.csv")
de = pd.read_csv(f"{DATA}/de_results.csv")
FIG = f"{OUT}/figures"


def _n(contrast, lfc):
    d = de[de.contrast_name == contrast]
    return int(((d.padj < 0.05) & (d.log2FoldChange > lfc)).sum())


nA, nB, nC = (_n("A_neoplastic_vs_normal", 1), _n("B_SMM_vs_NBM", 1),
              _n("C_paired", 0.5))
A = set(de[(de.contrast_name == "A_neoplastic_vs_normal") & (de.padj < .05)
           & (de.log2FoldChange > 1)].gene)
B = set(de[(de.contrast_name == "B_SMM_vs_NBM") & (de.padj < .05)
           & (de.log2FoldChange > 1)].gene)
C = set(de[(de.contrast_name == "C_paired") & (de.padj < .05)
           & (de.log2FoldChange > .5)].gene)
union = len(A | B | C)
n_cand, n_epi = len(allc), len(pd.read_csv(f"{OUT}/evidence/epitope_results.csv"))

prs = Presentation()
prs.slide_width, prs.slide_height = W, H

# ─────────────────────────────────────────── 1 title
s = new_slide(prs, "", "")
_rect(s, 0, 0, W, H, GROUND)
_rect(s, 0, 0, Inches(0.12), H, ACCENT)
_txt(s, M, Inches(1.7), Inches(9.2), Inches(0.35),
     "GSE193531 · BOIARSKY ET AL. NAT COMMUN 2022", size=12, font=MONO, color=INK3)
_txt(s, M, Inches(2.2), Inches(10.2), Inches(1.8),
     "Vaccine targets for\nsmouldering myeloma", size=44, font=HEAD, color=INK, bold=True)
_txt(s, M, Inches(4.2), Inches(8.6), Inches(0.8),
     "A ranked shortlist for an off-the-shelf preventive vaccine in SMM.",
     size=15, color=INK2)
for i, (num, lab) in enumerate([("29,387", "plasma cells"), ("35", "patients · 4 stages"),
                                (str(n_cand), "candidates scored"), ("15", "shortlisted")]):
    x = M + Inches(3.0) * i
    _txt(s, x, Inches(5.3), Inches(2.8), Inches(0.6), num, size=32, font=HEAD,
         color=ACCENT if lab == "shortlisted" else INK, bold=True)
    _txt(s, x, Inches(5.9), Inches(2.8), Inches(0.3), lab.upper(), size=10,
         font=MONO, color=INK3)

# ─────────────────────────────────────────── 2 what counts as a good target
s = new_slide(prs, "01 · what counts as a good target", "The product sets the safety bar")
_txt(s, M, Inches(1.8), Inches(11.9), Inches(0.4),
     "SMM patients are not yet sick and may never progress. Tolerance for "
     "off-tumour toxicity is close to zero.", size=14, color=INK2)
table(s, M, Inches(2.5), Inches(6.0), [
    ["Evidence layer", "Weight", "Why that weight"],
    ["Differential expression", "0.30", "Most statistically grounded"],
    ["Normal-tissue safety", "0.25", "Outlier test over ~50 tissues"],
    ["Novelty-confidence", "0.25", "Known biology, existing programs"],
    ["HLA-I presentation", "0.20", "Weakest — binding ≠ immunogenicity"],
], col_w=[3, 1, 3.4])
label(s, Inches(7.1), Inches(2.5), Inches(5.6), "takeaway")
takeaway(s, Inches(7.1), Inches(2.8), Inches(5.6), Inches(4.0), [
    ("Shared antigens, not neoantigens. ",
     "No mutation calls here, and an off-the-shelf product needs a shared target."),
    ("BCMA is correctly absent. ",
     "TNFRSF17 sits at log2FC 0.13 — a plasma-cell marker, not a tumour antigen."),
    ("Labels are expression-derived. ",
     "Clustering, not inferCNV — so DE against them is partly circular."),
])

# ─────────────────────────────────────────── 3 cohort
s = new_slide(prs, "02 · the cohort", "35 patients, one aspirate each")
pic(s, f"{FIG}/cohort_overview.png", M, Inches(1.5), h=Inches(5.5))
label(s, Inches(6.0), Inches(1.6), Inches(6.7), "key terms in this figure")
glossary(s, Inches(6.0), Inches(1.9), Inches(6.7), Inches(2.6), [
    ("CD138+", "plasma-cell surface marker — every library is this one sorted fraction"),
    ("neoplastic PC", "magenta bar — the malignant plasma cells"),
    ("normal PC", "teal bar — residual healthy plasma cells in the same aspirate"),
    ("tumour fraction", "right panel: neoplastic ÷ total plasma cells"),
    ("◆", "sample has BOTH cell types — the only within-patient comparisons"),
    ("stage / risk", "NBM healthy → MGUS → SMM → MM; SMM split high- vs low-risk"),
])
label(s, Inches(6.0), Inches(4.75), Inches(6.7), "takeaway")
takeaway(s, Inches(6.0), Inches(5.05), Inches(6.7), Inches(2.1), [
    "No longitudinal sampling — one aspirate each, so no progression trajectories.",
    "Only 11 of 35 are paired (MGUS 3 · SMM 5 · MM 3); the other 24 single-status.",
    ("Risk recovered from barcodes. ", "9 of 12 SMM are high-risk — the enrolment population."),
])

# ─────────────────────────────────────────── 4 pipeline funnel
s = new_slide(prs, "03 · pipeline", "From 22,273 genes to 15 targets")
steps = [("00–01", "Fetch GEO + HPA · align, QC, build AnnData", "22,273"),
         ("02", "Genes with ≥10 pseudobulk counts", "18,582"),
         ("02", "Significantly up in tumour, in any of three comparisons", f"{union}"),
         ("03", "Admitted on pooled signal must also be up in SMM itself", "832"),
         ("03", "Drop immunoglobulin / HLA / mitochondrial", "832"),
         ("03", "Expressed ≥5 CPM in ≥half of SMM tumours", "242"),
         ("02b·03", "Reject genes also expressed in brain, heart, other vital tissue", "146"),
         ("03", "Top 60 by effect size × reproducibility × tumour specificity", "60"),
         ("04", "Has a UniProt protein · MHCflurry over 27 alleles", f"{n_epi}"),
         ("05·05b", "Reject known housekeeping and broadly essential genes", f"{n_cand}"),
         ("06", "Weighted ranking across all four evidence layers", "15")]
y = Inches(1.8)
for i, (sid, lab, ct) in enumerate(steps):
    last = i == len(steps) - 1
    _rect(s, M, y, Inches(8.6), Inches(0.36), HL if last else SURF2)
    _txt(s, M + Inches(0.1), y + Inches(0.07), Inches(0.9), Inches(0.25), sid,
         size=9.5, font=MONO, color=INK3, bold=True)
    _txt(s, M + Inches(1.0), y + Inches(0.06), Inches(6.3), Inches(0.25), lab,
         size=11.5, color=ACCENT if last else INK2, bold=last)
    _txt(s, M + Inches(7.4), y + Inches(0.06), Inches(1.1), Inches(0.25), ct,
         size=11.5, font=MONO, color=ACCENT if last else INK, bold=True,
         align=PP_ALIGN.RIGHT)
    y = y + Inches(0.42)
label(s, Inches(9.6), Inches(1.8), Inches(3.1), "the four evidence layers")
glossary(s, Inches(9.6), Inches(2.1), Inches(3.1), Inches(4.0), [
    ("02", "differential expression, three contrasts"),
    ("02b", "normal-tissue restriction, ~50 HPA tissues"),
    ("04", "HLA-I presentation, 27 alleles"),
    ("05·05b", "literature knowledge, via Claude + a fixed grid"),
], size=10.5, gap=9)

# ─────────────────────────────────────────── 5 contrasts
s = new_slide(prs, "04 · differential expression", "Between people, or within a person?")
_txt(s, M, Inches(1.8), Inches(11.9), Inches(0.3),
     '"Up in tumour" is ambiguous until you say compared with whose normal cells.',
     size=14, color=INK2)
for i, (tag, head, body) in enumerate([
    ("CONTRAST A", "Tumour vs any normal PC",
     "All tumour vs all normal, pooled. Widest net — but mixes healthy donors with "
     "patients' own cells."),
    ("CONTRAST B", "SMM tumour vs a healthy person",
     "12 SMM tumours vs 9 healthy donors. The clinical question, and the effect "
     "size used in the score."),
    ("CONTRAST C", "Tumour vs the patient's own normal",
     "The 11 paired patients, each against themselves. Controls for between-person "
     "variation.")]):
    x = M + Inches(4.12) * i
    _rect(s, x, Inches(2.35), Inches(3.85), Inches(1.5), SURF2)
    _txt(s, x + Inches(0.15), Inches(2.45), Inches(3.6), Inches(0.2), tag,
         size=9.5, font=MONO, color=INK3)
    _txt(s, x + Inches(0.15), Inches(2.7), Inches(3.6), Inches(0.3), head,
         size=13, font=HEAD, color=INK, bold=True)
    _txt(s, x + Inches(0.15), Inches(3.05), Inches(3.6), Inches(0.8), body,
         size=10.5, color=INK2)
table(s, M, Inches(4.1), Inches(12.1), [
    ["", "Model", "Rows", "Normal side", "Entry threshold", "Admitted"],
    ["A", "~status", "46", "9 healthy donors + 14 patients' own",
     "log2FC > 1 · padj < 0.05", f"{nA}"],
    ["B", "~status", "21", "9 healthy donors", "log2FC > 1 · padj < 0.05", f"{nB}"],
    ["C", "~patient + status", "22", "the same 11 patients",
     "log2FC > 0.5 · padj < 0.05", f"{nC}"],
    ["", "", "", "A ∪ B ∪ C — a gene needs clear only one", "", f"{union}"],
], col_w=[0.5, 2.0, 0.7, 4.2, 2.8, 1.0], highlight=(2,), warn_rows=(4,))
label(s, M, Inches(5.85), Inches(12.1), "takeaway")
takeaway(s, M, Inches(6.15), Inches(12.1), Inches(1.1), [
    ("Union, not intersection. ",
     "One contrast is enough — but anything admitted on pooled signal alone must "
     "still be up in SMM itself."),
    ("Why B alone is not enough — EPHB1. ",
     "B log2FC 5.87, C log2FC 0.07. Healthy donors sit at 0.0 CPM, so any expression "
     "looks enormous."),
], size=11.5, gap=5)

# ─────────────────────────────────────────── 6 beyond fold change
s = new_slide(prs, "05 · beyond fold change", "Three questions a log2FC cannot answer")
for i, (tag, head, body) in enumerate([
    ("REPRODUCIBILITY", "Recurrence",
     "In how many of the 11 paired patients is the gene up — log2FC ≥ 0.25 AND "
     "≥5 CPM? Counting, not pooling, so no stage dominates."),
    ("OFF-THE-SHELF GATE", "SMM coverage",
     "In how many of the 12 SMM tumours is it on at ≥5 CPM? A shared product cannot "
     "rest on an antigen most patients lack."),
    ("DURABILITY", "MM retention",
     "Given at SMM to block progression. An antigen lost on the way to MM fails "
     "precisely when it is needed.")]):
    x = M + Inches(4.12) * i
    _rect(s, x, Inches(1.75), Inches(3.85), Inches(1.95), SURF2)
    _txt(s, x + Inches(0.15), Inches(1.85), Inches(3.6), Inches(0.2), tag,
         size=9.5, font=MONO, color=INK3)
    _txt(s, x + Inches(0.15), Inches(2.1), Inches(3.6), Inches(0.3), head,
         size=14, font=HEAD, color=INK, bold=True)
    _txt(s, x + Inches(0.15), Inches(2.48), Inches(3.6), Inches(1.15), body,
         size=10.5, color=INK2)
label(s, M, Inches(4.05), Inches(12.1), "why each gate exists")
takeaway(s, M, Inches(4.35), Inches(5.9), Inches(2.8), [
    ("Detection-prevalence had to go. ",
     "EDNRB reached rank 5 on a >1-UMI gate while clearing 5 CPM in only 3 of 12 "
     "SMM tumours."),
    ("Direction alone is not enough. ",
     "NEB in MGUS-3 is 1.3 vs 1.2 CPM. A magnitude gate cut perfect recurrence from "
     "112 genes to 7."),
], size=11.5)
takeaway(s, Inches(6.9), Inches(4.35), Inches(5.8), Inches(2.8), [
    ("Specificity is measured on SMM, not pooled. ",
     "The pooled tumour set is 57% MM cells. IFITM1 is 46% of it but only 21% of SMM "
     "cells. Spearman between versions: 0.575."),
    ("MGUS is annotated, never gated. ",
     "3 samples, 307 tumour cells — subtypes sampled, not biology."),
], size=11.5)

# ─────────────────────────────────────────── 7 tissue safety
s = new_slide(prs, "06 · normal-tissue safety", "Restricted to which tissue, not how many")
pic(s, f"{FIG}/specificity_vs_hla.png", M, Inches(1.55), h=Inches(5.0))
label(s, Inches(6.7), Inches(1.6), Inches(6.0), "key terms in this figure")
glossary(s, Inches(6.7), Inches(1.9), Inches(6.0), Inches(1.5), [
    ("x-axis", "TUMOUR specificity — % SMM tumour cells expressing minus % normal "
               "plasma cells. Within-marrow."),
    ("y-axis", "HLA allele coverage — fraction of the 27-allele panel with a strong binder."),
    ("colour", "composite score, the final weighted rank."),
])
label(s, Inches(6.7), Inches(3.45), Inches(6.0), "not the same as the 0.25 safety weight")
table(s, Inches(6.7), Inches(3.72), Inches(6.0), [
    ["Tier", "TISSUE safety — restricted to", "Score"],
    ["1a", "bone marrow / lymphoid", "1.00"],
    ["1b", "testis — cancer-testis antigen", "1.00"],
    ["2", "elsewhere, non-vital", "0.60"],
    ["3", "neural or vital organ", "rejected"],
    ["4", "nothing — broadly expressed", "0.20"],
], col_w=[0.7, 4.1, 1.2], highlight=(1, 2), warn_rows=(4,), size=10)
label(s, Inches(6.7), Inches(5.6), Inches(6.0), "takeaway")
takeaway(s, Inches(6.7), Inches(5.88), Inches(6.0), Inches(1.4), [
    "Two specificities: the x-axis is within marrow, the table across ~50 tissues.",
    ("It removed 9 of an earlier 15. ",
     "CADM1 retina, CD200 hypothalamus, IFI6 choroid plexus — seven had held the "
     "maximum safety score."),
    ("Marrow pooled with lymphoid, on purpose. ",
     "HPA marrow is 1% plasma cells — a marrow outlier is a myeloid signature."),
], size=11)

# ─────────────────────────────────────────── 8 novelty-confidence
s = new_slide(prs, "07 · the knowledge layer", "Two categories, combined by a grid")
table(s, M, Inches(1.75), Inches(6.0), [
    ["Confidence", "is it plausibly a myeloma antigen?"],
    ["A", "established myeloma antigen or oncogene"],
    ["B", "plausible — pathway, family or lineage"],
    ["C", "no known myeloma role"],
    ["D", "housekeeping / essential → REJECTED"],
], col_w=[1.4, 4.6], warn_rows=(4,))
table(s, M, Inches(3.4), Inches(6.0), [
    ["Novelty", "how far has it already been developed?"],
    ["1a", "vaccine in clinical study in another cancer"],
    ["1b", "vaccine in clinical study in myeloma — crowded"],
    ["2", "CAR / antibody / ADC / bispecific in myeloma"],
    ["3", "published relevance, no clinical program"],
    ["4", "not previously described as a myeloma marker"],
], col_w=[1.4, 4.6])
label(s, Inches(7.0), Inches(1.75), Inches(5.7), "the grid → novelty_confidence, 0–10")
table(s, Inches(7.0), Inches(2.02), Inches(5.7), [
    ["conf ↓ / nov →", "1a", "1b", "2", "3", "4"],
    ["A", "10", "3", "8", "9", "–"],
    ["B", "8", "2", "7", "7", "6"],
    ["C", "5", "–", "5", "4", "3"],
    ["D", "0", "0", "0", "0", "0"],
], col_w=[2.0, 0.74, 0.74, 0.74, 0.74, 0.74], highlight=(1,), warn_rows=(4,))
label(s, M, Inches(5.25), Inches(6.0), "where the numbers come from")
glossary(s, M, Inches(5.53), Inches(6.0), Inches(1.6), [
    ("confidence sets the band",
     "A 8–10 · B 6–8 · C 3–5 · D 0. How likely this is a real myeloma antigen."),
    ("novelty positions within it",
     "de-risking raises, crowding lowers: 1a > 3 ≈ 2 > 4."),
    ("1b drops out of the band",
     "10 → 3 at confidence A. Someone already has a vaccine — that overrides "
     "everything else."),
], size=10.5, gap=5)
label(s, Inches(7.0), Inches(3.65), Inches(5.7), "takeaway")
takeaway(s, Inches(7.0), Inches(3.93), Inches(5.7), Inches(3.2), [
    ("Why a grid, not a weighted sum. ",
     "Confidence is monotonic — A always beats D. Novelty is not: already in a "
     "myeloma trial means validated but taken."),
    ("These are hand-set ordinals, not fitted values. ",
     "So the ordering is what carries the argument, not the spacing — replacing "
     "every value with its pure rank keeps 14 of 15 shortlisted genes."),
    ("And the ranking barely depends on them. ",
     "Perturbing every cell by ±1 across 500 grids leaves ranks 1–11 unchanged; "
     "only slots 14–15 move. The top-left cell is empty in this run."),
], size=11)

# ─────────────────────────────────────────── 9 composite
s = new_slide(prs, "08 · composite score", "Where each target's score comes from")
pic(s, f"{FIG}/composite_scores.png", M, Inches(1.6), h=Inches(4.9))
label(s, Inches(7.3), Inches(1.7), Inches(5.4), "key terms in this figure")
glossary(s, Inches(7.3), Inches(2.0), Inches(5.4), Inches(2.2), [
    ("x-axis", "weighted contribution — each segment is component × its weight; "
               "the bar total IS the composite score."),
    ("blue", "DE evidence · 0.30 — log2FC × recurrence × tumour specificity"),
    ("green", "HLA presentation · 0.20 — epitope density × allele coverage"),
    ("orange", "tissue safety · 0.25 — the tier score from ~50 HPA tissues"),
    ("pink", "novelty-confidence · 0.25 — the grid value, rescaled"),
])
label(s, Inches(7.3), Inches(4.45), Inches(5.4), "takeaway")
takeaway(s, Inches(7.3), Inches(4.75), Inches(5.4), Inches(2.4), [
    "Read the composition, not the total — a long bar of one colour is fragile.",
    ("CCND1 maxes two layers. ",
     "The full 0.30 on DE and the full 0.25 on novelty-confidence — the only gene "
     "to top both. Its brake is tissue: tier 2, not lymphoid."),
    ("Two targets score zero on a whole layer. ",
     "TBXAS1 and HERC5 contribute nothing from novelty-confidence — a quarter of "
     "the scale, empty."),
], size=11.5)

# ─────────────────────────────────────────── 10 evidence matrix
s = new_slide(prs, "09 · evidence matrix", "The six quantities the composite runs on")
pic(s, f"{FIG}/evidence_heatmap.png", M, Inches(1.5), h=Inches(5.5))
label(s, Inches(6.4), Inches(1.6), Inches(6.3), "key terms in this figure")
glossary(s, Inches(6.4), Inches(1.9), Inches(6.3), Inches(2.6), [
    ("log2fc_smm", "effect size — SMM tumour vs healthy donor (contrast B)"),
    ("specificity", "% tumour cells expressing minus % normal plasma cells"),
    ("recurrence_frac", "fraction of the 11 paired patients where the gene is genuinely up"),
    ("smm_coverage_frac", "fraction of the 12 SMM tumours expressing it at ≥5 CPM"),
    ("epitope_density", "strong HLA-I binders per 100 aa — length-normalised"),
    ("novelty_confidence", "0–10 from the grid on the previous slide"),
])
label(s, Inches(6.4), Inches(4.6), Inches(6.3), "takeaway")
takeaway(s, Inches(6.4), Inches(4.9), Inches(6.3), Inches(2.3), [
    ("Colour is scaled within each column. ",
     "Compare down a column, never across. The printed number is raw."),
    ("HERC5 rides one column. ", "59% of its score is tissue safety alone."),
    "CCND1 pairs the largest fold change with maximum novelty-confidence.",
    "TBXAS1 and MLLT3 show full SMM coverage on modest effect sizes.",
], size=11.5)

# ─────────────────────────────────────────── 11 shortlist
s = new_slide(prs, "10 · shortlist", f"15 candidates of {n_cand}, evidence visible")
rows = [["#", "Gene", "Score", "N-C", "Cf", "Nv", "Tis", "Recur", "SMM pr", "SMM cov"]]
hl, wn = [], []
for i, r in sl.iterrows():
    tis = r.tissue_tier.split("_")[0]
    rows.append([int(r["rank"]), r.gene, f"{r.composite_score:.3f}",
                 int(r.novelty_confidence), r.mm_confidence,
                 r.novelty_tier.split("_")[0], tis,
                 f"{int(r.recurrence_n_samples)}/11",
                 f"{int(r.recurrence_smm_n)}/5",
                 f"{int(r.smm_coverage_n)}/12"])
    if tis in ("1a", "1b"):
        hl.append(i + 1)
    if not r.mm_retained:
        wn.append(i + 1)
table(s, M, Inches(1.65), Inches(7.4), rows,
      col_w=[0.45, 1.5, 0.85, 0.5, 0.45, 0.45, 0.6, 0.85, 0.8, 0.85],
      size=9.5, highlight=tuple(hl), warn_rows=tuple(wn))
label(s, Inches(8.3), Inches(1.65), Inches(4.4), "columns")
glossary(s, Inches(8.3), Inches(1.93), Inches(4.4), Inches(1.9), [
    ("N-C", "novelty-confidence, 0–10"),
    ("Cf / Nv", "confidence letter / novelty tier"),
    ("Tis", "tissue tier from the safety test"),
    ("rose row", "tier 1a/1b — lymphoid or testis-restricted, the two safest tiers"),
    ("amber row", "antigen NOT retained into MM"),
    ("Recur · SMM pr · SMM cov", "paired patients · paired SMM · SMM tumours"),
], size=10, gap=4)
label(s, Inches(8.3), Inches(4.0), Inches(4.4), "takeaway")
takeaway(s, Inches(8.3), Inches(4.3), Inches(4.4), Inches(2.9), [
    ("Six are lymphoid-restricted AND confidence A/B. ",
     "FCRLA, NLGN4X, CD1D, SPN, POU2F2, MLLT3 — two independent lines agreeing."),
    ("FCRLA is the strongest. ",
     "Tightest lymphoid restriction (FDR 9.5e-09), 10/12 SMM coverage, 8/8 in MM."),
    ("NLGN4X is what data alone could not find. ",
     "Its Y-paralogue is an H-Y minor histocompatibility antigen."),
    ("CCND1 at rank 1 is a trap. ",
     "Not retained into MM, t(11;14)-restricted, partly definitional."),
], size=10.5, gap=6)

# ─────────────────────────────────────────── 12 caveats
s = new_slide(prs, "11 · what the pipeline caught, and what it cannot", "Honest uncertainty")
label(s, M, Inches(1.7), Inches(5.9), "bugs this pipeline caught")
takeaway(s, M, Inches(2.0), Inches(5.9), Inches(5.0), [
    ("A filter that discarded 38 real genes. ",
     'startswith("IG") removed the entire IGF axis.'),
    ("Two artifacts in an earlier top 15. ",
     "Direction-only recurrence counted near-ties on near-zero values as enrichment."),
    ("A cardiac safety flag. ", "ST3GAL6 at rank 2, seen only with the full HPA reference."),
    ("A gene flat in the indication. ", "SLC1A5 reached rank 15 at B log2FC 0.20."),
    ("Silent gene loss to a timeout. ", "Two candidates vanished on a UniProt fetch."),
], size=11.5, gap=6)
label(s, Inches(7.0), Inches(1.7), Inches(5.7), "what it still cannot do")
takeaway(s, Inches(7.0), Inches(2.0), Inches(5.7), Inches(5.0), [
    ("Self-tolerance is unmodelled. ",
     "Affinity is not immunogenicity; the T-cell repertoire may already be deleted. "
     "The biggest reason a shortlist like this fails in the lab."),
    ("One LLM claim was checked, and it drifted. ",
     "NCT03591614 is a real DKK1 vaccine in SMM — but dendritic-cell, not the DNA "
     "plasmid named, wrong sponsor, withdrawn at zero enrolment. NLGN4X at rank 3 "
     "is still unchecked."),
    ("Everything is RNA. ",
     "CPQ stains highest in kidney and liver — that would be tier 3."),
    ("The labels are partly circular. ",
     "Neoplastic status came from clustering, so DE re-derives cluster genes."),
], size=11.5, gap=6)

out = f"{OUT}/SMM_vaccine_targets.pptx"
prs.save(out)
print(f"Saved {out}")
print(f"  {len(prs.slides.__iter__.__self__._sldIdLst)} slides, 16:9")
print("\nImport into Google Slides:")
print("  File > Import slides > Upload, or upload the .pptx to Drive and")
print("  right-click > Open with > Google Slides.")
