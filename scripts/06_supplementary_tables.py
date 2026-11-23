"""STEP 6 - Supplementary Tables S1 (primers) and S2 (exploratory mtDNA data) -> results/Supplementary_Tables_S1_S2.xlsx"""
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
from qpcr_lib import RAW, META, RES

# ------------------------------------------------------------------ Table S1: primers
pr = pd.read_excel(RAW / "sample_records" / "Primers.xlsx", header=0)
pr.columns = ["n", "primer", "sequence", "gene_name", "length"]
seq = dict(zip(pr.primer.str.strip(), pr.sequence.str.strip()))
ew = pd.read_csv(RES / "01_expression_wells.csv").groupby("gene").E_mean.mean()
dw = pd.read_csv(RES / "02_mtdna_wells.csv").groupby("assay").E_mean.first()
S1 = [  # gene, function, used on, forward, reverse, efficiency, note
    ("ACTB", "Beta-actin", "cDNA (reference); DNA (replicate run, nuclear reference)", "actB-F", "actB-R", f"{ew['ACTB']:.3f} (cDNA); 1.856 (DNA)", ""),
    ("RPL12", "Ribosomal protein L12", "cDNA (reference)", "rpl12-F", "rpl12-R", f"{ew['RPL12']:.3f}", ""),
    ("SRP54", "Signal recognition particle 54", "cDNA (rejected reference); DNA (nuclear reference)", "srp54-F", "srp54-R", f"{ew['SRP54']:.3f} (cDNA); {dw['SRP54']:.3f} (DNA)", "Intron-containing amplicon on genomic DNA [confirm]"),
    ("D2R", "Dopamine receptor 2", "cDNA", "d2rF", "d2rR", f"{ew['D2R']:.3f}", ""),
    ("CRY1", "Cryptochrome 1", "cDNA", "cry1AF or cry1BF", "cry1AR or cry1BR", f"{ew['CRY1']:.3f}", "Primer pair used (A or B) [confirm]"),
    ("mtMutS", "Mitochondrial mismatch repair", "cDNA", "mtMutSF", "mtMutSR", f"{ew['mtMutS']:.3f}", ""),
    ("COI", "Cytochrome c oxidase I", "cDNA", "[confirm]", "[confirm]", f"{ew['COI']:.3f}", "Primer pair used for expression [confirm]"),
    ("SOD", "Superoxide dismutase", "cDNA", "sodF1 or sodF2", "sodR", f"{ew['SOD']:.3f}", "Forward primer used [confirm]"),
    ("CAT", "Catalase", "cDNA", "catF", "catR", f"{ew['CAT']:.3f}", ""),
    ("GPX", "Glutathione peroxidase", "cDNA", "gpxF", "gpxR", f"{ew['GPX']:.3f}", ""),
    ("HSP70", "Heat-shock protein 70", "cDNA", "hsp70F", "hsp70R", f"{ew['HSP70']:.3f}", ""),
    ("HSP90", "Heat-shock protein 90", "cDNA", "hsp90F", "hsp90R", f"{ew['HSP90']:.3f}", ""),
    ("COI-small", "mtDNA, ~100 bp", "DNA (copy number; damage reference)", "COI_smlF", "[reverse primer missing from list]", f"{dw['COI-Sml']:.3f}", ""),
    ("COI-large", "mtDNA, 1,124 bp (COII-igr-COI region)", "DNA (damage)", "[confirm]", "[confirm]", f"{dw['COI-Lar']:.3f}", "Possibly COI_octoF/R [confirm]"),
]
rows1 = []
for g, fn, used, f, r, e, note in S1:
    fs = " / ".join(seq.get(x.strip(), "") for x in f.split(" or ")) if "[" not in f else ""
    rs = " / ".join(seq.get(x.strip(), "") for x in r.split(" or ")) if "[" not in r else ""
    rows1.append([g, fn, used, f, fs, r, rs, "[to add]", e, note])
H1 = ["Gene / assay", "Product", "Used on", "Forward primer", "Forward sequence (5'-3')", "Reverse primer",
      "Reverse sequence (5'-3')", "Amplicon (bp)", "Efficiency (amplification factor)", "Note"]

# ------------------------------------------------------------------ Table S2: exploratory mtDNA data
T = pd.read_csv(RES / "03_mtdna_per_colony.csv")
M = pd.read_csv(META / "samples.csv").set_index("sample")
x = T[T.analysis_group.isin(["mt3_lineageA_2014_YKS", "mt4_archived_DNA"])].copy()
x["group"] = x.analysis_group.map({"mt3_lineageA_2014_YKS": "sp.3, YKS, after gamete release (17 Jul 2014)", "mt4_archived_DNA": "Archived DNA (2007-2011)"})
x["site"] = x["sample"].map(M.site); x["subclade"] = x["sample"].map(M.subclade); x["date"] = x["sample"].map(M.date).fillna("unknown")
x["DNA"] = x["sample"].map(M.DNA_stock_ng_per_ul); x["note"] = x["sample"].map(M.notes).fillna("")
x = x.sort_values(["group", "site", "sample"])
HIST = {"A": "HC-A", "A'": "HC-A′", "B": "HC-B"}
rows2 = [[r.group, r["sample"], r.site, f"{r.subclade} ({HIST[r.lineage]})", r.date, r.DNA,
          "" if pd.isna(r.mtCN_fold_vs_B2014_before) else round(r.mtCN_fold_vs_B2014_before, 2),
          round(r.lesions_poisson_vs_B2014_before, 1), r.note] for _, r in x.iterrows()]
H2 = ["Group", "Colony", "Site", "Subclade (historical label)", "Collection date", "DNA stock (ng/ul)",
      "mtDNA copy number (fold vs 2014 coe.3 before brooding)", "Lesions per 10 kb (Poisson, vs 2014 coe.3 before brooding)", "Note"]

# ------------------------------------------------------------------ write workbook
wb = Workbook()
thin = Side(style="thin", color="999999")
def sheet(ws, caption, footer, header, rows, widths):
    ws.append([caption]); ws["A1"].font = Font(bold=True, size=11); ws["A1"].alignment = Alignment(wrap_text=True, vertical="top")
    ws.merge_cells(start_row=1, start_column=1, end_row=1, end_column=len(header)); ws.row_dimensions[1].height = 48
    ws.append([]); ws.append(header)
    for c in ws[3]:
        c.font = Font(bold=True); c.fill = PatternFill("solid", fgColor="E7E6E6")
        c.alignment = Alignment(wrap_text=True, vertical="top"); c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
    for r in rows:
        ws.append(r)
        for c in ws[ws.max_row]:
            c.alignment = Alignment(wrap_text=True, vertical="top"); c.border = Border(top=thin, bottom=thin, left=thin, right=thin)
    ws.append([]); ws.append([footer])
    ws.merge_cells(start_row=ws.max_row, start_column=1, end_row=ws.max_row, end_column=len(header))
    ws.cell(ws.max_row, 1).alignment = Alignment(wrap_text=True, vertical="top"); ws.row_dimensions[ws.max_row].height = 60
    for i, w in enumerate(widths, 1):
        ws.column_dimensions[get_column_letter(i)].width = w
    ws.freeze_panes = "A4"

ws = wb.active; ws.title = "Table S1"
sheet(ws, "Table S1. Primers used in this study, with amplification efficiencies (LinRegPCR amplicon means for cDNA plates; "
          "LinRegPCR-style estimates from raw curves for the SRP54 and COI-small DNA plates; LinRegPCR for COI-large).",
      "Items marked [confirm] or [to add] must be completed from laboratory records. Primers for cry2, cycle and clock appear in the "
      "laboratory primer list but no qPCR data for these genes were found and they are not part of this study.",
      H1, rows1, [12, 26, 30, 16, 26, 18, 26, 12, 22, 34])
ws2 = wb.create_sheet("Table S2")
sheet(ws2, "Table S2. Exploratory mtDNA data not used in the main tests. Copy number and lesions are expressed relative to the 2014 "
           "coe.3 colonies sampled before brooding (13 July 2014; all three sites pooled). Subclades follow Taninaka et al. (2021).",
      "Subclade, site, collection year and DNA storage are confounded in these groups, so no tests were run. High apparent lesion "
      "frequencies in archived DNA most likely reflect degradation during storage. Copy number is blank where SRP54 failed. "
      "Values are regenerated by scripts/06_supplementary_tables.py in Supplementary Data S1.",
      H2, rows2, [34, 10, 16, 8, 14, 12, 22, 24, 40])
out = RES / "Supplementary_Tables_S1_S2.xlsx"
wb.save(out)
print(f"STEP 6 done: {out.name} (Table S1: {len(rows1)} assays; Table S2: {len(rows2)} colonies)")
