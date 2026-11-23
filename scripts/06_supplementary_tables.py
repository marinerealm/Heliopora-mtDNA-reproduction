"""STEP 6 - Supplementary Tables S1 (primers) and S2 (exploratory mtDNA data) -> results/Supplementary_Tables_S1_S2.xlsx"""
import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
from qpcr_lib import RAW, META, RES

# ------------------------------------------------------------------ Table S1: primers
seq = {  # primer sequences, 5'-3' (from the laboratory primer list)
    "actB-F": "TCTGGTATGTGCAAAGCTGG", "actB-R": "CCATATCATCCCAGTTCGTGAC",
    "rpl12-F": "ATGCCTCCTAAGTTCGATCC", "rpl12-R": "GAACAGTCACTTTGAGTCCC",
    "srp54-F": "GTTACACAGAAATGGATCCTGT", "srp54-R": "CTCGCCTGAGATTCACAAGC",
    "d2rF": "TCTTTCGGTGAATGACTCTG", "d2rR": "AACGCCAGTAAGACAAGAC",
    "cry1BF": "GGTTCTATCTGCTGAAGGAG", "cry1BR": "CCAAGGTATGATACAGCAGG",
    "mtMutSF": "TGACAGAGACAAAGTAGCTCG", "mtMutSR": "CTAGCATTAGGTTCAGAGGG",
    "COI_smlF": "TAATTCTGCCGGGATTTGG", "COI_octoF": "CCATAACAGGACTAGCAGCATC", "COI_octoR": "ATCATAGCATAGACCATACC",
    "sodF1": "GATGCAGAAAGGCACGTTGG", "sodR": "TGCACCACAATTGAACGACCA",
    "catF": "GGTGAATCTGGTAGTGCTGA", "catR": "GTGATGAAGTCCCAAAACATGT",
    "gpxF": "TTTGGTAATCAGGAACCAGAT", "gpxR": "CTTBGCAAAATTCCACTTTAT",
    "hsp70F": "ATCATCGCGAATGACCAAGG", "hsp70R": "GCATCTGTGACCTTCTTGCC",
    "hsp90F": "GGAAATTTTCCTCAGGGAACTC", "hsp90R": "AACCCAACACCAAATTGACC",
}
ew = pd.read_csv(RES / "01_expression_wells.csv").groupby("gene").E_mean.mean()
dw = pd.read_csv(RES / "02_mtdna_wells.csv").groupby("assay").E_mean.first()
D2R_NOTE = ("D2R-like: partial class A G-protein-coupled receptor (TSA IABP01020958). Its closest matches are coral genes "
            "annotated as dopamine or octopamine receptor-like (automated predictions). No D2-type dopamine receptor has been "
            "described in corals; the only coral dopamine receptor sequences reported are D1-like (Isomura et al. 2013).")
COI_NOTE = "Same primer pair used for expression and COI-small amplification"
S1 = [  # gene, product, used on, forward, reverse, amplicon (bp), efficiency, note
    ("ACTB", "Beta-actin", "cDNA (reference); DNA (replicate run, nuclear reference)", "actB-F", "actB-R", 107, f"{ew['ACTB']:.3f} (cDNA); 1.856 (DNA)", ""),
    ("RPL12", "Ribosomal protein L12", "cDNA (reference)", "rpl12-F", "rpl12-R", 142, f"{ew['RPL12']:.3f}", ""),
    ("SRP54", "Signal recognition particle 54", "cDNA (rejected reference); DNA (nuclear reference)", "srp54-F", "srp54-R", 184, f"{ew['SRP54']:.3f} (cDNA); {dw['SRP54']:.3f} (DNA)", ""),
    ("D2R", "Dopamine receptor 2-like", "cDNA", "d2rF", "d2rR", 150, f"{ew['D2R']:.3f}", D2R_NOTE),
    ("CRY1", "Cryptochrome 1", "cDNA", "cry1BF", "cry1BR", 150, f"{ew['CRY1']:.3f}", ""),
    ("mtMutS", "Mitochondrial MutS homologue", "cDNA", "mtMutSF", "mtMutSR", 149, f"{ew['mtMutS']:.3f}", ""),
    ("COI", "Cytochrome c oxidase I", "cDNA", "COI_smlF", "COI_octoR", 96, f"{ew['COI']:.3f}", COI_NOTE),
    ("SOD", "Superoxide dismutase", "cDNA", "sodF1", "sodR", 134, f"{ew['SOD']:.3f}", ""),
    ("CAT", "Catalase", "cDNA", "catF", "catR", 209, f"{ew['CAT']:.3f}", ""),
    ("GPX", "Glutathione peroxidase", "cDNA", "gpxF", "gpxR", 183, f"{ew['GPX']:.3f}", ""),
    ("HSP70", "Heat-shock protein 70", "cDNA", "hsp70F", "hsp70R", 344, f"{ew['HSP70']:.3f}", ""),
    ("HSP90", "Heat-shock protein 90", "cDNA", "hsp90F", "hsp90R", 291, f"{ew['HSP90']:.3f}", ""),
    ("COI-small", "mtDNA, 96 bp", "DNA (copy number; damage reference)", "COI_smlF", "COI_octoR", 96, f"{dw['COI-Sml']:.3f}", COI_NOTE),
    ("COI-large", "mtDNA, 1,124 bp (COII-igr-COI region)", "DNA (damage)", "COI_octoF", "COI_octoR", 1124, f"{dw['COI-Lar']:.3f}", ""),
]
rows1 = [[g, fn, used, f, seq[f], r, seq[r], bp, e, note] for g, fn, used, f, r, bp, e, note in S1]
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
      "Amplification efficiencies were measured on the plates of this study. Values are regenerated by "
      "scripts/06_supplementary_tables.py in Supplementary Data S1.",
      H1, rows1, [12, 26, 30, 16, 26, 18, 26, 12, 22, 34])
ws2 = wb.create_sheet("Table S2")
sheet(ws2, "Table S2. Exploratory mtDNA data not used in the main tests. Copy number and lesions are expressed relative to the 2014 "
           "coe.3 colonies sampled before brooding (13 July 2014; all three sites pooled). Subclades follow Taninaka et al. (2021).",
      "Subclade, site, collection year and DNA storage are confounded in these groups, so no tests were run. High apparent lesion "
      "frequencies in archived DNA most likely reflect some degradation during storage. Copy number is blank where SRP54 failed. "
      "Values are regenerated by scripts/06_supplementary_tables.py in Supplementary Data S1.",
      H2, rows2, [34, 10, 16, 8, 14, 12, 22, 24, 40])
out = RES / "Supplementary_Tables_S1_S2.xlsx"
wb.save(out)
print(f"STEP 6 done: {out.name} (Table S1: {len(rows1)} assays; Table S2: {len(rows2)} colonies)")
