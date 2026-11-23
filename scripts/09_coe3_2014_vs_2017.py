"""STEP 9 - Subclade coe.3 sampled outside the brooding period in 2014 and 2017: are damage and copy number similar?

Why: coe.3 (historical HC-B) was sampled in both years in the Yonara Channel area:
  2014 - three sites (KR, YIS, YKS), before brooding on 13 July 2014 (three colonies per site);
  2017 - three coe.3 colonies of the subclade panel, Yonara, 13 June 2017 (about 400 m north of YIS).
All colonies of both years were amplified on the same plate for each assay (SRP54, short COI, long COI),
so their values can be compared directly. This step asks whether coe.3 colonies sampled outside the brooding
period gave similar estimates in the two years, and how much each measure varies among coe.3 groups.

Comparisons (2017 coe.3, n = 3, against):
  A. all 2014 coe.3 colonies sampled before brooding (n = 9; the reference group of Supplementary Table S2);
  B. YIS before brooding only (n = 3; the 2014 site nearest to the 2017 site).
Measures (per colony, from results/03_mtdna_per_colony.csv, log2 scale):
  damage            = -LS = log2 N0(short COI) - log2 N0(long COI); reported as lesions per 10 kb
                      (Poisson; difference x ln 2 x 10,000 / 1,124 bp)
  copy number       = mtCN = log2 N0(short COI) - log2 N0(SRP54); reported as fold difference
  short COI / input = COIs_per_input (mtDNA quantity per 100 ng DNA, no SRP54); reported as fold difference
  SRP54 / input     = COIs_per_input - mtCN (the nuclear reference); reported as fold difference
  SHC701 (SRP54 failed; 2014 YKS during brooding) is excluded from the SRP54-based measures; it is not in
  comparisons A or B, only in the group-mean summary.
Tests: Welch's t-test on the log2 values. Null hypothesis: equal means in 2017 coe.3 and the 2014 reference group.
  95% confidence intervals use the Welch-Satterthwaite df. No correction for multiple testing (exploratory).
Group-mean summary: the mean of each measure in each coe.3 group (six 2014 groups: site x before/during;
  plus 2017 coe.3) and the range of these means (largest minus smallest, as fold), with the colony SD.

Caution: year is confounded with sampling month (June 2017 v July 2014), exact location, sample storage time before
the 2018 qPCR (about 3.7 v 0.8 years) and possibly preservation. This is a consistency check, not an estimate of
change between years.

Outputs: 09_coe3_2014_vs_2017_tests.csv, 09_coe3_group_means.csv, 09_coe3_group_mean_spread.csv; a "Step 9" section in SUMMARY.md and three
sheets in all_results.xlsx (running this step again replaces them). Run after steps 1-8."""
import inspect
import math
import numpy as np
import pandas as pd
from scipy import stats
from qpcr_lib import RES, lesions_per_10kb

LONG_COI_BP = inspect.signature(lesions_per_10kb).parameters["length_bp"].default   # 1,124 bp, set in qpcr_lib
LESIONS_PER_LOG2 = math.log(2) * 10_000 / LONG_COI_BP

d = pd.read_csv(RES / "03_mtdna_per_colony.csv")
d = d[d.analysis_group.str.match(r"mt[12]_") & (d.lineage == "B")].copy()       # coe.3 colonies, main data sets
d["damage"] = -d["LS"]
d["srp54_per_input"] = d["COIs_per_input"] - d["mtCN"]
d.loc[d.exclude_copy_number.astype(bool), ["mtCN", "srp54_per_input"]] = np.nan
d["group"] = np.where(d.analysis_group.str.startswith("mt2"), "2017 coe.3 (Yonara, 13 Jun)",
                      "2014 " + d.site_code + " " + d.state2.map({"before": "before", "reproductive": "during"}))

MEASURES = [("damage", "Estimated mtDNA damage", "lesions per 10 kb"),
            ("mtCN", "Copy number (short COI / SRP54)", "fold"),
            ("COIs_per_input", "Short COI per DNA input", "fold"),
            ("srp54_per_input", "SRP54 per DNA input", "fold")]


def to_scale(v, col):
    """log2 difference -> reported scale (lesions per 10 kb for damage, fold otherwise)."""
    return v * LESIONS_PER_LOG2 if col == "damage" else 2.0 ** v


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    diff, se = a.mean() - b.mean(), math.sqrt(va + vb)
    df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
    tcrit = stats.t.ppf(0.975, df)
    p = float(stats.ttest_ind(a, b, equal_var=False).pvalue)
    return diff, diff - tcrit * se, diff + tcrit * se, df, p


y17 = d[d.group.str.startswith("2017")]
REFS = [("A: all 2014 coe.3 before brooding", d[d.group.str.endswith("before")]),
        ("B: 2014 YIS before brooding", d[d.group == "2014 YIS before"])]
rows = []
for ref_name, ref in REFS:
    for col, label, unit in MEASURES:
        a, b = y17[col].dropna(), ref[col].dropna()
        diff, lo, hi, df, p = welch(a, b)
        rows.append(dict(comparison=f"2017 coe.3 v {ref_name[3:]}", reference=ref_name[0], measure=label,
                         n_2017=len(a), n_2014=len(b), log2_difference=diff, df_welch=df, p_welch=p,
                         difference=to_scale(diff, col), ci95_low=to_scale(lo, col), ci95_high=to_scale(hi, col),
                         unit=("lesions per 10 kb (2017 minus 2014)" if col == "damage"
                               else "fold (2017 / 2014)")))
T = pd.DataFrame(rows)
T.to_csv(RES / "09_coe3_2014_vs_2017_tests.csv", index=False)

order = ["2014 KR before", "2014 KR during", "2014 YIS before", "2014 YIS during", "2014 YKS before",
         "2014 YKS during", "2017 coe.3 (Yonara, 13 Jun)"]
G = d.groupby("group")[[m[0] for m in MEASURES]].agg(["mean", "count"]).reindex(order)
G.columns = [f"{c}_{s}" for c, s in G.columns]
G = G.reset_index()
G.to_csv(RES / "09_coe3_group_means.csv", index=False)

spread = []
for col, label, unit in MEASURES:
    m = G[f"{col}_mean"]
    rng_log2 = m.max() - m.min()
    sd_within = d.groupby("group")[col].std().mean()          # mean within-group colony SD, log2
    spread.append(dict(measure=label, range_of_group_means_log2=rng_log2,
                       range_of_group_means=(f"{rng_log2 * LESIONS_PER_LOG2:.1f} lesions per 10 kb" if col == "damage"
                                             else f"{2 ** rng_log2:.2f}-fold"),
                       mean_within_group_colony_sd_log2=sd_within))
S = pd.DataFrame(spread)
S.to_csv(RES / "09_coe3_group_mean_spread.csv", index=False)

# ------------------------------------------------------------------ add to SUMMARY.md and all_results.xlsx
def fmt(r):
    if r.measure.startswith("Estimated"):
        return f"{r.difference:+.1f} ({r.ci95_low:+.1f} to {r.ci95_high:+.1f})"
    return f"{r.difference:.2f} ({r.ci95_low:.2f} to {r.ci95_high:.2f})"


lines = ["", "## Step 9 - coe.3 outside the brooding period, 2014 v 2017", "",
         "2017 coe.3 (n = 3) against 2014 coe.3 sampled before brooding. Welch's t-tests on log2 values; null "
         "hypothesis: equal means. Damage as lesions per 10 kb (2017 minus 2014); other measures as fold (2017 / 2014); "
         "95% CI in brackets. Year is confounded with month, exact location and sample storage time.", "",
         "| Reference | Measure | Difference (95% CI) | p |", "|---|---|---|---|"]
for _, r in T.iterrows():
    lines.append(f"| {r.reference} ({r.n_2014}) | {r.measure} | {fmt(r)} | {r.p_welch:.3f} |")
lines += ["", "Range of group means across the seven coe.3 groups (six 2014 site x period groups and 2017):", "",
          "| Measure | Range of group means | Mean colony SD within groups (log2) |", "|---|---|---|"]
for _, r in S.iterrows():
    lines.append(f"| {r.measure} | {r.range_of_group_means} | {r.mean_within_group_colony_sd_log2:.2f} |")
lines.append("")
summ = RES / "SUMMARY.md"
text = summ.read_text(encoding="utf-8")
start = text.find("\n## Step 9 - ")
if start >= 0:                                   # replace only this step's section
    nxt = text.find("\n## Step ", start + 1)
    text = text[:start] + (text[nxt:] if nxt >= 0 else "")
summ.write_text(text.rstrip("\n") + "\n" + "\n".join(lines), encoding="utf-8")
xl = RES / "all_results.xlsx"
if xl.exists():
    with pd.ExcelWriter(xl, mode="a", engine="openpyxl", if_sheet_exists="replace") as w:
        T.round(4).to_excel(w, sheet_name="09_coe3_2014_vs_2017", index=False)
        G.round(4).to_excel(w, sheet_name="09_coe3_group_means", index=False)
        S.round(4).to_excel(w, sheet_name="09_coe3_spread", index=False)

pd.set_option("display.width", 200)
print(T[["reference", "measure", "n_2017", "n_2014", "difference", "ci95_low", "ci95_high", "p_welch"]]
      .round(3).to_string(index=False))
print(S.round(3).to_string(index=False))
print("STEP 9 done: 09_coe3_2014_vs_2017_tests.csv, 09_coe3_group_means.csv, 09_coe3_group_mean_spread.csv")
