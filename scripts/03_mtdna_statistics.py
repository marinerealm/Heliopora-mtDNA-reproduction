"""STEP 3 - mtDNA copy number and damage: colony values and tests (Figure 2 is drawn by 07_figures.py).

Colony value = mean of duplicate wells on the log2 scale.
  log2 mtCN (SRP54-normalised) = log2N0(COI-small) - log2N0(SRP54)          [relative values only]
  log2 COI-small per input     = log2N0(COI-small)                          [SRP54-free check; equal DNA input]
  LS                           = log2N0(COI-large) - log2N0(COI-small)
  r                            = 2^(LS - mean LS of the reference group)
  lesions/10 kb (Poisson)      = -ln(r) * 10,000 / 1,124
  lesions/10 kb (Rothfuss)     = (1 - r) * 10,000 / 1,124
Reference group: before-brooding colonies of the same site (2014) or coe.3 (2017 panel).
Labels: lineage A = subclade sp.3, A' = sp.2, B = coe.3 (Taninaka et al. 2021); 'reproductive' = during brooding.

Null hypotheses (all two-sided):
  within site  : equal mean before vs during the reproductive period (Welch t)
  two-way      : no site, state or site x state effect (OLS, type II ANOVA)
  KR vs YIS    : equal change at KR and YIS (interaction contrast; same date and stage)
  subclades    : equal means for A (sp.3), A' (sp.2) and B (coe.3) (one-way ANOVA; Tukey HSD pairwise)
Copy number excludes colonies whose SRP54 failed (samples.csv: exclude_copy_number). Damage keeps them (no SRP54).
"""
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm
from qpcr_lib import META, RES, lesions_per_10kb

W = pd.read_csv(RES / "02_mtdna_wells.csv")
M = pd.read_csv(META / "samples.csv").set_index("sample")
M = M[M.in_mtDNA]

rows = {}
for col in ["log2N0", "log2N0_E2"]:
    s = W.groupby(["sample", "assay"])[col].mean().unstack()
    suf = "" if col == "log2N0" else "_E2"
    rows[f"mtCN{suf}"] = s["COI-Sml"] - s["SRP54"]
    rows[f"COIs_per_input{suf}"] = s["COI-Sml"]
    rows[f"LS{suf}"] = s["COI-Lar"] - s["COI-Sml"]
T = pd.DataFrame(rows).join(M[["analysis_group", "site_code", "lineage", "state2", "repro_state", "sex", "exclude_copy_number"]])

b14 = T[T.analysis_group.str.startswith("mt1")].copy()
y17 = T[T.analysis_group.str.startswith("mt2")].copy()

# lesions relative to the reference group
for suf in ["", "_E2"]:
    T[f"lesions_poisson{suf}"] = np.nan; T[f"lesions_rothfuss{suf}"] = np.nan
    for site in ["KR", "YIS", "YKS"]:
        idx = b14[b14.site_code == site].index
        ref = b14.loc[idx][b14.loc[idx].state2 == "before"][f"LS{suf}"].mean()
        p_, r_ = lesions_per_10kb(T.loc[idx, f"LS{suf}"] - ref)
        T.loc[idx, f"lesions_poisson{suf}"], T.loc[idx, f"lesions_rothfuss{suf}"] = p_, r_
    idx = y17.index
    ref = y17[y17.lineage == "B"][f"LS{suf}"].mean()
    p_, r_ = lesions_per_10kb(T.loc[idx, f"LS{suf}"] - ref)
    T.loc[idx, f"lesions_poisson{suf}"], T.loc[idx, f"lesions_rothfuss{suf}"] = p_, r_
# exploratory: everything relative to the 2014 coe.3 (lineage B) before-brooding colonies
ref_all = T[T.analysis_group.str.startswith("mt1") & (T.state2 == "before")].LS.mean()
T["lesions_poisson_vs_B2014_before"] = lesions_per_10kb(T.LS - ref_all)[0]
refcn = T[T.analysis_group.str.startswith("mt1") & (T.state2 == "before") & ~T.exclude_copy_number].mtCN.mean()
T["mtCN_fold_vs_B2014_before"] = np.where(T.exclude_copy_number, np.nan, 2 ** (T.mtCN - refcn))
T.round(4).to_csv(RES / "03_mtdna_per_colony.csv")
b14 = T[T.analysis_group.str.startswith("mt1")].copy()
y17 = T[T.analysis_group.str.startswith("mt2")].copy()

tests = []


def add(analysis, contrast, measure, n, diff=np.nan, lo=np.nan, hi=np.nan, p=np.nan, df=np.nan):
    tests.append(dict(analysis=analysis, contrast=contrast, measure=measure, n=n, diff=diff, CI_lo=lo, CI_hi=hi, p=p, df=df))


def welch(a, b):
    r = stats.ttest_ind(b, a, equal_var=False); ci = r.confidence_interval()
    return np.mean(b) - np.mean(a), ci.low, ci.high, r.pvalue, r.df


measures_cn = ["mtCN", "COIs_per_input", "mtCN_E2"]
measures_dmg = ["lesions_poisson", "lesions_rothfuss", "lesions_poisson_E2"]
# ---------------------------------------------------------------- within-site, B 2014
NEEDS_SRP54 = ["mtCN", "mtCN_E2"]          # only these use SRP54; other measures keep all colonies
for m in measures_cn + measures_dmg:
    d0 = b14[~b14.exclude_copy_number] if m in NEEDS_SRP54 else b14
    for site in ["KR", "YIS", "YKS"]:
        d = d0[d0.site_code == site]
        a, b_ = d[d.state2 == "before"][m].values, d[d.state2 == "reproductive"][m].values
        add("B2014 reproductive vs before", site, m, f"{len(a)} v {len(b_)}", *welch(a, b_))
# ---------------------------------------------------------------- two-way + KR vs YIS (copy number)
anova_out = []
for m in measures_cn + ["lesions_poisson"]:
    d0 = b14[~b14.exclude_copy_number] if m in NEEDS_SRP54 else b14
    fit = smf.ols(f"{m} ~ C(site_code, Treatment('KR')) * C(state2, Treatment('before'))", data=d0).fit()
    a2 = anova_lm(fit, typ=2); a2.insert(0, "measure", m); anova_out.append(a2)
    sub = d0[d0.site_code.isin(["KR", "YIS"])]
    fi = smf.ols(f"{m} ~ C(site_code) * C(state2, Treatment('before'))", data=sub).fit()
    k = [c for c in fi.params.index if ":" in c][0]
    add("B2014 KR vs YIS difference in change", "YIS change - KR change", m, "3/3/3/3", fi.params[k], *fi.conf_int().loc[k], fi.pvalues[k], fi.df_resid)
    c0 = d0[d0.state2 == "before"]
    if not m.startswith("lesions"):          # lesions are zero-centred on each site's before mean by construction
      add("B2014 baseline sites (13 Jul)", "KR = YIS = YKS", m, "3/3/3", p=stats.f_oneway(*[c0[c0.site_code == s][m] for s in ["KR", "YIS", "YKS"]]).pvalue)
pd.concat(anova_out).round(4).to_csv(RES / "03_two_way_anova.csv")
# ---------------------------------------------------------------- subclades (lineages A, A', B), Yonara 2017
for m in measures_cn + measures_dmg:
    d = y17[~y17.exclude_copy_number] if m in NEEDS_SRP54 else y17
    g = [d[d.lineage == l][m].values for l in ["A", "A'", "B"]]
    add("Yonara 2017 lineages", "A = A' = B (one-way ANOVA)", m, "3/3/3", p=stats.f_oneway(*g).pvalue, df=6)
    tk = stats.tukey_hsd(*g); ci = tk.confidence_interval()
    for i, j, lab in [(0, 2, "A - B"), (1, 2, "A' - B"), (0, 1, "A - A'")]:
        add("Yonara 2017 lineages (Tukey HSD)", lab, m, "3 v 3", g[i].mean() - g[j].mean(), ci.low[i, j], ci.high[i, j], tk.pvalue[i, j])
R = pd.DataFrame(tests)
R.round(4).to_csv(RES / "03_mtdna_tests.csv", index=False)

# ---------------------------------------------------------------- robustness: leave-one-colony-out (KR vs YIS copy number)
sub = b14[~b14.exclude_copy_number & b14.site_code.isin(["KR", "YIS"])]
loo = []
for s in sub.index:
    fi = smf.ols("mtCN ~ C(site_code) * C(state2, Treatment('before'))", data=sub.drop(s)).fit()
    k = [c for c in fi.params.index if ":" in c][0]
    loo.append(dict(dropped=s, diff=fi.params[k], p=fi.pvalues[k]))
pd.DataFrame(loo).round(4).to_csv(RES / "03_loo_KR_vs_YIS_copy_number.csv", index=False)
# robustness: damage using each COI-large duplicate alone
s = W.groupby(["sample", "assay"]).log2N0.mean().unstack()
single = []
for rep in [1, 2]:
    lar = W[(W.assay == "COI-Lar") & (W.rep == rep)].set_index("sample").log2N0
    ls = lar - s["COI-Sml"]
    for site in ["KR", "YIS", "YKS"]:
        d = b14[b14.site_code == site]
        c, t = [x for x in d.index if d.loc[x, "state2"] == "before"], [x for x in d.index if d.loc[x, "state2"] == "reproductive"]
        les = pd.Series(lesions_per_10kb(ls - ls[c].mean())[0], index=ls.index)
        single.append(dict(COI_large_well=rep, contrast=site, diff=les[t].mean() - les[c].mean(), p=stats.ttest_ind(les[t], les[c], equal_var=False).pvalue))
    g = [ls[y17[y17.lineage == l].index] for l in ["A", "A'", "B"]]
    single.append(dict(COI_large_well=rep, contrast="Yonara lineages ANOVA", p=stats.f_oneway(*g).pvalue))
pd.DataFrame(single).round(4).to_csv(RES / "03_damage_single_duplicate_check.csv", index=False)
# SRP54 vs COI-small precision across colonies
ok = T[~T.exclude_copy_number]
s2 = s.loc[ok.index]
pd.DataFrame({"SD_across_colonies_log2": [s2["SRP54"].std(), s2["COI-Sml"].std()]}, index=["SRP54", "COI-small"]).round(3).to_csv(RES / "03_assay_spread.csv")

# ---------------------------------------------------------------- group summary (incl. exploratory groups)
G = T.groupby(["analysis_group", "site_code", "lineage", "repro_state"]).agg(
    n=("LS", "size"), mtCN_fold_mean=("mtCN_fold_vs_B2014_before", "mean"), mtCN_n=("mtCN_fold_vs_B2014_before", "count"),
    lesions_vs_B2014_before_mean=("lesions_poisson_vs_B2014_before", "mean"),
    lesions_min=("lesions_poisson_vs_B2014_before", "min"), lesions_max=("lesions_poisson_vs_B2014_before", "max"))
G.round(3).to_csv(RES / "03_group_summary.csv")

# Figure 2 is drawn by 07_figures.py

print("STEP 3 done. Key tests:")
key = R[R.measure.isin(["mtCN", "COIs_per_input", "lesions_poisson"])]
print(key[["analysis", "contrast", "measure", "n", "diff", "p"]].round(3).to_string(index=False))
