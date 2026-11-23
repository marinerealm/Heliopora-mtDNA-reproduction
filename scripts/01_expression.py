"""STEP 1 - Gene expression, Kohama Reef, before (13 Jul 2014) vs during brooding (25 Jul 2014).

Inputs : data/raw/gene_expression/19032018_plate1_LinRegPCR.xls  (ACTB, RPL12, D2R, mtMutS, CRY1, COI)
         data/raw/gene_expression/19032018B_plate2_LinRegPCR.xls (SRP54, CAT, GPX, HSP70, HSP90, SOD, ACTB)
         data/metadata/samples.csv  (C1-C3 = SHC417-419, T1-T3 = SHC621-623, sex)
Outputs: results/01_* tables (Figure 2 is drawn by 07_figures.py)

Design: 3 colonies per group (independent colonies), technical duplicates.
Unit of analysis = colony. Duplicates are averaged on the log2 scale BEFORE any test.
Null hypothesis: equal mean normalised log2 expression before and during brooding.
Main test: two-sided Student t (df = 4); Benjamini-Hochberg across the 9 targets.
"""
import numpy as np
import pandas as pd
from scipy import stats
from qpcr_lib import RAW, META, RES, read_linregpcr, perm_p

TARGETS = ["SOD", "D2R", "mtMutS", "GPX", "HSP70", "HSP90", "CAT", "CRY1", "COI"]
REFS = ["ACTB", "RPL12"]
C, T = ["C1", "C2", "C3"], ["T1", "T2", "T3"]

# ---------------------------------------------------------------- 1. load wells
rows = []
for plate, f in [("plate1", "19032018_plate1_LinRegPCR.xls"), ("plate2", "19032018B_plate2_LinRegPCR.xls")]:
    o = read_linregpcr(RAW / "gene_expression" / f)
    for _, r in o.iterrows():
        lab, gene = str(r["name"]).split("_", 1)
        gene = {"d2r": "D2R", "MutS": "mtMutS", "cry1": "CRY1"}.get(gene, gene.upper())
        rows.append(dict(plate=plate, label=lab, gene=gene, Cq=float(r["Cq"]), E_ind=float(r["indiv PCR eff"]),
                         E_mean=float(r["mean_PCR_eff"]), thr=float(r["threshold"]),
                         baseline_error=" 2 " in f" {str(r['Quality_checks'])} "))
W = pd.DataFrame(rows)
W["rep"] = W.groupby(["plate", "label", "gene"]).cumcount() + 1
W["log2N0"] = np.log2(W.thr) - W.Cq * np.log2(W.E_mean)          # LinRegPCR amplicon-mean efficiency (Ruijter et al. 2009)
W["log2N0_ind"] = np.log2(W.thr) - W.Cq * np.log2(W.E_ind)
W["log2N0_E2"] = -W.Cq                                            # 100 % efficiency (sensitivity check)

# Technical outlier rule: a duplicate pair differing by > 1 cycle where a third measurement of the same cDNA
# (ACTB was run on both plates) identifies which well is aberrant. Only one well meets it:
W["outlier"] = (W.gene == "ACTB") & (W.plate == "plate1") & (W.label == "C1") & (W.rep == 1)
W.to_csv(RES / "01_expression_wells.csv", index=False)

dup = W.pivot_table(index=["gene", "plate", "label"], columns="rep", values="Cq")
dup["abs_diff"] = (dup[1] - dup[2]).abs()
dup.round(3).to_csv(RES / "01_expression_duplicate_agreement.csv")


def colony_table(col="log2N0", drop_outlier=True):
    x = W[~W.outlier] if drop_outlier else W
    x = x.assign(v=x[col] - x.groupby(["gene", "plate"])[col].transform("mean"))   # centre each gene@plate
    return x.groupby(["label", "gene"])["v"].mean().unstack()                         # ACTB = mean of both plates


S = colony_table()
meta = pd.read_csv(META / "samples.csv").dropna(subset=["expression_label"]).set_index("expression_label")

# ---------------------------------------------------------------- 2. reference gene stability
def genorm_M(X):
    return pd.Series({j: np.mean([np.std(X[j] - X[k], ddof=1) for k in X.columns if k != j]) for j in X.columns})

ref = pd.DataFrame({"geNorm_M": genorm_M(S[["ACTB", "RPL12", "SRP54"]])})
ref["raw_log2_T_minus_C"] = [S.loc[T, g].mean() - S.loc[C, g].mean() for g in ref.index]
ref["raw_p_welch"] = [stats.ttest_ind(S.loc[T, g], S.loc[C, g], equal_var=False).pvalue for g in ref.index]
ref["amplicon_efficiency"] = [W[W.gene == g].E_mean.mean() for g in ref.index]
ref["median_duplicate_dCq"] = [dup.loc[g]["abs_diff"].median() for g in ref.index]
ref["decision"] = ["use (both plates averaged)", "use", "reject: low efficiency, late Cq, noisy duplicates, drifts vs ACTB/RPL12"]
ref.round(3).to_csv(RES / "01_reference_genes.csv")


# ---------------------------------------------------------------- 3. relative expression
def analyse(S, refs, test="student"):
    NF = S[refs].mean(axis=1)
    out = []
    for g in TARGETS + [r for r in ["ACTB", "RPL12", "SRP54"] if r not in refs]:
        y = S[g] - NF
        c, t = y.loc[C].values, y.loc[T].values
        r = stats.ttest_ind(t, c, equal_var=(test == "student"))
        ci = r.confidence_interval()
        out.append(dict(gene=g, role="target" if g in TARGETS else "reference candidate", log2FC=t.mean() - c.mean(),
                        FC=2 ** (t.mean() - c.mean()), CI95_log2_lo=ci.low, CI95_log2_hi=ci.high, FC_CI_lo=2 ** ci.low,
                        FC_CI_hi=2 ** ci.high, t=r.statistic, df=r.df, p=r.pvalue, p_perm=perm_p(c, t)))
    o = pd.DataFrame(out)
    m = o.role == "target"
    o.loc[m, "q_BH"] = stats.false_discovery_control(o.loc[m, "p"])
    return o


main = analyse(S, REFS)
main.round(4).to_csv(RES / "01_expression_main.csv", index=False)
(S[TARGETS].sub(S[REFS].mean(axis=1), axis=0)).join(meta[["sample", "sex"]]).round(3).to_csv(RES / "01_expression_per_colony_log2NRQ.csv")

# sensitivity: every reference set, Welch, individual efficiency, outlier kept, 100 % efficiency
sens = []
for name, refs in {"ACTB": ["ACTB"], "RPL12": ["RPL12"], "SRP54": ["SRP54"], "ACTB+RPL12": REFS, "ACTB+RPL12+SRP54": ["ACTB", "RPL12", "SRP54"]}.items():
    sens.append(analyse(S, refs).assign(variant=f"refs={name}"))
sens.append(analyse(S, REFS, test="welch").assign(variant="Welch t"))
sens.append(analyse(colony_table("log2N0_ind"), REFS).assign(variant="individual-well efficiency"))
sens.append(analyse(colony_table(drop_outlier=False), REFS).assign(variant="ACTB outlier kept"))
sens.append(analyse(colony_table("log2N0_E2"), REFS).assign(variant="E = 2 (100 %)"))
pd.concat(sens).round(4).to_csv(RES / "01_expression_sensitivity.csv", index=False)

# leave-one-colony-out (influence diagnostic, NOT a reason to drop colonies)
NF = S[REFS].mean(axis=1)
loo = []
for g in TARGETS:
    y = S[g] - NF
    for drop in C + T:
        cc, tt = y.loc[[s for s in C if s != drop]], y.loc[[s for s in T if s != drop]]
        loo.append(dict(gene=g, dropped=drop, dropped_sample=meta.loc[drop, "sample"], log2FC=tt.mean() - cc.mean(),
                        p_welch=stats.ttest_ind(tt, cc, equal_var=False).pvalue))
pd.DataFrame(loo).round(4).to_csv(RES / "01_expression_leave_one_out.csv", index=False)

# power: colonies per group for 80 % power at the observed effect and SD (two-sided alpha 0.05)
from scipy.stats import nct, t as tdist
def n_needed(delta, sd):
    for n in range(2, 300):
        df = 2 * n - 2; tc = tdist.ppf(0.975, df); ncp = delta / (sd * np.sqrt(2 / n))
        if 1 - nct.cdf(tc, df, ncp) + nct.cdf(-tc, df, ncp) >= 0.8:
            return n
pw = []
for g in TARGETS:
    y = S[g] - NF
    sd = np.sqrt((y.loc[C].var() + y.loc[T].var()) / 2)
    d = abs(y.loc[T].mean() - y.loc[C].mean())
    pw.append(dict(gene=g, pooled_SD_log2=sd, abs_log2FC=d, n_per_group_80pct=n_needed(d, sd)))
pd.DataFrame(pw).round(3).to_csv(RES / "01_expression_power.csv", index=False)

# Figure 2 is drawn by 07_figures.py

print("STEP 1 done. Main result (Student t, df = 4; q = Benjamini-Hochberg over 9 targets):")
print(main[main.role == "target"][["gene", "FC", "FC_CI_lo", "FC_CI_hi", "p", "q_BH"]].round(3).to_string(index=False))
