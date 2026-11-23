"""STEP 4 - Independent replicate run (14 March 2018) on the six KR mtDNA colonies.

SHC413-415 (before, 13 Jul) and SHC608/612/614 (during brooding, 21 Jul); ACTB, SRP54, COI-small and COI-large in
triplicate on one plate, with no-template controls (NTC). This run tests whether the KR copy-number change seen on the
23 March SRP54 plate replicates, using TWO nuclear references.

Note: the curve export of this run lacks wells B4-B9 and F10-F12, so copy number uses the machine Ct values (all wells
present) with efficiencies estimated from the available curves. Damage uses curve-based Cq for COI-large and COI-small.
Null: equal mean before vs during (Welch t).
"""
import numpy as np
import pandas as pd
from scipy import stats
from qpcr_lib import RAW, RES, read_curves, read_machine_results, analyse_curves, lesions_per_10kb

D = RAW / "mtDNA_KR_replicate_20180314"
setup = pd.read_excel(D / "20180314_sample_setup.xls", "Sample Setup", header=None)
h = setup.index[setup.iloc[:, 0].astype(str).eq("Well")][0]
setup.columns = setup.iloc[h]; setup = setup.iloc[h + 1:]
setup = setup[setup["Target Name"].notna()][["Well", "Sample Name", "Target Name", "Task"]]
setup.columns = ["well", "sample", "target", "task"]

curves = analyse_curves(read_curves(D / "20180314_curves.xls")).merge(setup[["well", "sample", "task"]], on="well", how="left")
E = curves[curves.task != "NTC"].groupby("target").E_mean.first()
ntc = curves[curves.task == "NTC"][["well", "target", "Cq", "plateau"]]
ntc.round(2).to_csv(RES / "04_ntc.csv", index=False)

mach = read_machine_results(D / "20180314_machine_results.xls")
mach = mach[mach.Ct.notna()].copy()
summary = []
for label, eff in [("efficiency-corrected", True), ("E = 2", False)]:
    mach["q"] = -mach.Ct * (np.log2(mach.target.map(E)) if eff else 1.0)
    S = mach.groupby(["sample", "target"]).q.mean().unstack()
    S["mtCN_ACTB"] = S["COI-SMALL"] - S["ACTB"]
    S["mtCN_SRP54"] = S["COI-SMALL"] - S["SRP54"]
    S["ACTB_minus_SRP54"] = S["ACTB"] - S["SRP54"]
    grp = np.where(S.index.isin(["SHC413", "SHC414", "SHC415"]), "before", "during")
    for col in ["mtCN_ACTB", "mtCN_SRP54", "ACTB_minus_SRP54"]:
        c, d = S[col][grp == "before"], S[col][grp == "during"]
        summary.append(dict(quantity=col, method=label, log2_change=d.mean() - c.mean(), fold=2 ** (d.mean() - c.mean()),
                            p_welch=stats.ttest_ind(d, c, equal_var=False).pvalue))
    if eff:
        S.round(3).to_csv(RES / "04_replicate_per_colony.csv")
trip_sd = mach.groupby(["target", "sample"]).Ct.std().groupby("target").median()

# damage from curve-based quantities
X = curves[curves.task != "NTC"]
Q = X.groupby(["sample", "target"]).log2N0.mean().unstack()
LS = Q["COI-LARGE"] - Q["COI-SMALL"]
before = ["SHC413", "SHC414", "SHC415"]
les = pd.Series(lesions_per_10kb(LS - LS[before].mean())[0], index=LS.index).dropna()
dur = [s for s in les.index if s not in before]
summary.append(dict(quantity="lesions_poisson_per_10kb", method=f"curves (during n = {len(dur)})", log2_change=les[dur].mean() - les[before].mean(),
                    fold=np.nan, p_welch=stats.ttest_ind(les[dur], les[before], equal_var=False).pvalue))
out = pd.DataFrame(summary)
out.round(4).to_csv(RES / "04_replicate_summary.csv", index=False)

# same DNA on two dates: SRP54 spread on 23 Mar vs 14 Mar
W = pd.read_csv(RES / "02_mtdna_wells.csv")
q23 = W[W["sample"].isin(S.index)].groupby(["sample", "assay"]).log2N0.mean().unstack()
S_eff = pd.read_csv(RES / "04_replicate_per_colony.csv", index_col=0)
cmp_ = pd.DataFrame({"SRP54_23Mar": q23["SRP54"] - q23["SRP54"].mean(), "SRP54_14Mar": S_eff["SRP54"] - S_eff["SRP54"].mean()})
cmp_.round(3).to_csv(RES / "04_srp54_two_runs.csv")
print("STEP 4 done. Triplicate SD (median Ct):", trip_sd.round(3).to_dict())
print(out.round(3).to_string(index=False))
print("SRP54, same six DNAs - range on 23 Mar: %.2f log2; on 14 Mar: %.2f log2" %
      (np.ptp(cmp_.SRP54_23Mar), np.ptp(cmp_.SRP54_14Mar)))
