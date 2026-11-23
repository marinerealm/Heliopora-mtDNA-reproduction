"""STEP 2 - mtDNA assays: per-well Cq and efficiency.

Assays (total DNA, equal input):  SRP54 (nuclear, ~100 bp)   COI-small (mtDNA, ~100 bp)   COI-large (mtDNA, 1,124 bp)
  SRP54     : 23032018A plate  -> raw curves (Amplification Data) analysed here
  COI-small : 26032018 plate   -> raw curves analysed here (the first run, 23032018B, had 1-2.5-cycle duplicate
                                  disagreement and was repeated; it is used only to validate the estimator)
  COI-large : 23032018C plate  -> LinRegPCR 2017.1 output used directly

2a. VALIDATION: the LinRegPCR-style estimator (qpcr_lib.analyse_curves) is run on the two plates that also have real
    LinRegPCR output (23032018C COI-large, 23032018B COI-small) and compared with LinRegPCR.
2b. QUANTITIES for all three assays -> results/02_mtdna_wells.csv
2c. DIAGNOSTICS: duplicate agreement; sensitivity of SRP54 Ct to the machine threshold (two exports of one plate).
"""
import numpy as np
import pandas as pd
from qpcr_lib import RAW, RES, read_curves, read_machine_results, read_linregpcr, analyse_curves

D = RAW / "mtDNA"
RELABEL = {"SHC520": "SHC524"}   # plate label on SRP54 and COI-small plates; SHC520 is in no sample list

# ---------------------------------------------------------------- 2a. validation against LinRegPCR
val = []
for f, tgt, tag in [("23032018C_COIlarge_LinRegPCR.xls", "COI-Lar", "Lar"), ("23032018B_COIsmall_firstrun_LinRegPCR.xls", "COI-Sml", "Sml")]:
    R = analyse_curves(read_curves(D / f))
    R = R[R.target == tgt].copy(); R["k"] = R.groupby("well").cumcount()
    o = read_linregpcr(D / f); o = o[o.Amplicon == tag].copy()
    o["well"] = o.name.str.split("_").str[0]; o["k"] = o.groupby("well").cumcount()
    m = R.merge(o[["well", "k", "indiv PCR eff", "mean_PCR_eff", "Cq"]].rename(columns={"Cq": "Cq_LinReg"}), on=["well", "k"])
    d = m.Cq - m.Cq_LinReg
    val.append(dict(plate=f, n_wells=len(m), E_mean_estimator=m.E_mean.iloc[0], E_mean_LinRegPCR=m.mean_PCR_eff.iloc[0],
                    per_well_E_mean_abs_diff=np.mean(abs(m.E_ind - m["indiv PCR eff"])),
                    Cq_constant_offset=d.median(), Cq_SD_after_offset=d.std()))
pd.DataFrame(val).round(4).to_csv(RES / "02_estimator_validation.csv", index=False)

# ---------------------------------------------------------------- 2b. quantities
out = []
for f_curves, f_mach, assay in [("23032018A_SRP54_machine_and_curves.xls", "23032018A_SRP54_machine_and_curves.xls", "SRP54"),
                                ("26032018_COIsmall_curves.xls", "26032018_COIsmall_machine_results.xls", "COI-Sml")]:
    R = analyse_curves(read_curves(D / f_curves))
    mach = read_machine_results(D / f_mach)[["well", "sample", "Ct", "Tm1"]].rename(columns={"Ct": "Ct_machine"})
    m = R.merge(mach, on="well", how="left")
    m["sample"] = m["sample"].replace(RELABEL)
    m["assay"] = assay
    out.append(m[["assay", "well", "sample", "Cq", "E_ind", "E_mean", "thr", "log2N0", "Ct_machine", "Tm1"]])
o = read_linregpcr(D / "23032018C_COIlarge_LinRegPCR.xls"); o = o[o.Amplicon == "Lar"]
lar = pd.DataFrame({"assay": "COI-Lar", "well": None, "sample": o.name.str.split("_").str[0], "Cq": o.Cq.astype(float),
                    "E_ind": o["indiv PCR eff"], "E_mean": o.mean_PCR_eff, "thr": o.threshold})
lar["log2N0"] = np.log2(lar.thr) - lar.Cq * np.log2(lar.E_mean)
mL = read_machine_results(D / "23032018C_COIlarge_machine_results.xls")      # same well order ('Sample N' labels)
lar["Ct_machine"] = mL.Ct.values[:len(lar)]; lar["Tm1"] = mL.Tm1.values[:len(lar)]
W = pd.concat(out + [lar], ignore_index=True)
W["rep"] = W.groupby(["assay", "sample"]).cumcount() + 1
W["log2N0_E2"] = -W.Cq
W.to_csv(RES / "02_mtdna_wells.csv", index=False)

# ---------------------------------------------------------------- 2c. diagnostics
dup = W.pivot_table(index=["assay", "sample"], columns="rep", values="Cq")
dup["abs_diff"] = (dup[1] - dup[2]).abs()
summ = dup.groupby("assay").abs_diff.agg(median="median", p90=lambda s: s.quantile(.9), n_over_0_5=lambda s: (s > .5).sum())
summ["E_mean"] = W.groupby("assay").E_mean.first()
summ.round(3).to_csv(RES / "02_mtdna_duplicate_agreement.csv")

a = read_machine_results(D / "23032018A_SRP54_machine_and_curves.xls").set_index("well")
b = read_machine_results(D / "23032018A_SRP54_machine_alt_threshold.xls").set_index("well")
th = pd.DataFrame({"Ct_threshold_0.045": a.Ct, "Ct_threshold_1.95": b.Ct})
th["shift"] = th["Ct_threshold_1.95"] - th["Ct_threshold_0.045"]
th.round(3).to_csv(RES / "02_srp54_threshold_sensitivity.csv")
ok = th["Ct_threshold_0.045"] < 20
print("STEP 2 done.")
print(pd.DataFrame(val).round(3).to_string(index=False))
print(summ.round(3).to_string())
print(f"SRP54: moving the threshold 0.045 -> 1.95 shifts Ct by {th['shift'][ok].mean():.2f} on average, "
      f"but by {th['shift'][ok].min():.2f}-{th['shift'][ok].max():.2f} across wells (SD {th['shift'][ok].std():.2f}); "
      "well-behaved curves shift by a near-constant amount.")
