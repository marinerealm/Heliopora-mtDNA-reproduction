"""Shared helpers: file readers, a LinRegPCR-style curve analyser, and small statistics helpers.

Everything in the pipeline works on log2 starting quantities:
    log2 N0 = log2(threshold) - Cq * log2(E)
where E is the amplification factor per cycle (2.0 = 100 % efficiency).
The threshold term is a constant within an assay, so it cancels in every ratio we use.
"""
from pathlib import Path
import itertools
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
META = ROOT / "data" / "metadata"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)


# --------------------------------------------------------------------------- readers
def read_linregpcr(path):
    """LinRegPCR 2017.1 'Amplification Data_output' sheet -> one row per well, in plate order."""
    o = pd.read_excel(path, "Amplification Data_output", header=3)
    o = o[o["name"].notna() & o["name"].astype(str).str.contains("_")].reset_index(drop=True)
    return o


def read_machine_results(path):
    """StepOnePlus 'Results' sheet -> Well, Sample Name, Target Name, Ct, Tm1, threshold, baseline."""
    r = pd.read_excel(path, "Results", header=None)
    h = r.index[r.iloc[:, 0].astype(str).eq("Well")][0]
    r.columns = r.iloc[h]
    r = r.iloc[h + 1:]
    ct = [c for c in r.columns if str(c)[:1] == "C" and str(c)[1:] in ("т", "t", "T")][0]  # 'Cт' uses a Cyrillic letter
    r = r[r["Well"].astype(str).str.match(r"^[A-H]\d+$")]
    out = pd.DataFrame({"well": r["Well"].astype(str), "sample": r["Sample Name"].astype(str).str.strip(),
                        "target": r["Target Name"], "Ct": pd.to_numeric(r[ct], errors="coerce"),
                        "Tm1": pd.to_numeric(r["Tm1"], errors="coerce"), "threshold": r["Ct Threshold"]})
    return out.reset_index(drop=True)


def read_curves(path):
    """StepOnePlus 'Amplification Data' sheet -> long table (well, cycle, target, Rn). One curve per well."""
    a = pd.read_excel(path, "Amplification Data", header=None)
    h = a.index[a.iloc[:, 0].astype(str).eq("Well")][0]
    a = a.iloc[h + 1:, :5]
    a.columns = ["well", "cycle", "target", "Rn", "dRn"]
    a = a[a.target.notna()].copy()
    a["cycle"] = a.cycle.astype(int)
    a["Rn"] = a.Rn.astype(float)
    a["curve"] = (a.cycle == 1).cumsum()          # wells may be renamed to sample names; curve id = order of appearance
    return a


# --------------------------------------------------------------------------- LinRegPCR-style curve analysis
def _baseline_subtract(F):
    """Choose the baseline that makes the log-linear phase straightest (the LinRegPCR criterion:
    equal slopes in the lower and upper halves of the exponential phase)."""
    c = np.arange(1, len(F) + 1)
    ground = np.median(F[2:10])
    best = (np.inf, ground)
    span = max(F.max() - ground, 1e-6)
    for b in ground + np.linspace(-0.05, 0.05, 41) * span:
        G = F - b
        if G[G > 0].size < 8:
            continue
        top = G.max()
        idx = np.where((G > top * 2 ** -7) & (G < top * 2 ** -2))[0]
        idx = idx[idx < np.argmax(G)]
        if len(idx) < 5:
            continue
        h = len(idx) // 2
        s1 = np.polyfit(c[idx[:h + 1]], np.log10(G[idx[:h + 1]]), 1)[0]
        s2 = np.polyfit(c[idx[h:]], np.log10(G[idx[h:]]), 1)[0]
        if abs(s1 - s2) < best[0]:
            best = (abs(s1 - s2), b)
    return F - best[1]


def analyse_curves(a, win_log2=(-5.5, -2.7), min_pts=3):
    """Per-well efficiency and Cq from raw curves.
    Window of linearity per amplicon: from 2^-5.5 to 2^-2.7 of the median plateau (log2 scale),
    which reproduces LinRegPCR 2017.1 to within 0.012 in mean efficiency (see 02_mtdna_quantities.py).
    E_mean = mean of per-well efficiencies within 5 % of the amplicon median (LinRegPCR rule)."""
    curves = []
    for cid, g in a.groupby("curve", sort=True):
        g = g.sort_values("cycle")
        F = _baseline_subtract(g.Rn.values)
        curves.append(dict(curve=cid, well=g.well.iloc[0], target=g.target.iloc[0], F=F, plateau=np.max(F)))
    C = pd.DataFrame(curves)
    rows = []
    for amp, grp in C.groupby("target"):
        ok = grp[grp.plateau > 0.2 * grp.plateau.median()]
        plat = np.median(np.log2(ok.plateau))
        lo, hi = 2 ** (plat + win_log2[0]), 2 ** (plat + win_log2[1])
        thr = np.sqrt(lo * hi)
        for _, r in grp.iterrows():
            F = r.F
            cyc = np.arange(1, len(F) + 1)
            peak = np.argmax(F)
            idx = [i for i in range(peak) if lo <= F[i] <= hi and F[i] > 0]
            eff, cq = np.nan, np.nan
            if len(idx) >= min_pts:
                eff = 10 ** np.polyfit(cyc[idx], np.log10(F[idx]), 1)[0]
            above = np.where(F[:peak + 1] >= thr)[0]
            if len(above) and above[0] > 0:
                i = above[0]
                cq = cyc[i - 1] + (np.log10(thr) - np.log10(F[i - 1])) / (np.log10(F[i]) - np.log10(F[i - 1]))
            rows.append(dict(curve=r.curve, well=r.well, target=amp, E_ind=eff, n_window=len(idx), Cq=cq, thr=thr, plateau=r.plateau))
    R = pd.DataFrame(rows)

    def mean_eff(s):
        m = s.median()
        return s[(abs(s - m) / m) <= 0.05].mean()

    R["E_mean"] = R.groupby("target").E_ind.transform(mean_eff)
    R["log2N0"] = np.log2(R.thr) - R.Cq * np.log2(R.E_mean)
    return R.sort_values("curve").reset_index(drop=True)


# --------------------------------------------------------------------------- statistics helpers
def perm_p(a, b):
    """Exact two-sided permutation p for a difference in means (all relabelings). With 3 v 3 the minimum is 0.10."""
    v = np.r_[a, b]
    n = len(a)
    obs = abs(np.mean(b) - np.mean(a))
    k = t = 0
    for idx in itertools.combinations(range(len(v)), n):
        m = np.zeros(len(v), bool)
        m[list(idx)] = True
        k += abs(v[~m].mean() - v[m].mean()) >= obs - 1e-12
        t += 1
    return k / t


def lesions_per_10kb(log2_ratio_vs_ref, length_bp=1124):
    """SLR-qPCR lesion frequency from r = relative amplification of the long fragment (vs reference group).
    Returns (Poisson, Rothfuss):  -ln(r) * 1e4/L   and   (1 - r) * 1e4/L."""
    r = 2.0 ** np.asarray(log2_ratio_vs_ref, float)
    return -np.log(r) * 1e4 / length_bp, (1 - r) * 1e4 / length_bp
