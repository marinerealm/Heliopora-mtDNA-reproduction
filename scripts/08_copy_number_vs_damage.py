"""STEP 8 - Is mtDNA copy number (or the nuclear reference) associated with estimated mtDNA damage among colonies?

Question: within a sampling group, do colonies with more estimated mtDNA damage differ in
  (a) SRP54 quantity per DNA input (the nuclear reference used for copy number),
  (b) short-COI quantity per DNA input (mtDNA quantity without SRP54), or
  (c) SRP54-normalised copy number?

Colonies: the main data sets only - 2014 coe.3 colonies (KR, YIS, YKS; before and during brooding) and the
2017 subclade panel (coe.3, sp.3, sp.2). Exploratory and archived samples are not used. Colonies flagged
exclude_copy_number in samples.csv (SRP54 failed: SHC701) are left out of the comparisons that use SRP54.

Measures (per colony, from results/03_mtdna_per_colony.csv, log2 scale):
  damage        = log2 N0(short COI) - log2 N0(long COI) = -LS. Lesions per 10 kb are this value times a constant
                  (ln 2 x 10,000 / 1,124 bp), so correlations are identical; the figure shows lesions per 10 kb.
  SRP54/input   = log2 N0(SRP54) per DNA input = COIs_per_input - mtCN
  shortCOI/input= COIs_per_input
  mtCN          = log2 N0(short COI) - log2 N0(SRP54)
Group effects removed: every value is centred on its sampling-group mean (nine groups of three colonies:
site x period in 2014, subclade in 2017), so the tests ask whether colonies differ together WITHIN groups.

Tests. Null hypothesis: no linear association between the two centred measures.
  Pearson r with df = n - (number of groups) - 1 (the group means use up one degree of freedom each);
  within-group permutation test (values shuffled among colonies of the same group, 20,000 permutations, seed 1);
  Spearman rho of the centred values (descriptive); r after leaving out each group in turn (robustness).
Note: damage and short COI per input share the short-COI measurement; measurement error in it would push their
correlation towards positive values. SRP54 and long COI were run on different plates.

Outputs: 08_copy_number_vs_damage_tests.csv, 08_copy_number_vs_damage_per_colony.csv,
         FigS1_copy_number_vs_damage.svg/.png; a section is appended to SUMMARY.md and a sheet to all_results.xlsx.
Run after steps 1-7."""
import logging
import math
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)   # Arial falls back to Liberation Sans quietly
import matplotlib.pyplot as plt
import matplotlib.ticker
import inspect
from qpcr_lib import RES, lesions_per_10kb

LONG_COI_BP = inspect.signature(lesions_per_10kb).parameters["length_bp"].default   # 1,124 bp, set in qpcr_lib
LESIONS_PER_LOG2 = math.log(2) * 10_000 / LONG_COI_BP  # lesions per 10 kb for one log2 unit of damage
N_PERM, SEED = 20_000, 1

d = pd.read_csv(RES / "03_mtdna_per_colony.csv")
d = d[d.analysis_group.str.match(r"mt[12]_")].copy()
d["damage"] = -d["LS"]
d["srp54_per_input"] = d["COIs_per_input"] - d["mtCN"]
d["group"] = np.where(d.analysis_group.str.startswith("mt2"),
                      "2017 " + d.lineage.map({"B": "coe.3", "A": "sp.3", "A'": "sp.2"}),
                      "2014 " + d.site_code + " " + d.state2.map({"before": "before", "reproductive": "during"}))
excl = d.exclude_copy_number.astype(bool)
d.loc[excl, ["srp54_per_input", "mtCN"]] = np.nan          # SRP54 failed: not used

PAIRS = [("srp54_per_input", "SRP54 per DNA input"),
         ("COIs_per_input", "short COI per DNA input"),
         ("mtCN", "SRP54-normalised copy number")]


def centred(frame, col):
    return frame[col] - frame.groupby("group")[col].transform("mean")


def test(x):
    s = d[[x, "damage", "group", "sample"]].dropna().copy()
    # centre after dropping missing values, so each group mean uses the same colonies for both measures
    xc, yc = centred(s, x).to_numpy(), centred(s, "damage").to_numpy()
    n, k = len(s), s.group.nunique()
    r = float(np.corrcoef(xc, yc)[0, 1])
    df = n - k - 1
    p = float(2 * stats.t.sf(abs(r) * math.sqrt(df / (1 - r ** 2)), df))
    rng = np.random.default_rng(SEED)
    g = s.group.to_numpy(); idx = [np.where(g == u)[0] for u in np.unique(g)]
    hits = 0
    for _ in range(N_PERM):
        yp = yc.copy()
        for i in idx: yp[i] = rng.permutation(yp[i])
        hits += abs(np.corrcoef(xc, yp)[0, 1]) >= abs(r) - 1e-12
    p_perm = (hits + 1) / (N_PERM + 1)
    rho = float(stats.spearmanr(xc, yc)[0])
    loo = [float(np.corrcoef(xc[g != u], yc[g != u])[0, 1]) for u in np.unique(g)]
    return dict(n_colonies=n, n_groups=k, pearson_r=r, df=df, p=p, p_permutation=p_perm, spearman_rho=rho,
                r_leave_one_group_out_min=min(loo), r_leave_one_group_out_max=max(loo)), s.assign(x_c=xc, damage_c=yc)


rows, frames = [], {}
for col, label in PAIRS:
    res, s = test(col)
    rows.append(dict(measure=label, compared_with="estimated mtDNA damage", **res))
    frames[col] = s
T = pd.DataFrame(rows)
T.round(4).to_csv(RES / "08_copy_number_vs_damage_tests.csv", index=False)

pc = d[["sample", "group", "site_code", "lineage", "damage", "srp54_per_input", "COIs_per_input", "mtCN"]].copy()
for c in ["damage", "srp54_per_input", "COIs_per_input", "mtCN"]:
    pc[c + "_centred"] = pc[c] - pc.dropna(subset=[c]).groupby("group")[c].transform("mean")
pc["lesions_per_10kb_centred"] = pc["damage_centred"] * LESIONS_PER_LOG2
# damage centred on the colonies that also have SRP54 (used in the SRP54 tests; differs only in the group with SHC701)
ok = pc.dropna(subset=["srp54_per_input"])
pc["damage_centred_srp54_set"] = ok["damage"] - ok.groupby("group")["damage"].transform("mean")
pc.round(4).to_csv(RES / "08_copy_number_vs_damage_per_colony.csv", index=False)

# ------------------------------------------------------------------ Figure S1
plt.rcParams.update({"svg.fonttype": "none", "font.family": ["Arial", "Liberation Sans"], "font.size": 8,
                     "axes.linewidth": 0.6})
INK, MUTED = "#0b0b0b", "#52514e"
SITE_COL = {"KR": "#0072B2", "YIS": "#009E73", "YKS": "#E69F00"}
SUB_MK = {"coe.3": "o", "sp.3": "^", "sp.2": "s"}


def style(grp):
    if grp.startswith("2014"):
        _, site, period = grp.split()
        return dict(marker="o", mec=SITE_COL[site], mfc=SITE_COL[site] if period == "during" else "white")
    sub = grp.split()[1]
    return dict(marker=SUB_MK[sub], mec="#555555", mfc="#999999")


fig, axes = plt.subplots(1, 2, figsize=(174 / 25.4, 78 / 25.4))
for ax, (col, label), tag in zip(axes, PAIRS[:2], "AB"):
    s = frames[col]; res = T[T.measure == label].iloc[0]
    y = s.damage_c * LESIONS_PER_LOG2
    for grp, sub in s.assign(y=y).groupby("group"):
        st = style(grp)
        ax.plot(sub.x_c, sub.y, ls="none", ms=5, mew=0.9, **st)
    b = np.polyfit(s.x_c, y, 1); xs = np.linspace(s.x_c.min(), s.x_c.max(), 2)
    sig = res.p < .05
    ax.plot(xs, np.polyval(b, xs), color=INK if sig else MUTED, lw=0.8, ls="-" if sig else "--")
    ax.axhline(0, color="#e4e3df", lw=0.6, zorder=0); ax.axvline(0, color="#e4e3df", lw=0.6, zorder=0)
    ax.set_xlabel(f"{label}\n(log$_2$, difference from group mean)")
    ax.set_ylabel("Estimated mtDNA damage\n(lesions per 10 kb, difference from group mean)")
    ptxt = f"r = {res.pearson_r:.2f}, p = {res.p:.3f}" if res.p < .1 else f"r = {res.pearson_r:.2f}, p = {res.p:.2f}"
    ax.text(0.98, 0.97, f"{ptxt}\npermutation p = {res.p_permutation:.2f}\nn = {res.n_colonies} colonies, {res.n_groups} groups",
            transform=ax.transAxes, ha="right", va="top", fontsize=7, color=INK if sig else MUTED,
            fontweight="bold" if sig else "normal")
    lo, hi = y.min(), y.max(); pad = (hi - lo) * 0.06
    ax.set_ylim(lo - pad, hi + (hi - lo) * 0.32)   # headroom so the statistics box never covers a point
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(2))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:.0f}".replace("-", "\u2212")))
    ax.text(-0.2, 1.02, tag, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom")
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.5, width=0.6)
handles = [plt.Line2D([], [], ls="none", marker="o", mec=c, mfc="white", ms=5, label=f"{s} before") for s, c in SITE_COL.items()] + \
          [plt.Line2D([], [], ls="none", marker="o", mec=c, mfc=c, ms=5, label=f"{s} during") for s, c in SITE_COL.items()] + \
          [plt.Line2D([], [], ls="none", marker=m, mec="#555555", mfc="#999999", ms=5, label=f"2017 {s}") for s, m in SUB_MK.items()]
fig.legend(handles=handles, loc="lower center", ncol=9, fontsize=6.5, frameon=False, handletextpad=0.2, columnspacing=0.9)
fig.tight_layout(rect=(0, 0.07, 1, 1), w_pad=2.5)
fig.savefig(RES / "FigS1_copy_number_vs_damage.svg"); fig.savefig(RES / "FigS1_copy_number_vs_damage.png", dpi=300)

# ------------------------------------------------------------------ add to SUMMARY.md and all_results.xlsx
lines = ["", "## Step 8 - Copy number and the nuclear reference versus estimated mtDNA damage", "",
         "Within-group Pearson correlations (values centred on the mean of each sampling group; df = n - groups - 1; "
         "null hypothesis r = 0) and within-group permutation tests (20,000 permutations).", "",
         "| Measure | n | r | p | permutation p | r without one group |", "|---|---|---|---|---|---|"]
for _, r in T.iterrows():
    lines.append(f"| {r.measure} | {r.n_colonies} | {r.pearson_r:.2f} | {r.p:.3f} | {r.p_permutation:.3f} | "
                 f"{r.r_leave_one_group_out_min:.2f} to {r.r_leave_one_group_out_max:.2f} |")
lines += ["", "Figure S1: results/FigS1_copy_number_vs_damage.svg", ""]
summ = RES / "SUMMARY.md"
text = summ.read_text(encoding="utf-8")
start = text.find("\n## Step 8 - ")
if start >= 0:                                   # replace only this step's section (keeps later steps)
    nxt = text.find("\n## Step ", start + 1)
    text = text[:start] + (text[nxt:] if nxt >= 0 else "")
summ.write_text(text.rstrip("\n") + "\n" + "\n".join(lines), encoding="utf-8")
xl = RES / "all_results.xlsx"
if xl.exists():
    with pd.ExcelWriter(xl, mode="a", engine="openpyxl", if_sheet_exists="replace") as w:
        T.round(4).to_excel(w, sheet_name="08_copy_number_vs_damage", index=False)
        pc.round(4).to_excel(w, sheet_name="08_cn_vs_damage_per_colony", index=False)

print(T.round(3).to_string(index=False))
print("STEP 8 done: 08_copy_number_vs_damage_*.csv, FigS1_copy_number_vs_damage (.svg editable, .png preview)")
