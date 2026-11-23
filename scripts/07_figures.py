"""STEP 7 - Manuscript figures (brooding terminology; subclades sp.3, sp.2, coe.3 = historical HC-A, HC-A', HC-B), drawn from the result tables of steps 1, 3 and 4.
Fig1_expression: wider panel A, colonies dodged, sex marked (filled = female, open = male), SOD bracket,
    panel B with p and Benjamini-Hochberg q columns.
Fig2_mtDNA: dot plots with p above each comparison (as in the original Fig. 2), plus effect-size panels
    (difference, 95% CI) including the KR repeat run (14 Mar 2018) and the 2017 subclade contrasts.
Significance: * = unadjusted p < 0.05 (bold); grey = p >= 0.05. Run after steps 1-4."""
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from qpcr_lib import RES
OUT = RES

plt.rcParams.update({"svg.fonttype": "none", "font.family": ["Arial", "Liberation Sans"], "mathtext.fontset": "custom", "mathtext.rm": "Liberation Sans", "mathtext.it": "Liberation Sans:italic", "mathtext.bf": "Liberation Sans:bold", "mathtext.default": "regular", "font.size": 8, "axes.linewidth": 0.6})
CB, CR, GRID, INK, MUTED, NS = "#2a78d6", "#eb6834", "#e4e3df", "#0b0b0b", "#52514e", "#9a9893"
TARGETS = ["SOD", "D2R", "mtMutS", "GPX", "HSP70", "HSP90", "CAT", "CRY1", "COI"]


def fmt_p(p, star=True):
    s = f"{p:.3f}" if p < .1 else f"{p:.2f}"
    return s + ("*" if star and p < .05 else "")


def ptext(ax, x, y, p, prefix="p = ", **kw):
    ax.text(x, y, prefix + fmt_p(p), ha="center", va="bottom", fontsize=7,
            fontweight="bold" if p < .05 else "normal", color=INK if p < .05 else MUTED, **kw)


def bracket(ax, x1, x2, y, h, p, prefix="p = "):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color=INK if p < .05 else NS, lw=.7)
    ptext(ax, (x1 + x2) / 2, y + h * 1.3, p, prefix)


def forest(ax, rows, xlabel, xlim, extra=None):
    """rows: (label, diff, lo, hi, p, style) or ('--', header). style: main / rep (repeat run) / int (interaction)."""
    n = len(rows); y = n - 1; lx = xlim[0] - (xlim[1] - xlim[0]) * .04
    for k, r in enumerate(rows):
        if r[0] == "--":
            ax.text(lx, y, r[1], fontsize=7, fontweight="bold", va="center", ha="right", color=INK)
            if k: ax.axhline(y + .5, color=GRID, lw=.6, xmin=-.85, clip_on=False)
        else:
            lab, d, lo, hi, p, sty = r
            sig = p < .05; col = INK if sig else (MUTED if sty == "rep" else NS)
            ax.hlines(y, lo, hi, color=col, lw=1.3 if sig else 1.0)
            ax.scatter(d, y, s=30 if sig else 22, marker="D" if sty == "int" else "o", zorder=3, lw=1.1,
                       facecolor="white" if sty == "rep" else col, edgecolor=col)
            ax.text(lx, y, lab, va="center", ha="right", fontsize=7, color=INK, style="italic" if sty == "gene" else "normal")
            ax.text(1.04, y, fmt_p(p), va="center", ha="left", fontsize=7, transform=ax.get_yaxis_transform(),
                    fontweight="bold" if sig else "normal", color=INK if sig else MUTED)
            if extra is not None:
                ax.text(1.62, y, extra[lab], va="center", ha="left", fontsize=7, color=MUTED, transform=ax.get_yaxis_transform())
        y -= 1
    ax.text(1.04, n - .3, "p", ha="left", fontsize=7, style="italic", color=MUTED, transform=ax.get_yaxis_transform())
    if extra is not None:
        ax.text(1.62, n - .3, "q", ha="left", fontsize=7, style="italic", color=MUTED, transform=ax.get_yaxis_transform())
    ax.axvline(0, color=MUTED, lw=.6, ls=(0, (3, 2)))
    ax.set_xlim(*xlim); ax.set_ylim(-.7, n - .1); ax.set_yticks([])
    ax.set_xlabel(xlabel); ax.spines[["top", "right", "left"]].set_visible(False)


# ================================================================ Figure 1 v3
N = pd.read_csv(RES / "01_expression_per_colony_log2NRQ.csv").set_index("label")
main = pd.read_csv(RES / "01_expression_main.csv").set_index("gene").loc[TARGETS]
fig = plt.figure(figsize=(7.2, 3.7))
a = fig.add_axes([.105, .2, .54, .68]); b = fig.add_axes([.745, .2, .13, .68])
GAP = .22; DODGE = [-.075, 0, .075]
for i, g in enumerate(TARGETS):
    for labs, dx, col, mk in [(["C1", "C2", "C3"], -GAP, CB, "o"), (["T1", "T2", "T3"], GAP, CR, "^")]:
        v = N.loc[labs, g]; order = np.argsort(v.values)          # dodge in value order so close points separate
        for j, (val, s) in enumerate(zip(v.values[order], N.loc[labs, "sex"].values[order])):
            a.scatter(i + dx + DODGE[j], val, s=24, marker=mk, zorder=3, lw=1.0,
                      facecolor=col if s == "F" else "white", edgecolor=col)
        a.hlines(v.mean(), i + dx - .14, i + dx + .14, color=col, lw=1.6, zorder=2)
    if i: a.axvline(i - .5, color=GRID, lw=.5, zorder=0)
a.axhline(0, color=GRID, lw=.8, zorder=0)
a.set_xlim(-.5, len(TARGETS) - .5)
a.set_xticks(range(len(TARGETS)), TARGETS, style="italic")
a.set_ylabel("log$_2$ relative quantity\n(normalised to ACTB + RPL12, centred)")
top = N.loc[:, "SOD"].max()
bracket(a, -GAP, GAP, top + .15, .08, main.loc["SOD", "p"])
a.set_ylim(-1.8, 2.0)
h = [a.scatter([], [], color=CB, s=22, label="Before brooding (13 Jul)"),
     a.scatter([], [], marker="^", color=CR, s=24, label="During brooding (25 Jul)"),
     a.scatter([], [], color=MUTED, s=22, label="Female"),
     a.scatter([], [], facecolor="white", edgecolor=MUTED, lw=1.0, s=22, label="Male")]
a.legend(handles=h, frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=4, columnspacing=1.2, handletextpad=.3)
a.spines[["top", "right"]].set_visible(False)
rows = [(g, r.log2FC, r.CI95_log2_lo, r.CI95_log2_hi, r.p, "gene") for g, r in main.iterrows()]
forest(b, rows, "log$_2$ fold change,\nduring vs before (95% CI)", (-3.0, 2.2),
       extra={g: f"{q:.2f}" for g, q in main.q_BH.items()})
b.set_xticks([-2, -1, 0, 1, 2])
fig.text(.012, .95, "A", fontweight="bold", fontsize=9); fig.text(.662, .95, "B", fontweight="bold", fontsize=9)
fig.savefig(OUT / "Fig1_expression.svg"); fig.savefig(OUT / "Fig1_expression.png", dpi=300)

# ================================================================ Figure 2 v3
T = pd.read_csv(RES / "03_mtdna_per_colony.csv").set_index("sample")
R = pd.read_csv(RES / "03_mtdna_tests.csv")
Q = lambda an, ct, me: R[(R.analysis == an) & (R.contrast == ct) & (R.measure == me)].iloc[0]
b14 = T[T.analysis_group.str.startswith("mt1")]; y17 = T[T.analysis_group.str.startswith("mt2")]
WI, INT, TK, AN = ("B2014 reproductive vs before", "B2014 KR vs YIS difference in change",
                   "Yonara 2017 lineages (Tukey HSD)", "Yonara 2017 lineages")
LIN = "A = A' = B (one-way ANOVA)"

rep = pd.read_csv(RES / "04_replicate_per_colony.csv").set_index("sample")
bef, dur = ["SHC413", "SHC414", "SHC415"], ["SHC608", "SHC612", "SHC614"]
def welch(x, y):
    x, y = np.asarray(x), np.asarray(y); d = y.mean() - x.mean()
    vx, vy = x.var(ddof=1) / len(x), y.var(ddof=1) / len(y); se = np.sqrt(vx + vy)
    df = (vx + vy) ** 2 / (vx ** 2 / (len(x) - 1) + vy ** 2 / (len(y) - 1)); tq = stats.t.ppf(.975, df)
    return d, d - tq * se, d + tq * se, stats.ttest_ind(y, x, equal_var=False).pvalue
rep_rows = [(f"KR rerun, {ref}", *welch(rep.loc[bef, c], rep.loc[dur, c]), "rep") for ref, c in [("ACTB", "mtCN_ACTB"), ("SRP54", "mtCN_SRP54")]]
pd.DataFrame([dict(contrast=r[0], diff=r[1], CI_lo=r[2], CI_hi=r[3], p=r[4]) for r in rep_rows]).round(4).to_csv(
    OUT / "07_kr_replicate_welch_ci.csv", index=False)

def row(lab, an, ct, me, sty="main"):
    q = Q(an, ct, me); return (lab, q["diff"], q.CI_lo, q.CI_hi, q.p, sty)
cn_rows = [("--", "2014 sites"), row("KR", WI, "KR", "mtCN"), *rep_rows,
           row("YIS", WI, "YIS", "mtCN"), row("YKS", WI, "YKS", "mtCN"),
           row("YIS vs KR †", INT, "YIS change - KR change", "mtCN", "int"),
           ("--", "2017 subclades"), row("sp.3 vs coe.3", TK, "A - B", "mtCN"), row("sp.2 vs coe.3", TK, "A' - B", "mtCN")]
dm_rows = [("--", "2014 sites"), row("KR", WI, "KR", "lesions_poisson"), row("YIS", WI, "YIS", "lesions_poisson"),
           row("YKS", WI, "YKS", "lesions_poisson"), row("YIS vs KR", INT, "YIS change - KR change", "lesions_poisson", "int"),
           ("--", "2017 subclades"), row("sp.3 vs coe.3", TK, "A - B", "lesions_poisson"), row("sp.2 vs coe.3", TK, "A' - B", "lesions_poisson")]

rng = np.random.default_rng(1)
def strip(ax, x, vals, color, marker):
    vals = np.sort(vals)
    ax.scatter(x + np.array([-.06, 0, .06])[:len(vals)], vals, s=24, marker=marker, color=color, edgecolor="white", lw=.6, zorder=3)
    ax.hlines(np.mean(vals), x - .15, x + .15, color=color, lw=1.6, zorder=2)

fig = plt.figure(figsize=(7.2, 7.6))
ROWS = [(.575, .335), (.08, .335)]     # (bottom, height)
axs = []
for bot, hgt in ROWS:
    axs += [fig.add_axes([.105, bot, .275, hgt]), fig.add_axes([.53, bot, .17, hgt]), fig.add_axes([.875, bot, .12, hgt])]
A, B, C, D, E, F = axs
cnb = b14[~b14.exclude_copy_number.astype(bool)]
refm = cnb[cnb.state2 == "before"].mtCN.mean()
for i, site in enumerate(["KR", "YIS", "YKS"]):
    va, vd = [], []
    for k, (st, col, mk) in enumerate([("before", CB, "o"), ("reproductive", CR, "^")]):
        x = i + (k - .5) * .44
        v1 = (cnb[(cnb.site_code == site) & (cnb.state2 == st)].mtCN - refm).values
        v2 = b14[(b14.site_code == site) & (b14.state2 == st)].lesions_poisson.values
        strip(A, x, v1, col, mk); strip(D, x, v2, col, mk); va += list(v1); vd += list(v2)
    bracket(A, i - .22, i + .22, max(va) + .12, .07, Q(WI, site, "mtCN").p)
    bracket(D, i - .22, i + .22, max(vd) + 1.0, .6, Q(WI, site, "lesions_poisson").p)
lc, ld = [], []
for i, l in enumerate(["A", "A'", "B"]):
    v1 = (y17[y17.lineage == l].mtCN - y17[y17.lineage == "B"].mtCN.mean()).values
    v2 = y17[y17.lineage == l].lesions_poisson.values
    strip(C, i, v1, INK, "o"); strip(F, i, v2, INK, "o"); lc.append(v1.max()); ld.append(v2.max())
bracket(F, 0, 2, max(ld) + 1.4, .4, Q(TK, "A - B", "lesions_poisson").p)
A.set_ylim(-1.75, 1.75); D.set_ylim(-7, 21); C.set_ylim(-.9, 1.1); C.set_yticks([-.5, 0, .5, 1]); F.set_ylim(-3.2, 11)
for ax_, me in [(C, "mtCN"), (F, "lesions_poisson")]:
    ptext(ax_, 1, ax_.get_ylim()[1], Q(AN, LIN, me).p, prefix="ANOVA\np = ", clip_on=False)
A.set_xticks([0, 1, 2], ["KR\n21 Jul", "YIS\n21 Jul", "YKS\n28 Jul"])
D.set_xticks([0, 1, 2], ["KR", "YIS", "YKS"])
for ax_ in (C, F):
    ax_.set_xticks([0, 1, 2], ["sp.3", "sp.2", "coe.3"]); ax_.set_xlim(-.6, 2.6); ax_.tick_params(axis="x", labelsize=7)
C.set_xlabel("Subclade"); F.set_xlabel("Subclade\n(13 Jun 2017)")
A.set_ylabel("mtDNA copy number (log$_2$)\nvs before-brooding mean")
D.set_ylabel("mtDNA lesions per 10 kb\nvs before brooding, same site")
C.set_ylabel("Copy number (log$_2$)\nvs coe.3"); F.set_ylabel("Lesions per 10 kb\nvs coe.3")
A.scatter([], [], color=CB, s=22, label="Before brooding (13 Jul)")
A.scatter([], [], marker="^", color=CR, s=24, label="During brooding (21/28 Jul)")
A.legend(frameon=False, fontsize=7, loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, handletextpad=.3, columnspacing=1.0)
for ax_ in (A, C, D, F):
    ax_.axhline(0, color=GRID, lw=.8, zorder=0); ax_.spines[["top", "right"]].set_visible(False)
forest(B, cn_rows, "Difference in copy number\n(log$_2$, 95% CI)", (-3.0, 3.0))
forest(E, dm_rows, "Difference in lesions\nper 10 kb (95% CI)", (-20, 27))
for (bot, hgt), labs in zip(ROWS, ["ABC", "DEF"]):
    for x, lab in zip([.01, .405, .80], labs):
        fig.text(x, bot + hgt + .05, lab, fontweight="bold", fontsize=9)
fig.savefig(OUT / "Fig2_mtDNA.svg"); fig.savefig(OUT / "Fig2_mtDNA.png", dpi=300)
# name Arial first in every SVG text element (Liberation Sans, metric-compatible, is the fallback used to draw them)
for f in ("Fig1_expression.svg", "Fig2_mtDNA.svg"):
    fp = OUT / f
    fp.write_text(fp.read_text().replace("font-family: 'Liberation Sans'", "font-family: 'Arial', 'Liberation Sans'"))
print("STEP 7 done: Fig1_expression and Fig2_mtDNA (.svg editable, .png preview)")
