"""Figure 1: Heliopora sampling sites, with an aerial photograph of the Yonara Channel showing reef structure.

Panels
  A  Ryukyu Archipelago (Natural Earth land)                  -> locates the Yaeyama Islands
  B  Yaeyama Islands (MLIT National Land Numerical Info. N03) -> locates the Yonara Channel
  C  Yonara Channel on a GSI aerial photograph (seamless photo tiles) with the sampling sites;
     marker colour and shape = subclade sampled at that site

Not run by run_all.py: this script needs internet access (it downloads the map data and aerial tiles listed under
DATA_URLS and TILE_URL; downloads are cached, so re-running is fast) and the extra package geopandas.
    pip install geopandas pillow matplotlib        (or: mamba install -c conda-forge geopandas pillow matplotlib)
    python make_fig1_sites_aerial.py
Output folder: Fig1_sampling_sites/ next to where you run it (SVG, PDF, 600 dpi TIFF and PNG, README.txt).

To change what is shown, edit SITES and LINEAGES below.
"""
import math, time, json, urllib.request, io
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle, ConnectionPatch
from matplotlib.ticker import FuncFormatter, MultipleLocator
from shapely.geometry import box

OUT = Path.cwd() / "Fig1_sampling_sites"
CACHE = OUT / "map_data"
TILE_ZOOM = 16                      # ~2.2 m per pixel at 24° N
TILE_URL = "https://cyberjapandata.gsi.go.jp/xyz/seamlessphoto/{z}/{x}/{y}.jpg"
DATA_URLS = {
    "ne_10m_land.geojson": "https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_land.geojson",
    "N03_47381.json": "https://raw.githubusercontent.com/niiyz/JapanCityGeoJson/master/geojson/47/47381.json",   # Taketomi Town
    "N03_47207.json": "https://raw.githubusercontent.com/niiyz/JapanCityGeoJson/master/geojson/47/47207.json",   # Ishigaki City
}

# ------------------------------------------------------------------ what is plotted
def dms(d, m, s): return d + m / 60 + s / 3600
LINEAGES = {  # key: (legend label, colour, marker)
    "B":  ("coe.3 (HC-B)",  "#E69F00", "o"),
    "A":  ("sp.3 (HC-A)",   "#56B4E9", "^"),
    "A'": ("sp.2 (HC-A′)",  "#CC79A7", "s"),
}
SITES = [  # code, full name, lat, lon, subclades sampled there (from samples.csv), label offset (dx, dy in degrees)
    # positions: boat positions on the 13 July 2014 field sheets (RNAsampleList_v1.xlsx)
    ("KR",  "Kohama Reef",           dms(24, 21, 15.49), dms(123, 57, 32.69), ["B"],      (0.0011, 0.0006)),
    ("YKS", "Yonara Kohama Slope",   dms(24, 21, 11.87), dms(123, 57, 13.48), ["B", "A"], (-0.0013, -0.0011)),
    ("YIS", "Yonara Iriomote Slope", dms(24, 21, 47.93), dms(123, 56, 57.10), ["B"],      (-0.0012, 0.0008)),
    # 13 June 2017 subclade panel (coe.3, sp.3, sp.2): settlement-plate site on the western side of the Yonara Channel
    # (ヨナラ水道西着底板設置地点, 24°21.991' N, 123°57.035' E)
    ("2017", "Yonara Channel west (2017)", 24 + 21.991 / 60, 123 + 57.035 / 60, ["B", "A", "A'"], (0.0030, 0.0)),
]
EXT_A = (121.0, 131.5, 23.4, 31.6)
EXT_B = (123.62, 124.36, 24.17, 24.62)
EXT_C = (123.925, 123.990, 24.330, 24.381)   # north edge kept below 24°23' N, where the GSI photo coverage ends

LAND, LANDEDGE, SEA = "#d9d9d9", "#4d4d4d", "#eef4f8"
matplotlib.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                            "font.size": 7, "svg.fonttype": "none", "pdf.fonttype": 42, "axes.linewidth": 0.6})

# ------------------------------------------------------------------ downloads
def fetch(url, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Heliopora-site-map (academic figure)"})
            with urllib.request.urlopen(req, timeout=60) as r: return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 404: return None
            time.sleep(2 * (k + 1))
        except Exception:
            time.sleep(2 * (k + 1))
    raise SystemExit(f"Could not download {url}. Check the internet connection and run again.")

CACHE.mkdir(parents=True, exist_ok=True)
for name, url in DATA_URLS.items():
    f = CACHE / name
    if not f.exists() or f.stat().st_size == 0:
        print("downloading", name); f.write_bytes(fetch(url))

def lonlat_to_tile(lon, lat, z):
    n = 2 ** z
    x = (lon + 180) / 360 * n
    y = (1 - math.log(math.tan(math.radians(lat)) + 1 / math.cos(math.radians(lat))) / math.pi) / 2 * n
    return x, y
def tile_to_lat(y, z):
    return math.degrees(math.atan(math.sinh(math.pi * (1 - 2 * y / 2 ** z))))

def aerial_mosaic(ext, z):
    """download GSI photo tiles covering ext; return image resampled to a regular lon/lat grid and its extent"""
    lon0, lon1, lat0, lat1 = ext
    xa, ya = lonlat_to_tile(lon0, lat1, z); xb, yb = lonlat_to_tile(lon1, lat0, z)
    tx0, tx1, ty0, ty1 = int(xa), int(xb), int(ya), int(yb)
    tdir = CACHE / f"gsi_seamlessphoto_z{z}"; tdir.mkdir(exist_ok=True)
    W, H = (tx1 - tx0 + 1) * 256, (ty1 - ty0 + 1) * 256
    mosaic = Image.new("RGB", (W, H), (128, 128, 128))
    n_new = 0
    for tx in range(tx0, tx1 + 1):
        for ty in range(ty0, ty1 + 1):
            f = tdir / f"{tx}_{ty}.jpg"
            if not f.exists():
                data = fetch(TILE_URL.format(z=z, x=tx, y=ty)); n_new += 1
                f.write_bytes(data if data else b"")
                time.sleep(0.05)
            if f.stat().st_size:
                mosaic.paste(Image.open(f).convert("RGB"), ((tx - tx0) * 256, (ty - ty0) * 256))
    if n_new: print(f"downloaded {n_new} aerial photo tiles (zoom {z})")
    arr = np.asarray(mosaic)
    # pixel edges in lon (linear) and lat (Mercator): resample rows to a regular latitude grid
    lon_w = tx0 / 2 ** z * 360 - 180; lon_e = (tx1 + 1) / 2 ** z * 360 - 180
    lat_n = tile_to_lat(ty0, z); lat_s = tile_to_lat(ty1 + 1, z)
    lats_rows = np.array([tile_to_lat(ty0 + (i + 0.5) / 256, z) for i in range(H)])     # descending
    target = np.linspace(lat_n, lat_s, H)
    idx = np.clip(np.searchsorted(-lats_rows, -target), 0, H - 1)
    return arr[idx], (lon_w, lon_e, lat_s, lat_n)

# ------------------------------------------------------------------ drawing helpers
def deg_fmt(axis):
    def f(v, _):
        d = int(math.floor(v + 1e-9)); m = round((v - d) * 60)
        if m == 60: d, m = d + 1, 0
        hemi = "E" if axis == "x" else "N"
        return f"{d}°{m:02d}′{hemi}" if m else f"{d}°{hemi}"
    return FuncFormatter(f)

def setup(ax, ext, step_x, step_y, lat0, sea=True):
    x0, x1, y0, y1 = ext
    ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
    ax.set_aspect(1 / math.cos(math.radians(lat0)))
    if sea: ax.set_facecolor(SEA)
    ax.xaxis.set_major_locator(MultipleLocator(step_x)); ax.yaxis.set_major_locator(MultipleLocator(step_y))
    ax.xaxis.set_major_formatter(deg_fmt("x")); ax.yaxis.set_major_formatter(deg_fmt("y"))
    ax.tick_params(length=2, width=0.5, labelsize=6, pad=1.5)
    plt.setp(ax.get_yticklabels(), rotation=90, va="center")

HALO = [pe.withStroke(linewidth=2.2, foreground="black")]
def scalebar(ax, km, lat0, x_frac, y_frac, color="k", halo=None):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    dlon = km / (111.320 * math.cos(math.radians(lat0)))
    x = x0 + (x1 - x0) * x_frac; y = y0 + (y1 - y0) * y_frac
    ax.plot([x, x + dlon], [y, y], color=color, lw=1.4, solid_capstyle="butt", zorder=6, path_effects=halo)
    t = ax.text(x + dlon / 2, y + (y1 - y0) * 0.02, f"{km:g} km", ha="center", va="bottom", fontsize=6.5, color=color, zorder=6)
    if halo: t.set_path_effects(halo)

def north(ax, x_frac, y_frac, color="k", halo=None):
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    x = x0 + (x1 - x0) * x_frac; y = y0 + (y1 - y0) * y_frac; h = (y1 - y0) * 0.07
    ax.annotate("", xy=(x, y + h), xytext=(x, y), zorder=7,
                arrowprops=dict(arrowstyle="-|>,head_width=0.25,head_length=0.5", lw=1.0, color=color, path_effects=halo))
    t = ax.text(x, y + h * 1.1, "N", ha="center", va="bottom", fontsize=7.5, fontweight="bold", color=color, zorder=7)
    if halo: t.set_path_effects(halo)

def panel_label(ax, s):
    ax.text(0.02, 0.98, s, transform=ax.transAxes, ha="left", va="top", fontsize=9, fontweight="bold",
            bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="none"), zorder=8)

def frame(ax, ext):
    x0, x1, y0, y1 = ext
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, lw=0.9, ec="#c0392b", zorder=5))

# ------------------------------------------------------------------ data
ne = gpd.read_file(CACHE / "ne_10m_land.geojson")
n03 = gpd.GeoDataFrame(pd.concat([gpd.read_file(CACHE / "N03_47381.json"), gpd.read_file(CACHE / "N03_47207.json")]), crs="EPSG:4326")
n03["geometry"] = n03.geometry.make_valid()
photo, photo_ext = aerial_mosaic(EXT_C, TILE_ZOOM)

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(174 / 25.4, 128 / 25.4))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.25], height_ratios=[1.2, 1], left=0.055, right=0.995,
                      bottom=0.045, top=0.995, wspace=0.13, hspace=0.12)
axA, axB, axC = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])

setup(axA, EXT_A, 2, 2, 27.5)
ne.clip(box(*[EXT_A[i] for i in (0, 2, 1, 3)])).plot(ax=axA, color=LAND, edgecolor=LANDEDGE, lw=0.3, zorder=2)
frame(axA, EXT_B)
for t, x, y, kw in [("JAPAN", 130.3, 31.15, dict(fontsize=6.5, fontweight="bold")),
                     ("Okinawa Is.", 128.15, 26.55, dict(ha="left")), ("Amami-Oshima", 129.1, 28.75, {}),
                     ("Miyako Is.", 125.3, 25.25, {}), ("Yaeyama Is.", 123.95, 23.85, {}),
                     ("TAIWAN", 121.25, 23.75, dict(fontsize=6.5, fontweight="bold", ha="left")),
                     ("East China Sea", 124.3, 29.2, dict(style="italic", color="#3b6e8f")),
                     ("Pacific Ocean", 129.4, 25.1, dict(style="italic", color="#3b6e8f")),
                     ("Ryukyu Archipelago", 126.9, 27.25, dict(rotation=38, color="#555555"))]:
    axA.text(x, y, t, fontsize=kw.pop("fontsize", 6), ha=kw.pop("ha", "center"), va="center", zorder=6, **kw)
scalebar(axA, 200, 27.5, 0.62, 0.08)
panel_label(axA, "A")

setup(axB, EXT_B, 0.2, 0.1, 24.4)
n03.clip(box(*[EXT_B[i] for i in (0, 2, 1, 3)])).plot(ax=axB, color=LAND, edgecolor=LANDEDGE, lw=0.3, zorder=2)
frame(axB, EXT_C)
for t, x, y in [("Iriomote Is.", 123.80, 24.335), ("Ishigaki Is.", 124.235, 24.43), ("Kohama Is.", 123.975, 24.305),
                ("Taketomi Is.", 124.09, 24.355), ("Kuroshima", 124.065, 24.25)]:
    axB.text(x, y, t, fontsize=6, ha="center", va="center", zorder=6)
scalebar(axB, 10, 24.4, 0.05, 0.07)
panel_label(axB, "B")

setup(axC, EXT_C, 1 / 60, 1 / 60, 24.36, sea=False)
axC.set_anchor("S")                                       # align panel C with the bottom of panel B (shorter connector lines)
axC.imshow(photo, extent=photo_ext, origin="upper", interpolation="lanczos", zorder=1)
axC.set_xlim(EXT_C[0], EXT_C[1]); axC.set_ylim(EXT_C[2], EXT_C[3])
step = 0.0016                                            # spacing between markers of a site with >1 subclade
for code, name, lat, lon, subs, (dx, dy) in SITES:
    offs = [(i - (len(subs) - 1) / 2) * step for i in range(len(subs))]
    for o, s in zip(offs, subs):
        lab, col, mk = LINEAGES[s]
        axC.scatter(lon + o, lat, s=40, marker=mk, color=col, edgecolor="white", linewidth=0.9, zorder=7)
    t = axC.text(lon + dx, lat + dy, code, fontsize=8, fontweight="bold", color="white",
                 ha="left" if dx > 0 else "right", va="center", zorder=8)
    t.set_path_effects(HALO)
for t, x, y, kw in [("Iriomote Is.", 123.931, 24.352, {}), ("Kohama Is.", 123.978, 24.336, {}),
                    ("Yonara Channel", 123.9585, 24.3732, dict(style="italic", rotation=-62))]:
    tt = axC.text(x, y, t, fontsize=7, color="white", ha="center", va="center", zorder=6, **kw); tt.set_path_effects(HALO)
used = [k for k in LINEAGES if any(k in s[4] for s in SITES)]
handles = [plt.Line2D([], [], ls="none", marker=LINEAGES[k][2], ms=6.5, mfc=LINEAGES[k][1], mec="white", mew=0.9,
                      label=LINEAGES[k][0]) for k in used]
leg = axC.legend(handles=handles, title="Subclade sampled", title_fontsize=6.5, loc="upper left", bbox_to_anchor=(0.075, 0.995),
                 fontsize=6.5, frameon=True, framealpha=0.9, edgecolor="#999999", handletextpad=0.4, borderpad=0.5)
leg.get_frame().set_linewidth(0.5)
axC.text(0.015, 0.012, "Aerial photograph: GSI Japan", transform=axC.transAxes, ha="left", va="bottom",
         fontsize=5, color="white", zorder=8).set_path_effects([pe.withStroke(linewidth=1.5, foreground="black")])
scalebar(axC, 1, 24.36, 0.06, 0.06, color="white", halo=HALO)
north(axC, 0.93, 0.86, color="white", halo=HALO)
panel_label(axC, "C")


# connector lines: inset box in one panel -> corners of the enlarged panel
def connect(ax_from, ext, ax_to, pairs):
    """pairs: ((corner of box), (corner of target panel in axes fraction)); box corners as 'll','lr','ul','ur'"""
    x0, x1, y0, y1 = ext
    corner = {"ll": (x0, y0), "lr": (x1, y0), "ul": (x0, y1), "ur": (x1, y1)}
    for c_box, c_panel in pairs:
        fig.add_artist(ConnectionPatch(xyA=corner[c_box], coordsA=ax_from.transData, xyB=c_panel, coordsB=ax_to.transAxes,
                                       color="#c0392b", lw=0.6, ls=(0, (3, 2)), zorder=10, clip_on=False))
connect(axA, EXT_B, axB, [("ll", (0, 1)), ("lr", (1, 1))])        # A box -> top corners of B
connect(axB, EXT_C, axC, [("ur", (0, 1)), ("lr", (0, 0))])        # B box -> left corners of C
for ax in (axA, axB, axC): ax.set_xlabel(""); ax.set_ylabel("")
for ext in ("svg", "pdf"): fig.savefig(OUT / f"Fig1_sampling_sites.{ext}")
fig.savefig(OUT / "Fig1_sampling_sites.png", dpi=600)
fig.savefig(OUT / "Fig1_sampling_sites.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
(OUT / "README.txt").write_text(
    "Figure 1 - sampling sites\n\n"
    "Panel C background: GSI seamless aerial photographs (Geospatial Information Authority of Japan, "
    "https://maps.gsi.go.jp/development/ichiran.html), processed (mosaicked, resampled, annotated).\n"
    "Panels B, C coastlines: National Land Numerical Information, Administrative Areas N03 (1 Jan 2020), MLIT Japan.\n"
    "Panel A: Natural Earth 1:10m land (public domain).\n"
    "Site positions: KR, YKS, YIS = boat positions on the 13 July 2014 field sheets; "
    "2017 = settlement-plate site on the western side of the Yonara Channel (24°21.991' N, 123°57.035' E).\n")
print("done:", OUT)
