"""Robustness figures for the longitudinal CCW target-trial emulation (script 10).
Reads aggregated result CSVs only (no raw data, no R). Two panels:
  A. Subgroup forest of 60-day mortality RD (+ bootstrap CI), ordered to show
     the age/height equity gradient.
  B. Tornado of the cap / ceiling / grace / deviation-rule sensitivities,
     one bar per knob (others held at the primary specification).

Run (defaults to MIMIC; pass a site to override):
  uv run --with pandas --with matplotlib python figures/make_tte_robustness_figures.py MIMIC
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np
import sys, os

SITE = sys.argv[1] if len(sys.argv) > 1 else "MIMIC"
OUT = "figures"
os.makedirs(OUT, exist_ok=True)

def p(stub):
    return f"output/{SITE}_output/final/{stub}_{SITE}.csv"

OK = dict(orange="#E69F00", sky="#56B4E9", green="#009E73", blue="#0072B2",
          verm="#D55E00", purple="#CC79A7", yellow="#F0E442", grey="#999999")
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

# Primary specification (matches script 10 knobs: C_LOW=11, C_HIGH=16, GRACE=2,
# DAYW_CAP=5, rule="simple"). RD<0 = strain-limiting protective.
C_LOW, C_HIGH, GRACE, CAP = 11, 16, 2, 5
PROTECT = OK["green"]   # RD < 0  (strain-limiting reduces mortality)
HARM    = OK["verm"]    # RD > 0
def sgn_color(rd):
    return PROTECT if rd < 0 else HARM

overall = pd.read_csv(p("tte_ccw_overall")).iloc[0]
prim_rd = float(overall["rd"])

# =====================================================================
# FIG A — Subgroup forest (equity gradient)
# =====================================================================
sg = pd.read_csv(p("tte_ccw_subgroup"))
# Display order: gradient-revealing within each family.
order = [
    ("Age tertile",              ["Young", "Middle", "Old"]),
    ("Height tertile (w/in sex)", ["Short", "Middle", "Tall"]),
    ("Sex",                      ["Male", "Female"]),
    ("Race",                     ["WHITE", "BLACK", "OTHER"]),
]
fam_color = {"Age tertile": OK["green"], "Height tertile (w/in sex)": OK["blue"],
             "Sex": OK["orange"], "Race": OK["sky"]}
fam_label = {"Age tertile": "Age tertile", "Height tertile (w/in sex)": "Height tertile (within sex)",
             "Sex": "Sex", "Race": "Race"}

rows = []   # (y, label, rd, lo, hi, n, color, is_header)
y = 0.0
# Overall reference at top
rows.append((y, "OVERALL", float(overall["rd"]), float(overall["rd_lo"]),
             float(overall["rd_hi"]), int(overall["n_patients"]), OK["grey"], "overall"))
y -= 1.4
for fam, levels in order:
    rows.append((y, fam_label[fam], None, None, None, None, fam_color[fam], "header"))
    y -= 1.0
    sub = sg[sg.subgroup == fam].set_index("level")
    for lv in levels:
        if lv not in sub.index:
            continue
        r = sub.loc[lv]
        rows.append((y, f"   {lv}", float(r.rd), float(r.rd_lo), float(r.rd_hi),
                     int(r.n), fam_color[fam], "level"))
        y -= 1.0
    y -= 0.5

fig, ax = plt.subplots(figsize=(8.6, 0.42 * len([r for r in rows]) + 1.2))
for (yy, lab, rd, lo, hi, n, col, kind) in rows:
    if kind == "header":
        ax.text(-0.001, yy, lab, ha="right", va="center", fontsize=10.5,
                fontweight="bold", color=col, transform=ax.get_yaxis_transform())
        continue
    ec = sgn_color(rd) if kind != "overall" else "black"
    ax.plot([lo, hi], [yy, yy], "-", color=ec, lw=2.0, alpha=.85, zorder=2)
    marker = "D" if kind == "overall" else "o"
    ms = 10 if kind == "overall" else 7
    ax.scatter(rd, yy, s=ms**2, marker=marker, color=ec, zorder=3,
               edgecolors="black" if kind == "overall" else "white", linewidths=.8)
    ax.text(-0.001, yy, lab, ha="right", va="center", fontsize=10,
            fontweight="bold" if kind == "overall" else "normal",
            transform=ax.get_yaxis_transform())
    ax.text(1.002, yy, f"{rd*100:+.1f} ({lo*100:+.1f}, {hi*100:+.1f})  n={n}",
            ha="left", va="center", fontsize=8.5, color=OK["grey"],
            transform=ax.get_yaxis_transform())

ax.axvline(0, color="black", lw=1, ls="-", alpha=.6)
ax.set_yticks([])
ax.set_ylim(min(r[0] for r in rows) - 1, 1.0)
ax.set_xlabel("60-day mortality risk difference  (strain-limiting − permissive), percentage points",
              fontsize=10)
xmax = max(abs(sg.rd_lo.min()), abs(sg.rd_hi.max()), abs(prim_rd)) * 1.15
ax.set_xlim(-xmax, xmax)
xticks = ax.get_xticks()
ax.set_xticks(xticks)
ax.set_xticklabels([f"{t*100:+.0f}" for t in xticks])
ax.text(0.0, 1.0, "", transform=ax.transAxes)
ax.annotate("◀ strain-limiting protective", xy=(0.0, 1.0), xytext=(0.02, 1.02),
            xycoords="axes fraction", fontsize=9, color=PROTECT, ha="left", va="bottom")
ax.annotate("harm ▶", xy=(1.0, 1.0), xytext=(0.98, 1.02),
            xycoords="axes fraction", fontsize=9, color=HARM, ha="right", va="bottom")
fig.suptitle(f"Strain-limiting ventilation: subgroup mortality risk differences ({SITE})",
             fontsize=12.5, y=1.0)
fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.99))
fig.savefig(f"{OUT}/fig5_tte_subgroup_forest_{SITE}.png", dpi=150, bbox_inches="tight")

# =====================================================================
# FIG B — Tornado of design-choice sensitivities
# =====================================================================
wc = pd.read_csv(p("tte_ccw_sens_weightcap"))
cg = pd.read_csv(p("tte_ccw_sens_ceiling_grace"))
ru = pd.read_csv(p("tte_ccw_sens_rule"))

def rng(vals):
    vals = [float(v) for v in vals]
    return min(vals), max(vals)

bars = []  # (label, lo, hi, note)

# Weight cap: finite caps only (Inf is degenerate — noted separately)
wc_fin = wc[np.isfinite(wc.weight_cap)]
lo, hi = rng(wc_fin.rd)
bars.append(("IPCW weight cap\n(3 / 5 / 10)", lo, hi,
             f"∞ degenerate: {float(wc[~np.isfinite(wc.weight_cap)].rd.iloc[0])*100:+.1f}"
             if (~np.isfinite(wc.weight_cap)).any() else ""))

# Grace window: hold ceilings at primary, vary grace
g = cg[(cg.c_low == C_LOW) & (cg.c_high == C_HIGH)]
lo, hi = rng(g.rd)
bars.append((f"Grace window\n({'/'.join(str(int(x)) for x in sorted(g.grace.unique()))} d)", lo, hi, ""))

# Permissive ceiling c_high: hold c_low + grace at primary
ch = cg[(cg.c_low == C_LOW) & (cg.grace == GRACE)]
lo, hi = rng(ch.rd)
bars.append((f"Permissive ceiling\n(c_high {'/'.join(str(int(x)) for x in sorted(ch.c_high.unique()))}%)",
             lo, hi, ""))

# Strain-limiting ceiling c_low: hold c_high + grace at primary
cl = cg[(cg.c_high == C_HIGH) & (cg.grace == GRACE)]
lo, hi = rng(cl.rd)
bars.append((f"Strain-limit ceiling\n(c_low {'/'.join(str(int(x)) for x in sorted(cl.c_low.unique()))}%)",
             lo, hi, ""))

# Deviation rule: simple vs corrected
lo, hi = rng(ru.rd)
bars.append(("Deviation rule\n(simple / corrected)", lo, hi, ""))

# Classic tornado: widest range at top
bars.sort(key=lambda b: (b[2] - b[1]))   # ascending; plotted bottom→top so widest on top
fig, ax = plt.subplots(figsize=(9.2, 4.6))
for i, (lab, lo, hi, note) in enumerate(bars):
    span_protect = hi <= 0
    span_harm = lo >= 0
    c = PROTECT if span_protect else (HARM if span_harm else OK["grey"])
    ax.barh(i, hi - lo, left=lo, height=.62, color=c, alpha=.55, edgecolor=c, lw=1.2, zorder=2)
    ax.plot([lo, lo], [i - .31, i + .31], color=c, lw=1.6)
    ax.plot([hi, hi], [i - .31, i + .31], color=c, lw=1.6)
    ax.text(lo, i + .42, f"{lo*100:+.1f}", ha="center", va="bottom", fontsize=7.8, color=OK["grey"])
    ax.text(hi, i + .42, f"{hi*100:+.1f}", ha="center", va="bottom", fontsize=7.8, color=OK["grey"])
    if note:
        ax.text(0.995, i, note, transform=ax.get_yaxis_transform(), ha="right",
                va="center", fontsize=7.5, color=OK["grey"], style="italic")

ax.axvline(prim_rd, color="black", lw=1.4, ls="--", zorder=3)
ax.annotate(f"primary spec  {prim_rd*100:+.1f} pp", xy=(prim_rd, len(bars) - .55),
            xytext=(prim_rd, len(bars) - .15), ha="center", va="bottom", fontsize=8.5,
            fontweight="bold", bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="black", lw=.8),
            arrowprops=dict(arrowstyle="-", color="black", lw=1.0))
ax.axvline(0, color=OK["grey"], lw=.9, ls=":")
ax.set_yticks(range(len(bars)))
ax.set_yticklabels([b[0] for b in bars], fontsize=9)
ax.set_xlabel("60-day mortality risk difference (percentage points) as the design knob is varied", fontsize=10)
xs = [b[1] for b in bars] + [b[2] for b in bars] + [prim_rd, 0]
pad = (max(xs) - min(xs)) * 0.12 + 1e-4
ax.set_xlim(min(xs) - pad, max(xs) + pad)
xticks2 = ax.get_xticks()
ax.set_xticks(xticks2)
ax.set_xticklabels([f"{t*100:+.1f}" for t in xticks2])
ax.set_ylim(-.6, len(bars) + .35)
fig.suptitle(f"Design-choice robustness of the strain-limiting effect ({SITE})", fontsize=12.5, y=0.99)
ax.text(.5, -.22, "Each bar spans the RD range as one knob varies, others held at the primary spec.",
        transform=ax.transAxes, ha="center", fontsize=8.5, color=OK["grey"])
fig.tight_layout(rect=(0.0, 0.02, 1.0, 0.99))
fig.savefig(f"{OUT}/fig6_tte_tornado_{SITE}.png", dpi=150, bbox_inches="tight")

print("wrote:", f"fig5_tte_subgroup_forest_{SITE}.png,", f"fig6_tte_tornado_{SITE}.png",
      f"(site={SITE}, primary RD={prim_rd*100:+.2f} pp)")
