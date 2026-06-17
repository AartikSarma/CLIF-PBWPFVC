"""Working figures clarifying the strain / mechanical-power story.
Reads aggregated result CSVs (no raw data, no R). Run via:
  uv run --with pandas --with matplotlib python figures/make_strain_figures.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import pandas as pd
import os

OUT = "figures"
os.makedirs(OUT, exist_ok=True)
def p(site, stub):
    return f"output/{site}_output/final/{stub}_{site}.csv"

OK = dict(orange="#E69F00", sky="#56B4E9", green="#009E73", blue="#0072B2",
          verm="#D55E00", purple="#CC79A7", yellow="#F0E442", grey="#999999")
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

# =====================================================================
# FIG 1 — Strain is the structural (age-invariant) dose; pressure/power
#         steepen with age (the elastance ladder E^0 -> E^1 -> E^2)
# =====================================================================
rec = {s: pd.read_csv(p(s, "dose_recoil_absorption")).set_index("metric") for s in ("MIMIC", "UCSF")}
order = ["VT/PFVC", "DP", "MP/Crs"]
expo  = {"VT/PFVC": "E$^0$ (volume)", "DP": "E$^1$", "MP/Crs": "E$^2$"}
col   = {"VT/PFVC": OK["green"], "DP": OK["orange"], "MP/Crs": OK["verm"]}

fig, (a0, a1) = plt.subplots(1, 2, figsize=(11, 4.4))
m = rec["MIMIC"]
for met in order:
    a0.plot([0, 1], [m.loc[met, "slope_young"], m.loc[met, "slope_old"]],
            "-o", color=col[met], lw=2.4, ms=8, label=f"{met}  ({expo[met]})")
a0.axhline(0, color=OK["grey"], lw=.8, ls=":")
a0.set_xticks([0, 1]); a0.set_xticklabels(["younger\n(~10th pctile age)", "older\n(~90th pctile age)"])
a0.set_ylabel("per-SD mortality log-OR slope")
a0.set_title("MIMIC: dose-response slope vs age", fontsize=11)
a0.legend(frameon=False, fontsize=9.5, loc="upper left")
a0.text(.02, .02, "flat = age-invariant effect (strain)\nsteepening = age-contaminated (elastance)",
        transform=a0.transAxes, fontsize=8.5, color=OK["grey"])

x = range(len(order)); w = .38
for i, s in enumerate(("MIMIC", "UCSF")):
    vals = [abs(rec[s].loc[met, "slope_diff"]) for met in order]
    a1.bar([xx + (i - .5) * w for xx in x], vals, w,
           color=[col[m_] for m_ in order], alpha=1 - .45 * i,
           edgecolor="white", label=s)
for j, met in enumerate(order):
    a1.text(j, max(abs(rec["MIMIC"].loc[met, "slope_diff"]), abs(rec["UCSF"].loc[met, "slope_diff"])) + .01,
            expo[met], ha="center", fontsize=10)
a1.set_xticks(list(x)); a1.set_xticklabels(order)
a1.set_ylabel("|age-modification| (old − young slope)")
a1.set_title("Age-modification grows with elastance content\n(left bar MIMIC, right UCSF)", fontsize=10.5)
a1.text(.5, .92, "MP/Crs age-interaction p=0.003 (MIMIC)", transform=a1.transAxes,
        ha="center", fontsize=8.5, color=OK["verm"])
fig.suptitle("Strain (VT/PFVC) is the age-invariant dose; pressure/power are age-contaminated proxies",
             fontsize=12, y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/fig1_elastance_ladder.png", dpi=150, bbox_inches="tight")

# =====================================================================
# FIG 2 — Compliance-normalizing rescues the instrument:
#         raw MP beta_IV flips negative (invalid); MP/Crs is valid
# =====================================================================
iv = {s: pd.read_csv(p(s, "ivpolicy_exposure_diagnostic")).set_index("metric") for s in ("MIMIC", "UCSF")}
mets = ["DP", "MP", "MP/Crs", "Ers"]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, s in zip(axes, ("MIMIC", "UCSF")):
    d = iv[s]
    ypos = range(len(mets))
    for y, met in zip(ypos, mets):
        valid = bool(d.loc[met, "valid"])
        c = OK["green"] if valid else OK["verm"]
        biv, bols, F = d.loc[met, "beta_iv"], d.loc[met, "beta_ols"], d.loc[met, "first_stage_F"]
        ax.plot([bols, biv], [y, y], "-", color=c, lw=1.6, alpha=.6, zorder=1)
        ax.scatter(bols, y, s=55, facecolors="white", edgecolors=OK["grey"], zorder=2)
        ax.scatter(biv, y, s=90, color=c, zorder=3)
        ax.text(biv, y + .28, f"F={F:.0f}{'  ✓valid' if valid else '  ✗'}",
                ha="center", fontsize=8.5, color=c)
    ax.axvline(0, color=OK["grey"], lw=1, ls=":")
    ax.set_yticks(list(ypos)); ax.set_yticklabels(mets)
    ax.set_title(s, fontsize=11); ax.set_xlabel("effect on mortality (IV ●  vs  OLS ○)")
    ax.invert_yaxis()
axes[0].annotate("raw MP: IV flips\nNEGATIVE (invalid)", xy=(-0.04, 1), xytext=(0.18, 1.55),
                 fontsize=8.5, color=OK["verm"],
                 arrowprops=dict(arrowstyle="->", color=OK["verm"]))
axes[0].annotate("MP/Crs: IV positive,\nconcordant (valid)", xy=(0.86, 2), xytext=(0.30, 2.6),
                 fontsize=8.5, color=OK["green"],
                 arrowprops=dict(arrowstyle="->", color=OK["green"]))
fig.suptitle("Compliance-normalizing rescues the mechanical-power instrument (height IV)", fontsize=12, y=1.0)
fig.tight_layout(); fig.savefig(f"{OUT}/fig2_mpcrs_rescue.png", dpi=150, bbox_inches="tight")

# =====================================================================
# FIG 3 — Same protocol (VT/PBW 6-8), different delivered strain (VT/PFVC):
#         demographically shifted; ~70% of the variance is mis-sizing
# =====================================================================
sp = pd.read_csv(p("MIMIC", "dose_within_band_spread"))
vd = pd.read_csv(p("MIMIC", "dose_variance_decomposition"))
ax_order = ["Sex", "Race", "Age tertile", "Height tertile (w/in sex)"]
fig, (b0, b1) = plt.subplots(1, 2, figsize=(12, 4.4), gridspec_kw={"width_ratios": [3, 1]})
yt, ylab = [], []
y = 0
axc = {"Sex": OK["orange"], "Race": OK["sky"], "Age tertile": OK["green"],
       "Height tertile (w/in sex)": OK["blue"]}
for ax_name in ax_order:
    rows = sp[sp.subgroup == ax_name]
    for _, r in rows.iterrows():
        b0.plot([r.vtpfvc_p5, r.vtpfvc_p95], [y, y], color=axc[ax_name], lw=4, alpha=.5, solid_capstyle="round")
        b0.scatter(r.vtpfvc_median, y, s=45, color=axc[ax_name], zorder=3)
        yt.append(y); ylab.append(f"{r.level}")
        y += 1
    y += .6
b0.axvline(11, color=OK["verm"], lw=1.5, ls="--")
b0.text(11, y - .3, " 11% harm cut", color=OK["verm"], fontsize=9, va="top")
b0.set_yticks(yt); b0.set_yticklabels(ylab, fontsize=9)
b0.invert_yaxis(); b0.set_xlabel("Delivered strain  VT/PFVC (% predicted FVC)  — bar = p5–p95, dot = median")
b0.set_title("Same VT/PBW 6–8, ~2× spread in strain, shifted by body size", fontsize=10.5)

g = vd[vd.cohort == "guideline_6_8"].iloc[0]
shares = [g.share_missizing, g.share_protocol, g.share_cov]
labels = ["PBW/PFVC\nmis-sizing", "clinician's\nVT/PBW", "cov"]
cols = [OK["verm"], OK["sky"], OK["grey"]]
b1.bar([0], [shares[0]], .6, color=cols[0], label=labels[0])
b1.bar([0], [shares[1]], .6, bottom=shares[0], color=cols[1], label=labels[1])
b1.bar([0], [shares[2]], .6, bottom=shares[0] + shares[1], color=cols[2], label=labels[2])
b1.text(0, shares[0] / 2, f"{shares[0]*100:.0f}%", ha="center", va="center", color="white", fontweight="bold")
b1.text(0, shares[0] + shares[1] / 2, f"{shares[1]*100:.0f}%", ha="center", va="center", color="white")
b1.set_xticks([]); b1.set_ylim(0, 1); b1.set_ylabel("share of delivered-dose variance")
b1.set_title("What drives the\nstrain variation", fontsize=10)
b1.legend(frameon=False, fontsize=8, loc="upper right", bbox_to_anchor=(1.05, 1))
fig.suptitle("'Same treatment' (VT/PBW) delivers a different physiologic dose (strain), set mostly by body size",
             fontsize=12, y=1.02)
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_same_protocol_diff_strain.png", dpi=150, bbox_inches="tight")

# =====================================================================
# FIG 4 — Conceptual schematic (Dreyfuss): strain is the VILI mediator;
#         pressure/power are strain x elastance^k, so age (via E) decouples them
# =====================================================================
fig, ax = plt.subplots(figsize=(11, 5.2)); ax.axis("off"); ax.set_xlim(0, 12); ax.set_ylim(0, 7)
def box(x, y, w, h, text, fc, tc="black", fs=10):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                 fc=fc, ec="none"))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc, wrap=True)
def arrow(x0, y0, x1, y1, c="black"):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16, color=c, lw=1.6))

# mediator chain
box(0.3, 4.5, 1.9, 1.0, "Tidal volume\n(clinician sets\nVT/PBW)", OK["sky"], "white")
box(2.9, 4.5, 1.9, 1.0, "STRAIN\nVT / PFVC\n(volume÷size)", OK["green"], "white", 10.5)
box(5.5, 4.5, 1.9, 1.0, "biotrauma\n(VILI)", OK["yellow"])
box(8.1, 4.5, 1.9, 1.0, "mortality", OK["grey"], "white")
for x0, x1 in [(2.2, 2.9), (4.8, 5.5), (7.4, 8.1)]:
    arrow(x0, 5.0, x1, 5.0)
ax.text(6.1, 6.2, "Dreyfuss: VOLUME (strain), not pressure, is the VILI mediator → its effect is age-invariant",
        ha="center", fontsize=10.5, style="italic")
ax.text(3.85, 3.9, "age-invariant", ha="center", fontsize=8.5, color=OK["green"])

# the elastance equation + decoupling
box(0.3, 1.6, 4.4, 1.5,
    "stress = E$_{spec}$ × strain\nDP ∝ E$^1$ × strain\nMP/Crs ∝ E$^2$ × strain$^2$",
    "#f0f0f0", fs=11)
arrow(2.5, 1.6, 2.5, 0.95, OK["verm"])
box(0.3, 0.1, 4.4, 0.8, "pressure & power carry elastance E", OK["orange"], "white", 9.5)

box(6.0, 1.6, 5.6, 1.5,
    "Elastic recoil (E) changes with age.\nSo metrics that carry E (DP, MP/Crs)\nhave an age-CONTAMINATED effect;\nstrain (E$^0$) does not.",
    "#fde7e1", fs=10)
arrow(4.7, 2.35, 6.0, 2.35, OK["grey"])
ax.text(8.8, 0.5, "→ dose by STRAIN (PFVC-anchored), not by pressure/power targets",
        ha="center", fontsize=9.5, color=OK["green"], fontweight="bold")
fig.suptitle("Why strain is the right dosing target and pressure/power are age-contaminated", fontsize=12.5, y=0.98)
fig.savefig(f"{OUT}/fig4_dreyfuss_schematic.png", dpi=150, bbox_inches="tight")

print("wrote:", ", ".join(sorted(os.listdir(OUT))))
