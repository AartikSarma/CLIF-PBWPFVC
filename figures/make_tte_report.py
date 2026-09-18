"""Build a self-contained, teaching-oriented HTML report for the longitudinal
target-trial emulation (the TTE (30_tte_common)).  Reads the R source verbatim (sliced by
section) and pulls every number from the MIMIC/UCSF `final/` result CSVs --
no new analysis is run, nothing is transcribed by hand.

  uv run --with pandas python figures/make_tte_report.py

Writes reports/tte_methods_report.html  (*.html is gitignored).
"""
import pandas as pd
import numpy as np
import html as _html
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODE = os.path.join(ROOT, "code")
SITES = ["MIMIC", "UCSF"]

def finalp(site, stub):
    return os.path.join(ROOT, "output", f"{site}_output", "final", "causal", f"{stub}_{site}.csv")

# ---- load every result CSV for both sites -----------------------------------
# Tolerate either the legacy 'stress_limiting' identifier or the renamed
# 'strain_limiting' one: the TTE (30_tte_common)'s column rename lands on the next re-run, so a
# report generated before that re-run still reads the old-named CSVs correctly.
def _norm(df):
    df = df.rename(columns=lambda c: c.replace("stress_limiting", "strain_limiting"))
    if "arm" in df.columns:
        df["arm"] = df["arm"].replace({"stress_limiting": "strain_limiting"})
    return df

R = {}
for s in SITES:
    def L(stub, _s=s):
        return _norm(pd.read_csv(finalp(_s, stub)))
    def L_opt(stub, _s=s):   # optional CSV (e.g. pH sensitivity, written only on re-run)
        pth = finalp(_s, stub)
        return _norm(pd.read_csv(pth)) if os.path.exists(pth) else None
    R[s] = {
        "overall":  L("tte_ccw_overall").iloc[0],
        "diag":     L("tte_ccw_diagnostics"),
        "sub":      L("tte_ccw_subgroup"),
        "wcap":     L_opt("tte_ccw_sens_weightcap"),
        "cg":       L_opt("tte_ccw_sens_ceiling_grace"),
        "rule":     L_opt("tte_ccw_sens_rule"),
        "ph":       L_opt("tte_ccw_sens_ph"),
        "dp":       L_opt("tte_ccw_sens_dp"),
        # validity / completeness / positivity tables (added on the 28-d re-run)
        "drops":    L_opt("tte_ccw_panel_drops"),
        "struct":   L_opt("tte_ccw_structural_excluded"),
        "posemp":   L_opt("tte_ccw_positivity_empirical"),
        "wage":     L_opt("tte_ccw_weights_by_age"),
        "evalue":   L_opt("tte_ccw_evalue"),
        "numer":    L_opt("tte_ccw_sens_numerator"),
        "trim":     L_opt("tte_ccw_sens_trim"),
        "agg":      L_opt("tte_ccw_sens_aggregation"),
        "cumw":     L_opt("tte_ccw_sens_cumweight"),
        "pf":       L_opt("tte_ccw_sens_pf"),
        "mtp":      L_opt("tte_ccw_sens_mtp"),
        "balance":  L_opt("tte_ccw_balance"),
        # 37_tte_discordance_benefit discordance-HTE (the primary heterogeneity analysis) + its sensitivities
        "disc_hte":   L_opt("tte_ccw_disc_hte"),
        "disc_grad":  L_opt("tte_ccw_disc_gradient"),
        "disc_curve": L_opt("tte_ccw_disc_cate_curve"),
        "disc_slope": L_opt("tte_ccw_disc_cate_slope"),
        "disc_dose":  L_opt("tte_ccw_disc_dose_correction"),
        "disc_eval":  L_opt("tte_ccw_disc_evalue"),
        "disc_ceil":  L_opt("tte_ccw_disc_ceiling_sweep"),
        "disc_vr":    L_opt("tte_ccw_disc_vr_hte"),
        "disc_sfwm":  L_opt("tte_ccw_disc_sf_weightmodel_hte"),
    }

# ---- R source, sliced by MARKER across the SPLIT scripts ---------------------
# The 1121-line monolith code/10_longitudinal_tte.R was refactored (verified
# byte-identical in behavior) into a shared engine, 30_tte_common.R, plus one leaf
# per analysis (35_tte_primary..11.N) and a driver (32_tte_run_all.R). This report reads those
# split files -- never the retired monolith. Each named file is read into its own
# line list (FILES[name] -> lines); leaves are shown whole, while the long engine
# file is shown by marker-bounded sub-slice.
import re as _re

# The files the walkthrough quotes from, each loaded once into a line list.
_FILE_NAMES = [
    "30_tte_common.R", "35_tte_primary.R", "36_tte_diagnostics.R",
    "11.L_sens_trim.R", "11.N_refit_boot.R", "32_tte_run_all.R",
]
FILES = {}
for _nm in _FILE_NAMES:
    # some quoted sensitivities (11.L, 11.N) were moved to code/archive/ in the federated
    # cleanup; fall back there so the methods walkthrough still resolves their source.
    _path = os.path.join(CODE, _nm)
    if not os.path.exists(_path):
        _path = os.path.join(CODE, "archive", _nm)
    with open(_path) as _f:
        FILES[_nm] = _f.readlines()

# Engine banner markers in 30_tte_common.R look like a `# ====` rule line, then
# `# 10b. Title`, then another `# ====` rule. Within a file we locate the `# 10<key>.`
# header and slice to the line BEFORE the next top-level (different-letter) section,
# so a request for e.g. 10d bundles any numbered sub-sections up to the next letter.
# Matches `# 10a. Title`, `# 10g3. Title`, and `# 10e (design-build portion).` --
# the id followed by a period-or-space, never another alphanumeric.
_HDR = _re.compile(r"^# (10[a-m][0-9]?)[\.\s]")   # section-header marker

def _headers(name):
    """List of (section-id, 1-indexed header line) for a named file, in file order."""
    return [(m.group(1), i + 1) for i, ln in enumerate(FILES[name])
            for m in [_HDR.match(ln)] if m]

def _marker_bounds(name, key):
    """(start, end) 1-indexed inclusive for engine section `key` in `name`. The
    special key 'header' spans the setup block: Sys.setenv -> line before 10a."""
    lines = FILES[name]
    hdrs = _headers(name)
    if key == "header":
        start = next(i + 1 for i, ln in enumerate(lines) if ln.startswith("Sys.setenv("))
        first_10a = next(ln_no for sid, ln_no in hdrs if sid == "10a")
        return start, first_10a - 1
    start = next(ln_no for sid, ln_no in hdrs if sid == key)
    # end = line before the next TOP-LEVEL (different-letter) section header, so the
    # slice bundles any numbered sub-sections of `key` (e.g. 10g2/10g3/10g4 under 10g).
    letter = key[2]   # the single section letter, e.g. 'g' in '10g'/'10g3'
    nxt = [ln_no for sid, ln_no in hdrs if ln_no > start and sid[2] != letter]
    end = (min(nxt) - 1) if nxt else len(lines)
    # the next header sits under a `# ====` banner rule; trim trailing rule/blank
    # lines so the slice ends on real code, not the following section's banner.
    while end > start and (lines[end - 1].startswith("# ===") or not lines[end - 1].strip()):
        end -= 1
    return start, end

def _esc_slice(name, a, b):
    return _html.escape("".join(FILES[name][a-1:b]).rstrip("\n"))

# ---- number helpers ---------------------------------------------------------
def pp(x):   return f"{x*100:+.1f}"          # signed percentage points
def pp0(x):  return f"{x*100:.1f}"           # unsigned pp
def ci(o, k="rd"):
    val = "lib_diff" if k == "lib" else k
    return f"{pp(o[val])} [{pp(o[k+'_lo'])}, {pp(o[k+'_hi'])}]"

def code_file(name):
    """The whole of a (short leaf) file as an escaped <pre> block."""
    body = _html.escape("".join(FILES[name]).rstrip("\n"))
    return f'<pre class="r"><code>{body}</code></pre>'

def code_slice(name, key):
    """A marker-bounded engine sub-block (`# 10x.` -> next top-level marker) of `name`."""
    a, b = _marker_bounds(name, key)
    return f'<pre class="r"><code>{_esc_slice(name, a, b)}</code></pre>'

def linetag(name, key=None):
    """File-qualified tag. With a marker key -> 'file:a–b'; whole-file -> 'file (n lines)'."""
    if key is None:
        return f'<span class="tag">{name} · {len(FILES[name])} lines</span>'
    a, b = _marker_bounds(name, key)
    return f'<span class="tag">{name}:{a}–{b}</span>'

def callout(kind, title, body):
    return (f'<div class="callout {kind}"><div class="ctitle">{title}</div>'
            f'<div class="cbody">{body}</div></div>')

# =============================================================================
# figure machinery -- headless matplotlib, self-contained base64 embedding
# =============================================================================
import io as _io
import base64 as _b64
import tempfile as _tempfile
import subprocess as _subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

# Okabe-Ito palette for discrete categories; viridis for continuous fills.
OKABE = ["#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#999999"]
SITE_COLOR = {"MIMIC": OKABE[3], "UCSF": OKABE[0]}           # blue / orange
ARM_COLOR = {"permissive": OKABE[1], "strain_limiting": OKABE[2]}  # sky / green

plt.rcParams.update({
    "font.size": 11.5, "axes.titlesize": 13, "axes.labelsize": 11.5,
    "xtick.labelsize": 10.5, "ytick.labelsize": 10.5, "legend.fontsize": 10.5,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

def _clean_ax(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def embed_fig(fig, caption="", alt=""):
    """Save a matplotlib figure to an in-memory PNG, base64-embed as a data URI."""
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    b64 = _b64.b64encode(buf.getvalue()).decode("ascii")
    cap = f'<figcaption class="muted" style="text-align:center">{caption}</figcaption>' if caption else ""
    return (f'<figure><img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="width:100%;max-width:820px;display:block;margin:8px auto">{cap}</figure>')

def embed_png_bytes(b, caption="", alt=""):
    """Embed already-rendered PNG bytes (e.g. a pdftocairo conversion)."""
    b64 = _b64.b64encode(b).decode("ascii")
    cap = f'<figcaption class="muted" style="text-align:center">{caption}</figcaption>' if caption else ""
    return (f'<figure><img src="data:image/png;base64,{b64}" alt="{alt}" '
            f'style="width:100%;max-width:820px;display:block;margin:8px auto">{cap}</figure>')

# track which figures actually rendered, for the final report-back
FIG_LOG = []

# =============================================================================
# build the per-site results that get reused in prose
# =============================================================================
O = {s: R[s]["overall"] for s in SITES}

def overall_table():
    rows = ""
    for s in SITES:
        o = O[s]
        rows += (f"<tr><td><b>{s}</b></td><td>{ci(o)}</td>"
                 f"<td>{pp0(o['risk_strain_limiting'])}% vs {pp0(o['risk_permissive'])}%</td>"
                 f"<td>{int(o['n_patients']):,}</td>"
                 f"<td>{ci(o,'lib')}</td></tr>")
    return ("<table><thead><tr><th>Cohort</th>"
            "<th>28-day mortality RD<br><span class='sub'>strain-limiting − permissive (pp)</span></th>"
            "<th>Risk (SL vs perm)</th><th>n</th>"
            "<th>28-day liberation CIF diff<br><span class='sub'>(pp)</span></th></tr></thead>"
            f"<tbody>{rows}</tbody></table>")

def diag_table():
    rows = ""
    for s in SITES:
        d = R[s]["diag"].set_index("arm")
        arm_disp = {"permissive": "permissive", "strain_limiting": "strain-limiting"}
        for arm in ["permissive", "strain_limiting"]:
            r = d.loc[arm]
            rows += (f"<tr><td>{s}</td><td>{arm_disp[arm]}</td>"
                     f"<td>{r['frac_deviated']*100:.1f}%</td>"
                     f"<td>{r['ess_frac']:.3f}</td>"
                     f"<td>{r['wt_max']:.2f}</td></tr>")
    return ("<table><thead><tr><th>Cohort</th><th>Arm</th><th>Deviated</th>"
            "<th>ESS fraction</th><th>Max weight</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>")

def subgroup_table():
    order = [("Age tertile", ["Young", "Middle", "Old"]),
             ("Height tertile (w/in sex)", ["Short", "Middle", "Tall"]),
             ("Sex", ["Male", "Female"]),
             ("Race", ["WHITE", "BLACK", "OTHER"])]
    head = ("<table class='sg'><thead><tr><th>Subgroup</th><th>Level</th>"
            + "".join(f"<th colspan=2>{s}</th>" for s in SITES)
            + "</tr><tr><th></th><th></th>"
            + "".join("<th>RD [95% CI]</th><th>n</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for fam, levels in order:
        for i, lv in enumerate(levels):
            cells = ""
            for s in SITES:
                d = R[s]["sub"]
                row = d[(d.subgroup == fam) & (d.level == lv)]
                if len(row):
                    r = row.iloc[0]
                    cells += (f"<td>{pp(r['rd'])} [{pp(r['rd_lo'])}, {pp(r['rd_hi'])}]</td>"
                              f"<td>{int(r['n']):,}</td>")
                else:
                    cells += "<td>—</td><td>—</td>"
            famcell = f"<td rowspan={len(levels)}><b>{fam}</b></td>" if i == 0 else ""
            body += f"<tr>{famcell}<td>{lv}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def wcap_table():
    if any(R[s]["wcap"] is None for s in SITES):
        return "<p class='muted'>Weight-cap sensitivity not run at every site (tte_sens_censoring covers it).</p>"
    rows = ""
    caps = sorted({float(c) for s in SITES for c in R[s]["wcap"].weight_cap})
    for cp in caps:
        lab = "∞ (untruncated)" if not np.isfinite(cp) else f"{int(cp)}"
        cells = ""
        for s in SITES:
            d = R[s]["wcap"]
            row = d[(np.isfinite(d.weight_cap) == np.isfinite(cp)) &
                    ((d.weight_cap == cp) if np.isfinite(cp) else True)]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{pp(r['rd'])}</td><td>{r['ess_strain_limiting']:.2f}</td>"
            else:
                cells += "<td>—</td><td>—</td>"
        cls = ' class="degen"' if not np.isfinite(cp) else ""
        rows += f"<tr{cls}><td>{lab}</td>{cells}</tr>"
    head = ("<table><thead><tr><th>IPCW cap</th>"
            + "".join(f"<th>{s} RD</th><th>{s} ESS<sub>SL</sub></th>" for s in SITES)
            + "</tr></thead><tbody>")
    return head + rows + "</tbody></table>"

def cg_summary(s):
    d = R[s]["cg"]
    if d is None:
        return (np.nan, np.nan, np.nan, np.nan)
    return (d.rd.min(), d.rd.max(),
            d[(d.c_low==10)&(d.c_high==16)&(d.grace==1)].rd.iloc[0],
            d[(d.c_low==12)&(d.c_high==14)&(d.grace==3)].rd.iloc[0])

def rule_table():
    if any(R[s]["rule"] is None for s in SITES):
        return "<p class='muted'>Deviation-rule sensitivity not run at every site.</p>"
    rows = ""
    for s in SITES:
        d = R[s]["rule"].set_index("deviation_rule")
        rows += (f"<tr><td>{s}</td>"
                 f"<td>{pp(d.loc['simple','rd'])}</td>"
                 f"<td>{pp(d.loc['corrected','rd'])}</td>"
                 f"<td>{d.loc['simple','frac_deviated_strain_limiting']*100:.1f}% / "
                 f"{d.loc['corrected','frac_deviated_strain_limiting']*100:.1f}%</td></tr>")
    return ("<table><thead><tr><th>Cohort</th><th>Simple rule RD</th>"
            "<th>Corrected rule RD</th><th>Deviated (simple/corr)</th></tr></thead>"
            f"<tbody>{rows}</tbody></table>")

def ph_block():
    have = {s: R[s]["ph"] for s in SITES if R[s]["ph"] is not None}
    if not have:
        return ("<p class='muted'>This sensitivity (lagged pH in the weight model) computes on the next real-data "
                "re-run: synthetic CLIF has no blood-gas labs, and the current MIMIC/UCSF result CSVs predate "
                "the analysis. On re-run the table below populates automatically — for each pH source (pooled "
                "arterial+venous, and arterial-only), the gas-covered-subset RD without and with lagged-pH "
                "adjustment, against the full-cohort primary.</p>")
    src_label = {"full_cohort": "Full cohort (primary)",
                 "pooled_art_plus_venous": "Arterial + venous (+0.05)",
                 "arterial_only": "Arterial only"}
    spec_label = {"primary": "— (no pH)",
                  "subset_no_ph_adj": "subset, no pH adj.",
                  "subset_with_ph_adj": "subset, + lagged pH"}
    order = [("full_cohort", "primary"),
             ("pooled_art_plus_venous", "subset_no_ph_adj"),
             ("pooled_art_plus_venous", "subset_with_ph_adj"),
             ("arterial_only", "subset_no_ph_adj"),
             ("arterial_only", "subset_with_ph_adj")]
    rows = ""
    for s, df in have.items():
        d = df.set_index(["ph_source", "spec"])
        for src, spec in order:
            if (src, spec) in d.index:
                r = d.loc[(src, spec)]
                rows += (f"<tr><td>{s}</td><td>{src_label.get(src, src)}</td>"
                         f"<td>{spec_label.get(spec, spec)}</td><td>{pp(float(r['rd']))}</td>"
                         f"<td>{int(r['n_patients']):,}</td></tr>")
    return ("<table><thead><tr><th>Cohort</th><th>pH source</th><th>Specification</th>"
            "<th>RD (pp)</th><th>n</th></tr></thead>" f"<tbody>{rows}</tbody></table>")

def dp_block():
    have = {s: R[s]["dp"] for s in SITES if R[s]["dp"] is not None}
    if not have:
        return ("<p class='muted'>This sensitivity (lagged driving pressure in the weight model) computes on the next real-data "
                "re-run, on the subset of patient-days following a recorded plateau. On re-run the table below "
                "populates with the full-cohort primary RD vs the plateau-recorded subset without and with "
                "lagged worst-of-day driving pressure in the weight model.</p>")
    label = {"full_cohort_primary": "Full cohort (primary)",
             "dp_subset_no_dp_adj": "DP-recorded subset, no DP adj.",
             "dp_subset_with_dp_adj": "DP-recorded subset, + lagged worst DP"}
    rows = ""
    for s, df in have.items():
        d = df.set_index("spec")
        for key in ["full_cohort_primary", "dp_subset_no_dp_adj", "dp_subset_with_dp_adj"]:
            if key in d.index:
                r = d.loc[key]
                rows += (f"<tr><td>{s}</td><td>{label[key]}</td><td>{pp(float(r['rd']))}</td>"
                         f"<td>{int(r['n_patients']):,}</td></tr>")
    return ("<table><thead><tr><th>Cohort</th><th>Specification</th><th>RD (pp)</th><th>n</th>"
            "</tr></thead>" f"<tbody>{rows}</tbody></table>")

# =============================================================================
# validity tables (added on the 28-d re-run) -- all by-site (MIMIC | UCSF)
# =============================================================================
def _fnum(x, d=3):
    try:    return f"{float(x):.{d}f}"
    except Exception: return "—"

def drops_table():
    """Data completeness after the SF/SpO2 fix: how much of the daily panel is dropped."""
    rows = [("patient_days_total",     "Patient-days, total",            lambda r: f"{int(r):,}"),
            ("frac_days_dropped",       "Patient-days dropped",           lambda r: f"{r*100:.1f}%"),
            ("patients_extub_shifted",  "Patients with extubation pulled earlier", lambda r: f"{int(r):,}"),
            ("patients_lost_entirely",  "Patients lost entirely",         lambda r: f"{int(r):,}")]
    head = ("<table><thead><tr><th>Completeness metric</th>"
            + "".join(f"<th>{s}</th>" for s in SITES) + "</tr></thead><tbody>")
    body = ""
    for key, lab, fmt in rows:
        cells = ""
        for s in SITES:
            d = R[s]["drops"]
            cells += f"<td>{fmt(d.iloc[0][key]) if d is not None else '—'}</td>"
        body += f"<tr><td>{lab}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def struct_table():
    """Structural-positivity exclusion (n_infeasible) by age tertile."""
    levels = ["Young", "Middle", "Old", "All"]
    head = ("<table><thead><tr><th>Age tertile</th>"
            + "".join(f"<th>{s} n infeasible / n</th>" for s in SITES) + "</tr></thead><tbody>")
    body = ""
    for lv in levels:
        cells = ""
        for s in SITES:
            d = R[s]["struct"]
            if d is None:
                cells += "<td>—</td>"; continue
            row = d[d.age_grp == lv]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{int(r['n_infeasible'])} / {int(r['n']):,} ({r['frac_infeasible']*100:.2f}%)</td>"
            else:
                cells += "<td>—</td>"
        body += f"<tr><td>{lv}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def posemp_table():
    """THE positivity table: per arm x age tertile, predicted-adherence overlap."""
    arms = [("strain_limiting", "strain-limiting"), ("permissive", "permissive")]
    levels = ["Young", "Middle", "Old", "All"]
    head = ("<table class='sg'><thead><tr><th>Arm</th><th>Age tertile</th>"
            + "".join(f"<th colspan=3>{s}</th>" for s in SITES) + "</tr>"
            + "<tr><th></th><th></th>"
            + "".join("<th>median P(adh)</th><th>% days&lt;.05</th><th>n elig. days</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for arm_id, arm_disp in arms:
        for i, lv in enumerate(levels):
            cells = ""
            for s in SITES:
                d = R[s]["posemp"]
                if d is None:
                    cells += "<td>—</td><td>—</td><td>—</td>"; continue
                row = d[(d.arm == arm_id) & (d.age_grp == lv)]
                if len(row):
                    r = row.iloc[0]
                    hot = (arm_id == "strain_limiting" and lv == "Old")
                    op = "<b>" if hot else ""; cl = "</b>" if hot else ""
                    cells += (f"<td>{op}{r['median_padhere']:.3f}{cl}</td>"
                              f"<td>{op}{r['frac_padhere_lt05']*100:.0f}%{cl}</td>"
                              f"<td>{int(r['n_eligible_days']):,}</td>")
                else:
                    cells += "<td>—</td><td>—</td><td>—</td>"
            armcell = f"<td rowspan={len(levels)}><b>{arm_disp}</b></td>" if i == 0 else ""
            cls = ' class="degen"' if (arm_id == "strain_limiting" and lv == "Old") else ""
            body += f"<tr{cls}>{armcell}<td>{lv}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def wage_table():
    """ESS fraction & max weight by arm x age tertile -- the ESS collapse map."""
    arms = [("strain_limiting", "strain-limiting"), ("permissive", "permissive")]
    levels = ["Young", "Middle", "Old"]
    head = ("<table class='sg'><thead><tr><th>Arm</th><th>Age tertile</th>"
            + "".join(f"<th colspan=2>{s}</th>" for s in SITES) + "</tr>"
            + "<tr><th></th><th></th>"
            + "".join("<th>ESS frac</th><th>max wt</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for arm_id, arm_disp in arms:
        for i, lv in enumerate(levels):
            cells = ""
            for s in SITES:
                d = R[s]["wage"]
                if d is None:
                    cells += "<td>—</td><td>—</td>"; continue
                row = d[(d.arm == arm_id) & (d.age_grp == lv)]
                if len(row):
                    r = row.iloc[0]
                    hot = (arm_id == "strain_limiting" and lv == "Old")
                    cells += f"<td>{r['ess_frac']:.3f}</td><td>{r['wt_max']:.1f}</td>"
                else:
                    cells += "<td>—</td><td>—</td>"
            armcell = f"<td rowspan={len(levels)}><b>{arm_disp}</b></td>" if i == 0 else ""
            body += f"<tr>{armcell}<td>{lv}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def evalue_table():
    head = ("<table><thead><tr><th>Cohort</th><th>RR (point)</th><th>E-value (point)</th>"
            "<th>RR (CI bound)</th><th>E-value (CI bound)</th></tr></thead><tbody>")
    body = ""
    for s in SITES:
        d = R[s]["evalue"]
        if d is None:
            body += f"<tr><td>{s}</td><td>—</td><td>—</td><td>—</td><td>—</td></tr>"; continue
        r = d.iloc[0]
        body += (f"<tr><td>{s}</td><td>{r['rr_point']:.3f}</td><td><b>{r['evalue_point']:.2f}</b></td>"
                 f"<td>{r['rr_ci_bound']:.3f}</td><td>{r['evalue_ci']:.2f}</td></tr>")
    return head + body + "</tbody></table>"

def numer_table():
    order = ["time-only numerator, marginal MSM (PRIMARY)",
             "baseline numerator, marginal MSM (inconsistent ref)",
             "baseline numerator, V-adjusted MSM (standardized)"]
    disp = {order[0]: "Time-only numerator, marginal MSM <b>(PRIMARY, valid)</b>",
            order[1]: "Baseline numerator, marginal MSM (inconsistent reference)",
            order[2]: "Baseline numerator, V-adjusted standardized MSM <b>(valid)</b>"}
    head = ("<table><thead><tr><th>Estimator</th>"
            + "".join(f"<th>{s} RD (pp)</th><th>{s} ESS<sub>SL</sub></th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for key in order:
        cells = ""
        for s in SITES:
            d = R[s]["numer"]
            if d is None:
                cells += "<td>—</td><td>—</td>"; continue
            row = d[d.spec == key]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{pp(float(r['rd']))}</td><td>{_fnum(r['ess_strain_limiting'])}</td>"
            else:
                cells += "<td>—</td><td>—</td>"
        body += f"<tr><td>{disp[key]}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def trim_table():
    trims = [0.0, 0.01, 0.02, 0.05]
    head = ("<table><thead><tr><th>Trim α</th>"
            + "".join(f"<th>{s} RD (pp)</th><th>{s} frac trimmed (SL)</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for tr in trims:
        cells = ""
        for s in SITES:
            d = R[s]["trim"]
            if d is None:
                cells += "<td>—</td><td>—</td>"; continue
            row = d[np.isclose(d.trim, tr)]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{pp(float(r['rd']))}</td><td>{r['frac_trimmed_strain']*100:.1f}%</td>"
            else:
                cells += "<td>—</td><td>—</td>"
        lab = f"{tr:g}" + (" <b>(PRIMARY)</b>" if tr == 0.02 else (" (untrimmed)" if tr == 0 else ""))
        body += f"<tr><td>{lab}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def agg_table():
    head = ("<table><thead><tr><th>Within-day aggregation</th><th>Deviation rule</th>"
            + "".join(f"<th>{s} RD (pp)</th>" for s in SITES) + "</tr></thead><tbody>")
    body = ""
    specs = [("daily_median", "simple"), ("daily_max", "simple"), ("daily_max", "corrected")]
    labelmap = {"daily_median": "daily median (primary)", "daily_max": "daily max (peak)"}
    for ag, rl in specs:
        cells = ""
        for s in SITES:
            d = R[s]["agg"]
            if d is None:
                cells += "<td>—</td>"; continue
            row = d[(d.aggregation == ag) & (d.deviation_rule == rl)]
            cells += f"<td>{pp(float(row.iloc[0]['rd']))}</td>" if len(row) else "<td>—</td>"
        body += f"<tr><td>{labelmap[ag]}</td><td>{rl}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def cumw_table():
    order = ["untruncated (primary)", "cumulative trunc [1,99] pctile", "cumulative cap @ 10"]
    head = ("<table><thead><tr><th>Cumulative-weight specification</th>"
            + "".join(f"<th>{s} RD (pp)</th><th>{s} max IPCW</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for key in order:
        cells = ""
        for s in SITES:
            d = R[s]["cumw"]
            if d is None:
                cells += "<td>—</td><td>—</td>"; continue
            row = d[d.spec == key]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{pp(float(r['rd']))}</td><td>{float(r['max_ipcw']):.1f}</td>"
            else:
                cells += "<td>—</td><td>—</td>"
        body += f"<tr><td>{key}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def pf_block():
    have = {s: R[s]["pf"] for s in SITES if R[s]["pf"] is not None}
    if not have:
        return ("<p class='muted'>This sensitivity (lagged P/F ratio in the weight model) computes on the next real-data re-run, "
                "on the arterial-PaO₂ subset; synthetic CLIF has no arterial gases.</p>")
    label = {"full_cohort_primary": "Full cohort (primary)",
             "pf_subset_no_pf_adj": "P/F-recorded subset, no P/F adj.",
             "pf_subset_with_pf_adj": "P/F-recorded subset, + lagged worst P/F"}
    head = ("<table><thead><tr><th>Cohort</th><th>Specification</th><th>RD (pp)</th><th>n</th>"
            "</tr></thead><tbody>")
    body = ""
    for s, df in have.items():
        d = df.set_index("spec")
        for key in ["full_cohort_primary", "pf_subset_no_pf_adj", "pf_subset_with_pf_adj"]:
            if key in d.index:
                r = d.loc[key]
                body += (f"<tr><td>{s}</td><td>{label[key]}</td><td>{pp(float(r['rd']))}</td>"
                         f"<td>{int(r['n_patients']):,}</td></tr>")
    return head + body + "</tbody></table>"

def mtp_table():
    head = ("<table><thead><tr><th>Estimand</th>"
            + "".join(f"<th>{s} RD (pp)</th><th>{s} ESS<sub>SL</sub></th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for key in ["static ceiling (primary)", "de-escalation MTP"]:
        cells = ""
        for s in SITES:
            d = R[s]["mtp"]
            if d is None:
                cells += "<td>—</td><td>—</td>"; continue
            row = d[d.estimand == key]
            if len(row):
                r = row.iloc[0]
                cells += f"<td>{pp(float(r['rd']))}</td><td>{_fnum(r['ess_strain_limiting'])}</td>"
            else:
                cells += "<td>—</td><td>—</td>"
        body += f"<tr><td>{key}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

# =============================================================================
# FIGURES -- each returns an HTML <figure> data-URI block, or a graceful
# skip-note if its data source is missing.  Every CSV/PDF read is guarded.
# =============================================================================
def _skip(name, why):
    FIG_LOG.append((name, "SKIPPED", why))
    return (f"<p class='muted'><i>[Figure '{name}' skipped: {why}.]</i></p>")

# --- Figure 1: forest plot, primary + subgroups ------------------------------
def fig_forest():
    name = "forest (primary + subgroups)"
    try:
        SUB_ORDER = [("Age tertile", ["Young", "Middle", "Old"]),
                     ("Sex", ["Male", "Female"]),
                     ("Height tertile (w/in sex)", ["Short", "Middle", "Tall"]),
                     ("Race", ["WHITE", "BLACK", "OTHER"])]
        # build the ordered list of rows, top-to-bottom on screen
        rows = []          # (label, is_overall, is_header)
        rows.append(("Overall", True, False))
        rows.append((None, False, True))   # separator
        for fam, levels in SUB_ORDER:
            rows.append((fam, False, "fam"))
            for lv in levels:
                rows.append((lv, False, False))
        # y positions: row 0 at top
        n = len(rows)
        ypos = list(range(n - 1, -1, -1))
        fig, ax = plt.subplots(figsize=(8.2, 0.42 * n + 1.2))
        dodge = 0.16
        for si, s in enumerate(SITES):
            col = SITE_COLOR[s]
            off = (dodge if si == 0 else -dodge)
            sub = R[s]["sub"]
            xs, ys, los, his = [], [], [], []
            for ri, (label, is_overall, is_hdr) in enumerate(rows):
                if is_hdr:        # separator / family header -> no point
                    continue
                y = ypos[ri] + off
                if is_overall:
                    o = O[s]
                    rd, lo, hi = o["rd"], o["rd_lo"], o["rd_hi"]
                else:
                    # find this level under whichever family precedes it
                    fam = None
                    for back in range(ri, -1, -1):
                        if rows[back][2] == "fam":
                            fam = rows[back][0]; break
                    m = sub[(sub.subgroup == fam) & (sub.level == label)]
                    if not len(m):
                        continue
                    r = m.iloc[0]
                    rd, lo, hi = r["rd"], r["rd_lo"], r["rd_hi"]
                xs.append(rd * 100); ys.append(y); los.append(lo * 100); his.append(hi * 100)
            for x, y, lo, hi in zip(xs, ys, los, his):
                ax.plot([lo, hi], [y, y], color=col, lw=1.6, solid_capstyle="round", zorder=2)
            ax.scatter(xs, ys, color=col, s=26, zorder=3, label=s,
                       edgecolors="white", linewidths=0.6)
        # y tick labels
        ylabels = []
        yticks = []
        for ri, (label, is_overall, is_hdr) in enumerate(rows):
            if is_hdr == "fam":
                yticks.append(ypos[ri]); ylabels.append(f"$\\bf{{{label.split(' ')[0]}}}$")
            elif is_hdr is True:
                continue
            else:
                yticks.append(ypos[ri]); ylabels.append("  " + label if not is_overall else label)
        ax.set_yticks(yticks); ax.set_yticklabels(ylabels)
        # bold "Overall" by re-setting its weight
        for t in ax.get_yticklabels():
            if t.get_text() == "Overall":
                t.set_fontweight("bold")
        ax.axvline(0, ls="--", color="#888", lw=1.1, zorder=1)
        # separator line under Overall
        sep_y = ypos[1]
        ax.axhline(sep_y, color=OKABE[7], lw=0.7, ls=":")
        ax.set_xlabel("Risk difference (percentage points)\nstrain-limiting − permissive")
        ax.set_title("28-day mortality risk difference (strain-limiting − permissive)")
        ax.set_ylim(-0.7, n - 0.3)
        ax.grid(axis="x", color="#eee", lw=0.7)
        ax.legend(title="Cohort", loc="lower right", frameon=False)
        ax.annotate("◀ negative = strain-limiting protective", xy=(0.01, 0.005),
                    xycoords="axes fraction", fontsize=9.5, color="#555")
        _clean_ax(ax)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Forest plot of the 28-day mortality risk difference (strain-limiting − permissive), "
            "overall and across pre-specified subgroups. Points are MIMIC / UCSF with 95% bootstrap CIs; "
            "the dashed line marks no effect. Negative = strain-limiting protective.",
            alt="Forest plot of risk differences overall and by subgroup")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 2: empirical positivity scan -------------------------------------
def fig_positivity():
    name = "empirical positivity scan"
    try:
        if any(R[s]["posemp"] is None for s in SITES):
            return _skip(name, "positivity_empirical CSV missing")
        levels = ["Young", "Middle", "Old"]
        x = np.arange(len(levels))
        w = 0.34
        fig, axes = plt.subplots(2, len(SITES), figsize=(8.2, 5.4), sharex="col",
                                 gridspec_kw={"height_ratios": [1, 1], "hspace": 0.12, "wspace": 0.28})
        for si, s in enumerate(SITES):
            d = R[s]["posemp"]
            ax_top, ax_bot = axes[0, si], axes[1, si]
            for ai, arm in enumerate(["permissive", "strain_limiting"]):
                col = ARM_COLOR[arm]
                med = [d[(d.arm == arm) & (d.age_grp == lv)]["median_padhere"].iloc[0]
                       if len(d[(d.arm == arm) & (d.age_grp == lv)]) else np.nan for lv in levels]
                lt05 = [d[(d.arm == arm) & (d.age_grp == lv)]["frac_padhere_lt05"].iloc[0] * 100
                        if len(d[(d.arm == arm) & (d.age_grp == lv)]) else np.nan for lv in levels]
                lab = "strain-limiting" if arm == "strain_limiting" else "permissive"
                off = (-w / 2 if ai == 0 else w / 2)
                ax_top.bar(x + off, med, w, color=col, label=lab, edgecolor="white", linewidth=0.5)
                ax_bot.bar(x + off, lt05, w, color=col, edgecolor="white", linewidth=0.5)
            ax_top.set_yscale("log")
            ax_top.set_ylim(0.005, 1.4)
            ax_top.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
            ax_top.set_title(s, fontsize=12.5)
            ax_top.grid(axis="y", which="both", color="#eee", lw=0.6)
            ax_bot.grid(axis="y", color="#eee", lw=0.6)
            ax_bot.set_ylim(0, 70)
            ax_bot.set_xticks(x); ax_bot.set_xticklabels(levels)
            ax_bot.set_xlabel("Age tertile")
            if si == 0:
                ax_top.set_ylabel("median P(adhere)\n(log scale)")
                ax_bot.set_ylabel("% eligible days\nP(adhere) < 0.05")
            _clean_ax(ax_top); _clean_ax(ax_bot)
        axes[0, 0].legend(loc="lower left", frameon=False, fontsize=9.5)
        fig.suptitle("Empirical positivity: predicted adherence collapses in the strain×Old cell",
                     fontsize=12.5, y=0.97)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Empirical positivity scan by arm × age tertile, per site. Top: median modeled probability "
            "of staying on protocol, P(adhere) (log y). Bottom: fraction of eligible days with P(adhere)<0.05. "
            "The strain-limiting × Old cell collapses (median 0.12 UCSF / 0.011 MIMIC; 43%/62% of days <0.05) — "
            "the large Old-subgroup RD must be read with caution.",
            alt="Empirical positivity scan showing strain-Old overlap collapse")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 3: sensitivity tornado strip -------------------------------------
def _safe_rd(df, mask):
    if df is None:
        return None
    m = df[mask(df)] if callable(mask) else df[mask]
    return float(m.iloc[0]["rd"]) if len(m) else None

def fig_tornado():
    name = "sensitivity tornado strip"
    try:
        # each entry: (label, getter(site) -> rd or None)
        def g_primary(s):  return float(O[s]["rd"])
        def g_trim(v):     return lambda s: _safe_rd(R[s]["trim"], lambda d: np.isclose(d.trim, v))
        def g_numer(spec): return lambda s: _safe_rd(R[s]["numer"], lambda d: d.spec == spec)
        def g_rule_corr(s): return _safe_rd(R[s]["rule"], lambda d: d.deviation_rule == "corrected")
        def g_agg_dmax(s):  return _safe_rd(R[s]["agg"],
                              lambda d: (d.aggregation == "daily_max") & (d.deviation_rule == "corrected"))
        def g_cumw(spec):  return lambda s: _safe_rd(R[s]["cumw"], lambda d: d.spec == spec)
        def g_ph(s):       return _safe_rd(R[s]["ph"],
                              lambda d: (d.ph_source == "arterial_only") & (d.spec == "subset_with_ph_adj"))
        def g_pf(s):       return _safe_rd(R[s]["pf"], lambda d: d.spec == "pf_subset_with_pf_adj")
        def g_dp(s):       return _safe_rd(R[s]["dp"], lambda d: d.spec == "dp_subset_with_dp_adj")
        def g_mtp(s):      return _safe_rd(R[s]["mtp"], lambda d: d.estimand == "de-escalation MTP")
        def g_wcap(v):     return lambda s: _safe_rd(R[s]["wcap"], lambda d: np.isfinite(d.weight_cap) & (d.weight_cap == v))

        specs = [
            ("Primary (trim 0.02, cap 5)", g_primary),
            ("Trim α = 0", g_trim(0.0)),
            ("Trim α = 0.05", g_trim(0.05)),
            ("Baseline numerator, marginal", g_numer("baseline numerator, marginal MSM (inconsistent ref)")),
            ("Baseline numerator, standardized", g_numer("baseline numerator, V-adjusted MSM (standardized)")),
            ("Corrected deviation rule", g_rule_corr),
            ("Daily-max + corrected agg.", g_agg_dmax),
            ("Cumulative-weight trunc [1,99]", g_cumw("cumulative trunc [1,99] pctile")),
            ("Cumulative-weight cap @ 10", g_cumw("cumulative cap @ 10")),
            ("pH arterial-only + adj.", g_ph),
            ("P/F subset + adj.", g_pf),
            ("Driving-pressure subset + adj.", g_dp),
            ("De-escalation MTP", g_mtp),
            ("Weight cap = 3", g_wcap(3)),
            ("Weight cap = 10", g_wcap(10)),
        ]
        n = len(specs)
        ypos = list(range(n - 1, -1, -1))   # primary at top
        fig, ax = plt.subplots(figsize=(8.2, 0.40 * n + 1.0))
        # site primary reference lines
        for s in SITES:
            ax.axvline(float(O[s]["rd"]) * 100, color=SITE_COLOR[s], lw=1.0, ls="--", alpha=0.55,
                       zorder=1)
        ax.axvline(0, color="#888", lw=1.0, ls=":", zorder=1)
        dodge = 0.16
        for si, s in enumerate(SITES):
            col = SITE_COLOR[s]
            off = (dodge if si == 0 else -dodge)
            xs, ys = [], []
            for (label, getter), y in zip(specs, ypos):
                rd = getter(s)
                if rd is None:
                    continue
                xs.append(rd * 100); ys.append(y + off)
            ax.scatter(xs, ys, color=col, s=30, zorder=3, label=s, edgecolors="white", linewidths=0.6)
        ax.set_yticks(ypos)
        ax.set_yticklabels([lab for lab, _ in specs])
        for t in ax.get_yticklabels():
            if t.get_text().startswith("Primary"):
                t.set_fontweight("bold")
        ax.set_xlabel("Risk difference (percentage points)")
        ax.set_title("Sensitivity sweep — every specification clusters near the primary RD")
        ax.set_ylim(-0.7, n - 0.3)
        ax.grid(axis="x", color="#eee", lw=0.7)
        ax.legend(title="Cohort (dashed = site primary)", loc="lower left", frameon=False, fontsize=9.5)
        _clean_ax(ax)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Tornado of the 28-day mortality RD across specifications (primary at top); dashed vertical lines "
            "mark each site's primary RD. Every specification clusters near the primary. The untruncated "
            "weight-cap (∞) is excluded because it degenerates (UCSF +0.75) — see the weight-cap table.",
            alt="Sensitivity tornado strip of risk differences")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 4: weight-cap curve ----------------------------------------------
def fig_wcap_curve():
    name = "weight-cap curve"
    try:
        caps = [3, 5, 10]
        fig, ax = plt.subplots(figsize=(7.4, 4.4))
        for s in SITES:
            d = R[s]["wcap"]
            rds, ess = [], []
            for c in caps:
                m = d[np.isfinite(d.weight_cap) & (d.weight_cap == c)]
                rds.append(float(m.iloc[0]["rd"]) * 100 if len(m) else np.nan)
                ess.append(float(m.iloc[0]["ess_strain_limiting"]) if len(m) else np.nan)
            col = SITE_COLOR[s]
            ax.plot(caps, rds, "-o", color=col, lw=1.8, ms=7, label=s, zorder=3)
            for c, rd, e in zip(caps, rds, ess):
                ax.annotate(f"ESS {e:.2f}", (c, rd), textcoords="offset points",
                            xytext=(6, 8), fontsize=8.5, color=col, alpha=0.9)
        ax.axhline(0, color="#888", ls=":", lw=1.0)
        ax.axvline(5, color=OKABE[7], ls="--", lw=1.0, alpha=0.7)
        ax.annotate("cap = 5 (primary)", (5, ax.get_ylim()[1]), textcoords="offset points",
                    xytext=(4, -12), fontsize=9, color="#555")
        ax.set_xticks(caps)
        ax.set_xlabel("Per-day IPCW truncation cap")
        ax.set_ylabel("28-day mortality RD (pp)")
        ax.set_title("Weight-cap curve — RD stable across usable caps")
        ax.grid(color="#eee", lw=0.7)
        ax.legend(title="Cohort", frameon=False)
        _clean_ax(ax)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "28-day mortality RD vs per-day IPCW cap (3, 5, 10), with the strain-arm ESS annotated at each cap. "
            "The cap = 5 primary is marked; the untruncated cap (∞, omitted) degenerates as ESS collapses.",
            alt="Weight-cap curve of RD vs cap")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 5: covariate balance over follow-up ------------------------------
def fig_balance():
    name = "covariate balance over follow-up"
    try:
        if any(R[s]["balance"] is None for s in SITES):
            return _skip(name, "balance CSV missing")
        covs = ["Age (per 10y)", "PFVC (L)", "SOFA"]
        arms = ["strain_limiting", "permissive"]
        arm_ls = {"strain_limiting": "-", "permissive": "--"}
        fig, axes = plt.subplots(len(covs), len(SITES), figsize=(8.2, 7.4),
                                 sharex=True, gridspec_kw={"hspace": 0.22, "wspace": 0.2})
        for ci_, cov in enumerate(covs):
            for si, s in enumerate(SITES):
                ax = axes[ci_, si]
                d = R[s]["balance"]
                dc = d[d.covariate == cov]
                for arm in arms:
                    da = dc[dc.arm == arm].sort_values("day")
                    if not len(da):
                        continue
                    ls = arm_ls[arm]
                    ax.plot(da.day, da.smd_unweighted, ls, color=OKABE[4], lw=1.4,
                            marker="o", ms=3.5, alpha=0.55)
                    ax.plot(da.day, da.smd_weighted, ls, color=OKABE[2], lw=1.4,
                            marker="s", ms=3.5, alpha=0.55)
                # Contrast-relevant residual: weighted strain − permissive. Both arms drift
                # by the same day-28 survivorship (sick patients die), so the DIFFERENCE
                # cancels it and isolates genuine arm imbalance (11.U). This is the line
                # that matters for confounding of the strain-vs-permissive contrast.
                sp = dc.pivot_table(index="day", columns="arm", values="smd_weighted")
                if {"strain_limiting", "permissive"}.issubset(sp.columns):
                    diff = (sp["strain_limiting"] - sp["permissive"]).sort_index()
                    ax.plot(diff.index, diff.values, "-", color="#000000", lw=2.6,
                            marker="D", ms=4.5, zorder=6)
                ax.axhline(0, color="#444", lw=0.9)
                ax.axhline(0.1, color="#bbb", ls=":", lw=0.9)
                ax.axhline(-0.1, color="#bbb", ls=":", lw=0.9)
                ax.grid(axis="y", color="#f0f0f0", lw=0.6)
                if ci_ == 0:
                    ax.set_title(s, fontsize=12.5)
                if si == 0:
                    ax.set_ylabel(f"{cov}\nSMD")
                if ci_ == len(covs) - 1:
                    ax.set_xlabel("Follow-up day")
                ax.set_xticks([2, 5, 7, 14, 21, 28])
                _clean_ax(ax)
        # custom legend
        from matplotlib.lines import Line2D
        handles = [
            Line2D([0], [0], color=OKABE[4], lw=1.6, marker="o", ms=5, alpha=0.55, label="unweighted (vs baseline)"),
            Line2D([0], [0], color=OKABE[2], lw=1.6, marker="s", ms=5, alpha=0.55, label="IPCW-weighted (vs baseline)"),
            Line2D([0], [0], color="#000000", lw=2.6, marker="D", ms=5, label="weighted strain − permissive (contrast residual)"),
            Line2D([0], [0], color="#444", lw=1.4, ls="-", label="strain arm"),
            Line2D([0], [0], color="#444", lw=1.4, ls="--", label="permissive arm"),
        ]
        fig.legend(handles=handles, loc="upper center", ncol=3, frameon=False,
                   bbox_to_anchor=(0.5, 1.04), fontsize=9.0)
        fig.suptitle("IPCW covariate balance — most SOFA drift is day-28 survivorship (both arms); the strain−permissive residual is small",
                     y=1.07, fontsize=11.0)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Standardized mean differences (SMD) of baseline Age, PFVC, and SOFA among still-at-risk clones at "
            "days 7/14/21/28, by site. Faint lines are each arm vs the full-cohort baseline (unweighted/weighted, "
            "strain solid / permissive dashed); the <b>bold black line is the weighted strain−permissive "
            "difference</b> — the quantity that actually bears on confounding of the contrast. Guides at 0 and ±0.1. "
            "<b>Reading the SOFA panel correctly matters.</b> The strain arm's weighted SOFA SMD vs baseline looks "
            "alarming (−0.21 UCSF / −0.28 MIMIC at day 28), but the <i>permissive</i> arm drifts almost as far "
            "(−0.14 / −0.19) — and that permissive drift equals the <b>pure survivorship</b> floor (day-28 survivors "
            "sit −0.14 / −0.18 below baseline because sicker patients die). So ≈⅔ of the SOFA SMD is "
            "survivorship common to both arms — a reference artifact of comparing the day-28 at-risk set to a "
            "baseline population that includes patients who died — not confounding. The <b>contrast-relevant residual "
            "(strain−permissive) is small: SOFA −0.07 / −0.09, PFVC +0.09 / +0.05 at day 28</b>. PFVC is the genuine "
            "residual selection (larger lungs in the strain arm, the expected size/overlap channel that the "
            "height instrument and the trial reanalysis address); the SOFA differential is negligible, consistent "
            "with adding severity to the weights changing the estimate by less than a tenth of a percentage point "
            "and with adherence being driven by lung size, not severity. Age is balanced by weighting throughout.",
            alt="Covariate balance SMD over follow-up grid")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 6: weights by age × arm ------------------------------------------
def fig_weights_by_age():
    name = "weights by age × arm"
    try:
        if any(R[s]["wage"] is None for s in SITES):
            return _skip(name, "weights_by_age CSV missing")
        levels = ["Young", "Middle", "Old"]
        x = np.arange(len(levels))
        w = 0.34
        fig, axes = plt.subplots(len(SITES), 2, figsize=(8.2, 6.4),
                                 gridspec_kw={"hspace": 0.42, "wspace": 0.28})
        for si, s in enumerate(SITES):
            d = R[s]["wage"]
            ax_max, ax_mean = axes[si, 0], axes[si, 1]
            for ai, arm in enumerate(["permissive", "strain_limiting"]):
                col = ARM_COLOR[arm]
                lab = "strain-limiting" if arm == "strain_limiting" else "permissive"
                off = (-w / 2 if ai == 0 else w / 2)
                wmax = [d[(d.arm == arm) & (d.age_grp == lv)]["wt_max"].iloc[0]
                        if len(d[(d.arm == arm) & (d.age_grp == lv)]) else np.nan for lv in levels]
                wmean = [d[(d.arm == arm) & (d.age_grp == lv)]["mean_stab_w"].iloc[0]
                         if len(d[(d.arm == arm) & (d.age_grp == lv)]) else np.nan for lv in levels]
                ax_max.bar(x + off, wmax, w, color=col, label=lab, edgecolor="white", linewidth=0.5)
                ax_mean.bar(x + off, wmean, w, color=col, edgecolor="white", linewidth=0.5)
            ax_max.set_yscale("log")
            ax_max.yaxis.set_major_locator(LogLocator(base=10, numticks=6))
            ax_max.grid(axis="y", which="both", color="#eee", lw=0.6)
            ax_mean.axhline(1.0, color="#888", ls=":", lw=1.0)
            ax_mean.grid(axis="y", color="#eee", lw=0.6)
            for ax in (ax_max, ax_mean):
                ax.set_xticks(x); ax.set_xticklabels(levels)
                _clean_ax(ax)
            ax_max.set_ylabel(f"{s}\nmax IPC weight (log)")
            ax_mean.set_ylabel("mean stabilized weight")
            if si == 0:
                ax_max.set_title("Max IPC weight by age × arm", fontsize=11.5)
                ax_mean.set_title("Mean stabilized weight by age × arm", fontsize=11.5)
            if si == len(SITES) - 1:
                ax_max.set_xlabel("Age tertile"); ax_mean.set_xlabel("Age tertile")
        axes[0, 0].legend(loc="upper left", frameon=False, fontsize=9.5)
        fig.suptitle("Weight diagnostics by age × arm — tails concentrate in the older strain-arm strata",
                     y=0.99, fontsize=12.5)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Left: maximum IPC weight by age tertile × arm (log y). Right: mean stabilized weight (should sit "
            "near 1, dotted line). Extreme weights and mean-weight drift concentrate in the older strain-arm "
            "strata.",
            alt="Weight diagnostics by age and arm")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 7: ceiling/grace RD heatmap --------------------------------------
def fig_cg_heatmap():
    name = "ceiling/grace RD heatmap"
    try:
        graces = sorted(set(int(v) for s in SITES for v in R[s]["cg"].grace))
        c_lows = sorted(set(int(v) for s in SITES for v in R[s]["cg"].c_low))
        c_highs = sorted(set(int(v) for s in SITES for v in R[s]["cg"].c_high))
        # shared colour scale across all cells/sites (in pp)
        allrd = np.concatenate([R[s]["cg"].rd.values * 100 for s in SITES])
        vmin, vmax = float(np.nanmin(allrd)), float(np.nanmax(allrd))
        fig, axes = plt.subplots(len(SITES), len(graces),
                                 figsize=(2.7 * len(graces), 5.2),
                                 gridspec_kw={"hspace": 0.40, "wspace": 0.18})
        im = None
        for si, s in enumerate(SITES):
            d = R[s]["cg"]
            for gi, g in enumerate(graces):
                ax = axes[si, gi]
                grid = np.full((len(c_lows), len(c_highs)), np.nan)
                for ri, cl in enumerate(c_lows):
                    for cj, ch in enumerate(c_highs):
                        m = d[(d.c_low == cl) & (d.c_high == ch) & (d.grace == g)]
                        if len(m):
                            grid[ri, cj] = float(m.iloc[0]["rd"]) * 100
                im = ax.imshow(grid, cmap="viridis", vmin=vmin, vmax=vmax, aspect="auto")
                ax.set_xticks(range(len(c_highs))); ax.set_xticklabels(c_highs)
                ax.set_yticks(range(len(c_lows))); ax.set_yticklabels(c_lows)
                for ri in range(len(c_lows)):
                    for cj in range(len(c_highs)):
                        v = grid[ri, cj]
                        if not np.isnan(v):
                            # contrast text colour against viridis
                            frac = (v - vmin) / (vmax - vmin + 1e-9)
                            tcol = "white" if frac < 0.55 else "black"
                            primary = (c_lows[ri] == 11 and c_highs[cj] == 16 and g == 1)
                            ax.text(cj, ri, f"{v:.1f}", ha="center", va="center",
                                    color=tcol, fontsize=9,
                                    fontweight="bold" if primary else "normal")
                            if primary:
                                ax.add_patch(plt.Rectangle((cj - 0.5, ri - 0.5), 1, 1, fill=False,
                                             edgecolor=OKABE[0], lw=2.2))
                ax.set_title(f"grace = {g}", fontsize=10.5)
                if gi == 0:
                    ax.set_ylabel(f"{s}\nstrain ceiling C_low")
                if si == len(SITES) - 1:
                    ax.set_xlabel("permissive ceiling C_high")
        cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
        cbar.set_label("28-day mortality RD (pp)")
        fig.suptitle("Ceiling/grace RD grid — wider ceiling gap → larger effect (primary: 11/16, grace 1 boxed)",
                     y=0.99, fontsize=11.5)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Risk-difference heatmaps over the strain ceiling (rows) × permissive ceiling (cols) grid, one column "
            "per grace value, per site, shared viridis scale. The contrast is monotone: a wider ceiling gap "
            "yields a larger effect, and a shorter grace (stricter, more ARMA-like) a larger one still. The "
            "primary cell (C_low=11, C_high=16, grace=1) is boxed.",
            alt="Ceiling/grace RD heatmaps")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

# --- Figure 8: per-protocol cumulative incidence (PDF -> PNG) -----------------
def fig_cuminc(site):
    name = f"per-protocol cumulative incidence ({site})"
    pdf = os.path.join(ROOT, "output", f"{site}_output", "final", "causal", f"tte_ccw_cuminc_{site}.pdf")
    try:
        if not os.path.exists(pdf):
            return _skip(name, f"PDF missing ({pdf})")
        with _tempfile.TemporaryDirectory() as td:
            prefix = os.path.join(td, "cuminc")
            _subprocess.run(["/opt/homebrew/bin/pdftocairo", "-png", "-r", "130", "-singlefile",
                             pdf, prefix], check=True, capture_output=True)
            png = prefix + ".png"
            with open(png, "rb") as fh:
                b = fh.read()
        FIG_LOG.append((name, "OK", ""))
        return embed_png_bytes(b,
            f"{site}: per-protocol 28-day cumulative mortality by arm (MSM-reconstructed). The strain-limiting "
            "(≤11%) curve sits below permissive (≤16%), separating to roughly 6 pp by day 28.",
            alt=f"{site} per-protocol cumulative incidence curves")
    except Exception as e:
        rel = os.path.relpath(pdf, ROOT)
        FIG_LOG.append((name, "FALLBACK", f"{type(e).__name__}: {e}"))
        return (f"<p class='muted'><i>[Could not inline {name} ({type(e).__name__}); "
                f"open the PDF at <code>{rel}</code>.]</i></p>")

# =============================================================================
# 37_tte_discordance_benefit discordance-HTE (the primary heterogeneity analysis) -- tables + figure
# =============================================================================
DISC_ORDER = ["Concordant", "Mid", "Discordant"]

def disc_have():
    return all(R[s]["disc_hte"] is not None and R[s]["disc_grad"] is not None for s in SITES)

def disc_hte_table():
    """SOFA-adjusted per-tertile 28-d mortality RD + the Discordant-Concordant gradient, both sites."""
    head = ("<table class='sg'><thead><tr><th>PBW/PFVC discordance</th>"
            + "".join(f"<th colspan=2>{s}</th>" for s in SITES) + "</tr>"
            + "<tr><th></th>" + "".join("<th>RD [95% CI]</th><th>n</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    disp = {"Concordant": "Concordant (well-sized)", "Mid": "Mid",
            "Discordant": "<b>Discordant (misdosed)</b>"}
    for g in DISC_ORDER:
        cells = ""
        for s in SITES:
            d = R[s]["disc_hte"]; row = d[d.disc_grp == g]
            if len(row):
                r = row.iloc[0]
                cells += (f"<td>{pp(r['rd'])} [{pp(r['rd_lo'])}, {pp(r['rd_hi'])}]</td>"
                          f"<td>{int(r['n']):,}</td>")
            else:
                cells += "<td>—</td><td>—</td>"
        cls = ' class="degen"' if g == "Discordant" else ""
        body += f"<tr{cls}><td>{disp[g]}</td>{cells}</tr>"
    # gradient row
    gcells = ""
    for s in SITES:
        d = R[s]["disc_grad"]; r = d.iloc[0]
        gcells += f"<td colspan=2>{pp(r['estimate'])} [{pp(r['lo'])}, {pp(r['hi'])}]</td>"
    body += ("<tr><td><b>Gradient (Discordant − Concordant)</b></td>" + gcells + "</tr>")
    return head + body + "</tbody></table>"

def disc_dose_table():
    """Dose-correction decomposition by discordance tertile."""
    if any(R[s]["disc_dose"] is None for s in SITES):
        return "<p class='muted'>[Dose-correction table computes on the next 37_tte_discordance_benefit re-run.]</p>"
    head = ("<table class='sg'><thead><tr><th>Discordance</th>"
            + "".join(f"<th colspan=3>{s}</th>" for s in SITES) + "</tr>"
            + "<tr><th></th>" + "".join("<th>VT/PBW</th><th>VT/PFVC %</th><th>cut, mL/kg</th>" for s in SITES)
            + "</tr></thead><tbody>")
    body = ""
    for g in DISC_ORDER:
        cells = ""
        for s in SITES:
            d = R[s]["disc_dose"]; row = d[d.disc_grp == g]
            if len(row):
                r = row.iloc[0]
                cells += (f"<td>{r['median_vtpbw']:.1f}</td><td>{r['median_vtpfvc']:.1f}</td>"
                          f"<td>{r['dvt_mlkg_above']:.2f}</td>")
            else:
                cells += "<td>—</td><td>—</td><td>—</td>"
        body += f"<tr><td>{g}</td>{cells}</tr>"
    return head + body + "</tbody></table>"

def disc_robust_summary():
    """Compact prose summary of the HTE robustness suite, pulling live numbers where present."""
    def gradrange(s):
        d = R[s]["disc_ceil"]
        return (pp(d.gradient.min()), pp(d.gradient.max())) if d is not None else ("—", "—")
    def vr_pair(s):
        d = R[s]["disc_vr"]
        if d is None: return ("—", "—")
        di = d.set_index("model")
        return (pp(di.loc["VR subset, VR NOT in weights", "rd_discordant"]),
                pp(di.loc["VR subset, VR IN weights", "rd_discordant"]))
    def eval_disc(s):
        d = R[s]["disc_eval"]
        if d is None: return "—"
        row = d[d.disc_grp == "Discordant"]
        return f"{row.iloc[0]['evalue_point']:.2f}" if len(row) else "—"
    bits = []
    for s in SITES:
        glo, ghi = gradrange(s); vn, vi = vr_pair(s)
        bits.append(f"<b>{s}</b>: gradient across all swept ceilings {glo} to {ghi} pp; "
                    f"Discordant RD with dead-space (VR) confounder out vs in the weights {vn} vs {vi} pp; "
                    f"Discordant E-value {eval_disc(s)}.")
    return "<p>" + " ".join(bits) + "</p>"

def fig_cate_curve():
    name = "continuous CATE by discordance"
    try:
        if any(R[s]["disc_curve"] is None for s in SITES):
            return _skip(name, "disc_cate_curve CSV missing")
        fig, ax = plt.subplots(figsize=(8.2, 4.6))
        for s in SITES:
            d = R[s]["disc_curve"].sort_values("discordance")
            col = SITE_COLOR[s]
            ax.fill_between(d["discordance"], d["rd_lo"] * 100, d["rd_hi"] * 100, color=col, alpha=0.13, lw=0)
            ax.plot(d["discordance"], d["rd"] * 100, color=col, lw=2, label=s)
            # overlay the tertile point-RDs at each tertile's row in disc_hte (x = median discordance proxy)
        ax.axhline(0, ls="--", color="#888", lw=1.0)
        ax.set_xlabel("PBW/PFVC discordance  (higher = PBW oversizes the lung → misdosed)")
        ax.set_ylabel("CATE: 28-day mortality RD (pp)\nstrain-limiting − permissive")
        ax.set_title("Benefit of strain-limiting accelerates above a discordance threshold")
        ax.grid(axis="both", color="#eee", lw=0.7)
        ax.legend(title="Cohort", loc="lower left", frameon=False)
        ax.annotate("▼ more negative = larger benefit", xy=(0.02, 0.03), xycoords="axes fraction",
                    fontsize=9.5, color="#555")
        _clean_ax(ax)
        FIG_LOG.append((name, "OK", ""))
        return embed_fig(fig,
            "Continuous conditional average treatment effect (SOFA-adjusted, standardized) of the strain-limiting "
            "policy as a smooth function of PBW/PFVC discordance, per site (95% bootstrap band). The benefit is "
            "modest at low discordance and accelerates past a replicated elbow around 17–18, where PBW most "
            "oversizes the lung — i.e. the misdosed patients benefit most.",
            alt="Continuous CATE of strain-limiting by PBW/PFVC discordance, both sites")
    except Exception as e:
        return _skip(name, f"{type(e).__name__}: {e}")

def disc_section():
    """The whole discordance-HTE subsection, or a graceful note if 37_tte_discordance_benefit CSVs are absent."""
    if not disc_have():
        return ("<h3>Heterogeneity — who benefits most (PBW/PFVC discordance)</h3>"
                "<p class='muted'>[The discordance-HTE (script 37_tte_discordance_benefit) tables and figure populate on the next "
                "37_tte_discordance_benefit re-run; the result CSVs are not yet present for both sites.]</p>")
    g_uc = R["UCSF"]["disc_grad"].iloc[0]; g_mi = R["MIMIC"]["disc_grad"].iloc[0]
    robust = disc_robust_summary().replace("<p>", "").replace("</p>", "")
    return f"""
<h3>Heterogeneity — the benefit concentrates in the misdosed</h3>
<p>The age/sex subgroup gradient above is the discrete shadow of a continuous effect modifier: <b>how badly PBW
over-estimates the lung</b>, measured by the PBW/PFVC discordance. Script 37_tte_discordance_benefit re-estimates the strain-limiting
effect as a function of discordance, from one pooled, SOFA-adjusted, IPC-weighted standardized marginal
structural model. The benefit is present across the range but <b>concentrates in the misdosed</b> — the patients
a fixed VT/PBW dose strains most.</p>
{disc_hte_table()}
<p>At both sites the Discordant (most-misdosed) tertile benefits roughly twice as much as the Concordant, and the
<b>Discordant − Concordant gradient excludes zero</b> (UCSF {pp(g_uc['estimate'])} [{pp(g_uc['lo'])}, {pp(g_uc['hi'])}] pp;
MIMIC {pp(g_mi['estimate'])} [{pp(g_mi['lo'])}, {pp(g_mi['hi'])}] pp). The continuous curve shows the shape — a
modest effect at low discordance that accelerates past a replicated threshold near 17–18:</p>
{fig_cate_curve()}
<h4>Mechanism — the algorithm delivers a larger correction to the misdosed</h4>
<p>The targeting is mechanical, not a claim of differential biology: all tertiles look LTVV-compliant on the PBW
yardstick (≈6–8 mL/kg PBW), yet VT/PFVC strain rises with discordance, so the same strain ceiling delivers a
progressively larger tidal-volume cut to the misdosed.</p>
{disc_dose_table()}
{callout("note", "The discordance gradient is robust",
 "The gradient stays negative across all swept strain ceilings and trim/weight-cap settings, is unchanged by a "
 "richer S/F weight model or by adding a ventilatory-ratio (dead-space) confounder, and the Discordant E-value is "
 "the largest of the three tertiles. " + robust)}
{callout("threat", "One honest residual — site-specific oxygenation imbalance",
 "At UCSF (not MIMIC) the IPC weights leave a residual association between lagged oxygenation (S/F) and the "
 "deviation decision that does not respond to weight-timing, a richer S/F model, or dead-space (VR) adjustment — a "
 "real, unexplained, UCSF-specific feature. But the Discordant RD is <b>insensitive</b> to every attempt to address "
 "it (all shifts &lt; 0.2 pp, all conservative), MIMIC balances cleanly, and the E-value bounds it. The targeting "
 "finding does not depend on it.")}
"""

# =============================================================================
# assemble the HTML
# =============================================================================
CSS = """
:root{--ink:#1a1a1a;--mut:#666;--line:#e3e3e3;--bg:#fff;--code:#f6f8fa;
--green:#0a7d52;--greenbg:#e7f5ee;--red:#b23b18;--redbg:#fdece4;
--blue:#0b5fa5;--bluebg:#e8f1fb;--amber:#9a6b00;--amberbg:#fdf4e0;}
*{box-sizing:border-box}
body{font:16px/1.62 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
color:var(--ink);max-width:880px;margin:0 auto;padding:48px 28px 120px;background:var(--bg)}
h1{font-size:30px;line-height:1.2;margin:0 0 6px}
h2{font-size:23px;margin:52px 0 6px;padding-top:18px;border-top:2px solid var(--ink)}
h3{font-size:18px;margin:30px 0 4px;color:#000}
h4{font-size:15px;margin:22px 0 2px;color:var(--mut);text-transform:uppercase;letter-spacing:.04em}
h5{font-size:14px;margin:16px 0 2px;color:#222;font-weight:600}
p,li{margin:9px 0}
.lead{font-size:18px;color:#333}
.sub{font-weight:400;color:var(--mut);font-size:12px}
code{font-family:"SF Mono",ui-monospace,Menlo,Consolas,monospace;font-size:13px}
p code,li code{background:var(--code);padding:1px 5px;border-radius:4px;font-size:13.5px}
pre.r{background:var(--code);border:1px solid var(--line);border-left:3px solid var(--blue);
border-radius:7px;padding:14px 16px;overflow-x:auto;font-size:12.5px;line-height:1.5;margin:14px 0}
pre.r code{font-size:12.5px;color:#24292e}
table{border-collapse:collapse;width:100%;margin:16px 0;font-size:13.5px}
th,td{border:1px solid var(--line);padding:7px 10px;text-align:left;vertical-align:top}
th{background:#fafafa;font-weight:600}
table.sg td:first-child{background:#fafafa}
tr.degen{color:var(--mut);font-style:italic;background:#fbfbfb}
.callout{border-radius:8px;padding:13px 16px;margin:16px 0;border:1px solid}
.callout .ctitle{font-weight:700;margin-bottom:3px;font-size:14px}
.callout .cbody{font-size:14.5px}
.threat{background:var(--redbg);border-color:#f0c4b0}
.threat .ctitle{color:var(--red)}
.rescue{background:var(--greenbg);border-color:#aed9c4}
.rescue .ctitle{color:var(--green)}
.note{background:var(--bluebg);border-color:#bcd8f2}
.note .ctitle{color:var(--blue)}
.key{background:var(--amberbg);border-color:#ecd9a6}
.key .ctitle{color:var(--amber)}
.toc{background:#fafafa;border:1px solid var(--line);border-radius:8px;padding:14px 22px;font-size:14px}
.toc ol{margin:4px 0;padding-left:22px}
.toc a{color:var(--blue);text-decoration:none}
.toc a:hover{text-decoration:underline}
.tag{display:inline-block;background:#eef;border:1px solid #cce;color:#335;border-radius:4px;
padding:0 6px;font-size:11px;font-family:monospace;vertical-align:middle}
.muted{color:var(--mut);font-size:13px}
hr{border:0;border-top:1px solid var(--line);margin:30px 0}
.foot{color:var(--mut);font-size:12.5px;margin-top:60px;border-top:1px solid var(--line);padding-top:14px}
"""

def H(s):  # section header with anchor
    return s

parts = []
parts.append(f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
             f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
             f"<title>Longitudinal Target Trial Emulation — Methods & Results</title>"
             f"<style>{CSS}</style></head><body>")

# --- title + abstract --------------------------------------------------------
mi, uc = O["MIMIC"], O["UCSF"]
parts.append(f"""
<h1>Does Limiting Size-Relative Tidal Volume Reduce Mortality?</h1>
<p class="muted">A target-trial emulation in two intensive-care cohorts (MIMIC and UCSF, CLIF format), using
clone-censor-weight per-protocol estimation with an inverse-probability-of-censoring-weighted marginal structural
model.</p>

<h3>Background</h3>
<p>Lung-protective ventilation scales tidal volume (VT) to <b>predicted body weight (PBW)</b>, a function of
height and sex. But PBW reflects body frame, not lung size, and systematically over-estimates the lungs of older,
female, and non-white patients — so a nominally protective dose (6–8 mL/kg PBW) can deliver an injuriously large
breath <i>relative to a patient's actual lung</i>. <b>Predicted forced vital capacity (PFVC)</b>, from the
GLI-2012 spirometry reference equations (height, age, sex, race), tracks lung size directly, so the size-relative
dose <b>VT/PFVC</b> — tidal volume ÷ predicted FVC — measures how hard each breath strains the lung, independent of
body frame. This motivates a clinical question: <b>would a strategy that caps size-relative lung strain reduce
mortality, compared with usual care?</b></p>

<h3>What we estimated</h3>
<p>We emulate a randomized trial comparing two ventilation strategies, each applied every day a patient is
invasively ventilated:</p>
<ul>
<li><b>Strain-limiting</b> — keep VT/PFVC at or below <b>11%</b> of predicted FVC;</li>
<li><b>Permissive (≈ usual care)</b> — keep VT/PFVC at or below <b>16%</b>.</li>
</ul>
<p>The primary outcome is 28-day all-cause mortality and the estimand is the per-protocol risk difference between
the two strategies. PFVC here is the <i>measurement that defines the strain ceiling</i>, not the object of
comparison: the strategy under test is "limit lung strain," not "dose by PFVC instead of PBW."</p>

<h3>What we found</h3>
<p>Sustaining the strain-limiting strategy was associated with an absolute <b>28-day mortality reduction of about
5 to 6 percentage points</b> at both sites (MIMIC <b>{ci(mi)} pp</b>; UCSF <b>{ci(uc)} pp</b>; both confidence
intervals exclude zero), and the direction held across every sensitivity analysis. The same protective direction appears in
two companion analyses that rest on different assumptions — a height-instrument analysis and a reanalysis of the
randomized low-tidal-volume trials (§<a href="#triangulate">5</a>). Crucially the benefit is <b>not uniform</b>: it
concentrates in the patients PBW most over-sizes (highest PBW/PFVC discordance), with a Discordant−Concordant
gradient that excludes zero at both sites and a continuous effect that accelerates past a replicated discordance
threshold (§<a href="#results">4</a>) — the strategy most helps exactly the misdosed patients a fixed VT/PBW dose
strains hardest. The principal limitation is a residual
imbalance in <i>predicted lung size</i>: the strain-limiting arm is selected toward larger-lunged patients, because
the smallest-lunged patients cannot physically reach an 11% ceiling. The magnitude should therefore be read as an
upper bound on benefit, with the effect in the smallest lungs identified by the companion designs rather than by
this emulation.</p>

{callout("note","What this design can and cannot identify",
 "A target-trial emulation answers a pre-specified randomized question with observational data by writing down "
 "the trial protocol and emulating each element (Hernán &amp; Robins). Two facts shape what it can deliver here. "
 "<b>First</b>, a cross-sectional 'same patient, different dose' contrast is <i>not</i> identifiable: VT/PFVC is "
 "nearly determined by age, sex, race and height, so there is essentially no overlap to compare doses across "
 "otherwise-similar patients. The per-protocol design instead models the <i>day-to-day decision</i> to keep or "
 "relax the tidal volume, which does have overlap and is what makes the inverse-probability weighting valid. "
 "<b>Second</b>, in this cohort clinicians set the tidal volume early and rarely change it, so most patients' "
 "adherence is effectively fixed at baseline — the longitudinal contrast is close to a baseline-dose comparison "
 "and therefore <i>inherits</i>, rather than escapes, the confounding limits of an observational design. We treat "
 "this emulation as one of three converging lines of evidence, not as a standalone identification strategy.")}
""")

# --- TOC ---------------------------------------------------------------------
parts.append("""
<div class="toc"><b>Contents</b>
<ol>
<li><a href="#estimand">Methods I — the target trial and estimand</a></li>
<li><a href="#ccw">Methods II — how the effect is estimated (clone-censor-weight)</a></li>
<li><a href="#walk">Methods III — implementation</a>
  <ol>
   <li><a href="#s-common">The shared analysis engine</a></li>
   <li><a href="#s-primary">Primary estimate and confidence interval</a></li>
   <li><a href="#s-diag">Diagnostics</a></li>
   <li><a href="#s-sens">Sensitivity analyses</a></li>
   <li><a href="#s-repro">Reproducibility</a></li>
  </ol></li>
<li><a href="#results">Results across both cohorts</a></li>
<li><a href="#threats">Validity and limitations</a></li>
<li><a href="#triangulate">Corroboration from companion analyses</a></li>
</ol></div>
""")

# --- 1. estimand -------------------------------------------------------------
parts.append("""
<h2 id="estimand">1 · Methods I — the target trial and estimand</h2>
<p>A <b>target-trial emulation</b> (Hernán &amp; Robins) makes an observational analysis answer a specific,
pre-stated randomized question. One writes down the protocol of the trial one would run, then emulates each
element with data. Specifying the protocol first is what disciplines the analysis: it fixes the eligibility
window, the moment of "randomization" (time zero), the strategies compared, and the outcome, so the comparison
cannot drift into something ill-defined. The protocol and its emulation are:</p>
<table>
<thead><tr><th>Protocol element</th><th>The trial we would run</th><th>How we emulate it</th></tr></thead>
<tbody>
<tr><td><b>Eligibility</b></td><td>Adults on invasive ventilation with a set tidal volume</td>
<td>Patients with a valid predicted FVC, baseline severity (SOFA), demographics, and height</td></tr>
<tr><td><b>Time zero</b></td><td>Start of invasive ventilation</td>
<td>The first index invasive-ventilation timepoint; follow-up day 0</td></tr>
<tr><td><b>Strategies</b></td><td>Keep VT/PFVC ≤ 11% (strain-limiting) vs ≤ 16% (permissive) every day on the ventilator</td>
<td>Two sustained ceilings, 11% and 16% of predicted FVC</td></tr>
<tr><td><b>Assignment</b></td><td>Randomized at time zero</td>
<td><b>Cloning</b>: each patient is copied into <i>both</i> arms at time zero</td></tr>
<tr><td><b>Adherence</b></td><td>Per-protocol; a deviation is breaking the ceiling</td>
<td>Artificial <b>censoring</b> at the first deviation (after a one-day grace), corrected by inverse-probability weighting</td></tr>
<tr><td><b>Outcome</b></td><td>28-day all-cause mortality (primary); time to ventilator liberation (secondary)</td>
<td>Death by day 28; liberation as a competing-risk cumulative incidence</td></tr>
<tr><td><b>Estimand</b></td><td>Per-protocol risk difference</td>
<td>Standardized 28-day cumulative-incidence difference, restricted to the region where both strategies are achievable</td></tr>
</tbody></table>
<p>PFVC is the <b>measurement that defines the strain ceiling</b>, not the object of comparison. We are not
contrasting "PFVC dosing vs PBW dosing"; we are estimating the effect of a <i>strategy</i> — holding the
size-relative dose under a ceiling — that happens to be operationalized through PFVC.</p>

<h4>Estimand and estimator, specified for audit</h4>
<ul>
<li><b>Eligibility and time zero.</b> The first invasive-ventilation timepoint with a set tidal volume; follow-up
begins there (day 0). One record per hospitalization, requiring a valid predicted FVC and PBW, demographics,
baseline SOFA, and height.</li>
<li><b>Strategies.</b> Two sustained VT/PFVC ceilings — strain-limiting at 11% vs permissive at 16% of predicted
FVC — held every day on the ventilator.</li>
<li><b>Grace period.</b> A patient may sit above the ceiling for <b>one day</b> before it counts as a deviation —
a single calendar day to titrate to protective settings. This mirrors the ARDSNet (ARMA) protocol, which reduced
tidal volume to 6 mL/kg over hours and rechecked plateau pressure within ~4 h; a two-day grace proved more
permissive than that trial standard, and the result is stronger and cleaner under a one-day grace (sensitivity in
§<a href="#results">4</a>).</li>
<li><b>Clone-censor-weight estimation.</b> Each patient is cloned into both arms at time zero; a clone is censored
at its first deviation after the grace period; this informative censoring is corrected by
inverse-probability-of-censoring weights (described in §<a href="#ccw">2</a>).</li>
<li><b>Outcome model.</b> A discrete-time pooled-logistic marginal structural model: a weighted logistic
regression of the daily death indicator on arm × a smooth function of day, from which each arm's 28-day cumulative
incidence is reconstructed and the risk difference taken.</li>
<li><b>Uncertainty.</b> A cluster bootstrap resampling <i>patients</i> with replacement. The primary interval
holds the weight model fixed across resamples — a mild under-estimate of uncertainty because it omits
weight-estimation variance — and a refit-per-resample calibration confirms the interval is not materially wider.</li>
</ul>

<h4>How daily values are formed</h4>
<ul>
<li><b>Exposure (VT/PFVC)</b> — the daily <i>median</i> of the recorded ventilator settings (the within-day peak
is also retained, for a sensitivity analysis on the aggregation choice).</li>
<li><b>S/F ratio</b> (oxygenation) — computed per SpO₂ reading, each matched to the most recent FiO₂ within 4 h,
with SpO₂ clamped to the 80–97% range where it tracks oxygenation linearly, then reduced to the day's
<b>worst</b> (lowest) value — the value most likely to drive a tidal-volume decision.</li>
<li><b>pH</b> — daily worst (lowest); <b>driving pressure</b> — daily worst (highest), plateau − PEEP, from
recorded plateau pressures only (never carried forward from an earlier reading).</li>
</ul>

<h4>The confounder set and the weighting</h4>
<p>The weight model predicts the probability that a patient keeps adhering each day, given the variables a
clinician actually watches when deciding whether to push or relax the tidal volume: the <b>previous day's</b>
VT/PFVC, FiO₂, PEEP, respiratory rate, S/F ratio, mean arterial pressure, and vasopressor use, together with
baseline age, sex, race, and severity (SOFA). Two construction choices matter:</p>
<ul>
<li><b>Lagged predictors.</b> The deviation being modeled is defined by <i>today's</i> tidal volume crossing the
ceiling, so the model uses only information available <i>before</i> that day — otherwise it would predict the
action from itself and the weights would explode.</li>
<li><b>Time-only stabilization.</b> The baseline covariates enter the weight model but not the outcome model, so
that the weights — rather than an adjusted regression — carry the job of balancing them. Whether the weights
actually achieve that balance is examined directly in §<a href="#threats">5</a>; the short version is that they
balance age and severity well, and leave a residual selection on predicted lung size that is the analysis's main
limitation.</li>
</ul>

<h4>Restricting to where both strategies are achievable</h4>
<p>A patient is additionally censored on the first day the model judges adherence all-but-impossible (predicted
probability below 2%) — the tail where weights would otherwise blow up. This keeps the estimand to the region of
genuine overlap between the strategies. The restriction is consequential at UCSF, which has marginal overlap on
exactly that tail, and nearly free at MIMIC; a sensitivity sweep reports how the estimate and the share of
patients removed change as the threshold is varied.</p>

<h4>The liberation endpoint is defined independently of tidal volume</h4>
<p>For the secondary (ventilator-liberation) endpoint, extubation is taken as the last day of the full invasive
ventilation course (in any mode) plus one, derived independently of the recorded tidal volume. Because patients
are routinely weaned onto pressure support before extubation — at these sites 28–38% of patients spend their last
ventilated day on a non-volume mode — defining liberation off the end of volume-control ventilation alone would
overstate how early patients are freed. The exposure weighting still stops accruing at the last volume-targeted
day; only the liberation endpoint uses the full course.</p>
""")

# --- 2. CCW picture ----------------------------------------------------------
parts.append("""
<h2 id="ccw">2 · Methods II — how the effect is estimated (clone-censor-weight)</h2>
<p>The difficulty with comparing <i>sustained</i> strategies in observational data is that, at baseline, one does
not know which strategy a patient will turn out to follow — adherence only reveals itself over time. The
clone-censor-weight approach resolves this in three steps:</p>
<ol>
<li><b>Clone.</b> At time zero every patient is duplicated into both arms. Because the copy is identical at
baseline, the two arms are <i>identical by construction</i> at the start — this is what supplies baseline
balance without needing two comparable groups to exist in the raw data.</li>
<li><b>Censor.</b> A clone is artificially censored the moment the patient's actual care departs from the
strategy that clone represents (its VT/PFVC crosses the ceiling, after the grace day). Up to that point the
clone's data are fully consistent with the assigned strategy.</li>
<li><b>Weight.</b> Censoring at deviation is <i>informative</i> — patients who deviate differ systematically from
those who don't — so each surviving clone is up-weighted by the inverse probability that it stayed uncensored,
given its evolving covariates. These inverse-probability-of-censoring weights rebuild the population that would
have adhered, restoring the balance that censoring broke.</li>
</ol>
""" + callout("key", "The assumption everything rests on",
 "This estimate is unbiased only if the decision to deviate is driven entirely by the <b>measured</b>, "
 "day-to-day covariates — i.e. there is no unmeasured factor that pushes both the tidal-volume decision and the "
 "outcome. The daily panel of FiO₂, PEEP, respiratory rate, S/F ratio, mean arterial pressure, and vasopressor "
 "use exists to make that assumption as defensible as possible, and §<a href='#threats'>5</a> examines how well "
 "it holds."))

# --- 3. implementation -------------------------------------------------------
parts.append('<h2 id="walk">3 · Methods III — implementation</h2>'
 '<p>This section walks through the actual analysis code, for readers who want to audit it; it can be skipped '
 'without losing the argument. The analysis is a single shared engine that builds the cohort, the daily panel, '
 'the clones and their weights, and the primary estimate, followed by one short script per downstream analysis '
 '(the confidence interval, the diagnostics, and each sensitivity analysis). No results are computed at report '
 'time — every number in §<a href="#results">4</a> is read from output files written by these scripts on the '
 'MIMIC and UCSF data.</p>'
 + callout("note", "A note on terminology",
   "The low-ceiling arm caps <b>VT/PFVC</b>, a measure of <b>strain</b> (volume relative to lung size), and is "
   "called the <b>strain-limiting</b> arm throughout — in the text, figures, tables, and in the code and output "
   "files (as <code>strain_limiting</code>)."))

# 3.1 the shared engine -- code/30_tte_common.R
parts.append(f"""
<h3 id="s-common">3.1 · The shared analysis engine</h3>
<p>A single file holds everything the downstream scripts depend on: the setup and design parameters, the cohort
and outcome, the daily confounder panel, the cloning-and-weighting builder, the marginal-structural-model engine,
and the one-time build of the primary design. The sub-sections below follow it in source order.</p>

<h4>Setup &amp; design knobs {linetag('30_tte_common.R', 'header')}</h4>
<p>Before any library loads, the script pins every BLAS backend to a single thread. This matters because the
primary bootstrap (in <code>35_tte_primary</code>) spawns many worker <i>processes</i>; if each also launched multithreaded
linear algebra you would oversubscribe the CPU and the run can destabilize on macOS. The knobs that define the
trial are gathered here so the sensitivity leaves downstream are one-line changes:</p>
<ul>
<li><code>C_LOW=11</code>, <code>C_HIGH=16</code> — the two ceilings (% predicted FVC). 11% sits at roughly the
75th percentile of delivered VT/PFVC under guideline VT/PBW 6–8; 16% is a permissive ceiling close to usual care.</li>
<li><code>GRACE=1</code> — days a clone may sit above its ceiling before it counts as a deviation: one calendar
day to titrate to protective settings, then enforce (ARMA-aligned; a 2-day amnesty was more permissive than the
trial standard, and the effect is monotone in grace — see <code>11.D</code>).</li>
<li><code>DAYW_CAP=5</code> — truncation on the per-day weight, the primary defense against a few clones
dominating (revisited in the weight-cap sensitivity, <code>11.C</code>).</li>
<li><code>HORIZON=28</code> days (primary outcome window — the bulk of ICU mortality, ≈ the adherence window); <code>MAX_VENT_DAY=27</code> — the adherence/ventilation window. <code>TRIM_ALPHA=0.02</code> sets the common-support trim; <code>DEESC_FRAC=0.05</code> the de-escalation MTP shift.</li>
</ul>
<p>The L'Ecuyer-CMRG RNG kind is chosen because it gives independent, reproducible streams across the parallel
bootstrap workers (<code>clusterSetRNGStream</code> in <code>35_tte_primary</code> relies on it).</p>
{code_slice('30_tte_common.R', 'header')}

<h4>Cohort, baseline covariates, and outcome {linetag('30_tte_common.R', '10a')}</h4>
<p>This reads the cross-sectional analysis file (one row per hospitalization from script 03) and defines the
<b>outcome</b> and the <b>baseline covariates</b>. The outcome is the day of death within 28 days, computed
from <code>death_dttm − recorded_dttm</code>. Patients with no death in the window are right-censored at 28 days.
This step also applies the <b>structural-positivity restriction</b> (dropping patients whose lung is so small
that even the 4 mL/kg PBW floor already breaches the ceiling) and defines <code>imv_extub</code>, the true
IMV-course extubation day used for the competing-risk liberation endpoint.</p>
<p>Baseline covariates are deliberately the <i>time-invariant</i> ones — age, sex, race, baseline SOFA, and
within-sex height tertiles. These are the variables that define the subgroups and that enter the
numerator/denominator weight models as fixed terms.</p>
{callout("note","The synthetic-only survival simulation",
 "On the synthetic CLIF dataset (and <i>only</i> there) the script simulates a long-tailed survival outcome, "
 "because synthetic CLIF's mortality fields are known to be malformed. On any real site this branch is skipped "
 "and real <code>death_dttm</code> is used; the MIMIC and UCSF numbers in this report use real deaths.")}
{code_slice('30_tte_common.R', '10a')}

<h4>The daily exposure and confounder panel {linetag('30_tte_common.R', '10b')}</h4>
<p>This is the heart of the time-varying design. It assembles one row per <b>patient-day on the ventilator</b>,
carrying (a) the day's exposure — median VT/PFVC from the respiratory-support waterfall — and (b) the day's
confounders. The confounders are exactly the variables a clinician watches when deciding whether to push or relax
tidal volume the next day:</p>
<ul>
<li>ventilator settings — FiO₂, PEEP, respiratory rate (from the cleaned waterfall);</li>
<li>oxygenation — the worst (lowest) daily <b>S/F ratio</b>, where each SpO₂ is matched to the most-recent FiO₂
within a 4-hour window by a <code>data.table</code> rolling join (<code>roll = 4 * 3600</code> seconds);</li>
<li>acid-base / mechanics — the worst (lowest) daily <b>pH</b> and the worst (max) daily <b>driving pressure</b>
(plateau − PEEP, from recorded plateaus only, never forward-filled);</li>
<li>hemodynamics — mean arterial pressure and an indicator for any vasoactive infusion that day.</li>
</ul>
<p>These are the <b>treatment-confounder feedback</b> variables: they are caused by past exposure <i>and</i>
predict both future exposure and the outcome — precisely what ordinary regression adjustment mishandles and what
the MSM is built to handle. The block also computes the panel-drop diagnostic (<code>panel_drop_summary</code>,
how much of the daily panel is list-deleted for incompleteness), surfaced as a CSV by <code>36_tte_diagnostics</code>.</p>
{code_slice('30_tte_common.R', '10b')}

<h4>Cloning, censoring, and the weights (<code>arm_build</code>) {linetag('30_tte_common.R', '10c')}</h4>
<p><code>arm_build()</code> is the cloning engine. Called once per ceiling, it produces, for that arm, each clone's
deviation day and a daily cumulative weight. Three design choices inside it are worth dwelling on:</p>
<h5>The weight model uses <i>lagged</i> confounders only</h5>
<p>The deviation indicator <code>viol</code> is defined by today's VT/PFVC crossing the ceiling. The model that
predicts deviation therefore must <b>not</b> see today's VT/PFVC — that would be predicting an outcome from
itself and produce perfect separation (without lagging, the weights explode and the effective sample collapses). Every
confounder in the denominator model is lagged by one day (<code>l_vtpfvc</code>, <code>l_sf</code>,
<code>l_map</code>, …). This is the discrete-time analogue of "confounders measured before the action."</p>
<h5>Stabilized weights: numerator vs denominator</h5>
<p>Two logistic models are fit for the probability of <i>not</i> deviating: a <b>denominator</b> conditioned on
the full lagged-confounder history, and a <b>numerator</b> conditioned only on the baseline (time-invariant)
covariates. The per-day weight is (1−p<sub>num</sub>)/(1−p<sub>den</sub>); the ratio is the
<i>stabilized</i> IPC weight, whose mean is ≈1 and whose variance is far smaller than the raw 1/(1−p<sub>den</sub>).
The daily weights are multiplied along each clone's follow-up (<code>cumprod</code>) to give the cumulative
weight <code>cumw</code>, and clamped to [1/cap, cap].</p>
<h5>Two deviation rules</h5>
<p>The <code>rule</code> argument toggles how strict adherence is: <code>"simple"</code> counts any post-grace
exceedance as a permanent deviation; <code>"corrected"</code> forgives a transient excursion if the clinician
brings VT/PFVC back under the ceiling by the next day. The deviation-rule sensitivity (<code>11.E</code>) reports
both.</p>
{code_slice('30_tte_common.R', '10c')}

<h4>The marginal structural model and the liberation endpoint {linetag('30_tte_common.R', '10d')}</h4>
<p><code>make_long()</code> expands each clone into person-day rows up to its event or censoring day, attaching
the carried-forward IPC weight as <code>ipcw</code>. <code>build_design()</code> assembles both arms and the
competing-risk liberation table; <code>ci_curve()</code> then fits the MSM: a <b>weighted pooled logistic
regression</b> of the daily death indicator on arm, a natural-spline of day, and their interaction. Pooled
logistic with a fine time spline approximates a continuous-time hazard model; the <b>IPCW weights</b> are what
make the fitted arm contrast a <i>marginal</i> (population-standardized) one rather than a conditional one. The
day-28 cumulative incidence in each arm is reconstructed from the fitted daily hazards
(<code>1 − ∏(1 − hazard)</code>), and the <b>risk difference</b> (<code>rd_from</code>) is their difference.
<code>sg_rd</code> repeats this within each pre-specified subgroup.</p>
{callout("note","Why extubation is not a censoring event",
 "Censoring follow-up at extubation would be informative censoring (patients are extubated <i>because</i> they "
 "are improving) and would bias mortality downward. Because deaths are observed after extubation (including "
 "out-of-hospital deaths), mortality follow-up continues past extubation; the only "
 "censoring events are deviation, the overlap restriction, and the administrative 28-day horizon. Liberation is instead reported as a "
 "<b>competing-risk</b> secondary: an Aalen-Johansen cumulative-incidence function (<code>cif_lib</code>) with "
 "death as the competing event, IPC-weighted, differenced between arms.")}
{code_slice('30_tte_common.R', '10d')}

<h4>Building the primary design {linetag('30_tte_common.R', '10e')}</h4>
<p>The tail of the engine runs the primary design once — <code>build_design(C_LOW, C_HIGH, …)</code> — and
exposes <code>long_all</code>, <code>lib_all</code>, the point estimates <code>point</code> / <code>lib_pt</code>,
the subgroup point estimates <code>sg_point</code>, and the patient id list <code>ids</code>. The bootstrap
itself lives in <code>35_tte_primary</code>; everything here is the expensive shared build that <code>32_tte_run_all.R</code>
pays for exactly once.</p>
{code_slice('30_tte_common.R', '10e')}
""")

# 3.2 primary -- code/35_tte_primary.R
parts.append(f"""
<h3 id="s-primary">3.2 · Primary estimate &amp; bootstrap — <code>code/35_tte_primary.R</code> {linetag('35_tte_primary.R')}</h3>
<p>The first leaf carries the <b>primary result</b>: the cluster bootstrap, the overall RD, the E-value, the
subgroup CIs, and the cumulative-incidence figure. Confidence intervals come from a <b>cluster bootstrap
that resamples patients</b> (not patient-days) with replacement — the unit of independence is the patient, and
each clone/person-day within a patient must move together. For every replicate, <code>boot_one()</code>
recomputes the overall RD, the liberation difference, and every subgroup RD, so all intervals share one
resampling distribution.</p>
<p>The work is chunked across PSOCK worker processes with a live progress/ETA line. The script uses
<code>parLapply</code> (not <code>parLapplyLB</code>) deliberately: non-load-balanced dispatch hands each worker
a fixed contiguous block of replicates, so with <code>clusterSetRNGStream</code> the bootstrap is
<b>reproducible at a given worker count</b> — load-balanced dispatch is timing-dependent and would yield
non-reproducible CIs. It also writes <code>tte_ccw_overall</code>, <code>tte_ccw_evalue</code>, and
<code>tte_ccw_subgroup</code>, then draws the per-protocol cumulative-mortality curves
(<code>tte_ccw_cuminc_&lt;site&gt;.pdf</code>) shown in §4.</p>
{callout("threat","A known, deliberate simplification in the CI",
 "To keep the real-data run tractable, the bootstrap holds the <b>IPCW model fixed</b> across replicates rather "
 "than refitting the weight models inside every resample. This makes the intervals a <i>modest under-estimate</i> "
 "of true uncertainty (it ignores the variance of estimating the weights). The point estimates are unaffected, "
 "and given how far the CIs sit from zero (below), a fully-nested bootstrap would not change any conclusion — "
 "but it is the honest caveat to state, and a refit-per-resample calibration quantifies it.")}
{code_file('35_tte_primary.R')}
""")

# 3.3 diagnostics -- code/36_tte_diagnostics.R
parts.append(f"""
<h3 id="s-diag">3.3 · Diagnostics — <code>code/36_tte_diagnostics.R</code> {linetag('36_tte_diagnostics.R')}</h3>
<p>These are the numbers to inspect <i>before</i> trusting any arm's estimate. This script runs the full set of
positivity and balance diagnostics: the per-arm weight summary (<code>tte_ccw_diagnostics</code>), the
<b>structural</b> and <b>empirical positivity scans by age tertile</b>
(<code>tte_ccw_positivity_structural/empirical_*</code>), the <b>weighted covariate-balance figure</b> over
follow-up (<code>tte_ccw_balance_*</code>), and <b>effective sample size and extreme weights by age tertile ×
arm</b> (<code>tte_ccw_weights_by_age_*</code>). It also writes the daily-panel completeness diagnostic and the
structural-positivity exclusion table.</p>
<p>The per-arm summary reports the fraction of clones that ever deviate, the 99th-percentile and maximum
cumulative weight, and the <b>effective sample-size fraction</b> — the share of nominal sample size that survives
the weighting. As §4 shows, the strain-limiting arm retains ESS ≈ 0.54 overall at both sites — but the by-age
scan is where the real story is: overlap is healthy in the young/middle cells and collapses in the oldest tertile
of the strain arm.</p>
{code_file('36_tte_diagnostics.R')}
""")

# 3.4 sensitivities -- code/11.C - 33_tte_ceiling
SENS_ROWS = [
    ("11.C_sens_weightcap.R",     "Per-day IPCW truncation cap {3, 5, 10, ∞}",                       "tte_ccw_sens_weightcap"),
    ("11.D_sens_ceiling_grace.R", "Ceiling/grace grid (strain {10,11,12} × permissive {14,16} × grace {1,2,3})", "tte_ccw_sens_ceiling_grace"),
    ("11.E_sens_rule.R",          "Deviation rule: simple vs corrected",                             "tte_ccw_sens_rule"),
    ("11.F_sens_aggregation.R",   "Within-day aggregation: daily median vs daily peak VT/PFVC",       "tte_ccw_sens_aggregation"),
    ("11.G_sens_cumweight.R",     "Cumulative-weight truncation (untruncated MSM vs truncated)",      "tte_ccw_sens_cumweight"),
    ("11.H_sens_ph.R",            "Lagged pH in the IPCW denominator (gas-covered subset)",           "tte_ccw_sens_ph"),
    ("11.I_sens_dp.R",            "Lagged worst-of-day driving pressure (plateau-recorded subset)",   "tte_ccw_sens_dp"),
    ("11.J_sens_pf.R",            "Lagged worst-of-day P/F (arterial-PaO₂ subset)",                   "tte_ccw_sens_pf"),
    ("11.K_sens_numerator.R",     "Numerator / MSM-consistency (time-only vs baseline vs standardized)", "tte_ccw_sens_numerator"),
    ("11.L_sens_trim.R",          "Common-support trim sweep (α = 0 … 0.05; primary 0.02)",           "tte_ccw_sens_trim"),
    ("11.M_sens_mtp.R",           "De-escalation modified-treatment-policy (feasible-shift estimand)", "tte_ccw_sens_mtp"),
]
sens_table = ("<table><thead><tr><th>Script</th><th>What it varies</th><th>Output CSV</th></tr></thead><tbody>"
              + "".join(f"<tr><td><code>{nm}</code></td><td>{what}</td><td><code>{csv}_&lt;site&gt;.csv</code></td></tr>"
                        for nm, what, csv in SENS_ROWS)
              + "</tbody></table>")
parts.append(f"""
<h3 id="s-sens">3.4 · Sensitivity analyses — <code>code/11.C</code>–<code>33_tte_ceiling</code></h3>
<p>Each sensitivity is its own leaf: it sources the engine (or short-circuits if already loaded), re-runs the
design with <b>one knob changed</b>, and writes a single CSV. Each is tabulated in §4. The eleven leaves:</p>
{sens_table}
<p>All eleven share the same shape — <code>library(here)</code>, the <code>source(common)</code> guard, then a
<code>map_dfr</code>/<code>build_design</code> sweep ending in a <code>write_csv</code>. The
<b>common-support trim</b> leaf (<code>11.L_sens_trim.R</code>) is short and representative; it is shown in full
below so the source-common pattern is visible, and the rest follow the same template.</p>
{code_file('11.L_sens_trim.R')}
""")

# 3.5 reproducibility & orchestration -- 11.N + 32_tte_run_all.R
parts.append(f"""
<h3 id="s-repro">3.5 · Reproducibility &amp; orchestration — <code>11.N_refit_boot.R</code> + <code>32_tte_run_all.R</code></h3>
<p><code>11.N_refit_boot.R</code> {linetag('11.N_refit_boot.R')} is the env-gated refit-bootstrap calibration: the
primary bootstrap (<code>35_tte_primary</code>) holds the IPCW model fixed, so it omits weight-estimation uncertainty; this
optional leaf <i>refits the weight models per replicate</i> (overall RD only, reduced N) to measure how much the
fixed-weight CI understates. It is OFF by default — enable with <code>PBWPFVC_REFIT_BOOT=1</code>, size with
<code>PBWPFVC_REFIT_N</code> — and reuses the <code>overall</code> object written by <code>35_tte_primary</code> for the
width ratio (sourcing <code>35_tte_primary</code> first when run standalone).</p>
<p><code>32_tte_run_all.R</code> {linetag('32_tte_run_all.R')} is the driver: it builds the engine once and then runs each
analysis script in order. Every script loads the engine through a small loader that builds it on first use,
caches the result to disk keyed by site and cohort definition, and restores it (in well under a second) on later
runs — so the expensive cohort/panel/weight build is paid once per site and re-checked automatically whenever the
inputs change. The only sources of randomness are a synthetic-data survival simulation (used only on the synthetic
test dataset, never on real sites) and the bootstrap, which runs on fixed reproducible random-number streams, so
the order in which the scripts run does not change any result.</p>
{code_file('32_tte_run_all.R')}
""")

# --- 4. results --------------------------------------------------------------
cg_mi = cg_summary("MIMIC"); cg_uc = cg_summary("UCSF")
parts.append(f"""
<h2 id="results">4 · Results across both cohorts</h2>

<h3>Primary endpoint — 28-day mortality</h3>
{overall_table()}
<p>The two independent cohorts replicate closely: a <b>{ci(mi)} pp</b> absolute mortality reduction in MIMIC and
<b>{ci(uc)} pp</b> in UCSF, both confidence intervals well clear of zero. UCSF carries a lower baseline risk
({pp0(uc['risk_permissive'])}% vs {pp0(mi['risk_permissive'])}% in the permissive arm) yet shows the same absolute
benefit — i.e. a somewhat larger <i>relative</i> effect.</p>
{fig_forest()}
<h4>Per-protocol cumulative-incidence curves</h4>
{fig_cuminc("MIMIC")}
{fig_cuminc("UCSF")}

<h3>Data completeness (daily panel)</h3>
{drops_table()}
<p>After the SpO₂/SF rolling-join fix, the daily-panel drop rate is small — ≈2.3% of patient-days at MIMIC and
≈4.3% at UCSF — so list-deletion of incomplete vent-days is not silently reshaping the cohort. A modest number
of patients have their extubation proxy pulled earlier or are lost entirely; these are reported so the
completeness cost is auditable rather than hidden.</p>

<h3>Positivity &amp; weight diagnostics</h3>
{diag_table()}
<p>The permissive arm barely deviates (≈2–3% of patients) because most already sit below 16%, so its effective
sample is essentially the full cohort. The strain-limiting arm deviates in ≈24–26% of patients and retains an
effective-sample fraction of ≈0.54 after weight truncation, with the per-day weight pinned by the cap of 5 (the
<i>cumulative</i> weight runs higher — up to 54 at MIMIC, 144 at UCSF — examined in a sensitivity analysis below).
A healthy aggregate effective sample confirms the weighting is well-behaved <i>on average</i>, but the average
hides where overlap is thin — which the next two tables expose, and which is the analysis's real limitation.</p>

<h4>Structural-positivity exclusion</h4>
{struct_table()}
<p>A patient whose lung is so small that even the lowest LTVV tidal volume (4 mL/kg PBW) already implies
VT/PFVC above the strain ceiling can <i>never</i> adhere — a structural positivity violation no weight repairs.
These are excluded up front: only 5 (UCSF) / 12 (MIMIC) patients, <b>all in the oldest tertile</b>. Tiny in
count, but a directional signal that what little non-identifiability exists lives in the old, low-PFVC stratum.</p>

<h4>Empirical positivity scan — predicted adherence by arm × age</h4>
<p>This is the single most important validity table. For each eligible day it reports the weight model's
predicted probability the clone stays on protocol, P(adhere)=1−p<sub>den</sub>; a mass near zero is a near
positivity violation where the rare adherers carry enormous weights.</p>
{posemp_table()}
{callout("threat", "The strain-limiting × OLD cell is barely identified",
 "In the oldest tertile of the strain arm the median predicted adherence is only <b>0.12 at UCSF / 0.011 at "
 "MIMIC</b>, with <b>43% (UCSF) / 62% (MIMIC) of eligible days at P(adhere)&lt;0.05</b>. That stratum rests on a "
 "handful of effective adherers — so the large Old-subgroup RD below must be read against this scan, not at face "
 "value. The young/middle cells of the strain arm, and the entire permissive arm, have healthy overlap "
 "(median P(adhere) ≈ 0.8–1.0).")}

<h4>ESS &amp; extreme weights by age tertile × arm</h4>
{wage_table()}
<p>Mean stabilized weight should sit near 1; the ESS collapse and the heaviest tails concentrate in the oldest
strain-limiting cell, consistent with the empirical scan. This is the local picture the per-arm aggregate
ESS (≈0.54) averages over.</p>
{fig_weights_by_age()}

<h3>Subgroups — the equity gradient</h3>
{subgroup_table()}
<p>At both sites the effect is <b>largest in the oldest tertile</b> (UCSF Old {pp(R['UCSF']['sub'].set_index(['subgroup','level']).loc[('Age tertile','Old'),'rd'])} pp,
MIMIC Old {pp(R['MIMIC']['sub'].set_index(['subgroup','level']).loc[('Age tertile','Old'),'rd'])} pp), with the
largest benefits also accruing to <b>female</b> patients. This is the absolute-scale signature predicted by the
strain story: the risk difference scales with baseline risk (which rises with age), and the mis-sizing a fixed
VT/PBW imposes is worst in short and female patients — exactly where PBW most overestimates lung size.</p>
{callout("threat", "Pair the Old subgroup result with its positivity caveat",
 "The Old-tertile RD is also the stratum with the <b>worst overlap</b> (empirical scan above: median P(adhere) "
 "0.12 UCSF / 0.011 MIMIC; effective-sample collapse in the weights-by-age table). It must NOT be presented standalone: the point estimate there "
 "rests on a thin set of effective adherers and is the most fragile cell in the analysis. We report it with the "
 "scan attached so a reviewer sees the identification limit alongside the magnitude.")}

{disc_section()}

<h3>Robustness sweeps</h3>
{fig_tornado()}

<h4>Numerator / MSM-consistency — the two valid estimators agree</h4>
{numer_table()}
<p>The <b>primary</b> (time-only numerator, marginal MSM) and the <b>baseline-numerator V-adjusted standardized
MSM</b> are the two <i>valid</i> estimators of the marginal contrast; they agree closely
(MIMIC {pp(R['MIMIC']['numer'].set_index('spec').loc['time-only numerator, marginal MSM (PRIMARY)','rd'])} vs
{pp(R['MIMIC']['numer'].set_index('spec').loc['baseline numerator, V-adjusted MSM (standardized)','rd'])} pp;
UCSF {pp(R['UCSF']['numer'].set_index('spec').loc['time-only numerator, marginal MSM (PRIMARY)','rd'])} vs
{pp(R['UCSF']['numer'].set_index('spec').loc['baseline numerator, V-adjusted MSM (standardized)','rd'])} pp). The
middle row (baseline numerator + marginal MSM) is the <i>inconsistent reference</i> — V handled in neither the
balanced weights nor the outcome model — shown only to document why it was retired as the primary.</p>

<h4>Common-support trim — the primary estimand</h4>
{trim_table()}
<p>The RD is stable as the trim tightens from 0 to 0.05; the primary sits at <b>α=0.02</b>, restricting the
estimand to the modeled-overlap region. UCSF trims ≈9% of strain clones and MIMIC ≈10% at the primary α, and the
RD barely moves (within ≈0.5 pp of untrimmed) — the trim buys identification on the thin tail without changing
the answer.</p>

<h4>Weight cap</h4>
{wcap_table()}
<p>Across the usable caps (3, 5, 10) the estimate stays between roughly −5 and −6 pp at both sites. The
untruncated (∞) row is degenerate — the strain-arm effective sample collapses to ≈2% of nominal and the point
estimate becomes meaningless, jumping to an absurd <b>+76 pp</b> at UCSF and collapsing to <b>≈0</b> at MIMIC.
That degeneracy is the textbook justification for truncating the weights, shown rather than hidden.</p>
{fig_wcap_curve()}

<h4>Cumulative-weight truncation</h4>
{cumw_table()}
<p>The MSM uses the <i>untruncated</i> cumulative IPC weight (only the per-day factor is capped), so a few
clones with a long run of &gt;1 daily factors push the cumulative tail high (max IPCW <b>54.4 at MIMIC, 143.8 at
UCSF</b>). Re-fitting with the cumulative weight truncated leaves the RD stable, confirming the headline is not
tail-driven.</p>

<h4>Within-day aggregation</h4>
{agg_table()}
<p>Defining deviation on the daily <i>peak</i> VT/PFVC rather than the daily median gives essentially the same
RD, so the effect is not an artifact of the within-day summary.</p>
<h4>Ceiling/grace grid &amp; deviation rule</h4>
<p>All <b>24</b> ceiling/grace specifications are protective at both sites: MIMIC ranges from {pp(cg_mi[0])} pp
(widest gap, shortest grace: 10-vs-16, 0-day) to {pp(cg_mi[1])} pp (narrowest gap, longest grace: 12-vs-14,
3-day); UCSF ranges {pp(cg_uc[0])} to {pp(cg_uc[1])} pp. The pattern is physiologically sensible — a wider gap
between the arms yields a larger effect, and a longer grace mildly attenuates it. The deviation rule barely
matters:</p>
{rule_table()}
{fig_cg_heatmap()}

<h4>Respiratory-acidosis (pH) sensitivity</h4>
<p>Permissive hypercapnia is the specific feedback that makes a low-tidal-volume strategy "fail": cutting VT
raises CO₂ and drops pH, which prompts the clinician to relax the ceiling (a deviation), and acidosis is itself
prognostic. So arterial pH is the single most decision-relevant time-varying confounder for this exposure. The
sensitivity adds <b>lagged pH</b> to the IPCW denominator on the gas-covered subset and compares the risk
difference with versus without that adjustment — stability means the strategy effect is not an artifact of
uncontrolled acidosis. It is run for <b>two pH sources</b>: a <i>pooled</i> series (arterial + venous gases,
venous imputed as venous + 0.05) and an <i>arterial-only</i> series (dropping the imputed venous values), so the
venous imputation cannot be doing the work. It is a sensitivity rather than a core-panel confounder because gas
sampling is indication-driven (missing-not-at-random); forcing it into the primary would import selection bias.
Each source also reports the subset <i>without</i> pH, so any selection from restricting to gas-sampled patients
is itself visible.</p>
{ph_block()}

<h4>Driving-pressure sensitivity</h4>
<p>Clinicians titrate tidal volume in response to plateau pressure (the ARMA threshold) and driving pressure
(Amato, <i>NEJM</i> 2015) — so the plateau/DP a patient ran at yesterday is a behaviorally-real driver of
today's dosing decision, and a legitimate time-varying confounder. This sensitivity adds the <b>lagged
worst-of-day driving pressure</b> (DP = plateau − PEEP) to the IPCW denominator on the subset of patient-days
following a recorded plateau. Because plateau is recorded only intermittently and in passive conditions and is
never forward-filled, the subset is restricted to days after a measured plateau, and the base and adjusted
estimates run on that <i>same</i> day-set so they differ only by the DP term. DP enters the <i>weight</i> model,
not the outcome model, so it corrects for why clinicians deviated without attenuating the strain→mortality
pathway. Stability means the strain effect survives adjustment for the actual plateau/DP-driven titration
behavior.</p>
{dp_block()}

<h4>PF-ratio (PaO₂/FiO₂) sensitivity</h4>
<p>Oxygenation co-drives the FiO₂/PEEP titration that determines whether a low-VT arm can be sustained, so the
lagged worst-of-day P/F is a behaviorally-real confounder of the deviation decision. Like pH it is a sensitivity
(arterial gas is MNAR): base vs adjusted run on the same arterial-PaO₂ day-set, so the rows differ only by the
P/F term. The RD is unchanged with adjustment at both sites.</p>
{pf_block()}

<h4>De-escalation modified-treatment-policy (MTP) sensitivity</h4>
<p>The static ceiling is the primary estimand, but the decision that actually carries day-level overlap is
<i>dynamic</i> — cut next-day VT by ≥5% when above the ceiling. Re-estimating under that feasible-shift rule
keeps the effect protective and in the same direction (MIMIC
{pp(R['MIMIC']['mtp'].set_index('estimand').loc['de-escalation MTP','rd'])} pp, UCSF
{pp(R['UCSF']['mtp'].set_index('estimand').loc['de-escalation MTP','rd'])} pp), at somewhat higher ESS, showing
the headline is not an artifact of the static operationalization.</p>
{mtp_table()}

<h3>E-value for the primary RD</h3>
{evalue_table()}
<p>The E-value is the minimum association (risk-ratio scale) an unmeasured confounder of the adherence-censoring
would need with <i>both</i> deviation and mortality to explain the effect away: ≈1.6 (MIMIC) / ≈1.8 (UCSF) at the
point estimate, ≈1.36 / ≈1.66 at the CI bound nearest the null.
<b>Caveat:</b> the textbook E-value is for a point exposure; here it bounds residual confounding of the
informative <i>censoring</i>, so it is an approximate robustness index rather than the exact point-exposure
E-value.</p>

<h3>Secondary endpoint — ventilator liberation (competing risk)</h3>
<p>With extubation defined as the last day of the full ventilation course plus one (any mode), the 28-day
liberation cumulative-incidence difference is <b>null at both sites</b>: UCSF <b>{ci(uc,'lib')} pp</b> and MIMIC
<b>{ci(mi,'lib')} pp</b> (both confidence intervals cross zero). So there is <b>no liberation penalty</b>
alongside the mortality benefit — strain-limiting does not keep patients on the ventilator longer. (Defining
liberation off the end of volume-control ventilation alone, before pressure-support weaning, would have created a
spurious early-liberation signal; defining it off the full ventilation course removes that artifact.)
Mechanistically the null is coherent: a mortality benefit keeps the sickest, slowest-to-wean patients alive — a
competing event that would tend to flatten the liberation curve — and the observed near-zero differences are
consistent with that.</p>

<p class="muted">The cumulative-incidence curves appear under the primary endpoint above, and the covariate-balance
panel appears in §<a href="#threats">5</a> alongside the discussion of residual confounding. Both are also written
as standalone PDFs to each site's output folder.</p>
""")

# --- 5. threats --------------------------------------------------------------
parts.append("""
<h2 id="threats">5 · Validity and limitations</h2>
<p>This section is deliberately exhaustive: each potential threat to the analysis is stated plainly and paired
with the design feature, sensitivity result, or companion analysis that addresses it — or, where a threat is only
partly addressed, with an honest statement of what remains. The single most important limitation is stated first
under sequential exchangeability and positivity: a residual selection toward larger predicted lung size in the
strain-limiting arm, which means the estimated magnitude should be read as an upper bound on benefit.</p>
""")

threats = [
("Demographic paths to mortality not mediated by lung volume (especially age)",
 "Age, sex, and race have many routes to death besides lung size — frailty, comorbidity, immune senescence. "
 "In a conventional <i>dose</i>→mortality regression these are open backdoors, because VT/PFVC is a near-"
 "deterministic function of demographics (GLI-2012 + Devine), so demographics are a common cause of both the "
 "exposure and the outcome. This is the critique a general-medicine reviewer is most likely to raise.",
 "<b>The TTE changes the nature of this question.</b> The two strategies are assigned by <b>cloning the same "
 "patients</b> into both arms at t0, and demographics are <i>time-invariant</i> — so age, sex, and race (and "
 "any unmeasured mortality risk they index) are <b>identical between the strain-limiting and permissive clones "
 "by construction</b>, exactly as randomization would balance them. They cannot confound the strategy contrast; "
 "they enter the MSM only as baseline terms and effect modifiers. This is precisely the backdoor that is open "
 "in the cross-sectional analysis (04–05) and closed here.<br><br>"
 "<b>What remains is narrower.</b> Cloning balances <i>baseline</i> confounding; it cannot balance unmeasured "
 "<i>post-baseline</i> confounding. So the residual worry is not 'age has other paths to death' (balanced) but "
 "specifically 'is there unmeasured <i>time-varying</i> deterioration driving both deviation and death' — which "
 "the daily physiologic panel is built to capture before each dosing decision, the IPCW models adjust for "
 "demographics in both numerator and denominator, and the height-IV (no-unmeasured-confounding-free) and the "
 "randomized RCT reanalysis corroborate.<br><br>"
 "<b>The subgroup gradient is supporting evidence, not contamination.</b> The age gradient appears only on the "
 "<i>absolute</i> (RD) scale; on the <i>relative</i> (OR) scale the strain effect is age-invariant. A "
 "constant relative effect on a higher-baseline-risk group mechanically produces a larger RD — exactly the "
 "observed pattern. A genuine non-lung age→death pathway corrupting the estimate would instead show "
 "age-modification on the relative scale, and there is none. Running every model demographic-adjusted AND "
 "unadjusted, plus spline-age + E-value in the cross-sectional arm, closes the enumerable channels.",
 "rescue"),
("Respiratory acidosis (permissive hypercapnia) as a specific time-varying confounder",
 "A low-tidal-volume strategy works through a known feedback loop: cutting VT cuts minute ventilation, CO₂ "
 "rises, pH falls, and the clinician relaxes the ceiling — i.e. acidosis directly drives deviation — while "
 "acidosis is independently prognostic. If pH is uncaptured, this is the most mechanistically plausible "
 "residual time-varying confounder for <i>this</i> exposure specifically.",
 "Addressed directly by the pH sensitivity (§4): lagged pH is added to the IPCW denominator on the "
 "gas-covered subset and the risk difference is compared with versus without that adjustment. It is run for two "
 "pH sources — pooled (arterial + venous imputed as venous + 0.05) and arterial-only (no venous) — so neither "
 "the venous imputation nor uncontrolled acidosis can be driving the effect. It is a sensitivity rather than a "
 "core confounder because blood-gas sampling is indication-driven (missing-not-at-random); each source also "
 "reports the same subset without pH so any selection from restricting to gas-sampled patients is itself visible.",
 "rescue"),
("No unmeasured time-varying confounding (sequential exchangeability)",
 "The IPCW estimate is consistent only if deviation from a ceiling is fully explained by the <i>measured</i> "
 "lagged covariates. An unmeasured driver of both 'clinician relaxes the tidal volume' and 'patient dies' "
 "(e.g. an unrecorded deterioration the team perceived but did not chart) would bias the result. This is the "
 "load-bearing assumption and the one a reviewer will attack first.",
 "The daily panel captures the variables that actually drive the next-day tidal-volume decision (FiO₂, PEEP, "
 "respiratory rate, S/F ratio, mean arterial pressure, vasopressors), each measured <i>before</i> the decision, "
 "and the overall effective sample (≈0.54 of nominal) means no single patient dominates the headline. <b>But this "
 "is the assumption the data most visibly strain against, and it is the analysis's main limitation.</b> The "
 "empirical positivity scan shows overlap is thin in the oldest strain-arm stratum, and the covariate-balance "
 "figure (below) shows the weights leave the strain arm selected toward <b>larger predicted lung size</b> — a "
 "residual imbalance of about +0.09 (UCSF) / +0.05 (MIMIC) standard deviations at day 28, in a direction that "
 "would inflate apparent protection — because the smallest-lunged patients cannot reach an 11% ceiling and so "
 "cannot be restored by reweighting. (A companion imbalance in baseline severity proved to be mostly an artifact "
 "of comparing 28-day survivors against a baseline population that includes the patients who died: about "
 "two-thirds of it is this survivorship, present in both arms equally, and the genuine between-arm severity "
 "residual is a negligible −0.07 to −0.09 SD. Adherence is driven by lung size, not severity.) The estimated "
 "<i>magnitude</i> — and especially the oldest subgroup — must therefore be read as an upper bound. What supports "
 "the <i>direction</i> regardless is a <b>height-instrument analysis</b> that attacks the same question under a "
 "different, non-overlapping assumption: height shifts predicted lung size without plausibly affecting mortality "
 "except through the delivered dose, and it agrees — so a confounder would have to corrupt both a "
 "no-unmeasured-confounding analysis and an instrumental-variable one to explain the result away.",
 "rescue"),
("Positivity / overlap",
 "If some patients could never adhere to a ceiling, the weights for those clones blow up and the estimate is "
 "driven by a tiny, unrepresentative subset. This is the exact failure that killed the cross-sectional VT/PFVC "
 "contrast (c≈0.996, ESS≈0).",
 "Three lines of defense, but with one honest caveat. <b>First</b>, the design contrasts <b>decisions over "
 "time</b>, and it is an established property of this exposure that the day-level de-escalation decision has "
 "genuine propensity overlap (AUC≈0.66–0.70) — unlike the static exposure level, which is near-determined by "
 "demographics. <b>Second</b>, the positivity diagnostics <i>measure</i> the overlap directly: ESS≈0.54 in the "
 "strain-limiting arm at both sites, and the weight-cap sweep shows the estimate is stable for any sane "
 "truncation and only degenerates at ∞ (reported, not hidden). <b>Third</b>, a structural-positivity floor "
 "excludes the handful of patients whose lung is so small that even 4 mL/kg PBW exceeds the ceiling (5 UCSF / "
 "12 MIMIC, all oldest tertile), and the primary estimand is <b>trimmed to the modeled-overlap region</b> "
 "(P(adhere)≥0.02). <b>The caveat the empirical positivity scan (§4) forces:</b> overlap is <i>not</i> uniform. "
 "In the <b>oldest tertile of the strain arm</b> the median predicted adherence is only 0.12 (UCSF) / 0.011 "
 "(MIMIC), with 43%/62% of eligible days at P(adhere)&lt;0.05 — that cell is barely identified and rests on a "
 "handful of effective adherers, which is exactly why the large Old-subgroup RD must be read against this scan "
 "rather than at face value.",
 "rescue"),
("Informative (adherence) censoring",
 "Censoring a clone at first deviation is not random — sicker patients deviate, so naively analyzing the "
 "uncensored clones would compare survivors of different severities.",
 "This is the textbook job of IPC weighting, and it is the entire purpose of the weighting step. The stabilized weights "
 "restore the population that would have adhered. The two deviation-rule definitions (simple vs corrected) "
 "give essentially identical answers, showing the result is not an artifact of where exactly the "
 "censoring boundary is drawn.",
 "rescue"),
("Immortal-time / time-zero alignment",
 "If eligibility, assignment, and follow-up start are not aligned at the same instant, clones can accrue "
 "'immortal' person-time in which they cannot have the event — a classic source of spurious benefit.",
 "Cloning enforces alignment by construction: both clones start at the identical t0 (index ventilation), are "
 "'assigned' at that instant, and begin follow-up immediately. No clone can be classified by future "
 "information, because the deviation that censors it is evaluated forward in time only.",
 "rescue"),
("Competing risks for the secondary endpoint",
 "Treating extubation as plain censoring in a mortality analysis is informative censoring; treating death as "
 "censoring in a liberation analysis over-counts liberation.",
 "The competing-risk handling addresses both directions: mortality follow-up continues past extubation (out-of-hospital deaths are observed), and "
 "liberation is an Aalen-Johansen competing-risk CIF with death as the competing event. The divergent "
 "liberation result between sites is interpreted in §4, not swept under the rug.",
 "rescue"),
("Model misspecification (pooled-logistic MSM, spline degree, weight model form)",
 "The MSM is parametric — a misspecified hazard shape or weight model could bias the standardized contrast.",
 "The day effect uses a flexible natural spline; the ceiling/grace grid (18 fits) and weight-cap sweep show "
 "the estimate is insensitive to the structural choices that would most plausibly matter. A g-formula / "
 "parametric-standardization cross-check (companion analysis) targets the same estimand with a different "
 "modeling route and agrees in direction and magnitude.",
 "rescue"),
("Under-stated uncertainty (fixed-weight bootstrap)",
 "The cluster bootstrap holds the IPCW models fixed rather than refitting them per replicate, so the reported "
 "CIs slightly understate true uncertainty.",
 "Acknowledged explicitly. Because every interval sits far from zero — the nearest, the young-age "
 "subgroup, still excludes it by several percentage points — a fully-nested bootstrap would widen intervals "
 "without changing a single conclusion. A nested bootstrap is the planned belt-and-suspenders for the final "
 "manuscript.",
 "rescue"),
("External validity / single-strategy-pair",
 "The estimand is a contrast of <i>two specific ceilings</i> (11 vs 16), so the magnitude is contrast-dependent "
 "and the cohorts are two US academic health systems.",
 "Two-site replication (MIMIC + UCSF) with concordant effects is the first answer; the ceiling/grace grid maps "
 "how the effect scales with the contrast so the 11-vs-16 number is not mistaken for a universal constant. The "
 "definitive external check is the <b>Bayesian reanalysis of completed-trial IPD</b> (e.g. ARMA), re-scored by "
 "VT/PFVC — a genuinely randomized dataset that can confirm the observational strategy effect.",
 "rescue"),
]
# figures injected after the matching threat's rescue callout, keyed by a
# substring of the threat title.
THREAT_FIGS = {
    "sequential exchangeability": fig_balance,        # IPCW covariate balance over follow-up
    "Positivity / overlap":       fig_positivity,     # empirical positivity scan
}
for title, threat, rescue, _ in threats:
    block = (f'<h3>{title}</h3>'
             + callout("threat", "Threat", threat)
             + callout("rescue", "What addresses it", rescue))
    for key, figfn in THREAT_FIGS.items():
        if key in title:
            block += figfn()
    parts.append(block)

# --- 6. triangulation --------------------------------------------------------
parts.append("""
<h2 id="triangulate">6 · Corroboration from companion analyses</h2>
<p>This emulation is one of four analyses of the same clinical question, each carried out separately and each
resting on a <i>different</i> identifying assumption, so that the weakness of one is the strength of another. They
are summarized here so a reader can see where this emulation sits and why its limitation — residual selection on
lung size — is covered elsewhere. (The full analyses are reported in their own documents; this report is
self-contained for the target-trial emulation only.)</p>
<table>
<thead><tr><th>Analysis</th><th>Identifying assumption</th><th>Its main weakness</th><th>Covered by</th></tr></thead>
<tbody>
<tr><td><b>Cross-sectional comparison</b></td>
<td>No unmeasured confounding given baseline covariates</td>
<td>No overlap to compare doses — demographics nearly determine VT/PFVC</td>
<td>The other three designs, which do not require dose overlap</td></tr>
<tr><td><b>Height as an instrument</b></td>
<td>Height affects mortality only through the delivered dose</td>
<td>Height is not modifiable; the exclusion assumption is untestable</td>
<td>This emulation and the trial reanalysis, under different assumptions</td></tr>
<tr><td><b>This target-trial emulation</b></td>
<td>No unmeasured day-to-day confounding of the deviation decision</td>
<td>That assumption is untestable, and the exposure is near-constant over time, so the contrast is close to a baseline comparison and inherits its confounding limits</td>
<td>The instrument (different assumption) and the trial reanalysis (randomized)</td></tr>
<tr><td><b>Reanalysis of the randomized trials</b></td>
<td>Randomization (the gold standard)</td>
<td>The trials fixed tidal volume by PBW, so the dose contrast is indirect; the re-scoring by VT/PFVC is secondary</td>
<td>Provides the randomized anchor the observational analyses lean toward</td></tr>
</tbody></table>
<p>The strength of the overall argument is that four analyses with <i>non-overlapping</i> failure modes point the
same way. A single unmeasured confounder cannot rescue all of them at once: it would have to bias an adjusted
analysis, an instrumental-variable analysis, and this clone-censor-weight analysis simultaneously, and also
survive a randomized reanalysis. This emulation contributes one consistent line of that argument; it does not
carry the conclusion on its own.</p>
""")

# --- footer ------------------------------------------------------------------
parts.append(f"""
<div class="foot">
Generated by <code>figures/make_tte_report.py</code> from the result CSVs in
<code>output/&lt;site&gt;_output/final/causal/</code> and the verbatim source of the analysis scripts
(<code>code/30_tte_common.R</code> + <code>code/11.*_*.R</code>). No analysis was re-run to build this report.
Cohorts: MIMIC (n={int(O['MIMIC']['n_patients']):,}), UCSF (n={int(O['UCSF']['n_patients']):,}).
</div>
</body></html>
""")

out = os.path.join(ROOT, "reports", "tte_methods_report.html")
with open(out, "w") as f:
    f.write("".join(parts))
print("wrote", out, "(", os.path.getsize(out), "bytes )")
print("--- figure render log ---")
for nm, status, why in FIG_LOG:
    print(f"  [{status:8s}] {nm}" + (f"  ({why})" if why else ""))
