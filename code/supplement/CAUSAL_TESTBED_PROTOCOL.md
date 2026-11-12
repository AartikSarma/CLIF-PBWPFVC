# Simulation protocol: the planted-truth test of the causal designs

This is the protocol for `causal_testbed.R`. It follows the ADEMP layout for simulation
studies (Morris, White and Crowther, *Statistics in Medicine* 2019): aims,
data-generating mechanisms, estimand, methods, performance measures. The scenarios, the
decision rule and the performance measures were fixed before the first full run. They
were developed on single simulated sites and on runs of 8 to 12 replicates, which were
used to find coding errors and to check run time, not to choose scenarios by their
results. Section 11 lists what was changed after the first full run and why; those
changes were made with its results in hand and are not pre-specified.

## 1. Aims

The hypothesis is that, at a fixed VT/PBW, a high VT/PFVC causes organ injury. At a fixed
VT/PBW, VT/PFVC moves only with λ = log(PBW/PFVC), which age, sex and race explain almost
entirely (R² about 0.99). `causal_triangulation.R` therefore estimates one quantity from
fifteen designs that fail under different assumptions. The simulation asks:

1. Is each estimator correct when its model is true? (stage `matched`)
2. What does each design return under confounders it was and was not built for, when the
   data are not generated from its model? (stages `realistic`, `stress`)
3. Does a single rule fixed in advance, applied to all fifteen estimates, separate a strain
   effect from confounding, and how often is it wrong in each direction? (the decision rule)
4. What effect can a site of this size detect? (stage `power`)
5. When the effect does not accrue steadily, does a segment-by-segment estimate recover
   the time course, and what does it cost against a straight line? (stage `timecourse`)
6. Does pooling across sites keep its nominal error rates when sites differ? (stage `multisite`)

## 2. Data-generating mechanisms

**Kept from the site.** Each patient's cohort (ventilated or no support), age, sex, race,
height, PBW, λ, index VT/PBW, SOFA and SpO2/FiO2, and their hospital time before t0. No
marker value, measurement time or exit time is kept.

**Generated.** For patient *i* at time *t* (days from t0), the log marker is

    y_i(t) = level_i + m_c(t) + [u_i·t + a_i(t) + c_i(t) + s_i(t)] + AR_i(t) + e_i(t)

| Term | Meaning | Matched | Realistic base |
|---|---|---|---|
| `m_c(t)` | the cohort's mean path: its observed mean slope in each phase, joined at t0 and at liberation | yes | yes |
| `u_i` | a slope of the patient's own, plus `u_slope` × an unmeasured severity | severity off | `u_slope` 0.5 slope-SDs |
| `a_i(t)` | curvature: an acute change at t0 that recovers by day 7, with a patient-specific bend | off | on |
| `c_i(t)` | confounding: an index of the demographics times a slope, in every phase (`d_all`, scaled by `arm_scale` in the control), after t0 in both cohorts (`d_acute`), or on the ventilator only (`d_vent`) | per scenario | per scenario |
| `s_i(t)` | strain: β × dose × accrued ventilator time. Dose is log VT over the patient's true lung size, λ + log VT/PBW − (lung-size deviation) | linear, no deviation | linear, deviation SD 0.13 |
| `AR_i(t)` | a stationary autocorrelated process (correlation 0.7 per day) carrying `ar_share` of the residual variance | 0 | 0.5 |
| `e_i(t)` | measurement error, Gaussian or Student t with `error_df` degrees of freedom | Gaussian | t, 5 df |

**Measurement.** Two parts: whether a marker is followed at all in a phase, and then a
schedule of draws at roughly regular intervals. With `info_measure` > 0 the interval
shortens when the marker is worse than expected for that patient (0 matched, 0.3
realistic, 1 under stress).

**Exits.** On a 6-hour grid a ventilated patient may be liberated or die, and a control may
escalate or die, with hazards calibrated to the site. With `info_off` and `info_death` > 0,
liberation is slower, and escalation and death faster, when the marker is worse (0 matched;
0.3 and 0.5 realistic; 1 and 1.5 under stress). The analysis does not model exits.

**Calibration.** Residual SD and mean slope are read from the site by cohort and phase;
the slope SD, level SD, exit hazards, share of patients followed and draws per phase are
read from the site. The generator is then tuned (eight damped passes, at no effect) so that
the *analysis* reads the same residual SD, mean slope, slope SD and sampling off a simulated
site as off the real one. `causal_testbed_realism_<site>.csv` reports observed against
simulated for every summary, including two that are not tuned (the span of a phase and the
lag-1 autocorrelation of the residuals).

**Confounder size.** Unless stated, a confounder's slope is set so that the unadjusted
model sees exactly the λ gradient a strain effect β would give:
`d = β · var(λ) / cov(λ, index)`. Three scenarios depart from that: half, double, and the
age gradient estimated in the site's own controls.

**Scenarios.**

| Stage | Scenarios |
|---|---|
| `matched` (8) | null; strain; age in every phase; age, weaker in controls (30%); age from t0 in both cohorts; age on the ventilator only; sex on the ventilator only; strain plus age in every phase |
| `realistic` (14) | the same eight; a composite of age, sex, race and height on the ventilator; the demographics in GLI's own proportions on the ventilator; an age effect that starts at 65; age at half and at twice the size; age as observed in the controls |
| `stress` (22) | one change from the realistic base. With and without strain: strongly informative exits; strongly informative measurement; t errors with 3 df; autocorrelation 0.9 with 80% of the variance; VT set by predicted size and unmeasured severity; the same with the VT/PBW 6–8 gate. With strain only: a threshold in dose; an effect that saturates (3-day constant); an effect delayed by a day; an effect 0.3 times as large in patients with hospital time before intubation; positive pressure before the recorded intubation in 40%; 3 cm error in charted height; lung-size deviation SD 0.25; strain following GLI at age 25; strain following GLI without race. And strongly informative exits with age in every phase |
| `power` (5) | strain at 0.25, 0.5, 0.75, 1 and 1.5 times the planted effect |
| `timecourse` (4) | realistic base: null; strain accruing steadily; saturating (3-day constant); delayed by a day |
| `illness_death` (6) | realistic base, in which death already follows the marker: null; strain at 1 and 3 times the planted effect; strain at 3 times with a direct path from strain to death (log hazard 3 per log unit of dose); no strain, age raising death in both cohorts (log hazard 0.5 per SD); no strain, age steepening the marker in both cohorts (3 times the usual size) |

The planted effect is 0.15 log units per day per log unit of VT/PFVC in the direction of
injury (about 2% a day per SD of λ).

## 3. Estimand

The ventilation-specific change in the marker's slope per log unit of λ, at a fixed VT/PBW.
In every scenario the truth is computed from the planted quantities with no noise: the
least-squares gradient, on the λ the analysis uses and with log VT/PBW held, of each
ventilated patient's planted strain effect per day over the 7-day window. With a linear
effect on GLI's λ this is β. With a threshold, a saturating or delayed time course, or
another denominator, it is the linear summary a design can at best return, and it is
reported for each scenario (`truth`).

Two further truths are computed the same way. For the designs fitted without the VT/PBW
term, the gradient with log VT/PBW not held: the effect as practice delivers it, which
differs from the first when tidal volume follows λ. For the time course, the gradient of
the total effect accrued by days 2, 4 and 7; a segment's truth is the difference between
two of them over the segment's length.

## 4. Methods

The fifteen designs of `causal_common.R`, each a least-squares regression with an intercept
for every patient-phase and standard errors clustered on the patient: `unadjusted`,
`adjusted`, `fingerprint`, `control_did`, `within_onoff`, `within_triple`, `off_reversal`,
and `free_age`, `free_sex`, `free_race`, `free_height`, each also as a
ventilated-minus-control contrast (`_did`). A design needs 50 patients per arm with the
phases it reads; otherwise it is not estimable in that replicate and is counted
(`n_not_estimable`).

Every design is fitted twice: with the index log VT/PBW as a slope term on the ventilator
("held"), and without it ("dropped"). Each is scored against its own truth (section 3).

The time course (`fit_time_course()`) assumes no shape. The ventilator days are cut at
days 2 and 4; a design is fitted with the on-ventilator phase restricted to each segment
in turn, each segment with its own intercept; the cumulative effect by a cut is the sum
of segment gradient × segment length, with a standard error from the patients'
influences added across segments. It is run for `unadjusted`, `control_did`,
`free_age_did`, `within_onoff` and `off_reversal`, beside the straight line over the
whole window.

`causal_illness_death.R`'s analysis is run on the `illness_death` scenarios: the simulated
marker at t0 is the baseline, injury is the first measured crossing of the threshold on
the ventilator, death is followed over the 7 simulated days, and the λ gradient of each
transition's hazard is fitted by `unadjusted`, `control_did` and `free_age_did` (VT/PBW
held). A hazard gradient has no closed-form truth under this generator, so that stage
reports the mean estimate, the ratio of model to empirical standard error, and the
rejection rate; the scenarios without strain are where the rejection rate has a
target (5%).

## 5. Performance measures

Per marker, scenario and design, over the replicates in which the design was estimable:

| Measure | Definition | Monte Carlo SE |
|---|---|---|
| Bias | mean(estimate − truth) | SD(estimate − truth) / √n |
| Empirical SE | SD(estimate) | SE / √(2(n − 1)) |
| Model SE | √mean(SE²), and its ratio to the empirical SE | |
| Coverage | share of 95% intervals containing the truth | √(p(1 − p)/n) |
| Rejection rate | share with \|estimate / SE\| > 1.96 | √(p(1 − p)/n) |
| MDE | 2.80 × empirical SE (80% power, two-sided 5%) | |

With 1,000 replicates the Monte Carlo SE of a coverage of 95% is 0.7 percentage points and
of any rate at most 1.6. A design is labelled, per scenario: **biased** (bias above a
quarter of the planted effect and above three Monte Carlo SEs), **undercovers** (not biased,
coverage below 90% beyond three Monte Carlo SEs), **underpowered** (neither, MDE above the
planted effect), or **ok**.

## 6. The decision rule

One site's estimates support the hypothesis only if all of these hold, in the direction
of injury, at two-sided 5%, **with VT/PBW held and again with it dropped** (section 11):

1. **Between cohorts.** `control_did` is significant.
2. **Within patients.** `off_reversal` or `within_onoff` is significant, and neither is
   significant in the opposite direction.
3. **Not age, not sex.** `free_age_did` and `free_sex_did` are each significant.
4. **The rungs agree.** All four `free_*_did` estimates have the sign of injury and lie
   within a factor of two of `control_did`.

A design that is not estimable counts as not significant. The rule is scored as a whole:
the share of simulated sites it calls "supports", which is a false-positive rate where no
strain effect was planted and power where one was. It is compared with the rule applied
one way only (VT/PBW held; VT/PBW dropped) and with two single analyses a reader might
otherwise rely on: the unadjusted model alone and the control contrast alone.

The rule cannot separate strain from a path that acts only during ventilation through the
demographics in GLI's own proportions. That scenario (`gli_ventilation`) is included so that
its false-positive rate is reported, not hidden.

## 7. Multisite

Four sites per replicate. Each is the site's patients resampled with replacement, with its
own strain effect (between-site SD 30% of β), confounder size (log-normal, SD 0.3) and
residual SD (uniform 0.8 to 1.25). Each design is pooled by common-effect inverse variance
(the repository's convention for two or three sites) and by DerSimonian–Laird random
effects, under four truths (null, strain, strain plus age, age weaker in controls), and
scored against the mean effect β with the measures of section 5.

## 8. Verification

- `causal_tests.R` (no patient data): the estimator and its clustered standard error
  against least squares with a dummy for every patient-phase; exact recovery of the planted
  effect by all fifteen designs when there is no noise; exact removal of an age confounder
  by the designs that free age; the generator's own rules; the per-phase lines used for
  calibration against `lm`; the decision rule against its definition; the Poisson model
  for an outcome with zeros against `glm` with a dummy for every patient-phase; the time
  course (exact recovery segment by segment, the cumulative standard error, a saturating
  effect); the VT/PBW gradient.
- Stage `crosscheck`: `unadjusted`, `control_did`, `free_age_did` and `within_onoff` are
  refitted with a linear mixed model (`nlme::lme`, random level and slope per patient) on
  the same simulated sites, and the two estimates are compared.

## 9. Reproducibility

Replicate *r* of scenario *c* has seed `20261111 + 100000·c + r`, whatever the number of
parallel processes. Parallel work uses separate R processes (PSOCK), not forks, so the same
code runs on Windows, macOS and Linux. Packages are pinned in `uvr.lock`.

    uvr run code/supplement/causal_tests.R
    uvr run code/supplement/causal_testbed.R

Settings are environment variables listed in the script's header; in PowerShell set them
with `$env:PBWPFVC_CAUSAL_REPS = "1000"` before the command.

## 10. Limits

- The generator is calibrated to whichever site runs it. Until it is run at a site with
  real data, the measurement grid, noise and exit rates are those of synthetic CLIF data.
- Hazards are log-linear in how much worse the marker is than expected; other forms of
  informative exit are not covered.
- Whether a marker is followed is drawn independently in each phase.
- Confounders act on the slope, not on the level or on the chance of being measured.
- The two markers are simulated separately, so no joint verdict across markers is scored.
- The analysis under test is the longitudinal model alone. The joint model of
  `22_biotrauma_fit.R` is not part of this study.
- Only log-scale labs are simulated. The vasopressor outcomes (zero-inflated, charted
  continuously) and OSI (on the ventilator only) are fitted by `causal_triangulation.R`
  and unit-tested, but their error rates under this generator are unknown.
- In the `illness_death` stage death is simulated only on the ventilator and within 7
  days, so few controls die: the control contrast on the death transitions is often not
  estimable there. Injury is seen only when the marker is measured.
- The denominator comparison of `causal_triangulation.R` is not simulated beyond the two
  wrong-denominator stress scenarios.
- The separate-terms check of `causal_triangulation.R` (log PBW freed beside lambda with
  VT/PBW dropped, and rho) has unit tests but is not simulated, and is not in the rule.
- So do the age-mediation split (lambda's age part against the rest) and the aeration
  designs (index compliance over predicted size, with BMI). Both were built after the
  results from two sites had been seen, are not simulated, and are not in the rule.
- `causal_longitudinal.R` (within-patient changes in tidal volume) has unit tests but is
  not simulated: the generator does not produce a ventilator course. Its specification
  was fixed in the script header before any real-data run.

## 11. Changes after the first full run

The first full run (1,000 replicates, both markers) was made with the rule of section 6
applied with VT/PBW held only. Its results led to these changes:

1. **The rule must hold with VT/PBW held and dropped.** Where tidal volume was set by
   predicted size and unmeasured severity and no effect was planted, the rule with VT/PBW
   held called 20% (platelets) and 55% (creatinine) of sites "supports"; a re-run of the
   platelet stress stage with the term dropped gave 0.1%. VT/PBW is then a consequence of
   both λ and severity, and holding it lets severity into the λ gradient. Dropping it is
   not right everywhere either: it changes the estimand where tidal volume follows λ.
   Both versions are therefore fitted and reported, and the rule requires both.
2. **A truth for the designs without the VT/PBW term** (section 3), so that they are
   scored against what they estimate.
3. **The `illness_death` stage**, and the generator's two further routes to death
   (`death_direct`, `death_age`), added with `causal_illness_death.R`.
4. **The `timecourse` stage and estimator**, because a saturating effect cost the rule
   most of its power (26% to 35%) and markers are expected to differ in speed.

Because these were chosen with results in hand, the second run's error rates for the
amended rule are a description of it on these scenarios, not a pre-specified test.

### What the second run found (1,000 replicates, both markers)

These are results, recorded so the analyses are read accordingly; they were not
pre-specified.

- **The amended rule.** With no effect planted it said "supports" in at most 0.8% of sites
  in every scenario but one, including tidal volume set by indication (0.0%, against 20%
  and 55% with VT/PBW held only) and the VT/PBW gate as a collider (0.5% and 1.4%). The
  exception is a ventilation effect that follows the full PFVC equation without being
  strain (98.5% and 98.7%), which no design separates. Its power was 96% (platelets) and
  98% (creatinine) at the planted effect, 44% and 43% at half, and 8% and 4.5% at a quarter.
- **Low power** remains where the effect saturates over time (34%, 25%), under strong
  informative dropout (57%, 78%), when the true denominator has no race term (43%) and,
  most of all, when strain follows the age-free structural volume (PFVC at age 25: 1.3%,
  3.0%). A null under the rule does not exclude strain on that denominator.
- **Time course.** With no effect the segment and cumulative estimates of the control and
  age-freed contrasts rejected in at most 7.8% of sites. The first segment (days 0 to 2)
  was within about 5% of the truth under every time shape. For platelets the later
  segments ran about 15% too large in absolute value under a linear effect, with coverage
  kept at 93% to 95%; creatinine's did not. The before-and-after design cannot support the
  later segments.
- **Illness-death analysis.** On the hazard scale the plain control contrast
  (`control_did`) is biased: with age steepening the marker equally in both cohorts and
  no strain, its injury gradient was about -2 and rejected in 98% (platelets) and 76%
  (creatinine) of sites. The age-freed contrast (`free_age_did`) was not (4.6%, 5.5%).
  The unadjusted contrast of injury against death without injury is positive under that
  scenario too (92%, 100%), and negative under frailty. Read the illness-death designs
  from `free_age_did`. Death without injury has a positive gradient under strain alone,
  because worsening below the threshold also raises the hazard of death. The control
  arm has too few simulated deaths for the death transitions of the contrasts to be
  checked.
- **Multisite.** The common-effect pool under-covered when the effect varies between sites
  (coverage 44% to 72%); DerSimonian-Laird covered 80% to 88%.
