# Local CLIF project SOP

How a CLIF study moves from synthetic-data development in the cloud to a run on local
data, and where every file it reads or writes is allowed to live.

## Data zones

| Zone | What lives there | Where | In git? |
|------|------------------|-------|---------|
| 1. Core CLIF tables | The site's CLIF parquet files | `data_directory` (the shared data folder). Read in place, never copied. | No |
| 2. Intermediate | Row-level working files made by scripts | `intermediate_directory`, outside the repo. Default: `../clif_workspace/<repo>/intermediate` | No |
| 3. Results | Tables, figures and model outputs from a run | `results_directory`, outside the repo. Default: `../clif_workspace/<repo>/results` | No |
| 4. Release | Reviewed aggregates to share with the consortium | `output/final_no_phi/` in the repo | Yes, after review |

Logs and clifpy output go to `output_directory` (default `../clif_workspace/<repo>/logs`).

## The four standards

1. **Core tables stay on the server.** Scripts read them from `data_directory`. No copies,
   subsets or exports of whole core tables anywhere else.
2. **Intermediates default to outside the repo.** `load_config()` creates the workspace
   folders next to the repo and refuses any output path inside the repo or inside
   `data_directory`.
3. **Intermediates carry no PHI unless there is a written reason.** Do date arithmetic in
   memory and write offsets (hours or days from admission or another index time), not
   dates. Drop geography (zip codes, census tract or block, county, state and FIPS codes)
   and direct identifiers. Top-code ages at 90. `write_intermediate()` refuses a table
   with date, time, geography or identifier columns; if one is truly needed, pass it in
   `allow` and record the reason in the README.
4. **Final files go to the results folder.** `write_result()` writes there. Only
   `release_result()` copies a file into `output/final_no_phi/`, and only after it passes
   the small-cell check (`min_cell_size`, default 10). A person reviews every released
   file before it is committed.

## Set up a new study

1. Create the study repo from the consortium
   [CLIF-Project-Template](https://github.com/Common-Longitudinal-ICU-data-Format/CLIF-Project-Template).
2. Copy this `project-template/` folder over its root. `config/config_template.json` and
   `.gitignore` here are supersets of the template's.
3. Install the commit checks: `uv tool install pre-commit`, then `pre-commit install`.
4. In the cloud environment, copy `config/config_template.json` to `config/config.json`,
   set `data_directory` to `$CLIF_SYNTH_DATA_V3_0` (or `_V2_1`) and leave the three output
   folders empty to use the defaults.
5. Load paths only through the helpers:
   - Python: `from utils.clif_paths import load_config, clifpy_kwargs, write_intermediate, write_result`
   - R: `source(file.path("utils", "clif_paths.R"))`
6. Develop and test against synthetic data in both CLIF 2.1 and 3.0.

## Move to local data

1. Clone the repo in the local secure environment. Pull only from the project GitHub.
2. Create `config/config.json` from the template: set `site_name`, `clif_version`,
   `timezone` (the site's own zone), `data_directory` (the local CLIF folder) and
   `data_release`. Leave the output folders empty or point them at a project workspace
   that is outside the repo and outside the data folder.
3. Add the data and workspace folders to `~/.clif/phi-paths` so the clif-icu PHI guard
   blocks agent reads there.
4. Run with no agent present. Each run writes `run_manifest.json` (git commit, data
   release, package versions, counts) to the results folder.
5. Review results, then release aggregates with `release_result()`. Commit only
   `output/final_no_phi/`.

## Checklist

### Before the first commit
- [ ] Repo made from the consortium template, this folder copied over it
- [ ] `pre-commit install` run
- [ ] No local paths anywhere in code; every path comes from `config/config.json`

### Before moving to local data
- [ ] Scripts run end to end on synthetic CLIF 2.1 and 3.0
- [ ] Every intermediate goes through `write_intermediate()`, and every `allow` has a reason in the README
- [ ] Dates converted to offsets, geography dropped, ages top-coded before anything is written
- [ ] Every output goes through `write_result()`; nothing writes into the repo

### Local run
- [ ] `config/config.json` filled in (site, timezone, data path, data release), not committed
- [ ] Data and workspace folders listed in `~/.clif/phi-paths`
- [ ] Run done with no agent present; `run_manifest.json` written

### Before release
- [ ] Each file in `output/final_no_phi/` reviewed by a person
- [ ] No cells under `min_cell_size`, including counts that can be derived by subtraction
- [ ] No dates, geography or row-level data in released files
- [ ] `git status` shows only code and reviewed aggregates

## Imaging, waveform and note data

Keep these as pointers, not copies. CLIF tables carry the metadata (for example the CLIF 3.0
radiology and clinical notes tables), and a site resolver maps each pointer to the file on
whichever server holds it. In the cloud, the resolver points at small synthetic or public
samples; locally it points at the site's servers. Analysis code calls the resolver and never
hard-codes a server path. Derived features (for example a waveform summary per hour) are
intermediates and follow the same four standards.
