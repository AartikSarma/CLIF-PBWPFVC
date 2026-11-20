"""Data-zone helpers for local CLIF projects (see SOP.md).

Four rules, enforced here so scripts don't have to remember them:
  1. Core CLIF tables stay where they are: read them from data_directory, never copy them.
  2. Intermediate files go to intermediate_directory, which defaults to a folder OUTSIDE the repo.
  3. Intermediate files carry no dates, geography or direct identifiers unless a column is
     explicitly allowed (and the reason is written down in the README).
  4. Final files go to results_directory; only reviewed aggregates are copied to
     output/final_no_phi/ with release_result().

Usage (run from the project root):
    from utils.clif_paths import load_config, clifpy_kwargs, to_offsets, write_intermediate, write_result
    cfg = load_config()
    co = ClifOrchestrator(**clifpy_kwargs(cfg))
"""

from __future__ import annotations

import datetime as _dt
import json
import os
import platform
import re
import shutil
import subprocess
from pathlib import Path

CONFIG_PATH = Path("config") / "config.json"
RELEASE_DIR = Path("output") / "final_no_phi"

# Column names that suggest PHI. Matched case-insensitively against each column name.
PHI_NAME_PATTERNS = {
    "date or time": r"dttm|(^|_)date($|_)|birth|(^|_)dob($|_)",
    "geography": r"zip|census|county|fips|state_code|latitude|longitude|address|geocode",
    "direct identifier": r"mrn|ssn|first_name|last_name|patient_name|phone|email",
}
AGE_PATTERN = r"(^|_)age($|_)"
AGE_CAP = 90  # HIPAA Safe Harbor: ages over 89 are reported as 90+


class DataZoneError(RuntimeError):
    """Raised when a path or a table breaks the project's data-zone rules."""


# --------------------------------------------------------------------------- config

def _project_root() -> Path:
    return Path.cwd().resolve()


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def _expand(value: str) -> Path:
    return Path(os.path.expandvars(os.path.expanduser(value)))


def load_config(path: str | Path = CONFIG_PATH) -> dict:
    """Read config.json, fill in default folders and check the data-zone rules.

    Empty folder settings default to <repo parent>/clif_workspace/<repo name>/{logs,intermediate,results}.
    Paths may use environment variables, e.g. "$CLIF_SYNTH_DATA_V3_0".
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} not found: copy config/config_template.json to config/config.json")
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)

    root = _project_root()
    workspace = root.parent / "clif_workspace" / root.name
    data_dir = _expand(cfg["data_directory"]).resolve()
    if not data_dir.is_dir():
        raise DataZoneError(f"data_directory does not exist: {data_dir}")
    cfg["data_directory"] = str(data_dir)

    for key, sub in [("output_directory", "logs"),
                     ("intermediate_directory", "intermediate"),
                     ("results_directory", "results")]:
        value = cfg.get(key) or ""
        folder = (_expand(value) if value else workspace / sub).resolve()
        if _is_within(folder, root):
            raise DataZoneError(f"{key} ({folder}) is inside the project folder; point it outside the repo")
        if _is_within(folder, data_dir):
            raise DataZoneError(f"{key} ({folder}) is inside data_directory; core tables are read-only")
        folder.mkdir(parents=True, exist_ok=True)
        cfg[key] = str(folder)

    cfg.setdefault("min_cell_size", 10)
    cfg.setdefault("data_release", "unspecified")
    return cfg


def clifpy_kwargs(cfg: dict) -> dict:
    """Arguments for clifpy.ClifOrchestrator(**clifpy_kwargs(cfg))."""
    return {k: cfg[k] for k in ("data_directory", "filetype", "timezone", "output_directory")}


# --------------------------------------------------------------------------- PHI checks

def phi_columns(df, allow=()) -> dict:
    """Return {column: reason} for columns that look like PHI. Works on pandas or polars frames."""
    allow = set(allow)
    found = {}
    for col in df.columns:
        if col in allow:
            continue
        dtype = str(df[col].dtype).lower()
        if "datetime" in dtype or dtype.startswith("date") or "timestamp" in dtype:
            found[col] = "date or time type"
            continue
        for reason, pattern in PHI_NAME_PATTERNS.items():
            if re.search(pattern, col, flags=re.IGNORECASE):
                found[col] = reason
                break
        else:
            if re.search(AGE_PATTERN, col, flags=re.IGNORECASE):
                try:
                    if float(df[col].max()) > AGE_CAP:
                        found[col] = f"age above {AGE_CAP}; top-code it with top_code_age()"
                except (TypeError, ValueError):
                    pass
    return found


def to_offsets(df, index_col: str, unit: str = "h", drop_index: bool = True):
    """Replace every *_dttm column of a pandas frame with time since index_col.

    New columns are named <col>_offset_<unit> (e.g. recorded_dttm -> recorded_offset_h).
    The original timestamp columns are dropped, so the result can be written as an intermediate.
    """
    import pandas as pd

    per = {"h": pd.Timedelta(hours=1), "d": pd.Timedelta(days=1), "min": pd.Timedelta(minutes=1)}[unit]
    out = df.copy()
    index = pd.to_datetime(out[index_col])
    for col in [c for c in out.columns if c.endswith("_dttm") and c != index_col]:
        out[col.replace("_dttm", f"_offset_{unit}")] = (pd.to_datetime(out[col]) - index) / per
        out = out.drop(columns=col)
    return out.drop(columns=index_col) if drop_index else out


def top_code_age(df, col: str = "age_at_admission", cap: int = AGE_CAP):
    """Top-code age at 90 (HIPAA Safe Harbor) in a pandas frame."""
    out = df.copy()
    out[col] = out[col].clip(upper=cap)
    return out


# --------------------------------------------------------------------------- writers

def write_intermediate(df, name: str, cfg: dict, allow=()) -> Path:
    """Write a parquet file to intermediate_directory after checking for PHI-shaped columns.

    allow: columns you have a documented reason to keep (e.g. a hospitalization_id for joins
    is fine and is not flagged; a date column is not, unless listed here).
    """
    flagged = phi_columns(df, allow)
    if flagged:
        detail = ", ".join(f"{c} ({r})" for c, r in flagged.items())
        raise DataZoneError(f"Not writing {name}: possible PHI in {detail}. "
                            "Reduce it in memory (to_offsets, drop geography) or pass allow=[...] "
                            "and record the reason in the README.")
    path = Path(cfg["intermediate_directory"]) / f"{name}.parquet"
    if hasattr(df, "write_parquet"):  # polars
        df.write_parquet(path)
    else:
        df.to_parquet(path, index=False)
    return path


def write_result(df, name: str, cfg: dict) -> Path:
    """Write a final table (CSV) to results_directory. Review it before release_result()."""
    path = Path(cfg["results_directory"]) / f"{name}.csv"
    if hasattr(df, "write_csv"):  # polars
        df.write_csv(path)
    else:
        df.to_csv(path, index=False)
    return path


def small_cells(csv_path: Path, min_n: int) -> list[str]:
    """List count cells below min_n (columns named n, count, n_* or *_n)."""
    import pandas as pd

    df = pd.read_csv(csv_path)
    bad = []
    for col in df.columns:
        if re.fullmatch(r"n|count|n_.*|.*_n|.*_count", col, flags=re.IGNORECASE):
            vals = pd.to_numeric(df[col], errors="coerce")
            hits = vals[(vals > 0) & (vals < min_n)]
            bad += [f"{col} row {i}: {v:g}" for i, v in hits.items()]
    return bad


def release_result(name: str, cfg: dict) -> Path:
    """Copy a reviewed result from results_directory to output/final_no_phi/<site>_<name>.

    Refuses CSVs with count cells below min_cell_size. A person must still review every file
    (row-level data, IDs, exact dates, free text, small cells in figures) before it leaves the site.
    """
    src = Path(cfg["results_directory"]) / name
    if not src.exists():
        raise FileNotFoundError(src)
    if src.suffix == ".parquet":
        raise DataZoneError("Parquet files are not released; release aggregate CSVs or figures")
    if src.suffix == ".csv":
        bad = small_cells(src, int(cfg["min_cell_size"]))
        if bad:
            raise DataZoneError(f"{name} has cells below n={cfg['min_cell_size']}: " + "; ".join(bad[:10]))
    RELEASE_DIR.mkdir(parents=True, exist_ok=True)
    dst = RELEASE_DIR / f"{cfg['site_name']}_{src.name}"
    shutil.copy2(src, dst)
    return dst


def write_manifest(cfg: dict, counts: dict | None = None) -> Path:
    """Record what produced the results: commit, CLIF version, data release, packages, counts."""
    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                                check=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "status", "--porcelain"], capture_output=True,
                                    text=True).stdout.strip())
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = "unknown", None
    try:
        import clifpy
        clifpy_version = clifpy.__version__
    except ImportError:
        clifpy_version = None
    manifest = {
        "run_date": _dt.date.today().isoformat(),
        "git_commit": commit,
        "uncommitted_changes": dirty,
        "site_name": cfg.get("site_name"),
        "clif_version": cfg.get("clif_version"),
        "data_release": cfg.get("data_release"),
        "clifpy_version": clifpy_version,
        "python_version": platform.python_version(),
        "counts": counts or {},
    }
    path = Path(cfg["results_directory"]) / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
