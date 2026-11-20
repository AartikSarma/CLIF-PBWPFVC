## Configuration

1. Copy `config_template.json` to `config.json` (never committed; `.gitignore` and the
   commit checks keep it local).
2. Fill in `site_name`, `clif_version` (`2.1` or `3.0`), `timezone` (the site's own),
   `data_directory` (the folder holding the `clif_*` tables, read in place) and
   `data_release`.
3. Leave `output_directory`, `intermediate_directory` and `results_directory` empty to use
   `../clif_workspace/<repo>/logs`, `/intermediate` and `/results`, outside the project
   folder. A path inside the project or inside `data_directory` is refused.

`min_cell_size` (default 10) is the smallest count `release_result()` lets into
`output/final_no_phi/`. See [SOP.md](../SOP.md).
