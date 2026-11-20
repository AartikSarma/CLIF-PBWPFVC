# Intermediate zone I/O (SOP.md, standards 2 and 3). Sourced by utils/config.R.
#
# Every script reads and writes parquet with read_parquet() and write_parquet(). These two
# functions take arrow's place: for a file outside config$intermediate_directory (a core
# CLIF table, anything else) they are arrow's own; for a file inside it they keep dates,
# ages over 89 and other PHI out of what lands on disk, and put the dates back in memory
# when the file is read.
#
#   write   every date-time column becomes seconds from the hospitalization's admission
#           (hospitalization_id -> admission_dttm in the core clif_hospitalization table):
#           x_dttm is stored as x_offset_s, any other date-time column y as y_tsoffset_s.
#           Ages are top-coded at 90 (clif_paths.R: AGE_CAP). The template's
#           phi_columns() then checks what is left and the write stops if it flags a
#           column (a date, geography or a direct identifier).
#   read    the offsets are turned back into date-times from the same admission time,
#           to the second, under their original names, so the scripts' date arithmetic
#           (time from intubation, calendar day and hour, the index year) runs unchanged.
#           col_select works with the original names.
#
# The admission times come from the core table, read in place (standard 1), and are
# held in memory only. A table with a date-time column must carry hospitalization_id.

.arrow_read_parquet  <- arrow::read_parquet
.arrow_write_parquet <- arrow::write_parquet

.in_intermediate <- function(path)
  is.character(path) && length(path) == 1 && .is_within(path, config$intermediate_directory)

.admission_env <- new.env()
.admission_times <- function() {
  if (is.null(.admission_env$a)) {
    f <- file.path(config$tables_path, paste0("clif_hospitalization.", config$file_type))
    if (!file.exists(f)) stop("intermediate_io: no ", f, " to take admission times from")
    a <- switch(config$file_type,
                parquet = .arrow_read_parquet(f, col_select = c("hospitalization_id", "admission_dttm")),
                csv = readr::read_csv(f, col_select = c("hospitalization_id", "admission_dttm"), show_col_types = FALSE),
                fst = fst::read_fst(f, columns = c("hospitalization_id", "admission_dttm")))
    a <- as.data.frame(a)[!duplicated(a$hospitalization_id), ]
    if (!inherits(a$admission_dttm, "POSIXt")) stop("intermediate_io: clif_hospitalization$admission_dttm is not a date-time")
    .admission_env$a <- a
  }
  .admission_env$a
}
.anchor_for <- function(ids, what) {
  a <- .admission_times()
  t <- a$admission_dttm[match(as.character(ids), as.character(a$hospitalization_id))]
  if (anyNA(t[!is.na(ids)])) stop("intermediate_io: ", what, ": ", sum(is.na(t[!is.na(ids)])),
                                  " hospitalization_id(s) not in clif_hospitalization, so their dates cannot be offset")
  t
}

.stored_name <- function(col) if (grepl("_dttm$", col)) sub("_dttm$", "_offset_s", col) else paste0(col, "_tsoffset_s")
.original_name <- function(col) ifelse(grepl("_tsoffset_s$", col), sub("_tsoffset_s$", "", col),
                                       ifelse(grepl("_offset_s$", col), sub("_offset_s$", "_dttm", col), col))

write_parquet <- function(x, sink, ...) {
  if (!.in_intermediate(sink)) return(.arrow_write_parquet(x, sink, ...))
  what <- basename(sink)
  if (any(vapply(x, inherits, logical(1), "Date")))
    stop("intermediate_io: ", what, " has a Date column; keep date-times (POSIXct) or offsets")
  times <- which(vapply(x, inherits, logical(1), "POSIXt"))
  if (length(times)) {
    if (!"hospitalization_id" %in% names(x))
      stop("intermediate_io: ", what, " has date-time columns but no hospitalization_id to offset them from")
    anchor <- as.numeric(.anchor_for(x$hospitalization_id, what))
    for (j in times) {
      # in whole microseconds, truncated as arrow truncates a date-time when it writes
      # one, so the round trip is exact
      x[[j]] <- (trunc(as.numeric(x[[j]]) * 1e6) - trunc(anchor * 1e6)) / 1e6
      names(x)[j] <- .stored_name(names(x)[j])
    }
  }
  for (j in grep(AGE_PATTERN, names(x), ignore.case = TRUE))
    if (is.numeric(x[[j]])) x[[j]] <- pmin(x[[j]], AGE_CAP)
  flagged <- phi_columns(x)
  if (length(flagged))
    stop("intermediate_io: not writing ", what, ": possible PHI in ",
         paste0(names(flagged), " (", flagged, ")", collapse = ", "), ". Reduce it in memory first (SOP.md, standard 3).")
  .arrow_write_parquet(x, sink, ...)
}

read_parquet <- function(file, col_select = NULL, as_data_frame = TRUE, ...) {
  q <- rlang::enquo(col_select)
  if (!.in_intermediate(file)) return(.arrow_read_parquet(file, col_select = !!q, as_data_frame = as_data_frame, ...))
  if (!as_data_frame) stop("intermediate_io: read intermediate files as data frames (as_data_frame = TRUE)")
  stored <- names(arrow::open_dataset(file)$schema)
  original <- .original_name(stored)
  if (!rlang::quo_is_null(q)) {
    pick <- tidyselect::eval_select(q, stats::setNames(as.list(seq_along(original)), original))
    stored <- stored[pick]; original <- original[pick]
  }
  x <- .arrow_read_parquet(file, col_select = dplyr::all_of(stored), ...)
  offs <- which(original != stored)
  if (length(offs)) {
    if (!"hospitalization_id" %in% names(x)) {
      # the key was not selected: read it beside the offsets
      ids <- .arrow_read_parquet(file, col_select = "hospitalization_id")$hospitalization_id
    } else ids <- x$hospitalization_id
    anchor <- .anchor_for(ids, basename(file))
    # summed in whole microseconds and divided as arrow converts them, so each value is
    # the very double the date-time had before it was written
    tz <- attr(anchor, "tzone") %||% "UTC"
    for (j in offs) x[[j]] <- as.POSIXct((trunc(as.numeric(anchor) * 1e6) + round(x[[j]] * 1e6)) / 1e6,
                                         origin = "1970-01-01", tz = tz)
    names(x) <- original
  }
  x
}
