# =============================================================================
# The R session and the compiled libraries this machine runs the pipeline with
# =============================================================================
#
# For reporting a crash: R and package versions, the BLAS and LAPACK in use, and how
# many threads arrow and data.table will start. Prints no data and reads no file.
#
# Usage: uvr run code/tools/session_info.R
# =============================================================================

print(sessionInfo())
cat("\nLAPACK: ", La_library(), "\n", sep = "")
cat("\nThread settings in the environment:\n")
vars <- c("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "R_DATATABLE_NUM_THREADS")
print(Sys.getenv(vars, unset = "(not set)"))
cat("\narrow:\n")
print(arrow::arrow_info())
cat("\narrow CPU threads: ", arrow::cpu_count(), ", I/O threads: ", arrow::io_thread_count(), "\n", sep = "")
cat("\ndata.table:\n")
data.table::getDTthreads(verbose = TRUE)
