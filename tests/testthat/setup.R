withr::local_seed(20260513, .local_envir = teardown_env())

## A multithreaded BLAS (OpenBLAS on some Linux check machines) otherwise uses
## more CPU time than elapsed time in the characteristic-function fits.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  blas_threads <- RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1L)
  withr::defer(RhpcBLASctl::blas_set_num_threads(blas_threads), teardown_env())
}
dt_threads <- data.table::getDTthreads()
data.table::setDTthreads(1L)
withr::defer(data.table::setDTthreads(dt_threads), teardown_env())
