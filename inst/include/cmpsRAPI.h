#include <cmpsR.h>

#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>


SEXP COMPUTE_CROSS_CORR_(SEXP xx_in, SEXP yy_in, SEXP minoverlap_in) {
  static SEXP(*fun)(SEXP, SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun = (SEXP(*)(SEXP, SEXP, SEXP)) R_GetCCallable("cmpsR", "compute_cross_corr");
  return fun(xx_in, yy_in, minoverlap_in);
}

SEXP _NA_TRIM(SEXP seq_in) {
  static SEXP(*fun)(SEXP) = NULL;
  if (fun == NULL)
    fun = (SEXP(*)(SEXP)) R_GetCCallable("cmpsR", "na_trim_cmps");
  return fun(seq_in);  
}

SEXP LOCAL_MAX_(SEXP seq_in, SEXP MAX_MIN_in) {
  static SEXP(*fun)(SEXP, SEXP) = NULL;
  if (fun == NULL)
    fun = (SEXP(*)(SEXP, SEXP)) R_GetCCallable("cmpsR", "local_max_cmps");
  return fun(seq_in, MAX_MIN_in);  
}