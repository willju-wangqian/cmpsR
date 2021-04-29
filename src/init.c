#include "cmps.h"
#include <Rconfig.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const
R_CallMethodDef callMethods[] = {
        {"add",                 (DL_FUNC) &add_,                 2},
        {"compute_cross_corr",  (DL_FUNC) &COMPUTE_CROSS_CORR_,  3}, 
        {"na_trim_cmps",        (DL_FUNC) &NA_TRIM_,             1},
        {NULL,              NULL,                          0}
};

void R_init_CMPS(DllInfo *info)
{
//   SymbolShortcuts();
  R_registerRoutines(info,
                     NULL,
                     callMethods,
                     NULL,
                     NULL);

  R_useDynamicSymbols(info, TRUE);

  /* used by external packages linking to internal xts code from C */
  R_RegisterCCallable("CMPS","add",                     (DL_FUNC) &add_);
  R_RegisterCCallable("CMPS","compute_cross_corr",      (DL_FUNC) &COMPUTE_CROSS_CORR_);
  R_RegisterCCallable("CMPS","na_trim_cmps",            (DL_FUNC) &NA_TRIM_);
  
}
