
/* Symbol registration initialization.
 *
 * All functions called from R should be registered here (as well as declared
 * in tm.h). Note that we use a prefix (.fixup = "_C" defined in NAMESPACE).
 * Functions do _not_ have the C_ prefix in C but need the C_ prefix when
 * called (.Call()) from R.
 */
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdbool.h>

/* Include package header file */
#include "transitreg.h"

static const R_CallMethodDef callMethods[] = {
  {"treg_predict",             (DL_FUNC) &treg_predict,          9},
  {"treg_predict_pdfcdf",      (DL_FUNC) &treg_predict_pdfcdf,   6},
  {"treg_detect_cores",        (DL_FUNC) &treg_detect_cores,     0},
  {NULL, NULL, 0} // Termination entry
};

void R_init_TransitionModels(DllInfo *dll) {
    /* Registering .Call functions */
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    /* Disable dynamic symbol lookup for safety reasons */
    R_useDynamicSymbols(dll, FALSE);
}
