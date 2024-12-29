/* ------------- C routines for transition models --------------
 *
 * Functions called `C_*` are the ones called from R,
 * all other functions are only used internally (in C).
 *
 * Authors: Nikolaus Umlauf, Reto Stauffer
 * -------------------------------------------------------------
 */

#ifdef _OPENMP // needs to precede R.h
#include <omp.h>
#endif

#include <R.h>
#include <Rinternals.h>
#include "tm.h"

/* Helper function to check if OMP is ON (i.e.,
 * parallelization enabled
 *
 * @return Returns TRUE if OMP enabled, else FALSE.
 * */
SEXP tm_check_omp() {
    #ifdef OPENMP_ON
      return Rf_ScalarLogical(1);
    #else
      return Rf_ScalarLogical(0);
    #endif
}
