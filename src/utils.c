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
#include "transitreg.h"

/* Helper function to check if OMP is ON (i.e.,
 * parallelization enabled
 *
 * @return Returns an integer, 0 if OMP is not available,
 * else the number of available cores (procs).
 * */
SEXP treg_detect_cores() {
    int ncores = 0;
    #ifdef OPENMP_ON
        ncores = omp_get_num_procs();
    #endif
    return ScalarInteger(ncores);
}
