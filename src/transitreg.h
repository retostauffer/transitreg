
/* Strongly inspired by the great mgcv package! */

#include <stdbool.h>

// /* Most compilers with openMP support supply a pre-defined compiler macro
//  * _OPENMP. Following facilitates selective turning off (by testing value or
//  * defining multiple versions OPENMP_ON1, OPENMP_ON2...)
//  */
// #ifdef _OPENMP
// //#pragma message(" [dev] Defining OPENMP_ON 1")
// #define OPENMP_ON 1
// #else
// //#pragma message(" [dev] Defining OPENMP_ON 0")
// #define OPENMP_ON 0
// #endif


/* ... note also that there is no actual *need* to protect #pragmas with #ifdef
 * OPENMP_ON, since C ignores undefined pragmas, but failing to do so may
 * produce alot of compilation warnings if openMP is not supported In contrast
 * functions from omp.h must be protected, and there is non-avoidable use of
 * these in the mgcv code.
 */

//#define OMP_REPORT // define to have all routines using omp report on start and end.
#define OMP_REPORT

/* For safe memory handling from R */
#define CALLOC R_chk_calloc
#define FREE R_chk_free
#define REALLOC R_chk_realloc

/* BUT, this can mess up valgrinding for memory error checking - problems are.
 * sometimes missed because standard allocation is being circumvented. Then
 * errors can. corrupt R memory management without detection and trigger
 * nothing until R messes up internally because of corruption, which then makes
 * it look as if R is generating the problem. Hence better to reset for
 * checking. Also sizing errors in .C often generate no obvious valgrind error.
 */
//#define CALLOC calloc
//#define FREE free
//#define REALLOC realloc

/* ------------------------------------------------------------------------ */

/* Custom type; structured object with   
 * an integer vector with index and length of the resulting index vector.
 * NOTE: .index must be freed by the user! */
typedef struct {
    int* index;
    int length;
} integerVec;

/* Custom type: stuctured object with ...
 * a double vector and length.
 * NOTE: .values must be freed by the user! */
typedef struct {
    double* values;
    int length;
} doubleVec;

/* The custom struct is used by find_position.
 * Searches for 'x' (int) in integer vector 'y' with a max length of 'n'.
 *
 * Returns struct object with:
 *   -   .index:   integer vector. position of 'x' in 'y'.
 *   -   .length:  length of .index, or number of 'x' in 'y'.
 */
integerVec find_positions(int x, int* y, int n);

void eval_bins_quantile(double* res, double* tmp, int* positions, int count,
        double* bkptr, double* prob, int np, bool disc);

double interpolate_linear(double x1, double y1, double x2, double y2, double p);

doubleVec treg_calc_pdf(int* positions, int count, double* tpptr,
        double* bkptr, int nbins, int* y, int ny, bool disc);

doubleVec treg_calc_cdf(int* positions, int count, double* tpptr,
        double* bkptr, int nbins, int* y, int ny);

doubleVec treg_calc_quantile(int* positions, int count, double* tpptr,
        double* bkptr, double* prob, int np, bool disc);

doubleVec get_binmid(SEXP breaks);

double treg_calc_mode(int* positions, int count, double* pptr, double* binsptr);
double treg_calc_mean(int* positions, int count, double* tpptr, double* binsptr);

SEXP treg_predict(SEXP uidx, SEXP idx, SEXP tp, SEXP breaks, SEXP y, SEXP prob,
        SEXP type, SEXP ncores, SEXP elementwise, SEXP discrete);
SEXP treg_predict_pdfcdf(SEXP uidx, SEXP idx, SEXP tp, SEXP y,
        SEXP breaks, SEXP discrete, SEXP ncores);
SEXP treg_detect_cores();


