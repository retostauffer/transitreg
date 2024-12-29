
/* Strongly inspired by the great mgcv package! */

/* Most compilers with openMP support supply. a pre-defined compiler macro
 * _OPENMP. Following. facilitates selective turning off (by testing value or
 * defining multiple versions OPENMP_ON1, OPENMP_ON2...)
 */
#if defined _OPENMP
#define OPENMP_ON 1
#endif

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

void fun(double *y, double *H);
int* find_position(int x, int* y, int n, int* count);
double tm_calc_pdf(int* positions, int count, double* pptr);
double tm_calc_cdf(int* positions, int count, double* pptr);
double tm_calc_pmax(int* positions, int count, double* pptr);

SEXP tm_predict(SEXP uidx, SEXP idx, SEXP p, SEXP type);
SEXP tm_predict_pdfcdf(SEXP uidx, SEXP idx, SEXP p);


