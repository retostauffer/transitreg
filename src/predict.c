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
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>

#include <stdlib.h>
#include "tm.h" // Include header, required to have the OPENMP_ON macro available.

/* Helper function to find position of x in y
 *
 * TODO(R): Must update the docstring
 * @param x single integer (the 'needle').
 * @param y integer vector (the 'haystack').
 * @param n length of y.
 * @param count number of elements found (modified inside function).
 * 
 * @return returns an integer vector (must be freed outside function)
 * with the index position of x in y.
 */
tmWhich find_positions(int x, int* y, int n) {
    tmWhich which;
    which.index = (int*)malloc(n * sizeof(int));  // Allocate max possible size
    if (which.index == NULL) { error("Memory allocation failed for tmWhich.index."); }

    which.length = 0;
    for (int i = 0; i < n; i++) {
        if (y[i] == x) {
            which.index[which.length] = i; // Store the index
            which.length++;
        }
    }
    return which; // Return the struct with index position and length
}

/* Helper function for type = "pdf".
 *
 * Calculates (1 - p[count]) * prod(p[-count]) with p = pptr[positions]
 *
 * Note for future me: Using log-sums is slower as we need to take the
 * logarithm of each element in pptr.
 */
double tm_calc_pdf(int* positions, int count, double* pptr) {
    double res = 1.0; // Initialize with 1.0 for product
    for (int i = 0; i < (count - 1); i++) {
        // Calculates product over the first (count - 1) elements
        res *= pptr[positions[i]];
    }
    // Multiplies (1 - p[count]) * the product from above
    res *= (1.0 - pptr[positions[count - 1]]);
    return res;
}

/* Helper function for type = "cdf".
 *
 */
double tm_calc_cdf(int* positions, int count, double* pptr) {
    // Initialize with (1 - p[0])
    double res = 1.0 - pptr[positions[0]];
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        // Looping over all elements except first
        for (int i = 1; i < count; i++) {
            pprod *= pptr[positions[i - 1]]; // Multiply with previous element
            res += (1.0 - pptr[positions[i]]) * pprod;
        }
    }
    return res;
}

/* Helper function for type = "pmax"
 *
 */
double tm_calc_pmax(int* positions, int count, double* pptr) {
    // Initialize with (1 - p[0])
    double res = 1.0 - pptr[positions[0]];
    double max_res = res; // Current maximum
    int    pmax = 0;      // Position of highest value (max cdf)

    // Does the same as tm_calc_cdf except summing up, and only
    // keeps the index of the highest value.
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        for (int i = 1; i < count; i++) {
            pprod *= pptr[positions[i - 1]]; // Multiply with previous element
            res    = (1.0 - pptr[positions[i]]) * pprod;
            if (res > max_res) {
                max_res = res;
                pmax = i;
            }
        }
    }
    return pmax; // Return as is (zero based)
}

/* Calculating elementwise pdf
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param p probabilities, same length as idx vector.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Internally loops over all unique indices in `uidx`.
 * For each index (belonging to one observation) we check the position
 * of the elements in `idx`, and call the `predict_pdf_calc` function
 * which calculates and returns the pdf (double).
 * 
 * @return Returns SEXP double vector of length (length(uidx)).
 */
SEXP tm_predict(SEXP uidx, SEXP idx, SEXP p, SEXP type, SEXP ncores) {

    double *pptr    = REAL(p);
    int    *uidxptr = INTEGER(uidx);     // Unique indices in the dtaa
    int    *idxptr  = INTEGER(idx);      // Index vector
    int    nthreads = asInteger(ncores); // Number of threads for OMP
    int    n        = LENGTH(idx);
    int    un       = LENGTH(uidx);
    int    i;

    // Evaluate 'type' to define what to do
    const char* thetype = CHAR(STRING_ELT(type, 0));
    bool do_pdf = strcmp(thetype, "pdf") == 0;
    bool do_cdf = strcmp(thetype, "cdf") == 0;

    // Custom struct object to mimik "which()"
    tmWhich which;

    // Initialize results vector
    SEXP probs;
    PROTECT(probs = allocVector(REALSXP, un));
    double *probsptr = REAL(probs);


    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        if (do_pdf) {
            probsptr[i] = tm_calc_pdf(which.index, which.length, pptr);
        } else if (do_cdf) {
            probsptr[i] = tm_calc_cdf(which.index, which.length, pptr);
        } else {
            probsptr[i] = tm_calc_pmax(which.index, which.length, pptr);
        }
        free(which.index); // Free allocated memory
    }

    UNPROTECT(1); // Releasing protected objects
    return probs;
}


/* Calculating elementwise pdf and cdf (both at the same time)
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param p probabilities, same length as idx vector.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Does the same as c_tm_predict but calculates both PDF and
 * CDF simultanously, returning a named list. This is used in the
 * main `tm()` function, calculating both at the same time should
 * help to speed up the calculations.
 * 
 * @return Returns SEXP double vector of length (length(uidx)).
 */
SEXP tm_predict_pdfcdf(SEXP uidx, SEXP idx, SEXP p, SEXP ncores) {

    double *pptr    = REAL(p);
    int    *uidxptr = INTEGER(uidx);  // Unique indices in the dtaa
    int    *idxptr  = INTEGER(idx);   // Index vector
    int    nthreads = asInteger(ncores); // Number of threads for OMP
    int    n        = LENGTH(idx);
    int    un       = LENGTH(uidx);
    int    i;

    // Custom struct object to mimik "which()"
    tmWhich which;

    // Initialize results vector
    int nProtected = 0;

    SEXP pdf; PROTECT(pdf = allocVector(REALSXP, un)); ++nProtected;
    double *pdfptr = REAL(pdf);
    SEXP cdf; PROTECT(cdf = allocVector(REALSXP, un)); ++nProtected;
    double *cdfptr = REAL(cdf);


    /* Warning for future me: Do not use Rprintf inside omp -> segfault */
    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        pdfptr[i] = tm_calc_pdf(which.index, which.length, pptr);
        cdfptr[i] = tm_calc_cdf(which.index, which.length, pptr);
        free(which.index); // Free allocated memory
    }

    /* ----------------------------------- */
    /* Generating list */
    SEXP rval;
    PROTECT(rval = allocVector(VECSXP, 2)); ++nProtected;

    /* Adding data to list */
    SET_VECTOR_ELT(rval, 0, pdf);
    SET_VECTOR_ELT(rval, 1, cdf);

    SEXP names_rval;
    PROTECT(names_rval = allocVector(STRSXP, 2)); ++nProtected;

    SET_STRING_ELT(names_rval, 0, mkChar("pdf"));
    SET_STRING_ELT(names_rval, 1, mkChar("cdf"));

    setAttrib(rval, R_NamesSymbol, names_rval);

    UNPROTECT(nProtected);
    return rval;

}


