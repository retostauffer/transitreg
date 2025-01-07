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

/* VECTOR VERSION, RETURNS AN OBJECT OF CLASS doubleVec,
 * used for elementwise = FALSE */
doubleVec tm_calc_pdf(int* positions, int count, double* tpptr, double* binmidptr, double* y, int ny) {
    // Initialize return value/object
    doubleVec res;
    res.values = (double*)malloc(count * sizeof(double));  // Allocate vector result
    res.length = ny;

    if (res.values == NULL) { error("Memory allocation failed for doubleVec.values."); }

    // Set to true if 'binmidptr' is not provided (an NA). In this case we simply
    // iterate trough the entire tp vector and store the very last value.
    bool nobm = ISNAN(binmidptr[0]);

    // Start calculating
    double prod = 1.0; // Initialize product
    for (int i = 0; i < count; i++) {
        if (ISNAN(tpptr[positions[i]])) {
            error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        // In this case we are only interested in the final PDF (last bin)
        if (nobm) {
            res.values[0] = prod * (1.0 - tpptr[positions[i]]);
            printf(" ---- i = %d   count = %d   res[0] = %.5f\n", i, count, res.values[0]);
        } else {
            res.values[i] = prod * (1.0 - tpptr[positions[i]]);
        }
        // Updating vector product
        prod *= tpptr[positions[i]];
    }
    return res;
}

/* Helper function for type = "cdf".
 *
 */
doubleVec tm_calc_cdf(int* positions, int count, double* pptr) {
    // Initialize return value/object
    doubleVec res;
    res.values = (double*)malloc(count * sizeof(double));  // Allocate max possible size
    res.length = count;

    if (res.values == NULL) { error("Memory allocation failed for doubleVec.values."); }

    if (ISNAN(pptr[positions[0]])) {
        error("TODO(R): First element ISNAN, must be adressed in C");
        //return R_NaReal; 
    }

    // Else start calculation. As soon as we detect a missing
    // value, return NA as well.
    res.values[0] = 1.0 - pptr[positions[0]]; // Initialize with (1 - p[0])
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        // Looping over all elements except first
        for (int i = 1; i < count; i++) {
            if (ISNAN(pptr[positions[i - 1]])) {
                error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
                //return R_NaReal;
            }
            pprod *= pptr[positions[i - 1]]; // Multiply with previous element
            res.values[i] = res.values[i - 1] + (1.0 - pptr[positions[i]]) * pprod;
        }
    }
    return res;
}

/* Helper function for type = "quantile"
 */
doubleVec tm_calc_quantile(int* positions, int count, double* pptr, double* prob, int np) {

    int i, j;

    // Get max(prob), allows for early stopping in the loop
    double pmax = -1;
    for (i = 0; i < np; i++) {
        if (prob[i] > pmax) { pmax = prob[i]; }
    }

    // Initialize return value/object
    doubleVec res;
    res.values = (double*)malloc(np * sizeof(double));  // Allocate vector result
    res.length = np;

    // Temporary double vector to store calculated quantiles
    double* tmp;
    tmp = (double*)malloc(count * sizeof(double)); // Allocate max possible size

    if (res.values == NULL) { error("Memory allocation failed for doubleVec.values."); }

    if (ISNAN(pptr[positions[0]])) {
        error("TODO(R): First element ISNAN, must be adressed in C");
        //return R_NaReal; 
    }

    // TODO(R): Should not be needed, just initializing quantile with 0 for safety/testing
    for (i = 0; i < np; i++) { res.values[i] = -1; } 

    // Initialize with (1 - p[0])
    tmp[0] = 1.0 - pptr[positions[0]];

    // Looping over 'pptr' (counts) to calculate the quantiles; store in 'tmp'.
    // As soon as the calculated quantile tmp[i] is larger than pmax we can stop
    // as we will not need it.
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        for (int i = 1; i < count; i++) {
            if (ISNAN(pptr[positions[i - 1]])) {
                error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
                //return R_NaReal;
            }
            pprod *= pptr[positions[i - 1]]; // Multiply with previous element
            tmp[i] = tmp[i - 1] + (1.0 - pptr[positions[i]]) * pprod;

            // Break for loop as tmp[i] is already larger than pmax
            if (tmp[i] > pmax) { break; }
        }
    }

    // Assign correct quantile to each element in res.values.
    // i: Loops over the quantiles we are looking for
    // j: Loops over calculated quantiles
    j = 0;
    for (i = 0; i < np; i++) {
        for (j = j; j < count; j++) {
            if (prob[i] < tmp[j]) { res.values[i] = j; break; }
        }
    }
    free(tmp); // Freeing allocated memory

    return res;
}

/* Helper function for type = "pmax"
 */
double tm_calc_pmax(int* positions, int count, double* pptr) {
    // If the first value is NA: return NA immediately
    if (ISNAN(pptr[positions[0]])) { return R_NaReal; }

    // Initialize with (1 - p[0])
    double res = 1.0 - pptr[positions[0]];
    double max_res = res; // Current maximum
    int    pmax = 0;      // Position of highest value (max cdf)

    // Does the same as tm_calc_cdf except summing up, and only
    // keeps the index of the highest value. As soon as we
    // detect a missing value - return NA.
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        for (int i = 1; i < count; i++) {
            if (ISNAN(pptr[positions[i - 1]])) { return R_NaReal; }
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
 * @param tp transition probabilities, same length as idx vector.
 * @param binmid mid point of the corresponding bin, if used it must be of
 *        the same length as \code{tp} (see Details).
 * @param y for \code{type = "quantile"} y (numeric) specifies the probabilities
 *        at which the quantile has to be calculated. If \code{type = "cdf"} or
 *        \code{type = "pdf"}, \code{y} evaluates the threshold at which to
 *        evaluate CDF/PDF. See 'Details'.
 * @param y only used if type is 'quantile'; can be NULL (as not used) if
 *        the type is different. Else it is expected to be a vector of doubles
 *        with the same length as 'uidx'.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details This function has a series of different 'modes' for different purposes.
 * Allows to calculate the CDF, PDF, quantiles, as well as the expectation (pmax).
 * 
 * The three inputs which are always required are:
 *
 * * \code{uidx}: Unique 'distribution' index.
 * * \code{idx}: Integer vector specifying which transition probability
 *   (and binmid if used) corresponds to which \code{uidx}.
 * * \code{type}: What to return.
 * * \code{ncores}: Number of cores to be used when OpenMP is available.
 *
 * The use of the remaining arguments differ.
 *
 * Quantiles (\code{type = "quantile"})
 * ------------------------------------
 * \code{y} must be numeric, same length as \code{uidx}. Defines the probabilities
 * at which the different distributions have to be evaluated. In case
 * \code{elementwise} is set \code{true}, it specifies the probabilities at which
 * each distribution is evaluaed (thus, the length must not be equal to the
 * length of \code{uidx}).
 *
 *
 * @return TODO(R): Depends on mode.
 */
SEXP tm_predict(SEXP uidx, SEXP idx, SEXP tp, SEXP binmid, SEXP y,
                SEXP type, SEXP ncores, SEXP elementwise) {

    int    *uidxptr   = INTEGER(uidx);     // Unique indices in the dtaa
    int    *idxptr    = INTEGER(idx);      // Distribution index
    int    nthreads   = asInteger(ncores); // Number of threads for OMP
    int    n          = LENGTH(idx);
    int    un         = LENGTH(uidx);
    int    i, j, np;

    double *tpptr     = REAL(tp);          // Transition probabilities
    double *binmidptr = REAL(binmid);      // Bin mid point, not always needed
    double *yptr      = REAL(y);           // Where to evaluate the distribution

    bool   ewise      = asLogical(elementwise); // C boolean value

    // Evaluate 'type' to define what to do. Store a set of
    // boolean values to only do the string comparison once.
    const char* thetype = CHAR(STRING_ELT(type, 0));
    bool do_pdf  = strcmp(thetype, "pdf")  == 0;
    bool do_cdf  = strcmp(thetype, "cdf")  == 0;
    bool do_q    = strcmp(thetype, "quantile") == 0;
    // ... if none of them is true, it must be "do pmax"
    // bool do_pmax = strcmp(thetype, "pmax") == 0;

    // Boolean flat which is set true if 'binmid' is not provided.
    ////bool no_binmid = (LENGTH(binmid) == 1) & ISNAN(binmidptr[0]);
    ////printf(" --------------------- no_binmid : %d    %.5f\n", no_binmid, binmidptr[0]);

    // Allocating return vector.

    // If type is "pdf", "cdf", or "quantile" and ewise is false: the length of
    // the vector is the same as 'un', else 'un * ny' (number of distributions
    // times number of probabilites/thresholds at which each distribution is
    // evaluated).
    SEXP res;
    if (do_pdf | do_cdf | do_q) {
        np = (!ewise) ? 1 : LENGTH(y);
        PROTECT(res = allocVector(REALSXP, un * np));
    // Else it is type = "pmax", so length of res is equal to 'un'.
    } else {
        np = 1;
        PROTECT(res = allocVector(REALSXP, un));
    }

    // Pointer on results vector 'res'
    double *resptr = REAL(res);

    // If binmid is NA, but yptr is not: Error (this combination is not allowed)
    if (ISNAN(binmidptr[0]) & !ISNAN(yptr[0])) {
        error("Problem in C tm_predict(): binmid is NAN, but y is not.");
    }

    // If mode is not pdf, cdf, or quantile, elementwise must be false.
    // Else we throw an error here. That is for pmax where elementwise
    // makes no sense.
    if (!do_pdf & !do_cdf & !do_q & ewise) {
        printf("Using \"%s\" with elementwise = true not allowed.\n", CHAR(STRING_ELT(type, 0)));
        exit(99);
    }

    // Custom struct object to mimik "which()"
    tmWhich which;
    doubleVec tmp;

    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which, tmp)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        if (do_pdf | do_cdf | do_q) {
            // --- Calculating probability density
            if (do_pdf) {
                // Single PDF
                if (!ewise) {
                    tmp = tm_calc_pdf(which.index, which.length, tpptr, binmidptr, yptr, 1);
                    //tmp = tm_calc_pdf(which.index, which.length, tpptr, binmidptr, &yptr[i], 1);
                // Multiple PDFs (elementwise)
                } else {
                    tmp = tm_calc_pdf(which.index, which.length, tpptr, binmidptr, yptr, LENGTH(y));
                }
            // --- Calculating cumulative distribution
            } else if (do_cdf) {
                tmp = tm_calc_cdf(which.index, which.length, tpptr);
            } else {
                // Single quantile
                if (!ewise) {
                    tmp = tm_calc_quantile(which.index, which.length, tpptr, &yptr[i], 1);
                // Multiple quantiles (elementwise)
                } else {
                    tmp = tm_calc_quantile(which.index, which.length, tpptr, yptr, LENGTH(y));
                }
            }

            // Store results. If !elementwise, store last value
            if (!ewise) {
                resptr[i] = tmp.values[tmp.length - 1];
            // Else store the entire vector
            } else {
                for (j = 0; j < tmp.length; j++) { resptr[i * np + j] = tmp.values[j]; }
            }
            free(tmp.values); // Free allocated memory
        // Else it must be pmax
        } else {
            resptr[i] = tm_calc_pmax(which.index, which.length, tpptr);
        }
        free(which.index); // Free allocated memory
    }

    UNPROTECT(1); // Releasing protected objects
    return res;
}


/* Calculating pdf and cdf (both at the same time)
 *
 * This mode is used when estimating the Transition Model. Calculates
 * the PDF and CDF for each observation at the last bin (i.e., the bin
 * in which the observation falls into). Calls tm_calc_pdf and tm_calc_cdf
 * with a missing value on 'binmidptr' and 'ny = 1' which tells tm_calc_pdf/tm_calc_cdf
 * to use this specific "mode".
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param p probabilities, same length as idx vector.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Does something similar to tm_predict but calculates both PDF and
 * CDF simultanously for the very last bin in each distribution (observation),
 * returning a named list. This is used in the main `tm()` function,
 * calculating both at the same time should help to speed up the calculations.
 * 
 * @return Returns named list with two numeric vectors, each of which
 * has length(uidx) (vector with cdf and pdf).
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
    doubleVec tmppdf;
    doubleVec tmpcdf;

    // Initialize results vector
    int nProtected = 0;

    SEXP pdf; PROTECT(pdf = allocVector(REALSXP, un)); ++nProtected;
    double *pdfptr = REAL(pdf);
    SEXP cdf; PROTECT(cdf = allocVector(REALSXP, un)); ++nProtected;
    double *cdfptr = REAL(cdf);

    // Dummy value as we need a proper object when calling tm_calc_*() below
    double* na = malloc(sizeof(double)); // Single double pointer
    na[0] = 0.0 / 0.0; // Assign nan (missing value)

    /* Warning for future me: Do not use Rprintf inside omp -> segfault */
    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which, tmppdf, tmpcdf)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        tmppdf = tm_calc_pdf(which.index, which.length, pptr, na, na, 1);
        tmpcdf = tm_calc_cdf(which.index, which.length, pptr);
        // Store last value, that is the last bin provided for this distribution.
        pdfptr[i] = tmppdf.values[tmppdf.length - 1];
        cdfptr[i] = tmpcdf.values[tmpcdf.length - 1];
        // Free allocated memory
        free(which.index);
        free(tmppdf.values);
        free(tmpcdf.values);
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


