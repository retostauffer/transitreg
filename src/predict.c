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
doubleVec tm_calc_pdf(int* positions, int count, double* pptr) {
    // Initialize return value/object
    doubleVec res;
    res.values = (double*)malloc(count * sizeof(double));  // Allocate max possible size
    res.length = count;

    if (res.values == NULL) { error("Memory allocation failed for doubleVec.values."); }

    // Start calculating
    double prod = 1.0; // Initialize product
    for (int i = 0; i < count; i++) {
        if (ISNAN(pptr[positions[i]])) {
            error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        res.values[i] = prod * (1.0 - pptr[positions[i]]);
        // Updating vector product
        prod *= pptr[positions[i]];
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
    res.values = (double*)malloc(np * sizeof(double));  // Allocate max possible size
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
 * @param p probabilities, same length as idx vector.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param prob only used if type is 'quantile'; can be NULL (as not used) if
 *        the type is different. Else it is expected to be a vector of doubles
 *        with the same length as 'uidx'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Internally loops over all unique indices in `uidx`.
 * For each index (belonging to one observation) we check the position
 * of the elements in `idx`, and call the `predict_pdf_calc` function
 * which calculates and returns the pdf (double).
 * 
 * @return Returns SEXP double vector of length (length(uidx)).
 */
SEXP tm_predict(SEXP uidx, SEXP idx, SEXP p, SEXP type, SEXP prob,
                SEXP ncores, SEXP elementwise) {

    double *pptr    = REAL(p);
    int    *uidxptr = INTEGER(uidx);     // Unique indices in the dtaa
    int    *idxptr  = INTEGER(idx);      // Index vector
    int    nthreads = asInteger(ncores); // Number of threads for OMP
    int    n        = LENGTH(idx);
    int    un       = LENGTH(uidx);
    bool   ewise    = asLogical(elementwise); // C boolean value
    int    i, j, np;

    // This is only used if type == "quantile"; on the R side it is ensured
    // that 'prob' is a numeric vector of LENGTH(uidx) if type = 'quantile',
    // thus we do not check it here.
    double *probptr = REAL(prob);        // Probability used for quantiles


    // Evaluate 'type' to define what to do. Store a set of
    // boolean values to only do the string comparison once.
    const char* thetype = CHAR(STRING_ELT(type, 0));
    bool do_pdf  = strcmp(thetype, "pdf")  == 0;
    bool do_cdf  = strcmp(thetype, "cdf")  == 0;
    bool do_q    = strcmp(thetype, "quantile") == 0;
    // ... if none of them is true, it must be "do pmax"
    // bool do_pmax = strcmp(thetype, "pmax") == 0;

    // Allocating return vector.
    // If type is "pdf" or "df" and ewise is FALSE: the length of the
    //     vector is the same as 'un', else 'un * LENGTH(p);'.
    SEXP res;
    if (do_pdf | do_cdf) {
        np = (!ewise) ? 1 : LENGTH(p);
        PROTECT(res = allocVector(REALSXP, un * np));
    // If type is "quantile" and ewise is FALSE: the length of the
    //     vector is also the same length as 'un', else 'un' times LENGTH(prob);
    } else if (do_q) {
        np = (!ewise) ? 1 : LENGTH(prob);
        PROTECT(res = allocVector(REALSXP, un * np));
    // Else it is type = "pmax", so length of res is equal to 'un'.
    } else {
        np = 1;
        PROTECT(res = allocVector(REALSXP, un));
    }

    // Pointer on results vector 'res'
    double *resptr = REAL(res);

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
    #pragma omp parallel for num_threads(nthreads) private(which)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        if (do_pdf | do_cdf | do_q) {
            // --- Calculating probability density
            if (do_pdf) {
                tmp = tm_calc_pdf(which.index, which.length, pptr);
            // --- Calculating cumulative distribution
            } else if (do_cdf) {
                tmp = tm_calc_cdf(which.index, which.length, pptr);
            } else {
                // Single quantile
                if (!ewise) {
                    tmp = tm_calc_quantile(which.index, which.length, pptr, &probptr[i], 1);
                // Multiple quantiles (elementwise)
                } else {
                    tmp = tm_calc_quantile(which.index, which.length, pptr, probptr, LENGTH(prob));
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
            resptr[i] = tm_calc_pmax(which.index, which.length, pptr);
        }
        free(which.index); // Free allocated memory
    }

    UNPROTECT(1); // Releasing protected objects
    return res;
}


/* Calculating pdf and cdf (both at the same time)
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param p probabilities, same length as idx vector.
 * @param type character, either 'pdf', 'cdf', or 'pmax'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Does the same as tm_predict but calculates both PDF and
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
    doubleVec tmppdf;
    doubleVec tmpcdf;

    // Initialize results vector
    int nProtected = 0;

    SEXP pdf; PROTECT(pdf = allocVector(REALSXP, un)); ++nProtected;
    double *pdfptr = REAL(pdf);
    SEXP cdf; PROTECT(cdf = allocVector(REALSXP, un)); ++nProtected;
    double *cdfptr = REAL(cdf);


    /* Warning for future me: Do not use Rprintf inside omp -> segfault */
    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which, tmppdf, tmpcdf)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        tmppdf = tm_calc_pdf(which.index, which.length, pptr);
        tmpcdf = tm_calc_cdf(which.index, which.length, pptr);
        // Store last value
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


