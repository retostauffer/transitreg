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


/* Helper function for type = "mean".
 * Calculates elementwise pdf and multiplies it with center of the bin
 * ((lowerptr[i] + upperptr[i]) * 0.5) to get the weighted average.
 * Returns a (doubleVec of length 1) */
double tm_calc_mean(int* positions, int count, double* tpptr, double* lowerptr, double* upperptr) {

    // Initialize return value, initialize sum with 0
    double res = 0.0;

    // Calculate PDF for each bin given by the distribution for i = 0, ..., count - 1.
    // Store in double 'tmp', the required values will be extracted after this loop.
    double prod = 1.0; // Initialize product
    for (int i = 0; i < count; i++) {
        if (ISNAN(tpptr[positions[i]])) {
            error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        // Updating product of transition probabilities
        prod *= tpptr[positions[i]];
        // Updating weighted sum, simply "binmid * pdf(bin)" summed up
        // (lowerptr[i] + upperptr[i]) * 0.5 -> center of the bin
        res += (lowerptr[i] + upperptr[i]) * 0.5 * prod * (1.0 - tpptr[positions[i]]);
    }
    return res;
}


/* Helper function for type = "pdf" */
doubleVec tm_calc_pdf(int* positions, int count,
                      double* tpptr, double* lowerptr, double* upperptr,
                      double* y, int ny) {

    // Initialize return value/object.
    doubleVec res;
    res.values = (double*)malloc(ny * sizeof(double));  // Allocate vector result
    res.length = ny;

    // Temporary double vector to calculate PDF along i = 0, ..., count - 1
    double* tmp = malloc(count * sizeof(double)); // Single double pointer

    // Set to true if 'lowerptr[0]' or 'upperptr[0]' not given, this is used
    // when calling 'tm_predict_pdfcdf' as we do not need the bins there (i.e.,
    // 'lowerptr' and 'upperptr' are of length 1 and contain NA.
    // In this case we simply iterate trough the entire tp vector and store the
    // very last value.
    bool nobm = ISNAN(lowerptr[0]) | ISNAN(upperptr[0]);

    // Calculate PDF for each bin given by the distribution for i = 0, ..., count - 1.
    // Store in double 'tmp', the required values will be extracted after this loop.
    double prod = 1.0; // Initialize product
    for (int i = 0; i < count; i++) {
        if (ISNAN(tpptr[positions[i]])) {
            error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        // Updating product of transition probabilities
        prod *= tpptr[positions[i]];
        // Updating temporary PDF vector
        tmp[i] = prod * (1.0 - tpptr[positions[i]]);
    }

    // If nobm (no bin mids given) the result is simply the last value in tmp;
    // the PDF evaluated at the end of the entire distribution.
    if (nobm) {
        res.values[0] = tmp[count - 1];
    // Else we are calling eval_bins_pdf_cdf which evaluates
    // into which bin each of the thresholds/values in y fall into,
    // and assign the corresponding PDF to res.values | y.
    } else {
        eval_bins_pdf_cdf(res.values, tmp, positions, count, lowerptr, upperptr, y, ny);
    }

    // Free allocated memory, return result
    free(tmp);
    return res;
}

void eval_bins_pdf_cdf(double* res, double* tmp, int* positions, int count,
                       double* lowerptr, double* upperptr, double* y, int ny) {

    // Guesstimate/calculate bin width
    int i, j;

    i = 0;
    for (j = 0; j < ny; j++) {
        res[j] = NA_REAL; // Initial value
        for (i = i; i < count; i++) {
            // Current 'y' below current bin OR above highest bin (outside
            // upper support) we break the inner loop.
            if ((y[j] <= lowerptr[positions[i]]) | (y[j] > upperptr[positions[count - 1]])) {
                break;
            }
            // If y[j] in current bin: Store PDF, break loop and continue search
            if ((y[j] > lowerptr[positions[i]]) & (y[j] <= upperptr[positions[i]])) {
                res[j] = tmp[i];
                break;
            }
        }
    }
}

/* Helper function for type = "cdf" */
doubleVec tm_calc_cdf(int* positions, int count,
                      double* tpptr, double* lowerptr, double* upperptr,
                      double* y, int ny) {

    // Initialize return value/object.
    doubleVec res;
    res.values = (double*)malloc(ny * sizeof(double));  // Allocate vector result
    res.length = ny;

    // Temporary double vector to calculate PDF along i = 0, ..., count - 1
    double* tmp = malloc(count * sizeof(double)); // Single double pointer

    // Set to true if 'binmidptr' is not provided (an NA). In this case we simply
    // iterate trough the entire tp vector and store the very last value.
    bool nobm = ISNAN(lowerptr[0]) | ISNAN(upperptr[0]);

    // Calculate CDF for each bin given by the distribution for i = 0, ..., count - 1
    // Store in double 'tmp', the required values will be extracted after this loop.
    tmp[0] = 1.0 - tpptr[positions[0]]; // Initialize with (1 - p[0])
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        // Looping over all elements except first
        for (int i = 1; i < count; i++) {
            if (ISNAN(tpptr[positions[i - 1]])) {
                error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
                //return R_NaReal;
            }
            pprod *= tpptr[positions[i - 1]]; // Multiply with previous element
            // If 'nobm' is true (we have no binmidptr) we only need the PDF
            // of the final bin; overwrite res.values[0] in each iteration.
            tmp[i] = tmp[i - 1] + (1.0 - tpptr[positions[i]]) * pprod;
        }
    }

    // If nobm (no bin mids given) the result is simply the last value in tmp;
    // the PDF evaluated at the end of the entire distribution.
    if (nobm) {
        res.values[0] = tmp[count - 1];
    // Else we are calling eval_bins_pdf_cdf which evaluates
    // into which bin each of the thresholds/values in y fall into,
    // and assign the corresponding PDF to res.values | y.
    } else {
        eval_bins_pdf_cdf(res.values, tmp, positions, count, lowerptr, upperptr, y, ny);
    }

    // Free allocated memory, return result
    free(tmp);
    return res;
}


/* Linear interpolation for finer quantiles */
double interpolate_linear(double x1, double y1, double x2, double y2, double p) {
    if (ISNAN(x1) | ISNAN(y1) | ISNAN(x2) | ISNAN(y2)) {
        return NA_REAL;
    } else if (x1 == x2) {
        return x1;
    }
    // Performing interpolation
    return x1 + (p - y1) * (x2 - x1) / (y2 - y1);
}

/* Helper function for type = "quantile" */
doubleVec tm_calc_quantile(int* positions, int count,
                           double* tpptr, double* lowerptr, double* upperptr,
                           double* prob, int np, bool disc) {

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

    if (ISNAN(tpptr[positions[0]])) {
        error("TODO(R): First element ISNAN, must be adressed in C");
        //return R_NaReal; 
    }

    // Initialize the results vector with -nan; just for safety reasons, should not be needed
    for (i = 0; i < np; i++) { res.values[i] = NA_REAL; } 

    // Initialize with (1 - p[0])
    tmp[0] = 1.0 - tpptr[positions[0]];

    // Looping over 'tpptr' (counts) to calculate the quantiles; store in 'tmp'.
    // As soon as the calculated quantile tmp[i] is larger than pmax we can stop
    // as we will not need it.
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        for (int i = 1; i < count; i++) {
            if (ISNAN(tpptr[positions[i - 1]])) {
                error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
                //return R_NaReal;
            }
            pprod *= tpptr[positions[i - 1]]; // Multiply with previous element
            tmp[i] = tmp[i - 1] + (1.0 - tpptr[positions[i]]) * pprod;

            // Break for loop as tmp[i] is already larger than pmax
            if (tmp[i] > pmax) { break; }
        }
    }

    // Assign correct quantile to each element in res.values.
    // i: Loops over the quantiles we are looking for
    // j: Loops over calculated quantiles
    // Store bin mid (quantile we are looking for)
    j = disc ? 1 : 0; // Discrete distributions: start at j = 1, else j = 0
    for (i = 0; i < np; i++) {
        for (j = j; j < count; j++) {
            // Here we distinguish between discrete and non-discrete distributions.
            // If discrete, the quantile _must_ lie within a bin, and we take the
            // center of the bin.
            if (disc) {
                if ((prob[i] >= tmp[j - 1]) & (prob[i] < tmp[j])) {
                    res.values[i] = (lowerptr[positions[j]] + upperptr[positions[j]]) * 0.5;
                }
            // If not discrete we do the following:
            // (1) If j == 0 we set the quantile to the lower support of the distribution at first.
            // (2) If j == (count - 1), the last iteration, we store the upper end of the support.
            // (3) Else we check if we fall into the current bin. If so, we use linear interpolation
            //     within the bin (lower-upper); this overwrites rules (1) and (2).
            } else {
                // If j = 0 take lowest possible value
                if (j == 0) {
                    res.values[i] = lowerptr[positions[j]];
                // If j == (count - 1), last iteration, take highest possible value
                } else if (j == (count - 1)) {
                    res.values[i] = upperptr[positions[j]];
                } 
                // Check if we fall into this bin
                if ((prob[i] >= tmp[j - 1]) & (prob[i] < tmp[j])) {
                    // Perform linear interpolation between the two neighboring bin mids.
                    res.values[i] = interpolate_linear(lowerptr[positions[j]], tmp[j - 1],
                                                       upperptr[positions[j]], tmp[j],
                                                       prob[i]);
                    break; // Found what we were looking for, break inner loop
                }
            }
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
 * * \code{idx}: Integer vector specifying which transition probability.
 * * \code{lower}: Lower edge of the bin.
 * * \code{upper}: Upper edge of the bin.
 * * \code{y}: Where to evaluate the CDF, PDF, quantile. Unused for 
 *             type = "pmax" and type = "mean".
 * * \code{type}: What to return.
 * * \code{ncores}: Number of cores to be used when OpenMP is available.
 * * \code{elementwise}: Logical, should the 'type' be evaluated elementwise?
 *             Only used for PDF, CDF, quantile. If elementwise, y must have
 *             the same length as ui. Else each distribution is evaluated
 *             at all 'y's.
 * * \code{discrete}: Logical, if true the quantile will be fixed to the
 *             center of the bin, else we perform linear interpolation.
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
SEXP tm_predict(SEXP uidx, SEXP idx, SEXP tp, SEXP lower, SEXP upper, SEXP y,
                SEXP type, SEXP ncores, SEXP elementwise, SEXP discrete) {

    int    *uidxptr   = INTEGER(uidx);     // Unique indices in the dtaa
    int    *idxptr    = INTEGER(idx);      // Distribution index
    int    nthreads   = asInteger(ncores); // Number of threads for OMP
    int    n          = LENGTH(idx);
    int    un         = LENGTH(uidx);
    int    i, j, np;

    double *tpptr     = REAL(tp);          // Transition probabilities
    double *lowerptr  = REAL(lower);       // Lower bound of bin
    double *upperptr  = REAL(upper);       // Upper bound of bin
    double *yptr      = REAL(y);           // Where to evaluate the distribution

    bool   ewise      = asLogical(elementwise); // Elementwise calculation?
    int*   discptr    = LOGICAL(discrete);      // Discrete distribution? For quantile

    // Evaluate 'type' to define what to do. Store a set of
    // boolean values to only do the string comparison once.
    const char* thetype = CHAR(STRING_ELT(type, 0));
    bool do_pdf  = strcmp(thetype, "pdf")  == 0;
    bool do_cdf  = strcmp(thetype, "cdf")  == 0;
    bool do_q    = strcmp(thetype, "quantile") == 0;
    bool do_mean = strcmp(thetype, "mean") == 0;
    // ... if none of them is true, it must be "do pmax"
    // bool do_pmax = strcmp(thetype, "pmax") == 0;

    // Allocating return vector.

    // If type is "pdf", "cdf", or "quantile" and ewise is true: the length of
    // the vector is the same as 'un', else 'un * ny' (number of distributions
    // times number of probabilites/thresholds at which each distribution is
    // evaluated).
    SEXP res;
    if (do_pdf | do_cdf | do_q) {
        np = (ewise) ? 1 : LENGTH(y);
        PROTECT(res = allocVector(REALSXP, un * np));
    // Else it is type = "pmax", so length of res is equal to 'un'.
    } else {
        np = 1;
        PROTECT(res = allocVector(REALSXP, un));
    }

    // Pointer on results vector 'res'
    double *resptr = REAL(res);

    // If binmid is NA, but yptr is not: Error (this combination is not allowed)
    if ((do_cdf | do_pdf) & (ISNAN(lowerptr[0]) | ISNAN(upperptr[0])) & !ISNAN(yptr[0])) {
        error("Problem in C tm_predict(): lower or upper is NA, but y is not.");
    }

    // If mode is not pdf, cdf, or quantile, elementwise must be true.
    // Else we throw an error here. That is for pmax where elementwise
    // makes no sense.
    if (!do_pdf & !do_cdf & !do_q & !ewise) {
        printf("Using \"%s\" with elementwise = false not allowed.\n", CHAR(STRING_ELT(type, 0)));
        exit(99);
    }

    // Custom struct object to mimik "which()"
    tmWhich which;
    doubleVec tmp;

    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which, tmp, j)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);
        if (do_pdf | do_cdf | do_q) {
            // --- Calculating probability density
            if (do_pdf) {
                // Single PDF
                if (ewise) {
                    tmp = tm_calc_pdf(which.index, which.length, tpptr, lowerptr, upperptr, yptr, 1);
                // Multiple PDFs (elementwise)
                } else {
                    tmp = tm_calc_pdf(which.index, which.length, tpptr, lowerptr, upperptr, yptr, LENGTH(y));
                }
            // --- Calculating cumulative distribution
            } else if (do_cdf) {
                // Single CDF
                if (ewise) {
                    tmp = tm_calc_cdf(which.index, which.length, tpptr, lowerptr, upperptr, yptr, 1);
                // Multiple CDFs (elementwise)
                } else {
                    tmp = tm_calc_cdf(which.index, which.length, tpptr, lowerptr, upperptr, yptr, LENGTH(y));
                }
            } else {
                // Single quantile
                if (ewise) {
                    tmp = tm_calc_quantile(which.index, which.length, tpptr, lowerptr, upperptr, &yptr[i], 1, discptr[i] == 1);
                // Multiple quantiles (elementwise)
                } else {
                    tmp = tm_calc_quantile(which.index, which.length, tpptr, lowerptr, upperptr, yptr, LENGTH(y), discptr[i] == 1);
                }
            }

            // Store results. If !elementwise, store last value
            if (ewise) {
                resptr[i] = tmp.values[tmp.length - 1];
            // Else store the entire vector
            } else {
                for (j = 0; j < tmp.length; j++) { resptr[i * np + j] = tmp.values[j]; }
            }

            // Free allocated memory
            free(tmp.values);

        // -----------------------------------------------------------
        // These are the types where the result is always a single
        // value (one per distribution)
        // -----------------------------------------------------------
        // Calculate mean
        } else if (do_mean) {
            resptr[i] = tm_calc_mean(which.index, which.length, tpptr, lowerptr, upperptr);
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
SEXP tm_predict_pdfcdf(SEXP uidx, SEXP idx, SEXP tp, SEXP ncores) {

    double *tpptr   = REAL(tp);
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
    na[0] = NA_REAL; // Assign missing value

    /* Warning for future me: Do not use Rprintf inside omp -> segfault */
    #if OPENMP_ON
    #pragma omp parallel for num_threads(nthreads) private(which, tmppdf, tmpcdf)
    #endif
    for (i = 0; i < un; i++) {
        // Search position of uidxptr[i] in vector idxptr
        which = find_positions(uidxptr[i], idxptr, n);

        // Calculate pdf and cdf (for the last bin)
        // Input arguments are (in this order)
        //     positions, count, tpptr, lowerptr, upperptr, y, ny
        // Here lowerptr, upperptr, y, (and ny) are just dummy values!
        tmppdf = tm_calc_pdf(which.index, which.length, tpptr, na, na, na, 1);
        tmpcdf = tm_calc_cdf(which.index, which.length, tpptr, na, na, na, 1);

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


