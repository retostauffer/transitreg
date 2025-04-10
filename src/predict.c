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
#include "transitreg.h" // Include header

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
integerVec find_positions(int x, int* y, int n) {
    integerVec which;
    which.index = (int*)malloc(n * sizeof(int));  // Allocate max possible size
    if (which.index == NULL) { Rf_error("Memory allocation failed for integerVec.index."); }

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
 * ((bkptr[i + 1] + bkptr[i]) * 0.5) to get the weighted average.
 * Returns a (doubleVec of length 1) */
double treg_calc_mean(int* positions, int count, double* tpptr, double* bkptr) {

    // Initialize return value, initialize sum with 0
    double res = 0.0;

    // Calculate PDF for each bin given by the distribution for i = 0, ..., count - 1.
    // Store in double 'tmp', the required values will be extracted after this loop.
    double prod = 1.0; // Initialize product
    for (int i = 0; i < count; i++) {
        if (ISNAN(tpptr[positions[i]])) {
            Rf_error("TODO(R): Found missing value in tpptr, must be adressed in C");
            //return R_NaReal; 
        }
        // Updating weighted sum, simply "binmid * pdf(bin)" summed up
        // (bkptr[i + 1] + bkptr[i]) * 0.5 -> center of the bin
        res += (bkptr[i + 1] + bkptr[i]) * 0.5 * prod * (1.0 - tpptr[positions[i]]);
        // Updating product of transition probabilities
        prod *= tpptr[positions[i]];
    }
    return res;
}

///////* Calculate 'mid of bins' based on the breaks provided. */
//////doubleVec get_binmid(SEXP breaks) {
//////    int nb = LENGTH(breaks) - 1L;
//////    double* bkptr = REAL(breaks);
//////
//////    // Initialize/allocate return object
//////    doubleVec res;
//////    res.values = (double*)malloc(nb * sizeof(double));  // Allocate vector result
//////    res.length = nb;
//////
//////    // Calculating mid of bins
//////    for (int i = 0; i < nb; i++) { res.values[i] = (bkptr[i] + bkptr[i + 1]) / 2.0; }
//////
//////    return res;
//////}

/* Helper function for type = "pdf" */
doubleVec treg_calc_pdf(int* positions, int count, double* tpptr,
                        double* bkptr, double* binwidth, int nbins, int* y, int ny, bool disc) {

    // Temporary double vector to calculate PDF along i = 0, ..., count - 1
    double* tmp = malloc(count * sizeof(double)); // Single double pointer
    double tp_tmp;

    // Calculate PDF for each bin given by the distribution for i = 0, ..., count - 1.
    // Store in double 'tmp', the required values will be extracted after this loop.
    int i = 0;
    double prod = 1.0; // Initialize product
    for (i = 0; i < count; i++) {
        // Storing tpptr[positions[i]] once (used multiple times) for efficiency
        tp_tmp = tpptr[positions[i]];
        if (ISNAN(tp_tmp)) {
            Rf_error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        // Updating temporary PDF vector
        tmp[i] = prod * (1.0 - tp_tmp) / binwidth[i];
        // Updating product of transition probabilities
        prod *= tp_tmp;
    }

    // Initialize return value/object.
    doubleVec res;
    res.values = (double*)malloc(ny * sizeof(double));  // Allocate vector result
    res.length = ny;

    // Store CDF from tmp[y] to res.values
    for (int i = 0; i < ny; i++) {
        if ((y[i] < 0) || (y[i] >= nbins)) {
            res.values[i] = 0.0; // Below or above support
        } else {
            // discrete? Take pdf of the bin. Else interpolate
            // TODO(R): We could try to interpolate here
            res.values[i] = tmp[y[i]];
        }
    }

    // Free allocated memory, return result
    free(tmp);
    return res;
}

/* Helper function for type = "cdf" */
doubleVec treg_calc_cdf(int* positions, int count, double* tpptr,
                      double* bkptr,  int nbins, int* y, int ny) {

    // Temporary double vector to calculate PDF along i = 0, ..., count - 1
    double* tmp = malloc(count * sizeof(double)); // Single double pointer

    // Calculate CDF for each bin given by the distribution for i = 0, ..., count - 1
    // Store in double 'tmp', the required values will be extracted after this loop.
    tmp[0] = 1.0 - tpptr[positions[0]]; // Initialize with (1 - p[0])
    if (count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        // Looping over all elements except first
        for (int i = 1; i < count; i++) {
            if (ISNAN(tpptr[positions[i - 1]])) {
                Rf_error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
                //return R_NaReal;
            }
            pprod *= tpptr[positions[i - 1]]; // Multiply with previous element

            tmp[i] = tmp[i - 1] + (1.0 - tpptr[positions[i]]) * pprod; // Calculate cdf
        }
    }

    // Initialize return value/object.
    doubleVec res;
    res.values = (double*)malloc(ny * sizeof(double));  // Allocate vector result
    res.length = ny;

    // Store CDF from tmp[y] to res.values
    for (int i = 0; i < ny; i++) {
        if (y[i] < 0) {
            res.values[i] = 0.0;
        } else if (y[i] >= nbins) {
            res.values[i] = 1.0;
        } else {
            res.values[i] = tmp[y[i]];
        }
    }

    // Free allocated memory, return result
    free(tmp);
    return res;
}

/* Linear interpolation for finer quantiles */
double interpolate_linear(double x1, double y1, double x2, double y2, double p) {
    if (ISNAN(x1) || ISNAN(y1) || ISNAN(x2) || ISNAN(y2)) {
        return NA_REAL;
    } else if (x1 == x2) {
        return x1;
    }
    // Performing interpolation
    return x1 + (p - y1) * (x2 - x1) / (y2 - y1);
}

/* Helper function for type = "quantile" */
doubleVec treg_calc_quantile(int* positions, int count, double* tpptr,
                             double* bkptr, double* prob, int np, bool disc) {

    int i;

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
    if (res.values == NULL) { Rf_error("Memory allocation failed for doubleVec.values."); }

    if (ISNAN(tpptr[positions[0]])) {
        Rf_error("TODO(R): First element ISNAN, must be adressed in C"); //return R_NaReal; 
    }

    // Initialize with (1 - p[0])
    tmp[0] = 1.0 - tpptr[positions[0]];

    // If tmp[0] >= pmax: We have already found our solution, no need to calculate
    // the rest. Setting i = 0 (important for legacy mode).
    if (tmp[0] >= pmax) {
        i = 0; // Important
    // Else (tmp[0] < pmax) and we have count > 0 (i.e., more than only one
    // transition probability) calculate the remaining quantiles.
    } else if ( count > 0) {
        double pprod = 1.0; // Initialize with 1.0 for product
        for (i = 1; i < count; i++) {
            if (ISNAN(tpptr[positions[i - 1]])) {
                Rf_error("TODO(R): ISNAN must be implemented first (doubleVec; in cdf)");
            }

            // Update product of transition probabilities
            pprod *= tpptr[positions[i - 1]];

            tmp[i] = tmp[i - 1] + (1.0 - tpptr[positions[i]]) * pprod;

            // Break for loop as tmp[i] is already larger than pmax
            if (tmp[i] >= pmax) { break; }
        }
    }

    // Evaluate quantile
    eval_bins_quantile(res.values, tmp, positions, count, bkptr, prob, np, disc);

    free(tmp); // Freeing allocated memory

    return res;
}



void eval_bins_quantile(double* res, double* tmp, int* positions, int count,
                        double* bkptr, double* prob, int np, bool disc) {

    int i, j;
    double eps = sqrt(DBL_EPSILON);

    // Assign correct quantile to each element in res.values.
    // i: Loops over the quantiles we are looking for
    // j: Loops over calculated quantiles
    j = disc ? 1 : 0; // Discrete distributions: start at j = 1, else j = 0
    for (i = 0; i < np; i++) {
        for (; j < count; j++) {
            // Probability we are looking for too low?
            if (ISNAN(tmp[j])) { break; }

            // If prob[i] < 0 || > 1: Store NA
            if ((prob[i] < 0.0) || (prob[i] > 1.0)) {
                res[i] = R_NaReal;
                break;
            } else if (prob[i] < (tmp[0] + eps)) {
                // If prob[i] < tmp[0]: Store lowest value and break the loop
                if (disc) {
                    res[i] = disc ? (bkptr[1] + bkptr[0]) * 0.5 : bkptr[0];
                } else {
                    res[i] = interpolate_linear(bkptr[0], 0, bkptr[1], tmp[0], prob[i]);
                }
                break;
            }

            // Check if we fall into this bin
            if ((prob[i] >= (tmp[j - 1] - eps)) & (prob[i] < (tmp[j] + eps))) {
                // If discrete: Store center of the bin
                if (disc) {
                    res[i] = (bkptr[j + 1] + bkptr[j]) * 0.5;
                // Perform linear interpolation between the two neighboring bin mids.
                } else {
                    // TESTING //res[i] = ab_random(bkptr[j], bkptr[j + 1]);
                    res[i] = interpolate_linear(bkptr[j], tmp[j - 1], bkptr[j + 1], tmp[j], prob[i]);
                }
                break; // Found what we were looking for, break inner loop
            // If probability is == 1.0, take most upper break point
            } else if (prob[i] == 1.0) {
                res[i] = bkptr[count];
                break;
            }
        }

        // j == count?
        // This means we reached the end of the loop above but could not find a
        // bin we fall into. Outside support, return NA.
        if (j == count) {
            // Fill the results vector with "highest bin mid".
            res[i] = disc ? (bkptr[count - 1] + bkptr[count]) * 0.5 : R_NaReal; //bkptr[count];
        }
    }
    // void function, no return, we have updated/modified 'res'
}


/* Helper function for type = "mode"
 */
double treg_calc_mode(int* positions, int count, double* tpptr, double* bkptr) {

    // Temporary double vector to calculate PDF along i = 0, ..., count - 1
    double* tmp = malloc(count * sizeof(double)); // Single double pointer

    // If the first value is NA: return NA immediately
    if (ISNAN(tpptr[positions[0]])) { return R_NaReal; }

    // Calculate PDF for each bin given by the distribution for i = 0, ..., count - 1.
    // Store in double 'tmp', the required values will be extracted after this loop.
    int i = 0;
    double prod = 1.0; // Initialize product
    for (i = 0; i < count; i++) {
        if (ISNAN(tpptr[positions[i]])) {
            Rf_error("TODO(R): First element ISNAN, must be adressed in C");
            //return R_NaReal; 
        }
        // Updating temporary PDF vector
        tmp[i] = prod * (1.0 - tpptr[positions[i]]);
        // Updating product of transition probabilities
        prod *= tpptr[positions[i]];
    }
    // Divide PDF by width of the bin
    for (i = 0; i < count; i++) { tmp[i] = tmp[i] / (bkptr[i + 1] - bkptr[i]); }

    // Find bin with highest density
    i = 0;
    double pmax = tmp[0];
    int imode = 0;
    if (count > 0) {
        for (i = 1; i < count; i++) {
            if (tmp[i] > pmax) { pmax = tmp[i]; imode = i; }
        }
    }

    // Calculate mode; taking numeric bin midpoint
    double res = (bkptr[imode] + bkptr[imode + 1]) / 2.0;

    free(tmp);
    // Return highest value
    return res;
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
 * @param type character, either 'pdf', 'cdf', or 'mode'.
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details This function has a series of different 'modes' for different purposes.
 * Allows to calculate the CDF, PDF, quantiles, as well as the expectation (mode).
 * 
 * The three inputs which are always required are:
 *
 * * \code{uidx}: Unique 'distribution' index.
 * * \code{idx}: Integer vector specifying which transition probability.
 * * \code{lower}: Lower edge of the bin.
 * * \code{upper}: Upper edge of the bin.
 * * \code{y}: Where to evaluate the CDF, PDF, quantile. Unused for 
 *             type = "mode" and type = "mean".
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
SEXP treg_predict(SEXP uidx, SEXP idx, SEXP tp, SEXP breaks, SEXP censored, SEXP y,
                  SEXP prob, SEXP type, SEXP ncores, SEXP elementwise, SEXP discrete) {

    int    *uidxptr   = INTEGER(uidx);     // Unique indices in the dtaa
    int    *idxptr    = INTEGER(idx);      // Distribution index
    int    n          = LENGTH(idx);
    int    un         = LENGTH(uidx);
    int    i, j, np;

    // Number of threads for OMP, only used if _OPENMP support available.
    #if _OPENMP
    int    nthreads   = asInteger(ncores);
    #endif

    double* tpptr     = REAL(tp);            // Transition probabilities
    double* bkptr     = REAL(breaks);        // Breaks, points of intersection
    int     nbins     = LENGTH(breaks) - 1;  // Number of bins (3 breaks -> 2 bins -> bin "0" and bin "1"
    int*    yptr      = INTEGER(y);          // Where to evaluate the distribution; for cdf, pdf
    double* probptr   = REAL(prob);          // Probabilities to be evaluated; only for 'quantile'

    bool   ewise      = asLogical(elementwise); // Elementwise calculation?
    int*   discptr    = LOGICAL(discrete);      // Discrete distribution? For quantile

    // Evaluate 'type' to define what to do. Store a set of
    // boolean values to only do the string comparison once.
    const char* thetype = CHAR(STRING_ELT(type, 0));
    bool do_pdf  = strcmp(thetype, "pdf")  == 0;
    bool do_cdf  = strcmp(thetype, "cdf")  == 0;
    bool do_q    = strcmp(thetype, "quantile") == 0;
    bool do_mean = strcmp(thetype, "mean") == 0;
    // ... if none of them is true, it must be "do mode"
    // bool do_mode = strcmp(thetype, "modex") == 0;

    // Evalute 'censored' (str). Setting cens_left and/or cens_right to true if needed.
    const char* thecens = CHAR(STRING_ELT(censored, 0));
    bool cens_left  = (strcmp(thecens, "left")  == 0) || (strcmp(thecens, "both") == 0);
    bool cens_right = (strcmp(thecens, "right") == 0) || (strcmp(thecens, "both") == 0);

    // Calculating bin width only once (used to calculate pdf)
    double* binwidth = malloc(nbins * sizeof(double)); // Single double pointer
    for (i = 0; i < nbins; i++) { binwidth[i] = bkptr[i + 1] - bkptr[i]; }

    // Note: If censored (left or right) the last bin on these sides have a
    // width of 0.0; this width is used to calculate the PDF, where a width
    // of 0.0 would cause obvious issues. Thus, we replace the first/last
    // bin with 1.0 (so that we divide by 1; keep as is) in case needed.
    if (cens_left)  { binwidth[0]         = 1.0; }
    if (cens_right) { binwidth[nbins - 1] = 1.0; }

    // Allocating return vector.
    //
    // If type is "pdf", "cdf", or "quantile" and ewise is true: the length of
    // the vector is the same as 'un', else 'un * ny' (number of distributions
    // times number of probabilites/thresholds at which each distribution is
    // evaluated).
    SEXP res;
    if (do_pdf || do_cdf) {
        if (ewise & (LENGTH(y) != un)) { Rf_error("[C]: Length of 'y' must be equal to length of 'u'."); }
        np = (ewise) ? 1 : LENGTH(y);
        PROTECT(res = allocVector(REALSXP, un * np));
    } else if (do_q) {
        if (LENGTH(discrete) != un) { Rf_error("[C]: Length of 'discrete' must be equal to length of 'ui'."); }
        np = (ewise) ? 1 : LENGTH(prob);
        PROTECT(res = allocVector(REALSXP, un * np));
    // Else it is type = "mode", so length of res is equal to 'un'.
    } else {
        np = 1;
        PROTECT(res = allocVector(REALSXP, un));
    }

    // Pointer on results vector 'res'
    double *resptr = REAL(res);

    // If mode is not pdf, cdf, or quantile, elementwise must be true.
    // Else we throw an error here. That is for mode where elementwise
    // makes no sense.
    if (!do_pdf & !do_cdf & !do_q & !ewise) {
        Rf_error("Using \"%s\" with elementwise = false not allowed.\n", CHAR(STRING_ELT(type, 0)));
    }

    // Custom struct object to mimik "which()"
    integerVec which;
    doubleVec tmp;


    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) private(which, tmp, j)
    #endif
    for (i = 0; i < un; i++) {
        which = find_positions(uidxptr[i], idxptr, n);

        if (do_pdf || do_cdf || do_q) {
            // --- Calculating probability density
            if (do_pdf) {
                // Single PDF
                if (ewise) {
                    tmp = treg_calc_pdf(which.index, which.length, tpptr, bkptr,
                                        binwidth, nbins,
                                        &yptr[i], 1, discptr[i] == 1);
                // Multiple PDFs
                } else {
                    tmp = treg_calc_pdf(which.index, which.length, tpptr, bkptr,
                                        binwidth, nbins,
                                        yptr, LENGTH(y), discptr[i] == 1);
                }
            // --- Calculating cumulative distribution
            } else if (do_cdf) {
                // Single CDF
                if (ewise) {
                    tmp = treg_calc_cdf(which.index, which.length, tpptr, bkptr, nbins, &yptr[i], 1);
                // Multiple CDFs
                } else {
                    tmp = treg_calc_cdf(which.index, which.length, tpptr, bkptr, nbins, yptr, LENGTH(y));
                }
            } else {
                // Single quantile
                if (ewise) {
                    tmp = treg_calc_quantile(which.index, which.length, tpptr, bkptr, &probptr[i], 1, discptr[i] == 1);
                // Multiple quantiles (elementwise)
                } else {
                    tmp = treg_calc_quantile(which.index, which.length, tpptr, bkptr, probptr, LENGTH(prob), discptr[i] == 1);
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
            resptr[i] = treg_calc_mean(which.index, which.length, tpptr, bkptr);
        // Else it must be mode
        } else {
            resptr[i] = treg_calc_mode(which.index, which.length, tpptr, bkptr);
        }
        free(which.index); // Free allocated memory
    }

    /* free allocated memory */
    free(binwidth);

    UNPROTECT(1); // Releasing protected objects
    return res;
}


/* Calculating pdf and cdf (both at the same time)
 *
 * This mode is used when estimating the Transition Model. Calculates
 * the PDF and CDF for each observation at the last bin (i.e., the bin
 * in which the observation falls into). Calls treg_calc_pdf and treg_calc_cdf
 * with a missing value on 'binmidptr' and 'ny = 1' which tells treg_calc_pdf/treg_calc_cdf
 * to use this specific "mode".
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param tp probabilities, same length as idx vector.
 * @param breaks real, breaks for the different bins.
 * @param censored single character (uncensored, left, right, or both).
 * @param discrete single logical value, is this a discrete or continuous distribution?
 * @param ncores integer, number of cores to be used (ignored if OMP not available).
 *
 * @details Does something similar to treg_predict but calculates both PDF and
 * CDF simultanously for the very last bin in each distribution (observation),
 * returning a named list. This is used in the main `tm()` function,
 * calculating both at the same time should help to speed up the calculations.
 * 
 * @return Returns named list with two numeric vectors, each of which
 * has length(uidx) (vector with cdf and pdf).
 */
SEXP treg_predict_pdfcdf(SEXP uidx, SEXP idx, SEXP tp, SEXP y, SEXP breaks,
        SEXP censored, SEXP discrete, SEXP ncores) {

    double* tpptr    = REAL(tp);
    int*    uidxptr  = INTEGER(uidx);       // Unique indices in the dtaa
    int*    idxptr   = INTEGER(idx);        // Index vector
    int*    yptr     = INTEGER(y);          // Bin at which to evaluate the distribution
    double* bkptr    = REAL(breaks);        // Point intersection of bins
    int*    discptr  = LOGICAL(discrete);   // Discrete distribution?
    int     nbins    = LENGTH(breaks) - 1;  // Number of bins. 3 breaks -> 2 bins, bin "0" and bin "1"
    int     n        = LENGTH(idx);
    int     un       = LENGTH(uidx);
    int     i;

    // Evalute 'censored' (str). Setting cens_left and/or cens_right to true if needed.
    const char* thecens = CHAR(STRING_ELT(censored, 0));
    bool cens_left  = (strcmp(thecens, "left")  == 0) || (strcmp(thecens, "both") == 0);
    bool cens_right = (strcmp(thecens, "right") == 0) || (strcmp(thecens, "both") == 0);

    // Number of threads for OMP, only used if _OPENMP support available.
    #if _OPENMP
    int    nthreads   = asInteger(ncores);
    #endif

    // Calculating bin width only once (used to calculate pdf)
    double* binwidth = malloc(nbins * sizeof(double)); // Single double pointer
    for (i = 0; i < nbins; i++) { binwidth[i] = bkptr[i + 1] - bkptr[i]; }

    // Note: If censored (left or right) the last bin on these sides have a
    // width of 0.0; this width is used to calculate the PDF, where a width
    // of 0.0 would cause obvious issues. Thus, we replace the first/last
    // bin with 1.0 (so that we divide by 1; keep as is) in case needed.
    if (cens_left)  { binwidth[0]         = 1.0; }
    if (cens_right) { binwidth[nbins - 1] = 1.0; }

    // Custom struct object to mimik "which()"
    integerVec which;
    doubleVec tmppdf;
    doubleVec tmpcdf;

    // Initialize results vector
    int nProtected = 0;

    SEXP pdf; PROTECT(pdf = allocVector(REALSXP, un)); ++nProtected;
    double *pdfptr = REAL(pdf);
    SEXP cdf; PROTECT(cdf = allocVector(REALSXP, un)); ++nProtected;
    double *cdfptr = REAL(cdf);

    // Dummy value as we need a proper object when calling treg_calc_*() below
    double* na = malloc(sizeof(double)); // Single double pointer
    na[0] = NA_REAL; // Assign missing value

    /* Warning for future me: Do not use Rprintf inside omp -> segfault */
    #if _OPENMP
    #pragma omp parallel for num_threads(nthreads) private(which, tmppdf, tmpcdf)
    #endif
    for (i = 0; i < un; i++) {
        // Search position of uidxptr[i] in vector idxptr
        which = find_positions(uidxptr[i], idxptr, n);

        // Calculate pdf and cdf (for the last bin)
        //
        //     NOTE: This is the 'estimate model' mode and is only used while
        //     calling transitreg(). In this case we are only interested in the
        //     CDF and PDF of the very last value for the current distribution.
        //
        //     Important: To run this mode, 'y' must be NA (used to detect this
        //     mode) and 'ny' must be '1' (as we expect one single element in
        //     return). This is what the last two arguments to treg_calc_pdf()
        //     and treg_calc_cdf() do below [...., y = na, ny = 1)].
        //
        // Input arguments are (in this order)
        //     positions, count, tpptr, bkptr, y, ny
        tmppdf = treg_calc_pdf(which.index, which.length, tpptr,
                bkptr, binwidth, nbins, &yptr[i], 1, discptr[i] == 1);
        tmpcdf = treg_calc_cdf(which.index, which.length, tpptr, bkptr, nbins, &yptr[i], 1);

        // Store last value, that is the last bin provided for this distribution.
        pdfptr[i] = tmppdf.values[tmppdf.length - 1];
        cdfptr[i] = tmpcdf.values[tmpcdf.length - 1];

        // Free allocated memory
        free(which.index);
        free(tmppdf.values);
        free(tmpcdf.values);
    }

    /* free allocated memory */
    free(binwidth);

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


