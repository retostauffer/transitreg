
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdbool.h>

#include <stdlib.h>

/* Helper function to find position of x in y
 *
 * @param x single integer (the 'needle').
 * @param y integer vector (the 'haystack').
 * @param n length of y.
 * @param count number of elements found (modified inside function).
 * 
 * @return returns an integer vector (must be freed outside function)
 * with the index position of x in y.
 */
int* find_positions(int x, int* y, int n, int* count) {
    int* positions = (int*)malloc(n * sizeof(int)); // Allocate max possible size
    *count = 0;
    for (int i = 0; i < n; i++) {
        if (y[i] == x) {
            positions[*count] = i; // Store the index
            (*count)++;
        }
    }
    return positions; // Caller must free this memory
}

/* Helper function for type = "pdf".
 *
 * Calculates (1 - p[count]) * prod(p[-count]) with p = pptr[positions]
 */
double c_calc_pdf(int* positions, int count, double* pptr) {
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
double c_calc_cdf(int* positions, int count, double* pptr) {
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

/* Calculating elementwise pdf
 *
 * @param uidx integer vector with unique indices in data.
 * @param idx integer with indices, length of idx is sample size times breaks.
 * @param p probabilities, same length as idx vector.
 * @param type_pdf logical, if TRUE the PDF is calculated, else the CDF.
 *
 * @details Internally loops over all unique indices in `uidx`.
 * For each index (belonging to one observation) we check the position
 * of the elements in `idx`, and call the `predict_pdf_calc` function
 * which calculates and returns the pdf (double).
 * 
 * @return Returns SEXP double vector of length (length(uidx)).
 */
SEXP c_predict_pdf_cdf(SEXP uidx, SEXP idx, SEXP p, SEXP type_pdf) {

    int    do_pdf = LOGICAL(type_pdf)[0];
    double *pptr    = REAL(p);
    int    *uidxptr = INTEGER(uidx);  // Unique indices in the dtaa
    int    *idxptr  = INTEGER(idx);   // Index vector
    int    n = LENGTH(idx);
    int    un = LENGTH(uidx);
    int    i, count;

    // Initialize results vector
    SEXP probs;
    PROTECT(probs = allocVector(REALSXP, un));
    double *probsptr = REAL(probs);

    for (i = 0; i < un; i++) {
        int* positions = find_positions(uidxptr[i], idxptr, n, &count);
        if (do_pdf) {
            probsptr[i] = c_calc_pdf(positions, count, pptr);
        } else {
            probsptr[i] = c_calc_cdf(positions, count, pptr);
        }
        free(positions); // Free allocated memory
    }

    UNPROTECT(1); // Releasing protected objects
    return probs;
}


