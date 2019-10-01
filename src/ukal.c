/**
 * @file ukal.c
 * @brief Definitions of public uKal Kalman Filter library functions.
 *
 * @author Sargis Yonan
 * @date January 22, 2019
 **/

#include "ukal.h"

/**
 * @macro error_if_uninit
 * Helper macro for returning uninitialized object error.
 *
 * @param f A pointer to a filter or matrix object.
 *
 * @return Return uninitialized object error if the object is NULL.
 */
#define error_if_uninit(f)  do { \
                                if (!(f)) { \
                                    return filter_uninit_obj; \
                                } else { \
                                    break; \
                                } \
                            } while(1)

/**
 * @macro ret_iferr
 * Helper macro for returning error code if the code is not a success.
 *
 * @param err An error code.
 *
 * @return The error code passed in if the value passed in is not a success.
 */
#define ret_iferr( err )  do { \
                            if ((FilterError_t)(err) != filter_success) { \
                                return (FilterError_t)err; \
                            } else { \
                                break; \
                            } \
                          } while(1)

FilterError_t ukal_filter_create(Filter_t * const filter,
                                 const FilterType_t type,
                                 const Index_t n_states,
                                 const Index_t n_measurements,
                                 const Matrix_t * const Phi,
                                 const Matrix_t * const gamma,
                                 const Matrix_t * const x0,
                                 const Matrix_t * const Q,
                                 const Matrix_t * const P0,
                                 const Matrix_t * const H,
                                 const Matrix_t * const R) {
    /*
     * Check for invalid filter pointer.
     */
    error_if_uninit(filter);

    /*
     * The passed in number of states can not exceed the pre-allocated maximum
     * matrix dimensions for the uLAPack library.
     */
    if (n_states > UKAL_MAX_STATE_VECTOR_SIZE) {
        return filter_invalid_argument;
    }

    /*
     * The passed in number of observations can not exceed the pre-allocated 
     * maximum matrix dimensions for the uLAPack library.
     */
    if (n_measurements > UKAL_MAX_MEASUREMENT_VECTOR_SIZE) {
        return filter_invalid_argument;
    }

    /*
     * Check for invalid filter type.
     */
    switch (type) {
        /*
         * Set the filter type if it is a valid option.
         */
        case linear:
        case ekf:
        case sof:
            filter->type = type;
            break;

        default:
            return filter_invalid_argument;
    }

    /*
     * Set the number of states and the number of measurements in the filter.
     */
    filter->n_states = n_states;
    filter->n_measurements = n_measurements;

    /*
     * Initialize the filter's inner state space terms with the now known 
     * dimensions.
     */

    /*
     * Initialize x to be N x 1.
     */
    ret_iferr( ulapack_init(&filter->x, n_states, 1) );

    /*
     * Initialize the nonlinear state function if a nonlinear filter
     * is specified.
     */
    if (type != linear) {
        ret_iferr( ulapack_init(&filter->fx, n_states, 1) );
    }

    /*
     * Initialize Phi and Phi^T to be N x N.
     */
    ret_iferr( ulapack_init(&filter->Phi, n_states, n_states) );
    ret_iferr( ulapack_init(&filter->PhiT, n_states, n_states) );

    /*
     * Initialize gamma to N x J and gamma^T to J x N. J is not known until
     * we receive the gamma object.
     */
    if (gamma->n_rows != n_states) {
        return filter_invalid_argument;
    }
    ret_iferr( ulapack_init(&filter->gammaQgammaT, gamma->n_rows, gamma->n_rows) );

    /*
     * Initialize the covariance matrix to N x N.
     */
    ret_iferr( ulapack_init(&filter->P, n_states, n_states) );
    
    /*
     * Initialize the observation matrix and its transpose to M x N and N x M
     * respectively.
     */
    ret_iferr( ulapack_init(&filter->H, n_measurements, n_states) );
    ret_iferr( ulapack_init(&filter->HT, n_states, n_measurements) );

    /*
     * Preinitialize a place holder variable for H*x or for h(x) for
     * the nonlinear filters.
     */
    ret_iferr( ulapack_init(&filter->Hx, n_measurements, 1) );

    /*
     * Initialize the innovation term.
     */
    ret_iferr( ulapack_init(&filter->innovation, filter->n_measurements, 1) );

    /*
     * Initialize the sensor noise covariance matrix to M x M.
     */
    ret_iferr( ulapack_init(&filter->R, n_measurements, n_measurements) );

    /*
     * Initialize the Kalman gain matrix to N x M.
     */
    ret_iferr( ulapack_init(&filter->K, n_states, n_measurements) );

    /*
     * Initialize the identity matrix to N x N, and set it to I_{NxN}.
     */
    ret_iferr( ulapack_init(&filter->eye, n_states, n_states) );
    ret_iferr( ulapack_eye(&filter->eye) );

    /*
     * Set the internal terms to the given terms. Additionally set the transpose
     * of the terms for calculation efficiency down the road.
     */

    /*
     * Set the internal state propagation matrix to the given matrix.
     * Set the transpose of the internal state propagation matrix to the
     * transpose of the given matrix.
     */
    ret_iferr( ukal_set_phi(filter, Phi) );

    /*
     * Set the internal state propagation matrix to the given initial state.
     */
    ret_iferr( ukal_set_state(filter, x0) );

    /*
     * Set the internal state process noise matrix gamma. Set the internal
     * transpose to the transpose of the given gamma. Also set the process
     * noise covariance matrix, Q. The term gamma*Q*gamma^T is constant so
     * it will only be calculated once, and a direct getter for Q and gamma
     * will therefore not be accessible.
     */
    ret_iferr( ukal_set_process_noise(filter, gamma, Q) );

    /*
     * Set the internal state covariance/variance matrix to the given value.
     */
    ret_iferr( ukal_set_state_cov(filter, P0) );

    /*
     * Set the observation matrix to the given term.
     */
    ret_iferr( ukal_set_obs(filter, H) );

    /*
     * Set the sensor noise covariance matrix to the given term.
     */
    ret_iferr( ukal_set_obs_noise(filter, R) );

    return filter_success;
}

/**
 * @name ukal_set_fx
 * Set the value of f(x), the nonlinear observation dynamics matrix.
 *
 * @note Only call when a nonlinear filter type is specified.
 *
 * @param[in/out] filter The filter object to set f(x) to.
 * @param[in] fx The f(x) matrix.
 *
 * @return Filter status code.
 */
FilterError_t ukal_set_fx(Filter_t * const filter,
                          const Matrix_t * const fx) {
    error_if_uninit(filter);

    /*
     * This function is only intended for use with nonlinear filter types.
     */
    if (filter->type == linear) {
        return filter_error;
    }

    error_if_uninit(fx);

    /*
     * Check for illegal dimensions.
     */
    if (fx->n_rows != filter->n_states ||
        fx->n_cols != 1) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal term to the given term.
     */
    ret_iferr( ulapack_copy(fx, &filter->fx) );
    
    return filter_success; 
}


/**
 * @name ukal_set_hx
 * Set the value of h(x), the nonlinear observation dynamics matrix.
 *
 * @note Only call when a nonlinear filter type is specified.
 *
 * @param[in/out] filter The filter object to set h(x) to.
 * @param[in] hx The h(x) matrix.
 *
 * @return Filter status code.
 */
FilterError_t ukal_set_hx(Filter_t * const filter,
                          const Matrix_t * const hx) {
    error_if_uninit(filter);

    /*
     * This function is only intended for use with nonlinear filter types.
     */
    if (filter->type == linear) {
        return filter_error;
    }

    error_if_uninit(hx);

    /*
     * Check for illegal dimensions.
     */
    if (hx->n_rows != filter->n_measurements ||
        hx->n_cols != 1) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal term to the given term.
     */
    ret_iferr( ulapack_copy(hx, &filter->Hx) );
    
    return filter_success; 
}

FilterError_t ukal_set_phi(Filter_t * const filter, 
                            const Matrix_t * const Phi) {
    error_if_uninit(filter);
    error_if_uninit(Phi);

    /*
     * Check for illegal dimensions.
     */
    if (Phi->n_rows != filter->n_states ||
        Phi->n_cols != filter->n_states) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_transpose(Phi, &filter->PhiT) );
    ret_iferr( ulapack_copy(Phi, &filter->Phi) );   

    return filter_success;
}

FilterError_t ukal_set_state(Filter_t * const filter, 
                             const Matrix_t * const x) {
    error_if_uninit(filter);
    error_if_uninit(x);

    /*
     * Check for illegal dimensions.
     */
    if (x->n_rows != filter->n_states ||
        x->n_cols != 1) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_copy(x, &filter->x) );

    return filter_success;
}

FilterError_t ukal_set_process_noise(Filter_t * const filter, 
                                     const Matrix_t * const gamma,
                                     const Matrix_t * const Q) {
    error_if_uninit(filter);
    error_if_uninit(gamma);
    error_if_uninit(Q);

    /*
     * Temporary variables for matrix product to calculate the constant
     * term: gamma*Q*gamma^T.
     */
    Matrix_t gammaQ;
    Matrix_t gammaT;

    /*
     * Check for illegal dimensions.
     */
    if (gamma->n_rows != filter->n_states) {
        return filter_invalid_argument;
    }

    if (Q->n_rows != gamma->n_cols || 
        Q->n_cols != gamma->n_cols) {
        return filter_invalid_argument;
    }

    /*
     * Calculate gamma^T.
     */
    ret_iferr( ulapack_transpose(gamma, &gammaT) );
    /*
     * Calculate gamma * Q.
     */
    ret_iferr( ulapack_product(gamma, Q, &gammaQ) );

    /*
     * Calculate and internally store the constant term:
     * gamma * Q * gamma^T.
     */
    ret_iferr( ulapack_product(&gammaQ, &gammaT, &filter->gammaQgammaT) );

    return filter_success;
}

FilterError_t ukal_set_state_cov(Filter_t * const filter, 
                                 const Matrix_t * const P) {
    error_if_uninit(filter);
    error_if_uninit(P);

    /*
     * Check for illegal dimensions.
     */
    if (P->n_rows != filter->n_states ||
        P->n_cols != filter->n_states) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_copy(P, &filter->P) );

    return filter_success;
}

FilterError_t ukal_set_obs(Filter_t * const filter, 
                           const Matrix_t * const H) {
    error_if_uninit(filter);
    error_if_uninit(H);

    /*
     * Check for illegal dimensions.
     */
    if (H->n_rows != filter->n_measurements ||
        H->n_cols != filter->n_states) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_transpose(H, &filter->HT) );
    ret_iferr( ulapack_copy(H, &filter->H) );   

    return filter_success;
}

FilterError_t ukal_set_obs_noise(Filter_t * const filter, 
                                 const Matrix_t * const R) {
    error_if_uninit(filter);
    error_if_uninit(R);

    /*
     * Check for illegal dimensions.
     */
    if (R->n_rows != filter->n_measurements ||
        R->n_cols != filter->n_measurements) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_copy(R, &filter->R) );   

    return filter_success;
}

FilterError_t ukal_get_state_cov(const Filter_t * const filter, 
                                 Matrix_t * const P) {
    error_if_uninit(filter);
    error_if_uninit(P);

    /*
     * Check for illegal dimensions.
     */
    if (P->n_rows != filter->n_states ||
        P->n_cols != filter->n_states) {
        return filter_invalid_argument;
    }

    /*
     * Set the internal terms to the given terms.
     */
    ret_iferr( ulapack_copy(&filter->P, P) );

    return filter_success;
}

FilterError_t ukal_get_state(const Filter_t * const filter, 
                             Matrix_t * const x) {
    error_if_uninit(filter);
    error_if_uninit(x);
    
    /*
     * Check for illegal dimensions.
     */
    if (x->n_rows != filter->n_states ||
        x->n_cols != 1) {
        return filter_invalid_argument;
    }

    /*
     * Copy back the term requested.
     */
    ret_iferr( ulapack_copy(&filter->x, x) );

    return filter_success;
}

FilterError_t ukal_model_predict(Filter_t * const filter) {
    error_if_uninit(filter);

    /*
     * Temporary variable for Phi * x calculation.
     */
    Matrix_t Phix;

    /*
     * Temporary variables for Phi * P * Phi^T calculation.
     */
    Matrix_t PhiP;
    Matrix_t PhiPPhiT;
    
    ret_iferr( ulapack_init(&PhiP, filter->n_states, filter->n_states) );
    ret_iferr( ulapack_init(&PhiPPhiT, filter->n_states, filter->n_states) );

    /*
     * Calculate Phi * P.
     */
    ret_iferr( ulapack_product(&filter->Phi, &filter->P, &PhiP) );

    /*
     * Calculate Phi * P * Phi^T.
     */
    ret_iferr( ulapack_product(&PhiP, &filter->PhiT, &PhiPPhiT) );

    /*
     * Calculate and store Phi*P*Phi^T + gamma*Q*gamma^T.
     */
    ret_iferr( ulapack_add(&PhiPPhiT, &filter->gammaQgammaT, &filter->P) );

    switch(filter->type) {
        case linear:
            /*
             * Initialize the predicted state vector to size N x 1.
             */
            ret_iferr( ulapack_init(&Phix, filter->n_states, 1) );

            /*
             * Compute the total covariance matrix.
             */
            ret_iferr( ulapack_product(&filter->Phi, &filter->x, &Phix) );

            /*
             * Copy Phi * x into the internal state vector for the predicted 
             * state.
             */
            ret_iferr( ulapack_copy(&Phix, &filter->x) );

            break;

        case ekf:
        case sof:
            ret_iferr( ulapack_copy(&filter->fx, &filter->x) );

            break;

        default:
            /*
             * Unknown type.
             */
            return filter_error;
    }

    return filter_success;
}

/**
 * @name _ukal_update_kalman_gain
 * A static private helper function for updating the internal gain matrix.
 *
 * @note The filter object should have the most up to date predicted covariance
 *       matrix, P(-), at the point of calling this function.
 *
 * @param[in/out] filter The filter object to update the gain matrix.
 *
 * @return Filter status code.
 */
static FilterError_t _ukal_update_kalman_gain(Filter_t * const filter) {
    /*
     * Declare temporary variables for Kalman gain calculation.
     */
    Matrix_t PHT; /* P_{k+1}(-) * H_{k+1}^T */
    Matrix_t HPHT; /* H_{k+1} * P_{k+1}(-) * H_{k+1}^T */
    Matrix_t HPHTpR; /* H_{k+1} * P_{k+1}(-) * H_{k+1}^T + R_{k+1} */
    Matrix_t HPHTpR_inv; /* (H_{k+1} * P_{k+1}(-) * H_{k+1}^T + R_{k+1})^-1 */

    ret_iferr( ulapack_init(&PHT, filter->n_states, filter->n_measurements) );
    ret_iferr( ulapack_init(&HPHT, filter->n_measurements, filter->n_measurements) );
    ret_iferr( ulapack_init(&HPHTpR, filter->n_measurements, filter->n_measurements) );
    ret_iferr( ulapack_init(&HPHTpR_inv, filter->n_measurements, filter->n_measurements) );

    /*
     * Compute P_{k+1}(-) * H_{k+1}^T.
     */
    ret_iferr( ulapack_product(&filter->P, &filter->HT, &PHT) );
    /*
     * Compute H_{k+1} * P_{k+1}(-) * H_{k+1}^T.
     */
    ret_iferr( ulapack_product(&filter->H, &PHT, &HPHT) );

    /*
     * Compute H_{k+1} * P_{k+1}(-) * H_{k+1}^T + R_{k+1}.
     */
    ret_iferr( ulapack_add(&HPHT, &filter->R, &HPHTpR) );

    /*
     * Compute inv(H_{k+1} * P_{k+1}(-) * H_{k+1}^T + R_{k+1}).
     */
    ret_iferr( ulapack_inverse(&HPHTpR, &HPHTpR_inv) );

    /*
     * Compute the remainder of the Kalman gain matrix, and store it.
     */
    ret_iferr( ulapack_set(&filter->K, 0.0) );
    ret_iferr( ulapack_product(&PHT, &HPHTpR_inv, &filter->K) );

    return filter_success;
}

/**
 * @name _ukal_update_state_vector
 * A static private helper function for updating the internal state vector.
 *
 * @note The filter object should have the most up to date Kalman gain matrix.
 *
 * @param[in/out] filter The filter object to update the state vector within.
 * @param[in/out] y The most up to date measurement vector.
 *
 * @return Filter status code.
 */
static FilterError_t _ukal_update_state_vector(Filter_t * const filter,
                                               const Matrix_t * const y) {
    /*
     * Temporary variable for the weighed innovation term.
     */
    Matrix_t Kinn;

    /*
     * Compute the innovation term: (y - Hx) or (y - h(x)).
     * Assume H*x or h(x) is known before this point.
     */
    ret_iferr( ulapack_subtract(y, &filter->Hx, &filter->innovation) );

    /*
     * Weigh the measured and predicted observations via the Kalman gain matrix.
     * It is assumed that the gain matrix is known at this point.
     */
    ret_iferr( ulapack_init(&Kinn, filter->n_states, 1) );
    ret_iferr( ulapack_product(&filter->K, &filter->innovation, &Kinn) );

    /*
     * Add the predicted state vector to the weighted innovation.
     */
    ret_iferr( ulapack_add(&filter->x, &Kinn, &filter->x) );

    return filter_success;
}

/**
 * @name _ukal_update_cov
 * A static private helper function for updating the internal covariance matrix.
 *
 * @note The filter object should have the most up to date Kalman gain matrix.
 *
 * @param[in/out] filter The filter object to update the state vector within.
 *
 * @return Filter status code.
 */
static FilterError_t _ukal_update_cov(Filter_t * const filter) {
    /*
     * Temporary variable for K*H and (I - KH) * P.
     */
    Matrix_t temp;

    /*
     * Temporary variable for I - KH.
     */
    Matrix_t ImKH;

    /*
     * Initialize dimensions of temporary objects.
     */
    ret_iferr( ulapack_init(&temp, filter->n_states, filter->n_states) );
    ret_iferr( ulapack_init(&ImKH, filter->n_states, filter->n_states) );

    /*
     * Compute K*H.
     */
    ret_iferr( ulapack_product(&filter->K, &filter->H, &temp) );

    /*
     * Compute I - K*H.
     */
    ret_iferr( ulapack_subtract(&filter->eye, &temp, &ImKH) );

    /*
     * Calculate the covariance matrix.
     */
    ulapack_set(&temp, 0.0);
    ret_iferr( ulapack_product(&ImKH, &filter->P, &temp) );

    /*
     * Store the result in the internal filter.
     */
    ret_iferr( ulapack_copy(&temp, &filter->P) );

    return filter_success;
}

FilterError_t ukal_update(Filter_t * const filter,
                          const Matrix_t * const measurements) {

    error_if_uninit(filter);
    error_if_uninit(measurements);

    /*
     * Compute the new Kalman gain matrix.
     */
    ret_iferr( _ukal_update_kalman_gain(filter) );

    /*
     * For a linear filter, the product of H and x must be precomputed, as h(x) 
     * is assumed to be calculated before a call to this function.
     */
    if (filter->type == linear) {
        ret_iferr( ulapack_product(&filter->Hx, &filter->H, &filter->x) );
    }

    /*
     * Compute the new state vector estimate.
     */
    ret_iferr( _ukal_update_state_vector(filter, measurements) );

    /*
     * Compute the updated covariance matrix.
     */
    ret_iferr( _ukal_update_cov(filter) );

    return filter_success;
}
