/**
 * @file ukal.h 
 * @brief Declarations of public uKal Kalman Filter library functions.
 *
 * @author Sargis Yonan
 * @date January 22, 2019
 **/

/*
 * uKal is a discrete time implementation of the Kalman Filter.
 *
 * The definitions and declarations of the filter can be used with both linear 
 * and non-linear systems via an Extended Kalman Filter (EKF), or nth-Order 
 * Filter (nOF), like a second order Kalman filter (SOF).
 *
 * In the non-linear cases, the user must define their Jacobian (or higher order 
 * derivative) matrices, propagation function, and observation function in a
 * separate unit, and call the set f(x) and h(x) functions on their own.
 */

#ifndef _UKAL_H_
#define _UKAL_H_

#include "ukal_type.h"

#include "ulapack.h"

/**
 * @name ukal_filter_create
 * Create a Kalman filter object with all of the necessary state space system 
 * terms, initial state vector, and initial state covariance/variance matrix.
 *
 * @note The state space form of the system will abide by the following dynamics
 *       equations:
 *       The state vector, x, at filter iteration {k+1}, x_{k+1}, is equal to:
 *       Linear:
 *       x_{k+1} = \Phi_{k}*x_{k} + \gamma_{k}*w_{k}, or
 *       Nonlinear:
 *       x_{k+1} = f_{k}(x_{x}, w_{k}), where,
 *       w_{k} is a zero mean, stochastic, Wiener (Brownian motion) 
 *             process for the k^th iteration of the discrete filter.
 *             Note: Q_{k} = E[w_{k}*w_{k}^T];
 *
 *       \gamma_{k} is the k^th state equation process noise vector;
 *
 *       \Phi_{k} is the state propagation matrix. For linear systems, this
 *                matrix will likely stay constant, for non-linear systems, the 
 *                value can be modified via the set function.
 *
 *       The observation dynamics of the system will follow the following 
 *       equation:
 *       Linear:
 *       y_{k} = H_{k}x_{k} + v_{k}, or,
 *       Nonlinear:
 *       y_{k} = h_{k}(x_{k}) + v_{k}, where,
 *       y_{k} is the observation vector which will contain the kth sensor 
 *             observations;
 *
 *       H_{k} is the observation matrix which maps the state vector to the
 *             measurement states. For non-linear systems, the derivative of
 *             h_{k} can be set via its set function;
 *
 *       v_{k} is the zero mean sensor noise vector.
 *             Note: R_{k} = E[v_{k}*v{k}^T].
 *
 * @note For N states in the state vector, and M observations in the measurement 
 *       vector:
 *       size(x) = N x 1
 *       size(y) = M x 1
 *
 *       size(\Phi) = N x N
 *       size(\gamma) = N x J
 *       size(Q) = J x J
 *
 *       size(H) = M x N
 *       size(R) = M x M
 * @note For non-linear systems where:
 *                     x_{k+1} = f_{k}(x_{k}, w_{k}), and
 *                     y_{k} = h_{k}(x_{k}) + v_{k}
 *       The filter dynamics will propagate identically to the linear case, but
 *       and the observation matrix function, h, and the state propagation
 *       function, f, will not be updated accordingly. It is the responsibility
 *       of the user of this filter to handle the derivative matrices
 *       accordingly, but their values can be set via the set functions in this
 *       library.
 *
 * @note n_states must be less than or greater than UKAL_MAX_STATE_VECTOR_SIZE.
 * @note n_measurements must be less than or greater than 
 *       UKAL_MAX_MEASUREMENT_VECTOR_SIZE.
 *
 * @param[in/out] The filter to initialize, and set up.
 * @param type The filter type: linear or ekf.
 * @param n_state The number of states in the state vector, x.
 * @param n_measurements The number of observation states in the observation
 *        vector, y.
 * @param[in] Phi The state propagation matrix.
 * @param[in] gamma The state process noise matrix.
 * @param[in] x0 The initial state vector for the state space system.
 * @param[in] Q The process noise covariance matrix.
 * @param[in] P0 The initial state covariance/variance matrix.
 * @param[in] H The observation matrix.
 * @param[in] R The sensor noise covariance matrix.
 *
 * @return A FilterError status code indicating success, or failure.
 *         An invalid argument error can be returned if the dimensions
 *         of the given matrix and vector terms do not match the expected
 *         dimensions. An uninitialized object error can be returned if one
 *         or more matrix objects are not initialized, or the passed in filter
 *         object is not a valid pointer.
 */
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
                                 const Matrix_t * const R);

/*******************************
 * Nonlinear setter functions. *
 *******************************/

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
                          const Matrix_t * const fx);

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
                          const Matrix_t * const hx);

/*********************************
 * Filter propagation functions. *
 *********************************/

/**
 * @name ukal_model_predict
 * Predict the next state vector and state covariance matrix.
 *
 * @note The following equations are propagated:
 *       The predicted state vector, x(-):
 *       Linear case:
 *       x_{k+1}(-) = Phi_{k}*x_{k}
 *       Nonlinear case:
 *       x_{k+1}(-) = f(x_{k}, 0)
 *
 *       The predicted covariance matrix, P(-):
 *       P_{k+1}(-) = Phi_{k}*P_{k}*Phi_{k}^T + gamma_{k}*Q_{k}*gamma_{k}^T
 *
 * @param[in] filter The filter object associated with this model.
 *
 * @return Filter status code.
 */
FilterError_t ukal_model_predict(Filter_t * const filter);


/**
 * @name ukal_update
 * Update the state vector, covariance matrix, and Kalman gain with new 
 * measurement data.
 *
 * @note The following equations are propagated:
 *       
 *       The Kalman Gain Matrix, K:
 *       K_{k+1} = P_{k+1}(-)*H_{k+1}^T *
 *                (H_{k+1}P_{k+1}(-)H_{k+1}^T + R_{k+1})^-1
 *
 *       The update state vector, x(+):
 *       Linear case:
 *       x_{k+1}(+) = x_{k+1}(-) + K_{k+1} * 
 *                                 (y_{k+1} - H_{k+1}*x_{k+1}(-))
 *       Nonlinear case:
 *       x_{k+1}(+) = x_{k+1}(-) + K_{k+1} * 
 *                                 (y_{k+1} - h_{k+1}(x_{k+1}(-)))
 *
 *       The update covariance matrix, P(+):
 *       P_{k+1} = (I - K_{k+1} * H_{k+1}) * P_{k+1}(-)
 *
 * @param[in/out] filter The filter object associated with this model.
 * @param[in] obs_matrix The observation matrix for this filter.
 *
 * @return Filter status code.
 */
FilterError_t ukal_update(Filter_t * const filter,
                          const Matrix_t * const measurements);

/*********************
 * Getter functions. *
 *********************/

/**
 * @name ukal_get_state_cov
 * Getter for the current state covariance matrix in the filter.
 *
 * @note The size of P must match that of the covariance matrix in the filter.
 *
 * @param[in] filter The filter object to grab the covariance matrix from.
 * @param[out] P The covariance matrix to copy the filter's matrix back into.
 *
 * @return Filter status code.
 */
FilterError_t ukal_get_state_cov(const Filter_t * const filter, 
                                 Matrix_t * const P);

/**
 * @name ukal_get_state
 * Getter for the current state vector estimate.
 *
 * @note The size of x must match that of the state vector in the filter.
 *
 * @param[in] filter The filter object to grab the state vector from.
 * @param[out] x The vector to copy the filter's state back into.
 *
 * @return Filter status code.
 */
FilterError_t ukal_get_state(const Filter_t * const filter, 
                             Matrix_t * const x);

/*********************************************
 * State space system term setter functions. *
 *                                           *
 * Only used for on-the-fly term resetting.  *
 * The corresponding internal filter term    *
 * will be overwritten.                      *
 *********************************************/

FilterError_t ukal_set_prop(Filter_t * const filter, 
                            const Matrix_t * const Phi);

FilterError_t ukal_set_state(Filter_t * const filter, 
                             const Matrix_t * const x);

FilterError_t ukal_set_process_noise(Filter_t * const filter, 
                                     const Matrix_t * const gamma,
                                     const Matrix_t * const Q);

FilterError_t ukal_set_state_cov(Filter_t * const filter, 
                                 const Matrix_t * const P);

FilterError_t ukal_set_obs(Filter_t * const filter, 
                           const Matrix_t * const H);

FilterError_t ukal_set_obs_noise(Filter_t * const filter, 
                                 const Matrix_t * const R);

/*
 * End header guard.
 */
#endif