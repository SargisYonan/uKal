/**
 * @file ukal_type.h 
 * @brief Definitions of uKal data types.
 *
 * @author Sargis Yonan
 * @date January 22, 2019
 **/

#ifndef _UKAL_TYPE_H_
#define _UKAL_TYPE_H_

#include "ulapack.h"

/*
 * Definition of error types.
 */
typedef enum {
    /*@{*/
    filter_success = 1, /**< general success code. */
    filter_error = -1, /**< general error code. */
    filter_invalid_argument = -3, /**< bad argument given to function. */
    filter_uninit_obj = -4, /**< uninitialized object passed into function. */
    /*@}*/
} FilterError_t;

/*
 * Definitions for supported filter types.
 */
typedef enum {
    /*@{*/
    linear = 1, /**< Linear system case. */
    ekf = 2,  /**< Nonlinear system (EKF) case. */
    sof = 3, /**< Nonlinear second order filter (SOF) case. */
    /*@}*/
} FilterType_t;

/*
 * Definition of the statically allocated filter object type.
 */
typedef struct {
    /*@{*/
    Index_t n_measurements; /** < The number of states in the measurement model. */
    Index_t n_states; /** < The number of states in each model. */

    FilterType_t type; /** < Stores the filter type. */

    Matrix_t P; /** < The covariance matrix. */
    Matrix_t K; /** < The Kalman gain matrix. */
    Matrix_t eye; /** < Identity matrix for this filter. */

    /*
     * The model equation terms. 
     */
    Matrix_t x; /** < The state vector. */
    Matrix_t fx; /** < The nonlinear Jacobian result for f(x). */
    Matrix_t Phi; /** < The propagation matrix. */
    Matrix_t PhiT; /** < The propagation matrix transposed. */
    Matrix_t gammaQgammaT; /** < The product of the process noise vector,
                                 the process noise covariance matrix, and
                                 the transpose of the process noise. */
    /*
     * Measurement model.
     */
    Matrix_t H; /** < The measurement model matrix. */
    Matrix_t Hx; /** < The calculation place holder for H*x or h(x). */
    Matrix_t HT; /** < The transpose of the measurement model matrix. */
    Matrix_t innovation; /** < The innovation term: (y - Hx) or (y - h(x)). */
    Matrix_t R; /** < The sensor noise covariance. */
    /*@}*/
} Filter_t;

/*
 * uLAPack will be used as the backend matrix math library for the filter.
 * The library will be used in the static memory allocation mode.
 * The preprocessor options for the matrix library will be set before the 
 * library header is included.
 */
#ifndef ULAPACK_MATRIX_ENTRY_TYPE
#define ULAPACK_MATRIX_ENTRY_TYPE double
#endif

#ifndef ULAPACK_INDEX_TYPE
#define ULAPACK_INDEX_TYPE uint64_t
#endif

#ifndef ULAPACK_SIGNED_INDEX_TYPE
#define ULAPACK_SIGNED_INDEX_TYPE int64_t
#endif

#ifndef ULAPACK_USE_STATIC_ALLOC
#define ULAPACK_USE_STATIC_ALLOC
#endif

#if UKAL_MAX_MEASUREMENT_VECTOR_SIZE > UKAL_MAX_STATE_VECTOR_SIZE
  #ifndef ULAPACK_MAX_MATRIX_N_ROWS
  #define ULAPACK_MAX_MATRIX_N_ROWS UKAL_MAX_MEASUREMENT_VECTOR_SIZE
  #endif

  #ifndef ULAPACK_MAX_MATRIX_N_COLS
  #define ULAPACK_MAX_MATRIX_N_COLS UKAL_MAX_MEASUREMENT_VECTOR_SIZE
  #endif
#else
  #ifndef ULAPACK_MAX_MATRIX_N_ROWS
  #define ULAPACK_MAX_MATRIX_N_ROWS UKAL_MAX_STATE_VECTOR_SIZE
  #endif

  #ifndef ULAPACK_MAX_MATRIX_N_COLS
  #define ULAPACK_MAX_MATRIX_N_COLS UKAL_MAX_STATE_VECTOR_SIZE
  #endif
#endif

#ifndef ULAPACK_INITIALIZE_MEMORY
#define ULAPACK_INITIALIZE_MEMORY
#endif

#ifndef ULAPACK_INVERSE_VIA_LU
#define ULAPACK_INVERSE_VIA_LU
#endif

/*
 * The maximum number of states any one state vector will have.
 */
#ifndef UKAL_MAX_STATE_VECTOR_SIZE
    #error "Define the maximum number of possible states in any one uKal filter."
#endif

#ifndef UKAL_MAX_MEASUREMENT_VECTOR_SIZE
    #error "Define the maximum number of possible measurement states in any one uKal filter."
#endif

/*
 * End of header guard.
 */
#endif