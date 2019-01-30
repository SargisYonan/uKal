# μKal
The Micro Kalman Filter Library

[![Build Status](https://travis-ci.com/SargisYonan/uKal.svg?branch=master)](https://travis-ci.com/SargisYonan/uKal)

A library for Kalman filtering, state estimation, and sensor fusion on memory constrained microncontrollers and embedded systems. The library is built on top of [μLAPack](https://www.github.com/SargisYonan/ulapack), a micro linear algebra package optimized for memory constrained systems, and can run on any target that C can be compiled for.

μKal is a full discrete-time Kalman Filtering library. The library can filter linear, and non-linear, systems via an Extended Kalman Filter (EKF) or second-order filter (SOF).

## Features
All μKal functions are "safe" in that the matrix/vector operations are checked for initialization and mathematic dimensional legality before an operation takes place. The library contains various operations and manipulations that can be configured to your needs. The library can be configured to run safely in even the most constrained environments where memory allocation is a concern.

### Included Functionality

* Create a state space model filter
	- Set the filter type: `linear`, `ekf`, or `sof`
* Predict state values (Kalman Predictions)
* Filter Sensor measurements (Kalman Updates)
* Update functions - `ukal_set_fx` and `ukal_set_hx`
	- For non-linear systems, you may set the nonlinear state vector
	  and observation dynamics
* Safe getters and setters for the state vector, covariance matrix, and more!

### Kalman Filter Example Using an EKF
See [an example using μKal for an EKF](http://www.yonan.org/ukal/examples/ekf_example.html) to learn more about modeling a system, and using μKal to filter noisy sensor measurements.

#### Static Memory Allocation
The library is implemented on top of [μLAPack](https://www.github.com/SargisYonan/ulapack), and is configured to run only with statically allocated memory. This library is safe and ready to use on an embedded system where dynamic memory allocation is not feasible.

#### Choose Your Entry Container Type
The user of this library has the option of setting the matrix/vector element data type to a desired type via the `ULAPACK_MATRIX_ENTRY_TYPE` `#define` within `ukal.h`.

#### Clear Upon Initialization
All new matrix/vector objects made can be initialized to zeros if `ULAPACK_INITIALIZE_MEMORY` is defined at compile time.

#### Matrix Inversion Options
Define `ULAPACK_INVERSE_VIA_LU` to use LU decomposition for matrix inversions.

### Unit Tests
μLAPack was developed using test-driven practices. The unit tests and `Makefile` for building and running the unit tests are in the `tests/` directory. Included is an implementation of an Extended Kalman Filter to serve as an example.

### Error Codes
The following error codes can be returned from the library functions. See the in-file documentation for the function to check its return codes of type `FilterError_t`.

* `ukal_error` -  general error code
* `ukal_invalid_argument` - bad argument given to function
* `ukal_uninit_obj` - uninitialized object passed into function
* `ukal_success` - general success code

# Licensing & Support
μKal is free to use for personal and/or educational use without profit and for development purposes in your project(s) to verify if μKal is right for you. 

For developers and for derivative commercial use (with redistribution rights):
Fill out the [licensing inquiry form](https://goo.gl/forms/8QpSDgC3JthAGoTG2), and I will promptly reply with a response and proposed license agreement.

Refer to the [LICENSE](LICENSE) document for licensing details.

# TODO
* Clean up the EKF example.
* Push all unit tests after documenting.
* Set up CI build.
* Add examples to README.

Have a suggestion for a new feature/function? File an issue in this repository with your requests.
