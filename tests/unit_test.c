#include "ulapack.h"
#include "ukal.h"

#include <stdbool.h>
#include <stdio.h>

#include <math.h>

/*
 * Total unit test error counter.
 */
static unsigned int ut_error_counter = 0;

#define INFO(msg) \
    fprintf(stderr, "Error: %s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, "%s\n", msg)

#define ERROR(err) \
    fprintf(stderr, "Error code: %d\n\n", (int)err); \

#define ut_iserr(err, msg) do { \
                            if ((FilterError_t)err != filter_success) { \
                                INFO(msg); \
                                ERROR(err); \
                                ut_error_counter++; \
                                break; \
                            } else { \
                                break; \
                            } } while(1)

#define ut_incr_error_if(condition, msg) do { \
                                    if (condition) { \
                                        INFO(msg); \
                                        ut_error_counter++; \
                                        break; \
                                    } else { \
                                        break; \
                                    } } while(1)

#define ut_incr_error_ifnot(condition, msg) do { \
                                    if (! (condition) ) { \
                                        INFO(msg); \
                                        ut_error_counter++; \
                                        break; \
                                    } else { \
                                        break; \
                                    } } while(1)

static Index_t get_sensor_data(Matrix_t * const y) {
    static const MatrixEntry_t sensor_data[45][2] = {{44.9870, 31.7870},
                                                     {43.0940, 31.8840},
                                                     {41.1820, 32.0910},
                                                     {40.1250, 32.0920},
                                                     {39.0830, 32.1510},
                                                     {37.9010, 32.0940},
                                                     {36.0060, 31.8540},
                                                     {34.2810, 31.0730},
                                                     {33.0500, 29.4340},
                                                     {32.1820, 27.2490},
                                                     {31.9400, 25.2760},
                                                     {32.5190, 23.2850},
                                                     {33.1890, 21.4630},
                                                     {33.4130, 19.3840},
                                                     {32.6440, 17.5540},
                                                     {31.4740, 15.7790},
                                                     {29.5560, 14.8140},
                                                     {27.5360, 14.2560},
                                                     {25.3800, 13.9740},
                                                     {23.4250, 13.6860},
                                                     {21.3240, 13.4880},
                                                     {19.3340, 13.5470},
                                                     {17.1990, 13.4350},
                                                     {15.0330, 13.7050},
                                                     {12.8830, 13.9040},
                                                     {11.0570, 14.9960},
                                                     { 9.5555, 16.2590},
                                                     { 8.5324, 18.0680},
                                                     { 8.3078, 20.1080},
                                                     { 8.7145, 22.2830},
                                                     { 9.7376, 24.1000},
                                                     {11.4600, 25.6850},
                                                     {13.2850, 26.5540},
                                                     {14.7140, 26.8060},
                                                     {16.0280, 26.8710},
                                                     {17.3600, 26.9790},
                                                     {18.5490, 27.2470},
                                                     {19.8400, 27.6780},
                                                     {20.9680, 28.1090},
                                                     {22.7380, 28.1870},
                                                     {24.8510, 28.8490},
                                                     {26.9790, 28.9380},
                                                     {29.0090, 29.4170},
                                                     {30.9570, 29.2450},
                                                     {32.9430, 29.5060}};
    static Index_t iteration = 0;

    ulapack_edit_entry(y, 0, 0, sensor_data[iteration][0]);
    ulapack_edit_entry(y, 1, 0, sensor_data[iteration][1]);

    iteration++;
    return iteration;
}

static void get_state_jacobian(Matrix_t * const Phi, 
                               const Matrix_t * const x, 
                               const MatrixEntry_t dt) {
    MatrixEntry_t x3 = x->entry[2][0];
    MatrixEntry_t x4 = x->entry[3][0];

    ulapack_eye(Phi);

    ulapack_edit_entry(Phi, 0, 2, dt*cos(x4));
    ulapack_edit_entry(Phi, 1, 2, dt*sin(x4));

    ulapack_edit_entry(Phi, 0, 3, -1*dt*x3*sin(x4));
    ulapack_edit_entry(Phi, 1, 3, dt*x3*cos(x4));
}

static void get_predicted_state(Matrix_t * const fx,
                                const Matrix_t * const x,
                                const MatrixEntry_t dt) {
    MatrixEntry_t x3 = x->entry[2][0];
    MatrixEntry_t x4 = x->entry[3][0];

    ulapack_edit_entry(fx, 0, 0, x->entry[0][0] + dt * x3 * cos(x4));
    ulapack_edit_entry(fx, 1, 0, x->entry[1][0] + dt * x3 * sin(x4));
    ulapack_edit_entry(fx, 2, 0, x->entry[2][0]);
    ulapack_edit_entry(fx, 3, 0, x->entry[3][0]);
}

static void get_obs_jacobian(Matrix_t * const Hx, 
                             const Matrix_t * const x) {
    ulapack_edit_entry(Hx, 0, 0, x->entry[0][0]);
    ulapack_edit_entry(Hx, 1, 0, x->entry[1][0]);
}

int main(void) {
    /*
     * Filter object.
     */
    Filter_t filter;

    const Index_t n_states = 4;
    const Index_t n_measurements = 2;

    const MatrixEntry_t dt = 0.33333;

    const MatrixEntry_t stdx = 1.3940;
    const MatrixEntry_t stdy = stdx;

    const MatrixEntry_t varx = (stdx*stdx) / 3;
    const MatrixEntry_t vary = varx;
    const MatrixEntry_t varv = 2*sqrt(2)*(varx / (dt * dt));
    const MatrixEntry_t vartheta = sqrt(11);

    /*
     * Measurement vector.
     */
    Matrix_t y;
    ulapack_init(&y, n_measurements, 1);
    get_sensor_data(&y); get_sensor_data(&y);

    /*
     * State vector.
     */
    Matrix_t x;
    ulapack_init(&x, n_states, 1);
    ulapack_edit_entry(&x, 0, 0, y.entry[0][0]);
    ulapack_edit_entry(&x, 1, 0, y.entry[1][0]);
    ulapack_edit_entry(&x, 2, 0, 5.6865);
    ulapack_edit_entry(&x, 3, 0, 3.1416);

    /*
     * The state propagation matrix.
     */
    Matrix_t Phi;
    ulapack_init(&Phi, n_states, n_states);
    get_state_jacobian(&Phi, &x, dt);

    /*
     * The state process noise.
     */
    Matrix_t gamma;
    ulapack_init(&gamma, n_states, 2);
    ulapack_set(&gamma, 0.0);
    ulapack_edit_entry(&gamma, 2, 0, 1.0);
    ulapack_edit_entry(&gamma, 3, 1, 1.0);

    /*
     * The state process noise covariance matrix.
     */
    Matrix_t Q;
    ulapack_init(&Q, gamma.n_cols, gamma.n_cols);
    ulapack_set(&Q, 0.0);
    ulapack_edit_entry(&Q, 0, 0, 5*5*dt);
    ulapack_edit_entry(&Q, 1, 1, .5*.5*dt);

    /*
     * The state covariance matrix.
     */
    Matrix_t P;
    ulapack_init(&P, n_states, n_states);
    ulapack_edit_entry(&P, 0, 0, varx);
    ulapack_edit_entry(&P, 1, 1, vary);
    ulapack_edit_entry(&P, 2, 2, varv);
    ulapack_edit_entry(&P, 3, 3, vartheta);

    ulapack_edit_entry(&P, 0, 2, 2*varx / (dt * dt));
    ulapack_edit_entry(&P, 2, 0, 2*varx / (dt * dt));
    ulapack_edit_entry(&P, 1, 2, 2*varx / (dt * dt));
    ulapack_edit_entry(&P, 2, 1, 2*varx / (dt * dt));

    ulapack_edit_entry(&P, 0, 3, 2*varx);
    ulapack_edit_entry(&P, 3, 0, 2*varx);
    ulapack_edit_entry(&P, 1, 3, 2*varx);
    ulapack_edit_entry(&P, 3, 1, 2*varx);

    Matrix_t H;
    ulapack_init(&H, n_measurements, n_states);
    ulapack_set(&H, 0.0);
    ulapack_edit_entry(&H, 0, 0, 1.0);
    ulapack_edit_entry(&H, 1, 1, 1.0);

    Matrix_t R;
    ulapack_init(&R, n_measurements, n_measurements);
    ulapack_set(&R, 0.0);
    ulapack_edit_entry(&R, 0, 0, (stdx*stdx) / 3);
    ulapack_edit_entry(&R, 1, 1, (stdy*stdy) / 3);

    Matrix_t fx;
    ulapack_init(&fx, n_states, 1);
    get_predicted_state(&fx, &x, dt);

    Matrix_t Hx;
    ulapack_init(&Hx, n_measurements, 1);
    get_obs_jacobian(&Hx, &x);

    ut_iserr(ukal_filter_create(&filter, ekf, n_states, n_measurements,
                       &Phi, &gamma, &x, &Q,
                       &P,
                       &H, &R), "Cannot create filter.");

    for (Index_t itor = 0; itor < 45 - 2; itor++) {

        printf("k = %d\n", (int)itor);
        printf("*************************************\n");
        get_sensor_data(&y);

        get_state_jacobian(&Phi, &filter.x, dt);
        ut_iserr(ukal_set_phi(&filter, &Phi), "Cannot set prop matrix.");

        printf("Phi(-) = \n");
        ulapack_print(&filter.Phi, stdout);

        get_predicted_state(&fx, &filter.x, dt);
        ut_iserr(ukal_set_fx(&filter, &fx), "Cannot update f(x).");
        printf("f(x) = \n");
        ulapack_print(&filter.fx, stdout);

        ut_iserr(ukal_model_predict(&filter), "Cannot predict model.");
        printf("x(-) = \n");
        ulapack_print(&filter.x, stdout);
        printf("P(-) = \n");
        ulapack_print(&filter.P, stdout);

        get_obs_jacobian(&Hx, &filter.x);
        ut_iserr(ukal_set_hx(&filter, &Hx), "Cannot set h(x).");
        printf("h(x) = \n");
        ulapack_print(&filter.Hx, stdout);

        ut_iserr(ukal_update(&filter, &y), "Cannot update filter.");
        printf("K = \n");
        ulapack_print(&filter.K, stdout);
        printf("P(+) = \n");
        ulapack_print(&filter.P, stdout);
        printf("x(+) = \n");
        ulapack_print(&filter.x, stdout);
    }

    Matrix_t x_final_exp;
    ulapack_init(&x_final_exp, n_states, 1);
    /*
     * Parallel filter ran in MATLAB.
     */
    ulapack_edit_entry(&x_final_exp, 0, 0, 32.950418372726212);
    ulapack_edit_entry(&x_final_exp, 1, 0, 29.514796929099575);
    ulapack_edit_entry(&x_final_exp, 2, 0, 5.956417779623743);
    ulapack_edit_entry(&x_final_exp, 3, 0, 0.071855231860152);

    MatrixError_t isequal;

    ut_iserr (ulapack_is_equal(&filter.x, &x_final_exp, &isequal), "Cannot compare expected and actual states.");
    ut_iserr ( isequal, "Expected and actual state vectors do not match." );

    Matrix_t P_final_exp;
    ulapack_init(&P_final_exp, n_states, n_states);
    /*
     * Parallel filter ran in MATLAB.
     */
    ulapack_edit_entry(&P_final_exp, 0, 0, 0.517846656052970);
    ulapack_edit_entry(&P_final_exp, 1, 0, 0.005544392784906);
    ulapack_edit_entry(&P_final_exp, 2, 0, 1.033243899767675);
    ulapack_edit_entry(&P_final_exp, 3, 0, -0.012176599795471);

    ulapack_edit_entry(&P_final_exp, 0, 1, 0.005544392784906);
    ulapack_edit_entry(&P_final_exp, 1, 1, 0.457839474152853);
    ulapack_edit_entry(&P_final_exp, 2, 1, 0.101416749649471);
    ulapack_edit_entry(&P_final_exp, 3, 1, 0.122649098758324);

    ulapack_edit_entry(&P_final_exp, 0, 2,  1.033243899767674);
    ulapack_edit_entry(&P_final_exp, 1, 2,  0.101416749649471);
    ulapack_edit_entry(&P_final_exp, 2, 2, 12.491559373945359);
    ulapack_edit_entry(&P_final_exp, 3, 2,  0.004225656804946);

    ulapack_edit_entry(&P_final_exp, 0, 3, -0.012176599795471);
    ulapack_edit_entry(&P_final_exp, 1, 3,  0.122649098758324);
    ulapack_edit_entry(&P_final_exp, 2, 3,  0.004225656804946);
    ulapack_edit_entry(&P_final_exp, 3, 3,  0.151279207167239);

    ut_iserr (ulapack_is_equal(&filter.P, &P_final_exp, &isequal), "Cannot compare expected and cov mats.");
    ut_iserr ( isequal, "Expected and actual cov matrix do not match." );
        
    printf("Errors: %u\n", ut_error_counter);

    if (ut_error_counter > 0) {
        return -1 * ut_error_counter;
    }

    return 0;
}