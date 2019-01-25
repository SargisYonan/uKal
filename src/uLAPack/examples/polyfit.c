/**
 * @name polyfit.c
 * An example of using uLAPAck for polynomial regression.
 *
 * @note In this example, static memory allocation is used
 *       via the -D compiler options.
 * @note Build this example using: $ make polyfit
 * @note Run this example using: $ ./polyfit
 */

#include "ulapack.h"

#include <stdio.h>

/**
 * @name run_regression
 * A helper function to regress some data.
 *
 * @param[in] xdata The x coordinate data.
 * @param[in] ydata The y coordinate data.
 * @param The number of x and y points.
 * @param nth_degree The number of degrees to fit a polynomial to.
 * @param[out] polynomial_coefs The returning polynomial fit coefficients.
 */
static void run_regression(const double * xdata, 
						   const double * ydata, 
						   const uint64_t data_points, 
						   const uint64_t nth_degree,
						   Matrix_t * const polynomial_coefs) {
	/*
	 * Declare matrix objects.
	 */
	Matrix_t x;
    Matrix_t y;

    /*
     * Initialize matrix objects.
     */
    ulapack_init(&x, data_points, 1);
    ulapack_init(&y, data_points, 1);

    /*
     * Copy data points into vector objects.
     */
    for (uint64_t row_itor = 0; row_itor < data_points; row_itor++) {
        ulapack_edit_entry(&x, 
            row_itor, 0, 
            xdata[row_itor]);

        ulapack_edit_entry(&y, 
            row_itor, 0, 
            ydata[row_itor]);
    }

    /*
     * Run the regression.
     */
    ulapack_polyfit(&x, &y, nth_degree, polynomial_coefs);
}

int main(void) {

	double xdata[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	/*
	 * x^2 + noise
	 */
	double ydata[] = {0.679196, 3.215585, 
					  8.635037, 16.117271, 
					  25.174340, 35.784344, 
					  48.847389, 64.033688, 
					  81.458282, 101.281631};

	const uint64_t degree = 2;
	Matrix_t p;

    ulapack_init(&p, degree + 1, 1);

	run_regression(xdata, ydata, 10, 2, &p);

	printf("\nFit coefficients: \n");
	ulapack_print(&p, stdout);

}

/*
 * Corresponding MATLAB code for comparison.
 */

/*

clear all;
close all;
clc;

x = [1:1:10]';

%% Noisy x^2
y = [0.679196;
    3.215585;
    8.635037; 
    16.117271; 
    25.174340; 
    35.784344; 
    48.847389; 
    64.033688; 
    81.458282,
    101.281631];

p = polyfit(x, y, 2)';

ulap_p = [1.02164869318180340229673674912191927433013916015625;
         -0.09319982499999923675204627215862274169921875;
         -0.2981993499998480956492130644619464874267578125];

error = norm(p - ulap_p);

figure()
plot(x, y, 'o');
hold on;
plot(linspace(1, 10), polyval(p, linspace(1,10)))
hold on;
plot(linspace(1, 10), polyval(ulap_p, linspace(1,10)))

leg = legend('Actual Data', 'MATLAB', '$\mu$LAPack polyfit');
set(leg, 'FontSize', 12, 'Interpreter', 'Latex');

title('Comparing MATLAB and $\mu$LAPack Polynomial Fitting (Noisy $x^2$)', ...
'Interpreter', 'Latex', 'FontSize', 16)
xlabel('$x$', 'Interpreter', 'Latex', 'FontSize', 20)
ylabel('$y$', 'Interpreter', 'Latex', 'FontSize', 20)

grid on

fprintf('\nResidual Error = %f\n\n', error)

*/
