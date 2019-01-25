/**
 * @name pca.c
 * An example of using uLAPAck for principle component analysis.
 *
 * @note In this example, static memory allocation is used
 *       via the -D compiler options.
 * @note Build this example using: $ make pca
 * @note Run this example using: $ ./pca
 */

#include "ulapack.h"

#include <stdio.h>

int main(void) {

    double Adata[10][8] = { {0.4170, 0.4192, 0.8007, 0.0983, 0.9889, 0.0194, 0.1023, 0.9034},
                            {0.7203, 0.6852, 0.9683, 0.4211, 0.7482, 0.6788, 0.4141, 0.1375},
                            {0.0001, 0.2045, 0.3134, 0.9579, 0.2804, 0.2116, 0.6944, 0.1393},
                            {0.3023, 0.8781, 0.6923, 0.5332, 0.7893, 0.2655, 0.4142, 0.8074},
                            {0.1468, 0.0274, 0.8764, 0.6919, 0.1032, 0.4916, 0.0500, 0.3977},
                            {0.0923, 0.6705, 0.8946, 0.3155, 0.4479, 0.0534, 0.5359, 0.1654},
                            {0.1863, 0.4173, 0.0850, 0.6865, 0.9086, 0.5741, 0.6638, 0.9275},
                            {0.3456, 0.5587, 0.0391, 0.8346, 0.2936, 0.1467, 0.5149, 0.3478},
                            {0.3968, 0.1404, 0.1698, 0.0183, 0.2878, 0.5893, 0.9446, 0.7508},
                            {0.5388, 0.1981, 0.8781, 0.7501, 0.1300, 0.6998, 0.5866, 0.7260}};

    Matrix_t A;
    Matrix_t T;

    ulapack_init(&A, 10, 8);
    ulapack_init(&T, 10, 10);

    /*
     * Copy data points into vector objects.
     */
    for (uint64_t row_itor = 0; row_itor < 10; row_itor++) {
        for (uint64_t col_itor = 0; col_itor < 8; col_itor++) {
            ulapack_edit_entry(&A, 
            row_itor, col_itor, 
            Adata[row_itor][col_itor]);
        }
    }

    ulapack_pca(&A, &T);

    printf("\nA = \n");
    ulapack_print(&A, stdout);

    printf("\nT = \n");
    ulapack_print(&T, stdout);
}

/*
 * Corresponding MATLAB code for comparison.
 */

/*

clear all;
close all;
clc;

rng(1);

A = rand(10, 8);

[U,S,V] = svd(A * A');
T = U * S

*/