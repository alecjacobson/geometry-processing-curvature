#include "../include/internal_angles.h"
#include <igl/cotmatrix.h>


void internal_angles(
        const Eigen::MatrixXd &l_sqr,
        Eigen::MatrixXd &A) {

    int N = l_sqr.rows();
    A.resize(N, 3);

    for (int i = 0; i < N; i++) {
        //let's do a for loop with mod 3 here for simplicity.
        for (int j = 0; j < 3; j++) {

            double c = l_sqr(i, j % 3);
            double a = l_sqr(i, (j + 1) % 3);
            double b = l_sqr(i, (j + 2) % 3);

            double cos_ij = (a + b - c) / (2 * sqrt(a) * sqrt(b));

            A(i, j) = acos(cos_ij);
        }
    }

}
