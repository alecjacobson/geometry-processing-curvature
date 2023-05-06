#include "../include/angle_defect.h"
#include <igl/squared_edge_lengths.h>
#include "../include/internal_angles.h"

#define PI 3.14159

void angle_defect(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::VectorXd &D) {

    D = Eigen::VectorXd(V.rows());

    Eigen::MatrixXd L;
    igl::squared_edge_lengths(V, F, L);

    // Given sqr edge-lengths compute the internal
    // angles at each corner

    Eigen::MatrixXd A;
    internal_angles(L, A);

    D.setOnes();
    D = D * 2 * PI;

    for (int i = 0; i < F.rows(); i++) {
        for (int j = 0; j < 3; j++) {
            int idx = F(i, j);
            D(idx) = D(idx) - A(i, j);
        }
    }

}
