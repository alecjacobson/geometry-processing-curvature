#include "../include/angle_defect.h"
#include "../include/internal_angles.h"

#include <igl/squared_edge_lengths.h>

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    // Interior angles

    Eigen::MatrixXd l_sqr, A;
    igl::squared_edge_lengths(V, F, l_sqr);
    internal_angles(l_sqr, A);

    // Compute angle defect
    //      \Sigma_{i} = 2\pi - \sum_{f\in F(i)} InteriorAngle(i, f)

    D = Eigen::VectorXd::Ones(V.rows());
    D = 2 * M_PI * D;

    for (int i = 0; i < F.rows(); ++i) {
        for (int j = 0; j < 3; ++j) {
            auto v = F(i, j);
            D(v) = D(v) - A(i, j);
        }
    }
}