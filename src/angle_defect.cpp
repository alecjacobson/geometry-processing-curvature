#include "../include/angle_defect.h"
#include "internal_angles.h"
#include "igl/squared_edge_lengths.h"
#include "igl/massmatrix.h"

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    Eigen::MatrixXd l_sqr, A;
    Eigen::SparseMatrix<double> M;
    D.resize(V.rows());
    
    igl::squared_edge_lengths(V,F,l_sqr);
    igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
    internal_angles(l_sqr,A);
    
    for(int i = 0; i < F.rows(); i++){
        for(int k = 0; k < F.cols(); k++){
            D(F(i,k)) += A(i,k);
        }
    }
    D = ((2*M_PI - D.array())/M.diagonal().array()).matrix();
}
