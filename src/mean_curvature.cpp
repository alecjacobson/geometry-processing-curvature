#include "../include/mean_curvature.h"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  Eigen::MatrixXd Normals_per_vertex;
  igl::per_vertex_normals(V, F, Normals_per_vertex);

  Eigen::SparseMatrix<double> Laplacian;
  igl::cotmatrix(V, F , Laplacian);

  Eigen::SparseMatrix<double> Mass;
  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, Mass);

  Eigen::SparseMatrix<double> Mass_inv;
  igl::invert_diag(Mass, Mass_inv);
  
  // Approximation of the mean curvature normal 
  Eigen::MatrixXd H_approximation = Mass_inv * Laplacian * V;

  // Gather norms which approximate the mean curvature
  H = H_approximation.rowwise().norm();

  // We can use the normals to correct the sign of our approximation
  // by checking whether each row in H_approximation
  // agrees with consistently oriented per-vertex normals in Normals_per_vertex
  for(int vertexNumber = 0; vertexNumber < V.rows(); vertexNumber++){
    Eigen::VectorXd unsigned_Hn = H_approximation.row(vertexNumber);
    Eigen::VectorXd n = Normals_per_vertex.row(vertexNumber);
    double scalar_product = unsigned_Hn.dot(n);
    if (scalar_product > 0.0) // Flip if positive to match convention
    {
      H(vertexNumber) *= -1.0;
    }
  }
}
