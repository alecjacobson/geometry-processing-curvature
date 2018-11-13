#include "../include/mean_curvature.h"

#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <Eigen/Dense>
#include <igl/per_vertex_normals.h>

void mean_curvature(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & H)
{
  // Replace with your code
  H = Eigen::VectorXd::Zero(V.rows());

// calculate M and L
  	// calculate L, here, w is 1
  	Eigen::SparseMatrix<double> L;
  	igl::cotmatrix(V, F, L);

  	// M
  	Eigen::SparseMatrix<double> M;
  	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

  	// MH = LV
  	Eigen::MatrixXd b = L * V;

  	// svd
  	H = M.colPivHouseholderQr().solve(b);

  	// sign
  	Eigen::MatrixXd vnorms;
  	igl::per_vertex_normals(V, F, igl::PER_VERTEX_NORMALS_WEIGHTING_TYPE_AREA, vnorms);
  	int vnum = V.rows();
  	for (int i = 0; i< vnum; i++)
  	{
  		double sum = H.row(i).transpose() * vnorms.row(i);
  		if (sum > 0)
  		{
  			H.row(i) = -H.row(i);
  		}
  	}

}
