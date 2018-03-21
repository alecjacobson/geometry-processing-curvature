#include "../include/mean_curvature.h"
#include "igl/per_vertex_normals.h"
#include "igl/massmatrix.h"
#include "igl/cotmatrix.h"
#include "igl/invert_diag.h"

void mean_curvature(
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		Eigen::VectorXd & H) {

	Eigen::SparseMatrix<double> massMatrix;
	igl::massmatrix(V, F, 
			igl::MassMatrixType::MASSMATRIX_TYPE_DEFAULT, massMatrix);

	// Mass matrix is diagonal, so this should work.
	Eigen::SparseMatrix<double> inverseMassMatrix;
	igl::invert_diag(massMatrix, inverseMassMatrix);

	Eigen::SparseMatrix<double> laplacianMatrix;
	igl::cotmatrix(V, F, laplacianMatrix);

	Eigen::MatrixXd curvatureNormals = inverseMassMatrix * laplacianMatrix * V;

	// We can verify the sign of the curvature by comparing against the
	// normals.
	Eigen::MatrixXd normals;
	igl::per_vertex_normals(V, F, normals);

	int vertexCount = V.rows();
	H.resize(curvatureNormals.rows());
	for (int i = 0; i < vertexCount; i++) {

		// The curvature normal is supposed to face away from the surface
		// normal.
		Eigen::Vector3d curvatureNormal = curvatureNormals.row(i);
		Eigen::Vector3d normal = normals.row(i);
		double dot = normal.dot(curvatureNormal);

		H(i) = curvatureNormal.norm();
		if (dot > 0) {
			H(i) = -H(i);
		}
	}
}
