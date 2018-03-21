#include "../include/angle_defect.h"
#include "internal_angles.h"
#include "igl/squared_edge_lengths.h"

void angle_defect(const Eigen::MatrixXd & V, const Eigen::MatrixXi & F,
		Eigen::VectorXd & D) {

	Eigen::MatrixXd squaredLengths;
	igl::squared_edge_lengths(V, F, squaredLengths);

	Eigen::MatrixXd internalAngles;
	internal_angles(squaredLengths, internalAngles);

	int faceCount = F.rows();
	D = Eigen::VectorXd::Zero(V.rows());
	for (int i = 0; i < faceCount; i++) {

		// Iterate through each vertex of the face to contribute to the count
		// of internal angles.
		for (int v = 0; v < 3; v++) {

			int vertexIndex = F(i, v);
			D(vertexIndex) += internalAngles(i, v);
		}
	}

	D = -D + Eigen::VectorXd::Ones(V.rows()) * 2 * M_PI;

}
