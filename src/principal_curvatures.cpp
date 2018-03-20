#include "../include/principal_curvatures.h"
#include "igl/adjacency_matrix.h"
#include "igl/per_vertex_normals.h"
#include "igl/slice_mask.h"
#include <Eigen/Dense>
#include "igl/pinv.h"
#include <math.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  int n = V.rows();

  // Principal Curvatures (to be computed)
  K1 = Eigen::VectorXd::Zero(n);
  K2 = Eigen::VectorXd::Zero(n);
  // Principal Directions (to be computed)
  D1 = Eigen::MatrixXd::Zero(n, 3);
  D2 = Eigen::MatrixXd::Zero(n, 3);
 
  Eigen::SparseMatrix<double> adjacency_matrix;
  igl::adjacency_matrix(F, adjacency_matrix);
  // Ensuring positivity of Cotangent adjacency elements (adjacency is symmetric of course)
  Eigen::SparseMatrix<double> one_step = adjacency_matrix.cwiseAbs();
  // This will return all reachable areas in two hops (symmetric)
  Eigen::SparseMatrix<double> two_step = one_step * one_step;

  // Each row and column (by symmetry) which corresponds to a vertex number
  // will contain True for any vertices in the two ring of that vertex
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> twoRingNeighbors(n, n);
  twoRingNeighbors = (Eigen::MatrixXd::Identity(n , n) + Eigen::MatrixXd(one_step) + Eigen::MatrixXd(two_step)).array() > 0;

  // Normals will be used for alignment of height field
  Eigen::MatrixXd Normals_per_vertex;
  igl::per_vertex_normals(V, F, Normals_per_vertex);

  for(int i = 0; i < n; i++){
    // Obtain two ring of v
    Eigen::MatrixXd twoRing;
    igl::slice_mask(V, twoRingNeighbors.row(i), 1, twoRing);
    int numNeighbors = twoRing.rows();
    
    // Gather the positions of these k points relative to v
    Eigen::MatrixXd P = twoRing.rowwise() - V.row(i);

    // Perform PCA on P to find a locally best-fit plane at v, and the height offsets from this plane
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> PCA(P.transpose() * P);

    // The eigenvalues are sorted in increasing order

    // The “height” of each point in the least principal direction
    Eigen::Vector3d h = PCA.eigenvectors().col(0);
 
    // The coefficients for two most principal directions define a plane
    Eigen::Vector3d t1 = PCA.eigenvectors().col(1);
    Eigen::Vector3d t2 = PCA.eigenvectors().col(2);

    // Best align the height with the normal
    Eigen::Vector3d n = Normals_per_vertex.row(i);
    double scalar_product = h.dot(n);
    if (scalar_product < 0){
      // invert the height direction
      h = -1.0 * h;
    }

    // Project the differences onto the three components
    Eigen::VectorXd p1 = P * t1; 
    Eigen::VectorXd p2 = P * t2; 
    Eigen::VectorXd heights = P * h;

    // Find a best fit function from the (t₁, t₂) plane to the height = w parameterized by coeffients a
    // w₁ = a₁t₁ + a₂t₂ + a₃t₁² + a₄ t₁t₂ + a₅ t₂²
    // w₁ = a₁c₁ + a₂c₂ + a₃c₃ + a₄c₄ + a₅c₅ where C is a coefficient matrix
    // the system is therefore: w = Ca, where we wish to solve for a
    Eigen::MatrixXd C(numNeighbors, 5);
    C.col(0) = p1;
    C.col(1) = p2;
    C.col(2) = p1.array().square();
    C.col(3) = p1.array() * p2.array();
    C.col(4) = p2.array().square();
    
    // Moore–Penrose inverse
    Eigen::MatrixXd p_inverse;
    igl::pinv(C, p_inverse);
    // Solves the least squares fit for w = Ca
    Eigen::VectorXd a = p_inverse * heights;

    // Construct the shape operator from its fundamental forms 
    Eigen::Matrix2d L;
    double normalization = pow(a(0), 2) + pow(a(1), 2) + 1;
    double e = 2 * a(2) / normalization;
    double f = a(3) / normalization;
    double g = 2 * a(4) / normalization;
    L << e, f, f, g;

    Eigen::Matrix2d R;
    double E = 1 + pow(a(0), 2);
    double F = a(0) * a(1);
    double G = 1 + pow(a(1), 2);
    R << E, F, F, G;
    
    // The Shape Operator S is symmetric since L and R are symmetric
    Eigen::Matrix2d S = -1.0 * L * R.inverse();
    
    // Use self-adjoint decomposition of S to get principal curvatures and principal tangent directions
    // which will be real
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> PCA_S(S);

    Eigen::Vector2d eigen_values = PCA_S.eigenvalues();
    Eigen::Matrix2d eigen_vectors = PCA_S.eigenvectors();

    // Project the operator onto the eigen vectors to obtain the principal directions
    Eigen::Matrix2d principal_directions = S * eigen_vectors;

    // Assign the min and max directions respectively
    // Lift the principal tangent directions back to world coordinates
    K1(i) = eigen_values(0);
    D1.row(i) = eigen_vectors(0, 0) * t1 + eigen_vectors(1, 0) * t2;

    K2(i) = eigen_values(1);
    D2.row(i) = eigen_vectors(0, 1) * t1 + eigen_vectors(1, 1) * t2;
  }
}
