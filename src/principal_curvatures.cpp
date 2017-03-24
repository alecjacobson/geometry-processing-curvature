#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <igl/per_vertex_normals.h>
#include <Eigen/SVD>
#include <Eigen/Eigenvalues>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);

  /*
    Using neighbours of neighbours (P ϵℝ^{kx3} relative to the vᵢ) around each vertex V, define a quadratic surface as a height field above some 2D plane passing through V
    
    find plane by the SVD of (PᵀP), with S ϵℝ^{kx2} as the coefficients for most principal directions for each point in P, 
    and B ϵℝ^k be the height of each point in least pricipal direction

    Quadratic fn as a height-field through origin given by:
    w = a₁u + a₂v + a₃u² + a₄uv + a₅v²

    where we have k sets of u, v, and w. Can perform a least squares fitting (igl::pinv)
    
    the shape operator can be constructed as:

    S = - / e f \ /E F\⁻¹
          \ f g / \F G/

    where we have:
    E = 1 + a₁²
    F = a₁a₂
    G = 1 + a₂²
    e = 2a₃/sqrt(a₁² + 1 + a₂²)
    f =  a₄/sqrt(a₁² + 1 + a₂²)
    e = 2a₅/sqrt(a₁² + 1 + a₂²)    

    Adjacency matrix will give us the connectivity, and each row specifies which vertex is connected. If we want the neighbours neighbours, can square the matrix and find non-zero rows in that
   */

  Eigen::SparseMatrix<double> adj;
  igl::adjacency_matrix(F, adj);
  Eigen::SparseMatrix<double> adj2 = adj * adj; 

  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  int32_t n = V.rows();
  int32_t f = F.rows();

  for(int32_t i = 0; i < adj2.outerSize(); i++)
  {
    int32_t numNeighbours = adj2.innerVector(i).nonZeros();
    Eigen::MatrixXd P(numNeighbours, 3);

    int32_t r = 0;
    for(Eigen::SparseMatrix<double>::InnerIterator it(adj2, i); it; ++it)
    {
      P.row(r) = V.row(it.index());
      r++;
    }

    // P needs to relative to the "middle" point:
    P.rowwise() -= V.row(i);

    Eigen::Matrix3d PTP = P.transpose() * P;
    // regular eigensolve returns a complex matrix, but here we can take advantage of the "self-adjoint-ness" of the matrix
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigPTP(PTP);
    
    Eigen::Matrix3d evec = (eigPTP.eigenvectors());
    // w is the "least-principal", uv are "most-principal" directions
    Eigen::VectorXd u, v, w;
    Eigen::RowVectorXd uu, vv, ww;
    ww = evec.col(0).transpose();
    vv = evec.col(1).transpose();
    uu = evec.col(2).transpose();

    if(N.row(i) * ww.transpose() < 0)
      ww = -ww;

    u = P * uu.transpose();
    v = P * vv.transpose();
    w = P * ww.transpose();

    // want to solve for the a's so, Qa = w and can do the pseudo-inverse of Q to perform least squares fitting
    Eigen::MatrixXd Q(w.rows(), 5);
    Q.col(0) = u;
    Q.col(1) = v;
    Q.col(2) = u.array().pow(2);
    Q.col(3) = u.cwiseProduct(v);
    Q.col(4) = v.array().pow(2);

    Eigen::MatrixXd psinv;

    igl::pinv(Q, psinv);
    Eigen::VectorXd A = psinv * w;

    double denom = std::sqrt(A(0)*A(0) + 1 + A(1)*A(1));

    Eigen::Matrix2d efg1, EFG2;

    efg1(0,0) = 2*A(2)/denom;
    efg1(1,0) = efg1(0,1) = A(3)/denom;
    efg1(1,1) = 2*A(4)/denom;

    EFG2(0,0) = 1 + A(0)*A(0);
    EFG2(1,0) = EFG2(0,1) = A(0)*A(0);
    EFG2(1,1) = 1 + A(1)*A(1);

    Eigen::Matrix2d S = -efg1 * EFG2.inverse();

    // now do an eigen decomp on S to find k₁ and k₂, and tangent directions d₁ and d₂:
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigS(S);
    Eigen::Vector2d eval = eigS.eigenvalues();

    // vectors are in the uv basis, need to transform back into global coords:
    Eigen::Matrix<double, 2,3> uv;
    uv << uu, vv;
    Eigen::Matrix<double, 2,3> dat1 = (uv.transpose() * eigS.eigenvectors()).transpose();
    if(eval(0) > eval(1))
    {
      K1(i) = eval(0);
      K2(i) = eval(1);
    
      D1.row(i) = dat1.row(1);
      D2.row(i) = dat1.row(0);
    }
    else
    {
      K1(i) = eval(1);
      K2(i) = eval(0);
      
      D1.row(i) = dat1.row(0);
      D2.row(i) = dat1.row(1);
    }    
  }
}
