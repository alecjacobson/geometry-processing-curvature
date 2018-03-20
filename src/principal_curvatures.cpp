#include "../include/principal_curvatures.h"

#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/pinv.h>
#include <map>

// helper function used to compute the relative positions of the 2-ring of a given vertex
void ring(const Eigen::SparseMatrix<double> Adj, const Eigen::MatrixXd & V, const int v, std::map<int, Eigen::Vector3d> & neighbours) {

    // iterate through the non-zero entries of the relevant row and add
    // each adjacent vertex to the list of neighbours
    for (Eigen::SparseMatrix<double>::InnerIterator it(Adj, v); it; ++it) {

        // if v isn't in neighbours, then we haven't fully populated neighbours
        // with the 1-ring of our desired vertex. When we call ring below
        // to compute the 2-ring, all of the passed vertices will already be keys in
        // neighbours
        if (neighbours.find(v) == neighbours.end()) {
            neighbours[it.index()] = V.row(it.index()) - V.row(v);
        } else if (neighbours.find(it.index()) == neighbours.end()) {
            // if the current neighbour is already in neighbours, then it is also
            // a neighbour of our desired vertex so there's no need to recompute
            // the relative position. Otherwise, we need to compute the 2-ring relative 
            // position
            Eigen::Vector3d diff = V.row(it.index()) - V.row(v);
            neighbours[it.index()] = diff + neighbours[v]; 
        }
    }  
}

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2) {
    
    // construct an adjacency matrix for the mesh
    Eigen::SparseMatrix<double> Adj;
    igl::adjacency_matrix(F, Adj);
    
    // compute the per vertex normals for the mesh
    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    
    // initialize output matrices
    K1 = Eigen::VectorXd::Zero(V.rows());
    K2 = Eigen::VectorXd::Zero(V.rows());
    D1 = Eigen::MatrixXd::Zero(V.rows(),3);
    D2 = Eigen::MatrixXd::Zero(V.rows(),3);
    
    // loop through the vertices
    for (int v = 0; v < V.rows(); v++) {
    
        // get the relative positions for v's 1-ring
        std::map<int, Eigen::Vector3d> one_ring;
        ring(Adj, V, v, one_ring);
       
        // get the relative positions v's 2-ring
        std::map<int, Eigen::Vector3d> two_ring = one_ring;
        for (auto n = one_ring.begin(); n != one_ring.end(); n++) {
            ring(Adj, V, n->first, two_ring);
        }
        
        // the ring function ignores that v will be a neighbour of each of its neighbours.
        // This means two_ring will contain an entry for v, which it makes sense to discard here
        two_ring.erase(v);
        
        // store the relative positions for each of the vertices in the 2-ring
        Eigen::MatrixXd P(two_ring.size(), 3);
        int i = 0;
        for (auto n = two_ring.begin(); n != two_ring.end(); n++) {
            P.row(i) = n->second;
            i++;
        }
        
        // compute the principal component analysis of P using eigen decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(P.transpose()*P);

        // extract the coefficients of the two most principal directions
        Eigen::MatrixXd S = es.eigenvectors().block(0, 1, 3, 2);

        // compute B (the height of each point in the least principal direction), remembering
        // to account for the orientation of the vertex normal
        Eigen::Vector3d w = es.eigenvectors().col(0);
        
        if (N.row(v).dot(w) < 0.0) {
            w *= -1;
        }
        
        Eigen::MatrixXd B = P*w;
        
        // build the least squares fitting problem for the quadratic function
        // as a height-field surface
        Eigen::MatrixXd A(P.rows(), 5);
        Eigen::MatrixXd uv = P*S;
        Eigen::MatrixXd uu = uv.col(0).asDiagonal();
        Eigen::MatrixXd vv = uv.col(1).asDiagonal();

        A.col(0) = uv.col(0);
        A.col(1) = uv.col(1);
        A.col(2) = (uu*uu).diagonal();
        A.col(3) = uv.rowwise().prod();
        A.col(4) = (vv*vv).diagonal();
        
        // solve the least squares fitting problem
        Eigen::MatrixXd Ainv;
        igl::pinv(A, Ainv);
        Eigen::VectorXd coeffs = Ainv*B;
        
        // construct the second and first fundamental forms
        double E = 1 + coeffs(0)*coeffs(0);
        double F = coeffs(0)*coeffs(1);
        double G = 1 + coeffs(1)*coeffs(1);
        
        double e = (2*coeffs(2))/(std::sqrt(coeffs(0)*coeffs(0) + 1 + coeffs(1)*coeffs(1)));
        double f = (coeffs(3))/(std::sqrt(coeffs(0)*coeffs(0) + 1 + coeffs(1)*coeffs(1)));
        double g = (2*coeffs(4))/(std::sqrt(coeffs(0)*coeffs(0) + 1 + coeffs(1)*coeffs(1)));
        
        Eigen::Matrix2d first;
        first << E, F,
                 F, G;
        
        Eigen::Matrix2d second;
        second << e, f,
                  f, g;
        
        // build the shape operator from the first and second fundamental forms
        Eigen::Matrix2d S_op = -second*first.inverse();
        
        // compute the principal component analysis of S_op using eigen decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> ses(S_op);
        Eigen::Vector2d curvs = ses.eigenvalues();
        
        // principal directions
        Eigen::Matrix2d vecs = ses.eigenvectors();
        Eigen::MatrixXd dirs = S*vecs;
        
        // assign max and min appropriately
        int max = curvs(0) > curvs(1) ? 0 : 1;
        int min = 1 - max;
        
        
        K1(v) = curvs(min);
        K2(v) = curvs(max);
        D1.row(v) = dirs.col(min);
        D2.row(v) = dirs.col(max);
                
    }
}
