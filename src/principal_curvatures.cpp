#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <Eigen/Eigenvalues> 
#include <igl/sort.h>
#include <Eigen/Core>
#include <igl/pinv.h>

void principal_curvatures(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::MatrixXd & D1,
  Eigen::MatrixXd & D2,
  Eigen::VectorXd & K1,
  Eigen::VectorXd & K2)
{
  // Replace with your code
  K1 = Eigen::VectorXd::Zero(V.rows());
  K2 = Eigen::VectorXd::Zero(V.rows());
  D1 = Eigen::MatrixXd::Zero(V.rows(),3);
  D2 = Eigen::MatrixXd::Zero(V.rows(),3);


    // We need to loop through each vertex and find the nearest vertices connected to that vertex (and the nearest vertices to each of those).
    // This can be done using igl's adjacency matrix function, which basically treats vertices as the nodes
    // of a big graph.

    // First, get the adjacency matrix
    Eigen::SparseMatrix<double> A;
    igl::adjacency_matrix(F, A);

    // Loop through all vertices
    for (int ii = 0; ii < V.rows(); ii++)
    {
        // For this row of the A matrix, every non-zero should be a connected vertex, so get that first list of vertices
        // Then, do the same thing within that loop to get the second ring of vertices.

        // Inner ring
        std::vector<int> adj_vertices; // container for the list of first ring vertices

        for (int jj = 0; jj < A.cols(); jj++)
        {
            if (A.coeff(ii, jj) != 0)
            {
                adj_vertices.push_back(jj);

                // Outer ring
                for (int kk = 0; kk < A.cols(); kk++)
                {
                    if (A.coeff(jj, kk) != 0)
                        adj_vertices.push_back(kk);
                }
            }
        }

        // Now we want to remove all the duplicates. Using C++'s vectors makes this easy.
        std::sort(adj_vertices.begin(), adj_vertices.end());
        adj_vertices.erase(std::unique(adj_vertices.begin(), adj_vertices.end()), adj_vertices.end());

        // Now find the relative vertices and put in the matrix P as per readme
        Eigen::MatrixXd P (adj_vertices.size(), 3);
        for (int jj = 0; jj < adj_vertices.size(); jj++)
            for (int kk = 0; kk < 3; kk++)
                P(jj, kk) = V(adj_vertices[jj], kk) - V(ii, kk);

        // Now, construct P^T*P, on which we need to do eigen decomposition
        Eigen::MatrixXd PTP = P.transpose()*P;

        // Inspired from examples on the Eigen documentation. Tried EigenSolver but that gave complex numbers.
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(PTP);
        Eigen::MatrixXd eigvecs = es.eigenvectors();
        Eigen::VectorXd eigvals = es.eigenvalues();

        // Sort eigen values so we can pick the two largest ones
        Eigen::VectorXi sort_ind;
        Eigen::VectorXd eigvals_sorted;
        igl::sort(eigvals.cwiseAbs(), 1, false, eigvals_sorted, sort_ind);

        // Strongest eigen vector
        Eigen::VectorXd u_vec = eigvecs.col(sort_ind(0)); 

        // Next strongest eigen vector
        Eigen::VectorXd v_vec = eigvecs.col(sort_ind(1));

        // Weakest eigen vector
        Eigen::VectorXd w_vec = eigvecs.col(sort_ind(2));

        // Compute the product with P to get the coefficients
        Eigen::VectorXd u = P*u_vec;
        Eigen::VectorXd v = P*v_vec;
        Eigen::VectorXd w = P*w_vec;

        // Now we can set up and solve the least squares system
        // Construct the system matrix
        Eigen::MatrixXd A2(P.rows(),5);
        A2.col(0) = u;
        A2.col(1) = v;
        A2.col(2) = u.cwiseProduct(u);
        A2.col(3) = u.cwiseProduct(v);
        A2.col(4) = v.cwiseProduct(v);

        // Solve
        Eigen::MatrixXd X2;
        igl::pinv(A2, X2);

        // Multiply by the right hand side
        Eigen::VectorXd a = X2*w;

        // Build the parts of the shape operator
        Eigen::MatrixXd S(2,2), Sa(2,2), Sb(2,2);
        Sa(0,0) = 2.0*a(2)/std::sqrt(a(0)*a(0) + a(1)*a(1) + 1.0);
        Sa(0,1) = a(3)/std::sqrt(a(0)*a(0) + a(1)*a(1) + 1.0);
        Sa(1,0) = a(3)/std::sqrt(a(0)*a(0) + a(1)*a(1) + 1.0);
        Sa(1,1) = 2.0*a(4)/std::sqrt(a(0)*a(0) + a(1)*a(1) + 1.0);
        Sb(0,0) = 1.0 + a(0)*a(0);
        Sb(0,1) = a(0)*a(1);
        Sb(1,0) = a(0)*a(1);
        Sb(1,1) = 1.0 + a(1)*a(1);

        // Compute S
        S = -Sa*Sb.inverse();

        // Now do eigen decomposition again to find max and min curvatures and corresponding directions
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es2(S);
        Eigen::MatrixXd eigvecs_S = es2.eigenvectors();
        Eigen::VectorXd eigvals_S = es2.eigenvalues();

        // Sort eigenvalues so we know which ones correspond to max and min curvature
        Eigen::VectorXi sort_ind_S;
        Eigen::VectorXd eigvals_sorted_S;
        igl::sort(eigvals_S.cwiseAbs(), 1, false, eigvals_sorted_S, sort_ind_S);

        // Store max and min curvatures for this vertex
        K1(ii) = eigvals_sorted_S(0);
        K2(ii) = eigvals_sorted_S(1);

        // Convert corresponding eigen vectors to 3D domain and store. First create the transformation matrix
        // based on the eigen vectors of the original system
        Eigen::MatrixXd to_3D (3,2);
        to_3D.col(0) = u_vec;
        to_3D.col(1) = v_vec;

        D1.row(ii) = to_3D*eigvecs_S.col(sort_ind_S(0));
        D2.row(ii) = to_3D*eigvecs_S.col(sort_ind_S(1));

    }


    return;

}
