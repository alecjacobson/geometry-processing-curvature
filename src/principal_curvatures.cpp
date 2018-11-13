#include "../include/principal_curvatures.h"
#include <igl/adjacency_matrix.h>
#include <igl/per_vertex_normals.h>
#include <set>
#include <Eigen/Eigenvalues>
#include <igl/sort.h>
#include <igl/slice.h>
#include <igl/pinv.h>
#include <math.h>
#include <iostream>

using namespace std;

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

  Eigen::SparseMatrix<double> A;
  igl::adjacency_matrix(F,A);

  vector<vector<int>> adj_list;
  for (int i = 0; i < V.rows(); i++) {
      vector<int> adj;
      // Getting the first ring of the ith vertex
      for (int j = 0; j < V.rows(); j++) {
          if (A.coeff(i, j) != 0 && i != j) {
              adj.push_back(j);
          }
      }
      adj_list.push_back(adj);
  }

  // Used to detect whether the w direction is consistent with normal
  Eigen::MatrixXd N;
  igl::per_vertex_normals(V, F, N);

  // Calculate the required quantities for each vertex
  for (int i = 0; i < V.rows(); i++) {
      set<int> sample;
      for (int j = 0; j < adj_list[i].size(); j++) {
          // Getting the first ring of the ith vertex
          int first_pt = adj_list[i][j];
          sample.insert(adj_list[i][j]);
          // Get the second ring of the ith vertex
          for (int k = 0; k < adj_list[first_pt].size(); k++) {
              if (adj_list[first_pt][k] != i) {
                  sample.insert(adj_list[first_pt][k]);
              }
          }
      }
      Eigen::MatrixXd P(sample.size(), 3);
      set<int>::iterator it; int k = 0;
      for (it = sample.begin(); it != sample.end(); it++) {
          P.row(k) = V.row(*it) - V.row(i);
          k++;
      }

      // Eigen decomposition on P^T * P
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(P.transpose() * P);
      Eigen::MatrixXd eigenval(es.eigenvalues().rows(), 1);
      eigenval.col(0) = es.eigenvalues();
      Eigen::MatrixXd Y, IX;
      igl::sort(eigenval, 1, false, Y, IX);

      // Sort the eigenvalue eigenvector pairs
      Eigen::MatrixXd sortedVec(3, 3);
      sortedVec.col(0) = es.eigenvectors().col(IX(0));
      sortedVec.col(1) = es.eigenvectors().col(IX(1));
      sortedVec.col(2) = es.eigenvectors().col(IX(2));

      // Get the coefficients
      Eigen::MatrixXd coeff = P * sortedVec;
      Eigen::MatrixXd ls_matrix(P.rows(), 5);

      ls_matrix.col(0) = coeff.col(0).array();
      ls_matrix.col(1) = coeff.col(1).array();
      ls_matrix.col(2) = coeff.col(0).array().pow(2);
      ls_matrix.col(3) = coeff.col(0).array() * coeff.col(1).array();
      ls_matrix.col(4) = coeff.col(1).array().pow(2);

      // Check whether the w direction is consistent with the normal
      Eigen::VectorXd rhs = coeff.col(2);
      if (N.row(i) * sortedVec.col(2) > 0) {
          rhs = coeff.col(2);
      } else {
          rhs = -coeff.col(2);
      }

      // Solve the least quare fitting problem
      Eigen::MatrixXd pseudoinv;
      igl::pinv(ls_matrix, pseudoinv);
      Eigen::VectorXd a = pseudoinv * rhs;

      // Calculate the first and second fundamental forms
      double E = 1 + pow(a(0), 2);
      double Fn = a(0) * a(1);
      double G = 1 + pow(a(1), 2);
      double denom = sqrt(pow(a(0), 2) + 1 + pow(a(1), 2) );
      double e = 2.*a(2) / denom;
      double f =    a(3) / denom;
      double g = 2.*a(4) / denom;
      Eigen::MatrixXd secnd_fund_form(2, 2);
      secnd_fund_form(0,0) = e;
      secnd_fund_form(0,1) = f;
      secnd_fund_form(1,0) = f;
      secnd_fund_form(1,1) = g;
      Eigen::MatrixXd first_fund_form(2, 2);
      first_fund_form(0,0) = E;
      first_fund_form(0,1) = Fn;
      first_fund_form(1,0) = Fn;
      first_fund_form(1,1) = G;
      Eigen::MatrixXd S = (-secnd_fund_form * first_fund_form.inverse()).transpose();

      // Perform eigen decomposition on S
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ss(S);
      Eigen::MatrixXd sseigenval(ss.eigenvalues().rows(), 1);;
      sseigenval.col(0) = ss.eigenvalues();
      Eigen::MatrixXd sY, sIX;
      igl::sort(sseigenval, 1, false, sY, sIX);
      Eigen::MatrixXd sseigenvec(2,2);
      sseigenvec.col(0) = ss.eigenvectors().col(0);
      sseigenvec.col(1) = ss.eigenvectors().col(1);

      K1(i) = ss.eigenvalues()[sIX(0)];
      D1.row(i) = (sortedVec.block(0, 0, 3, 2) * ss.eigenvectors().col(sIX(0))).transpose().normalized();
      K2(i) = ss.eigenvalues()[sIX(1)];
      D2.row(i) = (sortedVec.block(0, 0, 3, 2) * ss.eigenvectors().col(sIX(1))).transpose().normalized();
  }
}
