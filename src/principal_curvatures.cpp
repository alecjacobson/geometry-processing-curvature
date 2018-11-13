#include "../include/principal_curvatures.h"
#include <Eigen/Core>
#include <igl/adjacency_matrix.h>
#include <igl/pinv.h>
#include <igl/cat.h>

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
    Eigen::SparseMatrix<int> adj;
  igl::adjacency_matrix(F, adj);
  for (int i=0; i<V.rows(); i++){
  	adj.insert(i,i)=1;
  }
  adj = adj*adj;

  for (int i = 0; i < V.rows(); ++i)
  {
    Eigen::MatrixXd P(adj.col(i).nonZeros(), 3);
    int p = 0;
    for (Eigen::SparseMatrix<int>::InnerIterator it(adj, i); it; ++it)
    {
      P.row(p) = V.row(it.index()) - V.row(i);
      p++;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXd UU = svd.matrixV();
    Eigen::MatrixXd UV = P*UU.leftCols(2);
    Eigen::MatrixXd W = P*UU.col(2);

    Eigen::MatrixXd A(UV.rows(), 5);
    for (int j = 0; j < UV.rows(); ++j)
        A.row(j) << UV(j, 0), UV(j, 1), UV(j, 0)*UV(j, 0), UV(j, 0)*UV(j, 1), UV(j, 1)*UV(j, 1);

    Eigen::MatrixXd Ainv; 
	igl::pinv(A, Ainv);
    Eigen::VectorXd a = Ainv*W;

    double E=1+a(0)*a(0);
    double F1=a(0)*a(1);
    double G=1+a(1);
    double sqt=sqrt(a(0)*a(0)+1+a(1)*a(1));
    double e=2*a(2)/sqt;
    double f=a(3)/sqt;
    double g=2*a(4)/sqt;
    Eigen::MatrixXd S1(2,2),S2(2,2),S;
    S1 << e,f,
          f,g;
    S2 << E,F1,
          F1,G;
    S=-S1*S2.inverse();

	Eigen::JacobiSVD<Eigen::MatrixXd> svd1(S, Eigen::ComputeFullU);
    Eigen::MatrixXd svdU = svd1.matrixU();
	Eigen::VectorXd svdA=svd1.singularValues();
    K1(i)=svdA(0);
    K2(i)=svdA(1);
    D1.row(i)=UU.leftCols(2)*svdU.col(0);
	D2.row(i)=UU.leftCols(2)*svdU.col(1);
  }
  /*Eigen::SparseMatrix<int> adj;
  igl::adjacency_matrix(F,adj);
  for (int i=0; i<V.rows(); i++){
  	adj.insert(i,i)=1;
  }
  adj=adj*adj;
  
  int p;
  for (int i=0; i<V.rows(); i++){
  	Eigen::MatrixXd P(adj.col(i).nonZeros(),3);
  	p=0;
  	for (Eigen::SparseMatrix<int>::InnerIterator it(adj,i); it; ++it){
  	  	P.row(p)=V.row(it.index())-V.row(i);
  	  	p++;
	}
	printf("a\n");
	Eigen::MatrixXd PTP=P.transpose()*P;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(PTP, Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::MatrixXd UU = svd.matrixV();
    Eigen::MatrixXd UV=P*UU.leftCols(2);
    Eigen::MatrixXd W=P*UU.rightCols(1);
	printf("b\n");
    //Eigen::VectorXd U2=(U.array()*U.array()).matrix();
    //Eigen::VectorXd UV=(U.array()*V.array()).matrix();
    //Eigen::VectorXd V2=(V.array()*V.array()).matrix();
    Eigen::MatrixXd A(UV.rows(),5),Ainv;
	printf("c");
	for (int j=0; j<UV.rows(); j++){
	  A.row(j) << UV(j,0),UV(j,1),UV(j,0)*UV(j,0),UV(j,0)*UV(j,1),UV(j,1)*UV(j,1);
	}
    //A << U,V,U2,UV,V2;
	printf("d\n");
    igl::pinv(A,Ainv);
	printf("e\n");
    Eigen::VectorXd a=Ainv*W;
    double E=1+a(0)*a(0);
    double F1=a(0)*a(1);
    double G=1+a(1);
    double sqt=sqrt(a(0)*a(0)+1+a(1)*a(1));
    double e=2*a(2)/sqt;
    double f=a(3)/sqt;
    double g=2*a(4)/sqt;
    Eigen::MatrixXd S1,S2,S;
    S1 << e,f,
          f,g;
    S2 << E,F1,
          F1,G;
    S=-S1*S2.inverse();
	Eigen::JacobiSVD<Eigen::MatrixXd> svd1(S, Eigen::ComputeFullU);
    Eigen::MatrixXd svdU = svd1.matrixU();
	Eigen::VectorXd svdA=svd1.singularValues();
    K1(i)=svdA(0);
    K2(i)=svdA(1);
    D1.row(i)=UU.leftCols(2)*svdU.col(0);
	D2.row(i)=UU.leftCols(2)*svdU.col(1);
  }*/
  
}
