#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include "igl/squared_edge_lengths.h"
#include <math.h>
#include <iostream>
using namespace std;

void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
  D = Eigen::VectorXd::Zero(V.rows());
    Eigen::MatrixXd L, IntAngs;
    L.resizeLike(V);
    igl::squared_edge_lengths(V,F, L);
    
    internal_angles(L,IntAngs);
    int numF = F.rows();
    int numV = F.maxCoeff() + 1;
        
    double cosAng;
    /*for (int curV = 0; curV < numV; curV ++) {
        D(curV) = 2.0 * M_PI;
    }*/

    for (int curF = 0; curF < numF; curF ++) {
        
        for (int fNo = 0; fNo < 3; fNo ++){
        //Apply cosine law
            D(F(curF,fNo)) += IntAngs(curF,fNo);
        }
        
    }
    D.array() = 2.0*M_PI - D.array();
    //cout << D << "\n";
}
