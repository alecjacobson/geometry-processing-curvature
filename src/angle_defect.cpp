#include "../include/angle_defect.h"
#include "../include/internal_angles.h"
#include "igl/squared_edge_lengths.h"
#include <math.h>


void angle_defect(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  Eigen::VectorXd & D)
{
    //Reshape D
    D = Eigen::VectorXd::Zero(V.rows());
    
    //L is the squared edge lengths
    //IntAngs stores the interior angles
    Eigen::MatrixXd L, IntAngs;
    L.resizeLike(V);
    
    //Compute squared edge lengths
    igl::squared_edge_lengths(V,F, L);
    
    //Compute interior angles
    internal_angles(L,IntAngs);
    
    //Extract number of faces and vertices
    int numF = F.rows();
    int numV = F.maxCoeff() + 1;
        
    double cosAng;

    //Loop over the faces
    for (int curF = 0; curF < numF; curF ++) {
        
        //Loop over the vertices
        for (int fNo = 0; fNo < 3; fNo ++){
            //Sum the interior angles
            D(F(curF,fNo)) += IntAngs(curF,fNo);
        }
        
    }
    
    //Compute the angle defect
    D.array() = 2.0*M_PI - D.array();

}
