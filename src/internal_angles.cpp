#include "../include/internal_angles.h"
#include <math.h>


void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
    
    //Compute number of faces
    int numF = l_sqr.rows();
    
    A.resize(numF,3);
    Eigen::MatrixXd vals(3,2);
    
    //Makes it easier to reference other vertices
    vals(0,0) = 1;
    vals(0,1) = 2;
    vals(1,0) = 0;
    vals(1,1) = 2;
    vals(2,0) = 0;
    vals(2,1) = 1;
    
    double cosAng;    
    for (int curF = 0; curF < numF; curF ++) {
        
        for (int fNo = 0; fNo < 3; fNo ++){
            //Apply cosine law
            cosAng =l_sqr(curF,vals(fNo,0)) + l_sqr(curF,vals(fNo,1)) - l_sqr(curF,fNo);
            
            cosAng = cosAng / (2.0 * sqrt(l_sqr(curF,vals(fNo,0)) * l_sqr(curF,vals(fNo,1))));
            
            //Acos to get back the angle
            A(curF,fNo) = acos(cosAng);
            
        }

    }
    
}
