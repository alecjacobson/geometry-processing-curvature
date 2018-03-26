#include "../include/internal_angles.h"
#include <math.h>
#include <iostream>
using namespace std;

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
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
        /*cout << "Lengths 0: " << l_sqr(curF,0)*10000 << "\n";
        cout << "Lengths 1: " << l_sqr(curF,1)*10000 << "\n";
        cout << "Lengths 2: " << l_sqr(curF,2)*10000 << "\n";
         */
        for (int fNo = 0; fNo < 3; fNo ++){
            //Apply cosine law
            cosAng =l_sqr(curF,vals(fNo,0)) + l_sqr(curF,vals(fNo,1)) - l_sqr(curF,fNo);
            
            cosAng = cosAng / (2.0 * sqrt(l_sqr(curF,vals(fNo,0)) * l_sqr(curF,vals(fNo,1))));
            
            A(curF,fNo) = acos(cosAng);
            
        }
        /*
        cout << (A(curF,0)  + A(curF,1)  + A(curF,2))*180/M_PI << "\n";
        cout << A(curF,0) *180/M_PI << "\n";
        cout << A(curF,1) *180/M_PI << "\n";
        cout << A(curF,2) *180/M_PI << "\n";*/
    }
    
}
