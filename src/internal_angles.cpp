#include "../include/internal_angles.h"

void internal_angles(
  const Eigen::MatrixXd & l_sqr,
  Eigen::MatrixXd & A)
{
  int32_t r = l_sqr.rows();

  A.resize(r, 3);

  // internal angles for the triangles using cosine law, where l_sqr are the squared side lengths 
  for(int32_t i = 0; i < r; i++)
  {
    for(int32_t j = 0; j < 3; j++)
    {
      double a = l_sqr(i, (j + 0) % 3);
      double b = l_sqr(i, (j + 1) % 3);
      double c = l_sqr(i, (j + 2) % 3);

      if(b == 0 || c == 0)
	printf("ERR");
      A(i, j) = std::acos( -(a - b - c)/(2*std::sqrt(b*c)) );
    }
  }
  
}
