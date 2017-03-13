#include <igl/read_triangle_mesh.h>
#include <igl/parula.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Core>


int main(int argc, char *argv[])
{
  // Scale for the color axis
  double scale = 1.0;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  // Load input meshes
  igl::read_triangle_mesh(
    (argc>1?argv[1]:"../shared/data/cactus.obj"),V,F);
  igl::viewer::Viewer viewer;
  std::cout<<R"(
S,s      Stretch, squish color axis range
)";
  const auto update = [&]()
  {
    Eigen::VectorXd Z = V.col(1);
    Eigen::MatrixXd C;
    igl::parula(Z,-scale,scale,C);
    viewer.data.set_colors(C);
  };
  viewer.callback_key_pressed = 
    [&](igl::viewer::Viewer &, unsigned int key, int mod)
  {
    switch(key)
    {
      case 'S':
      case 's':
        scale *= key=='S' ? 2.0 : 0.5;
        std::cout<<"Color axis range: ["<<-scale<<","<<scale<<"]"<<std::endl;
        break;
      default:
        return false;
    }
    update();
    return true;
  };

  viewer.data.set_mesh(V,F);
  update();
  viewer.core.show_lines = false;
  viewer.data.face_based = true;
  viewer.launch();
  return EXIT_SUCCESS;
}
