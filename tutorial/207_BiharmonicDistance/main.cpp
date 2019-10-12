#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/exact_geodesic.h>
#include <igl/colormap.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/PI.h>
#include <iostream>
#include "tutorial_shared_path.h"
#include <imgui/imgui.h>
#include <igl/read_triangle_mesh.h>

// biharmonic distance
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/eigs.h>


template <typename DerivedV,typename DerivedF,typename DerivedD>
void biharmonic_distance(
  const Eigen::MatrixBase<DerivedV> &V,
  const Eigen::MatrixBase<DerivedF> &F,
  const int i,
  const double p,
  Eigen::PlainObjectBase<DerivedD> &D)
{
  Eigen::SparseMatrix<double> L,M;
  Eigen::MatrixXd eU; //eigen vectors 
  Eigen::VectorXd eS; //eigen values 

  igl::cotmatrix(V,F,L);
  igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_BARYCENTRIC,M);
  const size_t k = 10;

  if(!eigs(L,M,k+1,igl::EIGS_TYPE_SM,eU,eS))
  {
    std::cout<<"eigen decomp failed."<<std::endl;
  }

  Eigen::MatrixXd eSD = eS.asDiagonal();
  eSD = eSD.block(1,1,eSD.rows()-1,eSD.rows()-1);
  eU = eU.middleCols(1,eU.cols()-1);
  
  //biharmonic embedding
  Eigen::MatrixXd B = eU * eSD.array().abs().matrix().inverse();
  //create query array from source point
  D = (B.row(i).replicate(B.rows(),1) - B)
                          .array().square()
                          .matrix().rowwise().sum()
                          .cwiseSqrt();
  
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  igl::opengl::glfw::Viewer viewer;

  //very trusting...
  igl::read_triangle_mesh(argv[1], V, F);

  enum DistanceType { geodisic=0, biharmonic };
  enum ColorType { banded=0, gradient };
  static DistanceType distanceType = geodisic;
  static ColorType colorType = banded;

  const auto update_distance = [&](const int vid, const int distanceType)
  {
    Eigen::VectorXi VS,FS,VT,FT;
    // The selected vertex is the source
    VS.resize(1);
    VS << vid;
    // All vertices are the targets
    VT.setLinSpaced(V.rows(),0,V.rows()-1);
    Eigen::VectorXd d;

    switch(distanceType) {
      case geodisic:
        std::cout<<"Computing geodesic distance to vertex "<<vid<<"..."<<std::endl;
        igl::exact_geodesic(V,F,VS,FS,VT,FT,d);
        break;
      case biharmonic:
        std::cout<<"Computing biharmonic distance to vertex "<<vid<<"..."<<std::endl;
        biharmonic_distance(V,F,vid,2.0,d);
        break;
      default:
        std::cout << "invalid option...";
        break;
    }
    
    Eigen::MatrixXd C;
    const double strip_size = 0.01;
    // compute per vertex color
    switch(colorType) {
      case banded:
        d = (d/strip_size*igl::PI).array().sin().abs().eval();
        igl::colormap(igl::COLOR_MAP_TYPE_INFERNO,d,false,C);
        break;  
      case gradient:
        igl::colormap(igl::COLOR_MAP_TYPE_JET,d,true,C);
        break;

    }

    // Plot the mesh
    viewer.data().clear();
    viewer.data().add_points(V.row(vid),Eigen::RowVector3d(1,0,0));
    viewer.data().set_mesh(V, F);
    viewer.data().set_colors(C);
  };

  // Plot a distance when a vertex is picked
  viewer.callback_mouse_down =
  [&](igl::opengl::glfw::Viewer& viewer, int, int)->bool
  {
    int fid;
    Eigen::Vector3f bc;
    // Cast a ray in the view direction starting from the mouse position
    double x = viewer.current_mouse_x;
    double y = viewer.core().viewport(3) - viewer.current_mouse_y;
    if(igl::unproject_onto_mesh(
      Eigen::Vector2f(x,y),
      viewer.core().view,
      viewer.core().proj,
      viewer.core().viewport,
      V,
      F,
      fid,
      bc))
    {
      int max;
      bc.maxCoeff(&max);
      int vid = F(fid,max);
      update_distance(vid,distanceType);
      return true;
    }
    return false;
  };
  viewer.data().set_mesh(V,F);

  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Add content to the default menu window
  menu.callback_draw_viewer_menu = [&]()
  {
    // Expose an enumeration type
    ImGui::Combo("Distance Type", (int *)(&distanceType), "Geodisic\0Biharmonic\0\0");
    ImGui::Combo("Color Method", (int *)(&colorType), "Banded\0Gradient\0\0");

  };

  cout << "Click on mesh to define new source.\n" << std::endl;
  return viewer.launch();
}
