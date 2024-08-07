
#include "figure.h"
#include "uniform_mesh_cut.h"
#include "cut_mesh_algorithm0.h"
#include <cmath>
#include <iostream>


using namespace HEM;
using Mesh = UniformMeshCut<2>;
using Node = Mesh::Node;
using Edge = Mesh::Edge;
using Cell = Mesh::Cell;
using HalfEdge = Mesh::HalfEdge;
using Point = Mesh::Point;
using Vector = Mesh::Vector;

using CutMeshAlg = CutMeshAlgorithm<Mesh>;
using Interface = typename CutMeshAlg::Interface;
using Intersection = typename CutMeshAlg::Intersection;
using InterfacePoint = typename CutMeshAlg::InterfacePoint;


int test111()
{
  double a = 0.0;
  double b = 0.0;
  double hx = 0.2;
  double hy = 0.2;
  int nx = 5;
  int ny = 4;
  std::shared_ptr<Mesh> meshptr = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  Mesh & mesh = *meshptr;

  CutMeshAlg cutalg(meshptr);

  std::vector<Point> points;
  points.push_back({0.12, 0.13});
  points.push_back({0.43, 0.41});
  points.push_back({0.97, 0.52});
  points.push_back({0.82, 0.57});
  points.push_back({0.62, 0.78});


  std::vector<bool> is_fixed_points(5, false);
  is_fixed_points[0] = true;
  Interface interface(points, is_fixed_points, meshptr, true);

  cutalg.cut_by_loop_interface(interface);

  std::cout << " draw mesh ..." << std::endl;

  Figure fig0("out0", mesh.get_box());
  fig0.draw_mesh(mesh, true);
  fig0.draw_halfedge(mesh, true);
  fig0.draw_node(mesh, true);
  //for(uint32_t i = 0; i < points.size(); i++)
  //{
  //  fig0.draw_line(points[i], points[(i+1)%points.size()], 0.002);
  //}

  return 0;
}

int main(int, char ** argv[[maybe_unused]])
{
  //int NNN = std::stoi(argv[1]);
  //int test_time = std::stoi(argv[2]);
  //test(NNN, test_time);
  test111();
  //test(1);
  return 0;
}


