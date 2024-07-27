
#include "figure.h"
#include "geometry.h"
#include "entity.h"
#include "uniform_mesh.h"
#include "uniform_mesh_cut.h"
#include "cut_mesh_algorithm0.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <numeric>


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
  points.push_back({0.0, 0.0});
  points.push_back({0.0, 0.4});

  std::vector<bool> is_fixed_points(2, false);
  Interface interface(points, is_fixed_points, meshptr, false);

  std::vector<std::vector<Intersection> > intersections;
  cutalg.find_intersections_of_interface(interface, intersections);
  for (auto & inter : intersections[0])
  {
    std::cout << inter.point.x << " " << inter.point.y << " " << inter.type << std::endl; 
  }

  Figure fig0("out0", mesh.get_box());
  fig0.draw_mesh(mesh, true);
  fig0.draw_halfedge(mesh, true);
  fig0.draw_node(mesh, true);

  return 0;
}

int main(int, char ** argv)
{
  //int NNN = std::stoi(argv[1]);
  //int test_time = std::stoi(argv[2]);
  //test(NNN, test_time);
  test111();
  //test(1);
  return 0;
}


