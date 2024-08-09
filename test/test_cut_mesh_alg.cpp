
#include "figure.h"
#include "uniform_mesh_cut.h"
#include "cut_mesh_algorithm0.h"
#include <cmath>
#include <iostream>
#include <chrono>

using namespace std::chrono;

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
  double hx = 0.005;
  double hy = 0.005;
  int nx = 400;
  int ny = 800;

  auto start0 = high_resolution_clock::now();

  std::shared_ptr<Mesh> meshptr = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  auto start1 = high_resolution_clock::now();
  Mesh & mesh = *meshptr;

  CutMeshAlg cutalg(meshptr);

  std::vector<Point> points;
  points.push_back({0.0000000001, 0.13});
  points.push_back({0.53, 0.41});
  points.push_back({0.37, 0.62});
  points.push_back({0.82, 0.57});
  points.push_back({0.62, 0.28});


  std::vector<bool> is_fixed_points(5, false);
  is_fixed_points[0] = true;
  Interface interface(points, is_fixed_points, meshptr, true);

  auto start3 = high_resolution_clock::now();
  cutalg.cut_by_loop_interface(interface);
  auto start4 = high_resolution_clock::now();

  std::cout << " draw mesh ..." << std::endl;

  Figure fig0("out0", mesh.get_box());
  fig0.draw_mesh<Mesh, uint8_t>(mesh, false, "is_in_the_interface");
  //fig0.draw_halfedge(mesh, true);
  //fig0.draw_node(mesh, true);
  for(uint32_t i = 0; i < points.size(); i++)
  {
    fig0.draw_line(points[i], points[(i+1)%points.size()], 0.0002);
  }

  auto stop = high_resolution_clock::now();
  auto duration0 = duration_cast<microseconds>(start1 - start0);
  auto duration1 = duration_cast<microseconds>(start3 - start1);
  auto duration2 = duration_cast<microseconds>(start4 - start3);
  auto duration3 = duration_cast<microseconds>(stop - start4);
  std::cout << "time to create mesh: " << duration0.count()/ 1000000.0 << " seconds" << std::endl; 
  std::cout << "time to create cutalg: " << duration1.count()/ 1000000.0 << " seconds" << std::endl;
  std::cout << "time to cut: " << duration2.count()/ 1000000.0 << " seconds" << std::endl;
  std::cout << "time to draw: " << duration3.count()/ 1000000.0 << " seconds" << std::endl;
  std::cout << "total time: " << (duration0.count() + duration1.count() + duration2.count() + duration3.count())/ 1000000.0 << " seconds" << std::endl;
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


