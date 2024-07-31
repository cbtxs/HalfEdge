
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


int test_adj(int c_i, int n_j)
{
  double a = 0.0;
  double b = 0.0;
  double hx = 0.2;
  double hy = 0.2;
  int nx = 5;
  int ny = 4;
  std::shared_ptr<Mesh> meshptr = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  Mesh & mesh = *meshptr;

  auto & cell = *mesh.get_cell();
  auto & node = *mesh.get_node();

  auto & c_idx = *mesh.get_cell_indices();
  auto & n_idx = *mesh.get_node_indices();

  auto & c = cell[c_i];
  auto & n = node[n_j];


  std::cout << "The adj data of cell : " << c_idx[c.index()] << std::endl;
  auto c_adjs = c.adj_cells();
  auto e_adjs = c.adj_edges();
  auto n_adjs = c.adj_nodes();
  for (auto & c_adj : c_adjs)
    std::cout << "adj cell : " << c_idx[c_adj.index()] << std::endl;
  std::cout << " " << std::endl;
  for (auto & e_adj : e_adjs)
    std::cout << "adj edge : " << e_adj.index() << std::endl;
  std::cout << " " << std::endl;
  for (auto & n_adj : n_adjs)
    std::cout << "adj node : " << n_idx[n_adj.index()] << std::endl;

  std::cout << " " << std::endl;
  std::cout << "******************" << std::endl;
  std::cout << " " << std::endl;

  std::cout << "The adj data of node : " << n_idx[n.index()] << std::endl;
  auto n_adjs2 = n.adj_nodes();
  auto e_adjs2 = n.adj_edges();
  auto c_adjs2 = n.adj_cells();
  for (auto & n_adj : n_adjs2)
    std::cout << "adj node : " << n_idx[n_adj.index()] << std::endl;
  std::cout << " " << std::endl;
  for (auto & e_adj : e_adjs2)
    std::cout << "adj edge : " << e_adj.index() << std::endl;
  std::cout << " " << std::endl;
  for (auto & c_adj : c_adjs2)
    std::cout << "adj cell : " << c_idx[c_adj.index()] << std::endl;



  Figure fig0("out0", mesh.get_box());
  fig0.draw_mesh(mesh, true);
  //fig0.draw_halfedge(mesh, true);
  fig0.draw_node(mesh, true);
  fig0.draw_edge(mesh, true);

  return 0;
}

int main(int, char ** argv)
{
  int i = std::stoi(argv[1]);
  int j = std::stoi(argv[2]);

  //int test_time = std::stoi(argv[2]);
  //test(NNN, test_time);
  test_adj(i, j);
  //test(1);
  return 0;
}


