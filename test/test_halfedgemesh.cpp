
#include <memory>
#include "entity.h"
#include "halfedge_mesh.h"
#include "figure.h"
#include "uniform_mesh.h"

using namespace HEM;
using Mesh = HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>;
using UMesh = UniformMesh;

void check_mesh(Mesh & m)
{
  uint32_t NC = m.number_of_cells();
  uint32_t NE = m.number_of_edges();
  uint32_t NN = m.number_of_nodes();
  uint32_t NHE = m.number_of_halfedges();
  std::cout << " number of nodes : " << NN << std::endl;
  std::cout << " number of cells : " << NC << std::endl;
  std::cout << " number of edges : " << NE << std::endl;
  std::cout << " number of halfedges : " << NHE << std::endl;

  auto & node = *(m.get_node());
  auto & edge = *(m.get_edge());
  auto & cell = *(m.get_cell());
  auto & halfedge = *(m.get_halfedge());

  for(auto & n : node)
  {
    assert(n.halfedge()->node()==&n);
  }

  for(auto & e : edge)
  {
    assert(e.halfedge()->edge()==&e);
    assert(e.halfedge()->opposite()->edge()==&e);
  }

  for(auto & c : cell)
  {
    for(HalfEdge * h = c.halfedge(); h != c.halfedge()->previous(); h=h->next())
      assert(h->cell()==&c);
    for(HalfEdge * h = c.halfedge(); h != c.halfedge()->next(); h=h->previous())
      assert(h->cell()==&c);
  }

  uint32_t NBE = 0;
  for(auto & h : halfedge)
  {
    assert(h.previous()->next()==&h);
    if(h.opposite() != &h)
      assert(h.previous()->node()==h.opposite()->node());
    else
      NBE++;
  }
}

void print(Mesh & m)
{
  m.update();
  uint32_t NC = m.number_of_cells();
  uint32_t NE = m.number_of_edges();
  uint32_t NN = m.number_of_nodes();
  uint32_t NHE = m.number_of_halfedges();
  std::cout << " number of nodes : " << NN << std::endl;
  std::cout << " number of cells : " << NC << std::endl;
  std::cout << " number of edges : " << NE << std::endl;
  std::cout << " number of halfedges : " << NHE << std::endl;

  auto & node = *(m.get_node());
  auto & edge = *(m.get_edge());
  auto & cell = *(m.get_cell());
  auto & halfedge = *(m.get_halfedge());

  auto & node_indices = *(m.get_node_indices());
  auto & edge_indices = *(m.get_edge_indices());
  auto & cell_indices = *(m.get_cell_indices());
  auto & halfedge_indices = *(m.get_halfedge_indices());

  std::cout << "node : " << std::endl;
  std::cout << " indices   halfedge    x       y" << std::endl;
  for(auto & n : node)
  {
    std::cout << "    ";
    std::cout << node_indices[n.index()] << "          " << halfedge_indices[n.halfedge()->index()];
    std::cout << "       " << n.coordinate().x << "       " << n.coordinate().y <<std::endl;
  }
  std::cout << " \n " << std::endl;

  std::cout << "edge : " << std::endl;
  std::cout << " indices   halfedge   start node   end node" << std::endl;
  for(auto & e : edge)
  {
    std::cout << "    " << edge_indices[e.index()] 
      << "          " << halfedge_indices[e.halfedge()->index()] 
      << "          " << node_indices[e.halfedge()->node()->index()]
      << "           " << node_indices[e.halfedge()->previous()->node()->index()]<<std::endl;
  }
  std::cout << " \n " << std::endl;

  std::cout << "cell : " << std::endl;
  std::cout << " indices   halfedge" << std::endl;
  for(auto & c : cell)
    std::cout << "    " << cell_indices[c.index()] << "         " << halfedge_indices[c.halfedge()->index()] <<std::endl;
  std::cout << " \n " << std::endl;

  std::cout << "cell2node : " << std::endl;
  std::cout << " indices       NV        node" << std::endl;
  for(auto & c : cell)
  {
    uint32_t N = c.get_top();
    std::cout << "    " << cell_indices[c.index()] << "          " << N<< "      ";
    for(uint8_t i = 0; i < N; i++)
      std::cout << node_indices[c.cell2node[i]->index()] << "   ";
    std::cout << " " << std::endl;
  }
  std::cout << " \n " << std::endl;

  std::cout << "cell2edge : " << std::endl;
  std::cout << " indices       NV        edge" << std::endl;
  for(auto & c : cell)
  {
    uint32_t N = c.get_top();
    std::cout << "    " << cell_indices[c.index()] << "          " << N<< "      ";
    for(uint8_t i = 0; i < N; i++)
      std::cout << edge_indices[c.cell2edge[i]->index()] << "   ";
    std::cout << " " << std::endl;
  }
  std::cout << " \n " << std::endl;

  std::cout << "halfedge : " << std::endl;
  std::cout << " indices   node   edge   cell   next   prev   oppo" << std::endl;
  for(auto & h : halfedge)
  {
    std::cout << "    " << halfedge_indices[h.index()] << "       " << 
      node_indices[h.node()->index()] << "      " << 
      edge_indices[h.edge()->index()] << "      "  << 
      cell_indices[h.cell()->index()] << "      " << 
      halfedge_indices[h.next()->index()]<< "      " <<
      halfedge_indices[h.previous()->index()]<< "      " <<
      halfedge_indices[h.opposite()->index()]<< " " << std::endl;
  }
}

void test_shared_ptr()
{
  std::shared_ptr<uint32_t> a = std::make_shared<uint32_t>(10);
  std::cout << a.use_count() << std::endl;
  std::vector<std::shared_ptr<uint32_t>> b;
  b.push_back(a);
  std::cout << a.use_count() << std::endl;
  b.clear();
  std::cout << a.use_count() << std::endl;

}

void test_simple_mesh()
{
  double node[8] = {0, 0, 1, 0, 1, 2, 0, 2}; 
  uint32_t cell[6] = {0, 1, 2, 0, 2, 3};

  Mesh mesh(node, cell, 4, 2, 3);
  check_mesh(mesh);
  print(mesh);
  Figure fig("out", mesh.get_box());
  fig.draw_mesh(mesh, true);
  fig.draw_halfedge(mesh, true);
}

void test_uniform_mesh()
{
  UniformMesh mesh(0, 0, 0.1, 0.1, 10, 5);
  check_mesh(mesh);
  print(mesh);
  Figure fig("out", mesh.get_box());
  fig.draw_mesh(mesh, true);
  fig.draw_halfedge(mesh, true);
}

void test_splite_halfedge()
{
  UniformMesh mesh(0, 0, 0.1, 0.1, 10, 5);
  auto & h = (*(mesh.get_halfedge()))[147];
  mesh.splite_halfedge(&h);
  auto & h0 = (*(mesh.get_halfedge()))[150];
  auto & h1 = (*(mesh.get_halfedge()))[201];
  auto & h3 = (*(mesh.get_halfedge()))[190];
  mesh.splite_cell(h0.cell(), &h0, &h1);
  mesh.splite_halfedge(&h3);

  std::cout <<  h0.cell()->area() << std::endl;
  std::cout <<  h.cell()->area() << std::endl;

  check_mesh(mesh);
  print(mesh);
  Figure fig("out", mesh.get_box());
  fig.draw_mesh(mesh, true);
  fig.draw_halfedge(mesh, true);
  fig.draw_node(mesh, true);
  std::cout <<  h1.cell()->area() << std::endl;
  std::cout <<  h3.cell()->area() << std::endl;
}

void test_cut_mesh()
{
  double a = 0, b = 0, c = 1, d = 1;
  uint32_t nx = 10, ny = 10;
  double hx = (c-a)/nx, hy = (d-b)/ny;
  CutUniformMesh mesh(0, 0, hx, hy, nx, ny);
  double n[10] = {0.45001, 0.25001, 0.45, 0.4, 0.75, 0.4, 0.75, 0.2, 0.45, 0.2};
  std::vector<uint32_t> idx0 = {0, 1, 2, 3, 4, 0};

  std::vector<std::vector<uint32_t>> idx;
  idx.push_back(idx0);

  std::vector<bool> fix = {false, true, true, true, true, false};
  mesh.cut_by_interface(n ,fix, idx);
  check_mesh(mesh);
  //print(mesh);
  Figure fig("out", mesh.get_box());
  fig.draw_mesh(mesh, false);
  //fig.draw_halfedge(mesh, true);
  //fig.draw_node(mesh, true);
}

void test_bird_mesh()
{
  double a = 0, b = 0, c = 4, d = 1;
  uint32_t nx = 80, ny = 20;
  double hx = (c-a)/nx, hy = (d-b)/ny;
  CutUniformMesh mesh(0, 0, hx, hy, nx, ny);
  double n[20] =  {3.117, 0.49, 2.906, 0.523, 2.25, 0.32,
         2.062, 0.807, 1.98, 0.333, 1.781, 0.847,
         1.753, 0.377, 1.003, 0.312, 2.438, 0.177, 2.875, 0.477};
  std::vector<uint32_t> idx0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
  std::vector<std::vector<uint32_t>> idx;
  idx.push_back(idx0);

  std::vector<bool> fix(10, true);
  mesh.cut_by_interface(n ,fix, idx);
  check_mesh(mesh);
  Figure fig("bird", mesh.get_box());
  fig.draw_mesh(mesh, false);
}

int main()
{
  //test_uniform_mesh();
  test_cut_mesh();
  //test_simple_mesh();
  //test_splite_halfedge();
  test_bird_mesh();
}
