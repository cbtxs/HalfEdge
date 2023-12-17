
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
    c.get_top();
    std::cout << "    " << cell_indices[c.index()] << "          " << (uint32_t)c.N<< "      ";
    for(uint8_t i = 0; i < c.N; i++)
      std::cout << node_indices[c.cell2node[i]->index()] << "   ";
    std::cout << " " << std::endl;
  }
  std::cout << " \n " << std::endl;

  std::cout << "cell2edge : " << std::endl;
  std::cout << " indices       NV        edge" << std::endl;
  for(auto & c : cell)
  {
    c.get_top();
    std::cout << "    " << cell_indices[c.index()] << "          " << (uint32_t)c.N<< "      ";
    for(uint8_t i = 0; i < c.N; i++)
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

void test_cut_mesh()
{
  CutUniformMesh mesh(0, 0, 0.1, 0.1, 10, 5);
  double n[4] = {0.11, 0.21, 0.82, 0.41};
  std::vector<uint32_t> idx0 = {0, 1};
  std::vector<std::vector<uint32_t>> idx;
  idx.push_back(idx0);

  std::vector<bool> fix = {false, false};
  mesh.cut_by_interface(n ,fix, idx);
  check_mesh(mesh);
  print(mesh);
  Figure fig("out", mesh.get_box());
  fig.draw_mesh(mesh, true);
  fig.draw_halfedge(mesh, true);
}

int main()
{
  //test_uniform_mesh();
  test_cut_mesh();
  //test_simple_mesh();
}
