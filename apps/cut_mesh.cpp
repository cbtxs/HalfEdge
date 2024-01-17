#ifndef CUT_MESH_APP
#define CUT_MESH_APP


#include "geometry.h"
#include "entity.h"
#include "uniform_mesh.h"
#include "cut_mesh_algorithm.h"
#include <cmath>
#include <fstream>


using namespace HEM;
using Mesh = CutMesh<UniformMesh>;
using Node = Mesh::Node;
using HalfEdge = Mesh::HalfEdge;
using CutMeshAlg = CutMeshAlgorithm<UniformMesh>;
using Interface = typename CutMeshAlg::Interface;

extern "C"
{

void get_node(std::shared_ptr<Mesh> meshptr, double * point_out)
{
  auto & nindex = *(meshptr->get_node_indices());
  std::function<bool(Node &)> fun = [&nindex, &point_out](Node & n)->bool 
  { 
    uint32_t idx = nindex[n.index()];
    point_out[idx*2] = n.coordinate().x;
    point_out[idx*2+1] = n.coordinate().y;
    return true;
  };
  meshptr->for_each_entity(fun);
}

void get_inner_cell(std::shared_ptr<Mesh> meshptr, int * inner_cell)
{
  auto & is_in_the_interface = *(meshptr->get_cell_data<uint8_t>("is_in_the_interface"));
  auto & cindex = *(meshptr->get_cell_indices());
  std::function<bool(Cell &)> fun = [&cindex, &inner_cell, is_in_the_interface](Cell & c)->bool 
  { 
    uint32_t cidx = cindex[c.index()];
    inner_cell[cidx] = is_in_the_interface[c.index()];
    return true;
  };
  meshptr->for_each_entity(fun);
}

void get_halfedge(std::shared_ptr<Mesh> meshptr, int * halfedge_out)
{
  meshptr->update();
  auto & hindices = *(meshptr->get_halfedge_indices());
  auto & nindices = *(meshptr->get_node_indices());
  auto & eindices = *(meshptr->get_edge_indices());
  auto & cindices = *(meshptr->get_cell_indices());
  std::function<bool(HalfEdge & h)> fun = 
    [&halfedge_out, &hindices, &nindices, &eindices, &cindices](HalfEdge & h)->bool 
  { 
    uint32_t idx = hindices[h.index()];
    std::cout << idx << std::endl;
    halfedge_out[idx*6] = nindices[h.node()->index()];
    halfedge_out[idx*6+1] = cindices[h.cell()->index()];
    halfedge_out[idx*6+2] = hindices[h.next()->index()];
    halfedge_out[idx*6+3] = hindices[h.previous()->index()];
    halfedge_out[idx*6+4] = hindices[h.opposite()->index()];
    halfedge_out[idx*6+5] = eindices[h.edge()->index()];
    return true;
  };
  meshptr->for_each_entity(fun);
}

void generate_interface(double * point, 
                        bool * is_fixed_point, 
                        int * segment, 
                        int NP, 
                        int NS, 
                        Interface & iface)
{
  std::vector<Point> & points = iface.points;
  for(int i = 0; i < NP*2; i+=2)
    points.push_back(Point(point[i], point[i+1]));

  auto & is_fixed_points = iface.is_fixed_points;
  for(int i = 0; i < NP; i++)
    is_fixed_points.push_back(is_fixed_point[i]);

  auto & segments = iface.segments;
  for(int i = 0; i < NS; i++)
    segments.push_back(segment[i]);
}

struct MeshParameter
{
  double a, b, c, d;
  int nx, ny;
};

struct InterfaceParameter
{
  double * point; 
  bool * is_fixed_point;
  int * segment;
  int NP;
  int NS;
};

struct OutParameter
{
  double * point_out;
  int * halfedge_out;
  int * inner_cell;
  int * idx0;
  int * idx1;
  int NC0; 
  int NC1;
};

void get_cut_mesh(double a, double b, double c, double d, int nx, int ny, 
              double * point,
              bool * is_fixed_point,
              int * segment,
              int NP,
              int NS,
              double * point_out,
              int * halfedge_out, 
              int * N)
{
  double hx = (c-a)/nx, hy = (d-b)/ny;

  std::shared_ptr<Mesh> meshptr = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  CutMeshAlg cut(meshptr);

  Interface iface;
  generate_interface(point, is_fixed_point, segment, NP, NS, iface);

  cut.cut_by_loop_interface(iface);
  meshptr->update();
  get_node(meshptr, point_out);
  get_halfedge(meshptr, halfedge_out);
  N[0] = meshptr->number_of_nodes()*2;
  N[1] = meshptr->number_of_halfedges()*6;
}

void get_cut_mesh2(MeshParameter mp, 
                   InterfaceParameter i0, 
                   InterfaceParameter i1,
                   OutParameter out)
{
  double a = mp.a, b = mp.b, c = mp.c, d = mp.d;
  int nx = mp.nx, ny = mp.ny;
  double hx = (c-a)/nx, hy = (d-b)/ny;

  double * point0 = i0.point; 
  bool * is_fixed_point0 = i0.is_fixed_point; 
  int * segment0 = i0.segment; 
  int NP0 = i0.NP; 
  int NS0 = i0.NS;

  double * point1 = i1.point; 
  bool * is_fixed_point1 = i1.is_fixed_point; 
  int * segment1 = i1.segment; 
  int NP1 = i1.NP; 
  int NS1 = i1.NS;

  double * point_out = out.point_out;
  int * halfedge_out = out.halfedge_out;
  int * inner_cell = out.inner_cell;
  int * idx0 = out.idx0;
  int * idx1 = out.idx1;

  /**
   * meshptr0 : cut 两次的网格
   * meshprt1 : 被第 0 个界面 cut 的网格
   * meshprt2 : 被第 1 个界面 cut 的网格
   */
  std::cout << "asdasdasd" << std::endl;
  std::shared_ptr<Mesh> meshptr0 = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  std::shared_ptr<Mesh> meshptr2 = std::make_shared<Mesh>(*meshptr0);
  CutMeshAlg cut0(meshptr0);
  CutMeshAlg cut1(meshptr2);


  Interface iface0, iface1;
  generate_interface(point0, is_fixed_point0, segment0, NP0, NS0, iface0);
  generate_interface(point1, is_fixed_point1, segment1, NP1, NS1, iface1);

  cut0.cut_by_loop_interface(iface0);
  std::shared_ptr<Mesh> meshptr1 = std::make_shared<Mesh>(*meshptr0);

  cut0.cut_by_loop_interface(iface1); 
  cut1.cut_by_loop_interface(iface1);

  meshptr0->update();
  meshptr1->update();
  meshptr2->update();

  auto & cindex0 = *(meshptr0->get_cell_indices());
  auto & cindex1 = *(meshptr1->get_cell_indices());
  auto & cindex2 = *(meshptr2->get_cell_indices());

  std::function<bool(Cell &)> fun = [&](Cell & c)->bool 
  { 
    uint32_t cidx = cindex0[c.index()];
    Point p = c.inner_point();
    idx0[cidx] = cindex1[meshptr1->find_point(p)->index()];
    idx1[cidx] = cindex2[meshptr2->find_point(p)->index()];
    return true;
  };
  meshptr0->for_each_entity(fun);

  out.NC0 = meshptr0->number_of_cells();
  out.NC1 = meshptr2->number_of_cells();

  get_node(meshptr0, point_out);
  get_halfedge(meshptr0, halfedge_out);
  get_inner_cell(meshptr2, inner_cell);
}

}


#endif /* CUT_MESH_APP */ 
