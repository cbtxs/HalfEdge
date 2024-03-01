#ifndef CUT_MESH_APP
#define CUT_MESH_APP

//#include "figure.h"
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
  double * point_out1;
  int * halfedge_out1;
  int * inner_cell1;

  double * point_out2;
  int * halfedge_out2;

  int * idx0; /** 密网格单元在 0 号界面生成的网格单元中的编号 */
  int * idx1; /** 密网格单元在 1 号界面生成的网格单元中的编号 */

  int * N; /** NN1, NHE1, NC1, NN2, NHE2, NC2 */
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

  double * point_out1    = out.point_out1;    /** 坐标 */
  int    * halfedge_out1 = out.halfedge_out1; /** 半边 */
  int    * inner_cell1   = out.inner_cell1;   /** 内部单元 */

  double * point_out2    = out.point_out2;    /** 坐标 */
  int    * halfedge_out2 = out.halfedge_out2; /** 半边 */

  int    * idx0          = out.idx0;          /** 密网格单元在网格 0 中单元的编号 */
  int    * idx1          = out.idx1;          /** 密网格单元在网格 1 中单元的编号 */

  /**
   * meshptr2: cut 两次的网格
   * meshprt0 : 被第 0 个界面 cut 的网格
   * meshprt1 : 被第 1 个界面 cut 的网格
   */
  std::shared_ptr<Mesh> meshptr1 = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  std::shared_ptr<Mesh> meshptr2 = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  CutMeshAlg cut0(meshptr2);
  CutMeshAlg cut1(meshptr1);

  Interface iface0, iface1;
  generate_interface(point0, is_fixed_point0, segment0, NP0, NS0, iface0);
  generate_interface(point1, is_fixed_point1, segment1, NP1, NS1, iface1);

  cut0.cut_by_loop_interface(iface0);
  std::shared_ptr<Mesh> meshptr0 = std::make_shared<Mesh>(*meshptr2); 

  Interface iface2 = iface1;
  cut0.cut_by_loop_interface(iface1); 

  //auto & mesh = *meshptr2;
  //Figure fig("out", mesh.get_box());
  //fig.draw_mesh(mesh, true);
  //fig.draw_halfedge(mesh, true);
  //fig.draw_node(mesh, true);

  cut1.cut_by_loop_interface(iface2);

  meshptr0->update();
  meshptr1->update();
  meshptr2->update();

  /** 单元编号 */
  auto & cindex0 = *(meshptr0->get_cell_indices());
  auto & cindex1 = *(meshptr1->get_cell_indices());
  auto & cindex2 = *(meshptr2->get_cell_indices());
  std::function<bool(Cell &)> fun = [&](Cell & c)->bool 
  { 
    uint32_t cidx = cindex2[c.index()];
    Point p = c.inner_point();
    idx0[cidx] = cindex0[meshptr0->find_point(p, false)->index()];
    idx1[cidx] = cindex1[meshptr1->find_point(p, false)->index()];
    return true;
  };
  meshptr2->for_each_entity(fun);

  out.N[0] = meshptr1->number_of_nodes()*2;
  out.N[1] = meshptr1->number_of_halfedges()*6;
  out.N[2] = meshptr1->number_of_cells();

  out.N[3] = meshptr2->number_of_nodes()*2;
  out.N[4] = meshptr2->number_of_halfedges()*6;
  out.N[5] = meshptr2->number_of_cells();

  get_node(meshptr1, point_out1);
  get_halfedge(meshptr1, halfedge_out1);
  get_inner_cell(meshptr1, inner_cell1);

  get_node(meshptr2, point_out2);
  get_halfedge(meshptr2, halfedge_out2);
}

}

int main()
{
  MeshParameter mp{0.0, 0.0, 1.0, 1.0, 10, 10};
  double a = 0.421315;
  double b = 0.66124;
  double c = 0.71121301;
  double d = 0.84413012;

  //double point0[8] = {0.22131245, 0.21252151, 0.213515125, 0.6566125,  
  //      0.712341251235, 0.65, 0.713515125, 0.23};
  //double point1[8] = {0.22131245, 0.21252151, 0.213515125, 0.6566125,  
  //      0.712341251235, 0.65, 0.713515125, 0.23};
  double point0[8] = {a, c, a, d, b, d, b, c};
  double point1[8] = {a, c, a, d, b, d, b, c};
  for(uint32_t i = 1; i < 8; i+=2)
    point1[i] -= 0.002;

  bool is_fixed_point[4] = {1, 1, 1, 1};
  int segment[5] = {3, 2, 1, 0, 3};

  InterfaceParameter i0{point0, is_fixed_point, segment, 4, 5};
  InterfaceParameter i1{point1, is_fixed_point, segment, 4, 5};

  const static int N = 1000;
  double point_out0[N*2] = {0};
  int halfedge_out0[N*6] = {0};
  double point_out1[N*2] = {0};
  int halfedge_out1[N*6] = {0};
  int inner_cell[N] = {0};
  int idx0[N] = {0};
  int idx1[N] = {0};
  int NN[6] = {0};

  OutParameter out{point_out0, halfedge_out0, inner_cell, point_out1, halfedge_out1, idx0, idx1, NN};
  get_cut_mesh2(mp, i0, i1, out);
  for(int i = 0; i < NN[5]; i++)
    std::cout << i << " " << idx0[i] << " " << idx1[i] << std::endl;

}


#endif /* CUT_MESH_APP */ 
