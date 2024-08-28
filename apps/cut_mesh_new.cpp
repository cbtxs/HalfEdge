#ifndef CUT_MESH_APP
#define CUT_MESH_APP

//#include "figure.h"

#include "uniform_mesh_cut.h"
#include "cut_mesh_algorithm0.h"
#include <cmath>
#include <memory>
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

extern "C"
{

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
    assert(c.area()>0.0);
    for(HalfEdge * h = c.halfedge(); h != c.halfedge()->previous(); h=h->next())
      assert(h->cell()==&c);
    for(HalfEdge * h = c.halfedge(); h != c.halfedge()->next(); h=h->previous())
      assert(h->cell()==&c);
  }
}



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
                        std::vector<Point> & points,
                        std::vector<bool> & is_fixed_points,
                        bool & is_loop_interface)
{
  for(int i = 0; i < NP*2; i+=2)
    points.push_back(Point(point[i], point[i+1]));

  for(int i = 0; i < NP; i++)
    is_fixed_points.push_back(is_fixed_point[i]);

  is_loop_interface = segment[0] == segment[NS-1];
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

  std::vector<Point> points;
  std::vector<bool> is_fixed_points;
  bool is_loop_interface;
  generate_interface(point, is_fixed_point, segment, NP, NS, points, is_fixed_points, is_loop_interface);
  Interface iface(points, is_fixed_points, meshptr, is_loop_interface);

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
   * meshprt0 : 被第 0 个界面 cut 的网格
   * meshprt1 : 被第 1 个界面 cut 的网格
   * meshptr2 : cut 两次的网格
   * 1. 
   */
  std::shared_ptr<Mesh> meshptr0 = std::make_shared<Mesh>(a, b, hx, hy, nx, ny);
  std::shared_ptr<Mesh> meshptr1 = std::make_shared<Mesh>(*meshptr0);
  CutMeshAlg cut0(meshptr0);
  CutMeshAlg cut1(meshptr1);

  std::vector<Point> points0, points1;
  std::vector<bool> is_fixed_points0, is_fixed_points1;
  bool is_loop_interface0, is_loop_interface1;
  generate_interface(point0, is_fixed_point0, segment0, NP0, NS0, points0,
      is_fixed_points0, is_loop_interface0);
  generate_interface(point1, is_fixed_point1, segment1, NP1, NS1, points1,
      is_fixed_points1, is_loop_interface1);

  /** iface0 的背景网格是 meshptr2 */
  Interface iface0(points0, is_fixed_points0, meshptr0, is_loop_interface0);
  Interface iface1(points1, is_fixed_points1, meshptr1, is_loop_interface1);

  cut0.cut_by_loop_interface(iface0);
  cut1.cut_by_loop_interface(iface1);

  std::shared_ptr<Mesh> meshptr2 = std::make_shared<Mesh>(*meshptr0); 
  CutMeshAlg cut2(meshptr2);

  Interface iface2(points1, is_fixed_points1, meshptr2, is_loop_interface1);
  cut2.cut_by_loop_interface(iface2); 

  //auto & mesh0 = *meshptr0;
  //Figure fig0("out0", mesh0.get_box());
  //fig0.draw_mesh(mesh0, true);
  //fig0.draw_halfedge(mesh0, true);
  //fig0.draw_node(mesh0, true);

  //auto & mesh1 = *meshptr1;
  //Figure fig1("out1", mesh1.get_box());
  //fig1.draw_mesh(mesh1, true);
  //fig1.draw_halfedge(mesh1, true);
  //fig1.draw_node(mesh1, true);

  //auto & mesh2 = *meshptr2;
  //Figure fig2("out2", mesh2.get_box());
  //fig2.draw_mesh(mesh2, true);
  //fig2.draw_halfedge(mesh2, true);
  //fig2.draw_node(mesh2, true);

  meshptr0->update();
  meshptr1->update();
  meshptr2->update();

  //check_mesh(*meshptr0);
  //check_mesh(*meshptr1);
  //check_mesh(*meshptr2);
  
  /** 单元编号 */
  auto & cindex0 = *(meshptr0->get_cell_indices());
  auto & cindex1 = *(meshptr1->get_cell_indices());
  auto & cindex2 = *(meshptr2->get_cell_indices());
  std::function<bool(Cell &)> fun = [&](Cell & c)->bool 
  { 
    uint32_t cidx = cindex2[c.index()];
    Point p = c.inner_point();
    //idx0[cidx] = cindex0[meshptr0->find_point(p)->index()];
    //idx1[cidx] = cindex1[meshptr1->find_point(p)->index()];
    Cell * c0 = meshptr0->find_point(p);
    Cell * c1 = meshptr1->find_point(p);
    idx0[cidx] = cindex0[c0->index()];
    idx1[cidx] = cindex1[c1->index()];
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

int test111()
{
  MeshParameter mp{-0.0, 0.0, 1, 1, 10, 10};
  double point0[] = {
       0.6959059882504989, 0.7400197040177318, 0.6908278512800098,
       0.7596326245946713, 0.6837915623240461, 0.7786311710226634,
       0.6748693232289165, 0.796820392506216 , 0.6641526882414552,
       0.8140136430189582, 0.6517516245385582, 0.8300344965444442,
       0.6377933838151373, 0.8447185574458237, 0.6224211965095325,
       0.8579151473875395, 0.6057928020653924, 0.8694888514989499,
       0.5880788303115269, 0.879320907914148 , 0.5694610505689641,
       0.8873104264294158, 0.5501305064517441, 0.8933754237732406,
       0.5302855555009154, 0.897453664865622 , 0.5101298337677426,
       0.8995033014342103, 0.4898701662322575, 0.8995033014342103,
       0.4697144444990847, 0.897453664865622 , 0.4498694935482559,
       0.8933754237732406, 0.4305389494310359, 0.8873104264294158,
       0.4119211696884731, 0.879320907914148 , 0.3942071979346075,
       0.8694888514989499, 0.3775788034904675, 0.8579151473875395,
       0.3622066161848627, 0.8447185574458237, 0.3482483754614418,
       0.8300344965444442, 0.3358473117585448, 0.8140136430189582,
       0.3251306767710835, 0.796820392506216 , 0.3162084376759539,
       0.7786311710226634, 0.3091721487199902, 0.7596326245946713,
       0.3040940117495011, 0.7400197040177318, 0.301026135321621 ,
       0.7199936643974861, 0.3               , 0.6997599999999997,
       0.301026135321621 , 0.6795263356025134, 0.3040940117495011,
       0.6595002959822678, 0.3091721487199902, 0.6398873754053281,
       0.3162084376759539, 0.620888828977336 , 0.3251306767710836,
       0.6026996074937835, 0.3358473117585447, 0.5855063569810413,
       0.3482483754614417, 0.5694855034555553, 0.3622066161848626,
       0.5548014425541757, 0.3775788034904675, 0.5416048526124599,
       0.3942071979346075, 0.5300311485010496, 0.4119211696884731,
       0.5201990920858515, 0.4305389494310359, 0.5122095735705836,
       0.4498694935482559, 0.5061445762267588, 0.4697144444990846,
       0.5020663351343775, 0.4898701662322573, 0.5000166985657891,
       0.5101298337677425, 0.5000166985657891, 0.5302855555009152,
       0.5020663351343775, 0.5501305064517442, 0.5061445762267589,
       0.5694610505689641, 0.5122095735705836, 0.5880788303115267,
       0.5201990920858514, 0.6057928020653924, 0.5300311485010495,
       0.6224211965095325, 0.5416048526124599, 0.6377933838151373,
       0.5548014425541756, 0.6517516245385582, 0.5694855034555552,
       0.6641526882414552, 0.5855063569810413, 0.6748693232289165,
       0.6026996074937836, 0.683791562324046 , 0.620888828977336 ,
       0.6908278512800098, 0.639887375405328 , 0.6959059882504989,
       0.6595002959822677, 0.698973864678379 , 0.6795263356025132,
       0.7               , 0.6997599999999995, 0.698973864678379 ,
       0.7199936643974862};

  double point1[124] = {
       0.6959059882504989, 0.7399997040177317, 0.6908278512800098,
       0.7596126245946713, 0.6837915623240461, 0.7786111710226634,
       0.6748693232289165, 0.7968003925062159, 0.6641526882414552,
       0.8139936430189582, 0.6517516245385582, 0.8300144965444441,
       0.6377933838151373, 0.8446985574458237, 0.6224211965095325,
       0.8578951473875395, 0.6057928020653924, 0.8694688514989499,
       0.5880788303115269, 0.8793009079141479, 0.5694610505689641,
       0.8872904264294158, 0.5501305064517441, 0.8933554237732406,
       0.5302855555009154, 0.8974336648656219, 0.5101298337677426,
       0.8994833014342103, 0.4898701662322575, 0.8994833014342103,
       0.4697144444990847, 0.8974336648656219, 0.4498694935482559,
       0.8933554237732406, 0.4305389494310359, 0.8872904264294158,
       0.4119211696884731, 0.8793009079141479, 0.3942071979346075,
       0.8694688514989499, 0.3775788034904675, 0.8578951473875395,
       0.3622066161848627, 0.8446985574458237, 0.3482483754614418,
       0.8300144965444441, 0.3358473117585448, 0.8139936430189582,
       0.3251306767710835, 0.7968003925062159, 0.3162084376759539,
       0.7786111710226634, 0.3091721487199902, 0.7596126245946713,
       0.3040940117495011, 0.7399997040177317, 0.301026135321621 ,
       0.7199736643974861, 0.3               , 0.6997399999999997,
       0.301026135321621 , 0.6795063356025134, 0.3040940117495011,
       0.6594802959822678, 0.3091721487199902, 0.6398673754053281,
       0.3162084376759539, 0.620868828977336 , 0.3251306767710836,
       0.6026796074937835, 0.3358473117585447, 0.5854863569810412,
       0.3482483754614417, 0.5694655034555552, 0.3622066161848626,
       0.5547814425541757, 0.3775788034904675, 0.5415848526124599,
       0.3942071979346075, 0.5300111485010496, 0.4119211696884731,
       0.5201790920858514, 0.4305389494310359, 0.5121895735705836,
       0.4498694935482559, 0.5061245762267588, 0.4697144444990846,
       0.5020463351343775, 0.4898701662322573, 0.4999966985657891,
       0.5101298337677425, 0.4999966985657891, 0.5302855555009152,
       0.5020463351343775, 0.5501305064517442, 0.5061245762267589,
       0.5694610505689641, 0.5121895735705836, 0.5880788303115267,
       0.5201790920858513, 0.6057928020653924, 0.5300111485010495,
       0.6224211965095325, 0.5415848526124599, 0.6377933838151373,
       0.5547814425541756, 0.6517516245385582, 0.5694655034555551,
       0.6641526882414552, 0.5854863569810412, 0.6748693232289165,
       0.6026796074937836, 0.683791562324046 , 0.620868828977336 ,
       0.6908278512800098, 0.639867375405328 , 0.6959059882504989,
       0.6594802959822677, 0.698973864678379 , 0.6795063356025132,
       0.7               , 0.6997399999999995, 0.698973864678379 ,
       0.7199736643974862};

  const int NP = 62;
  bool is_fixed_point[NP] = {0};
  int segment[NP+1] = {0};
  std::iota(segment, segment+NP+1, 0);
  segment[NP] = 0;

  InterfaceParameter i0{point0, is_fixed_point, segment, NP, NP+1};
  InterfaceParameter i1{point1, is_fixed_point, segment, NP, NP+1};

  const static int N = 100000;
  double * point_out0 = new double[N*2];
  int * halfedge_out0 = new int[N*6];
  double * point_out1 = new double[N*2];
  int * halfedge_out1 = new int[N*6];
  int inner_cell[N] = {0};
  int idx0[N] = {0};
  int idx1[N] = {0};
  int NN[6] = {0};

  OutParameter out{point_out0, halfedge_out0, inner_cell, 
    point_out1, halfedge_out1, idx0, idx1, NN};
  get_cut_mesh2(mp, i0, i1, out);

  delete[] point_out0;
  delete[] halfedge_out0;
  delete[] point_out1;
  delete[] halfedge_out1;
  return 0;
}


int main(int, char ** )
{
  test111();
  return 0;
}


}

#endif /* CUT_MESH_APP */ 
