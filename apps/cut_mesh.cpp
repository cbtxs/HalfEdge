#ifndef CUT_MESH_APP
#define CUT_MESH_APP

//#include "figure.h"
#include "geometry.h"
#include "entity.h"
#include "uniform_mesh.h"
#include "cut_mesh_algorithm.h"
#include <cmath>
#include <fstream>
#include <sstream>


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
  std::shared_ptr<Mesh> meshptr2 = std::make_shared<Mesh>(*meshptr1);
  CutMeshAlg cut0(meshptr2);
  CutMeshAlg cut1(meshptr1);

  std::cout << std::to_string(a) << " " << std::to_string(b) << " " << std::to_string(c)  << " " << std::to_string(d) << "\n ";
  std::cout << " nx :" << std::to_string(nx) << " ny : " << std::to_string(ny) << " \n";
  std::cout << " NP0 : " << std::to_string(NP0) << " NS0 " << std::to_string(NS0) << "\n"; 
  std::cout << " NP1 : " << std::to_string(NP1) << " NS1 " << std::to_string(NS1) << "\n"; 
  std::cout << "p0 : " << " \n";
  for (int i = 0; i < NP0*2;  i+=2) 
  {
    std::cout << " x : " << std::to_string(point0[i]);
    std::cout << " y : " << std::to_string(point0[i+1]) <<" \n" ;
  }

  std::cout << "p1 : " << " \n";
  for (int i = 0; i < NP0*2;  i+=2) 
  {
    std::cout << " x : " << std::to_string(point1[i]);
    std::cout << " y : " << std::to_string(point1[i+1]) <<" \n" ;
  }

  std::cout << " segment0 : " << std::endl;
  for (int i = 0; i < NS0;  i+=1) 
    std::cout << std::to_string(segment0[i]);
  std::cout << std::endl;

  std::cout << " segment1 : " << std::endl;
  for (int i = 0; i < NS0;  i+=1) 
    std::cout << std::to_string(segment1[i]);
  std::cout << std::endl;

  std::cout << " is_fixed_point0 : " << std::endl;
  for (int i = 0; i < NP0;  i+=1) 
    std::cout << std::to_string(is_fixed_point0[i]);
  std::cout << std::endl;

  std::cout << " is_fixed_point1 : " << std::endl;
  for (int i = 0; i < NP0;  i+=1) 
    std::cout << std::to_string(is_fixed_point1[i]);
  std::cout << std::endl;

  Interface iface0, iface1;
  generate_interface(point0, is_fixed_point0, segment0, NP0, NS0, iface0);
  generate_interface(point1, is_fixed_point1, segment1, NP1, NS1, iface1);

  cut0.cut_by_loop_interface(iface0);
  std::shared_ptr<Mesh> meshptr0 = std::make_shared<Mesh>(*meshptr2); 

  Interface iface2 = iface1;
  cut0.cut_by_loop_interface(iface1); 

  cut1.cut_by_loop_interface(iface2);


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

int test(int NNN, int test_time = 1)
{
  for (int i = 0; i < test_time; i++) 
  {
    double s0 = -0.05 + static_cast<double>(std::rand()) / RAND_MAX * 0.1;
    double s1 = -0.05 + static_cast<double>(std::rand()) / RAND_MAX * 0.1;
    MeshParameter mp{-0.0, -0.0, 0.04, 0.08, 10, 20};
    std::cout << "s0 : "<< s0 << " s1 : " << s1 << std::endl;

    double point0[8] = {0.01365, 0.05272157, 0.01365, 0.05907157, 0.02635, 
                        0.05907157, 0.02635, 0.05272157};
    double point1[8] = {0.01365, 0.05158477, 0.01365, 0.05793477, 0.02635,    
                        0.05793477, 0.02635, 0.05158477};

    bool is_fixed_point[4] = {1, 1, 1, 1};
    int segment[5] = {3, 2, 1, 0, 3};

    InterfaceParameter i0{point0, is_fixed_point, segment, 4, 5};
    InterfaceParameter i1{point1, is_fixed_point, segment, 4, 5};

    const static int N = 100000;
    double * point_out0 = new double[N*2];
    int * halfedge_out0 = new int[N*6];
    double * point_out1 = new double[N*2];
    int * halfedge_out1 = new int[N*6];
    int inner_cell[N] = {0};
    int idx0[N] = {0};
    int idx1[N] = {0};
    int NN[6] = {0};

    OutParameter out{point_out0, halfedge_out0, inner_cell, point_out1, halfedge_out1, idx0, idx1, NN};
    get_cut_mesh2(mp, i0, i1, out);
  }
  return 0;
}

int test111(char **argv)
{
  int NNN = std::stoi(argv[1]);
  clock_t t0 = clock();
  MeshParameter mp{-0.0, -0.0, 1.0, 1.0, NNN, NNN};
  //double a = 0.4;
  //double b = 0.6;
  //double c = 0.7;
  //double d = 0.9;

  //double point0[8] = {0.22131245, 0.21252151, 0.213515125, 0.6566125,  
  //      0.712341251235, 0.65, 0.713515125, 0.23};
  //double point1[8] = {0.22131245, 0.21252151, 0.213515125, 0.6566125,  
  //      0.712341251235, 0.65, 0.713515125, 0.23};

  //double point0[8] = {a, c, a, d, b, d, b, c};
  //double point1[8] = {a, c, a, d, b, d, b, c};
  //for(uint32_t i = 1; i < 8; i+=2)
  //  point1[i] -= 0.002;

  //bool is_fixed_point[4] = {1, 1, 1, 1};
  //int segment[5] = {3, 2, 1, 0, 3};

  double point0[] ={ 
 0.69908052,  0.68915589, 0.69793314, 0.69867875, 0.69633052, 0.70813564,
 0.69427633,  0.71750482, 0.6917753 , 0.72676474, 0.68883317, 0.7358941 ,
 0.68545673,  0.7448719 , 0.68165373, 0.7536775 , 0.67743292, 0.76229063,
 0.672804  ,  0.77069149, 0.66777764, 0.77886076, 0.66236538, 0.78677964,
 0.65657967,  0.79442993, 0.65043383, 0.80179402, 0.64394198, 0.80885498,
 0.63711906,  0.81559658, 0.62998077, 0.82200329, 0.62254351, 0.82806039,
 0.6148244 ,  0.83375395, 0.60684119, 0.83907087, 0.59861224, 0.84399893,
 0.59015648,  0.84852677, 0.58149336, 0.85264401, 0.5726428 , 0.85634115,
 0.56362516,  0.8596097 , 0.55446118, 0.86244215, 0.54517193, 0.86483197,
 0.53577879,  0.86677367, 0.52630336, 0.86826279, 0.51676743, 0.86929589,
 0.50719294,  0.86987061, 0.4976019 , 0.86998562, 0.48801637, 0.86964066,
 0.47845841,  0.86883652, 0.46894999, 0.86757504, 0.45951299, 0.86585914,
 0.45016912,  0.86369275, 0.44093985, 0.86108087, 0.43184642, 0.85802949,
 0.42290975,  0.85454564, 0.4141504 , 0.85063733, 0.40558849, 0.84631355,
 0.39724374,  0.84158424, 0.38913533, 0.83646028, 0.38128191, 0.83095346,
 0.37370155,  0.82507643, 0.36641167, 0.81884273, 0.35942906, 0.81226669,
 0.35276976,  0.80536342, 0.3464491 , 0.79814882, 0.34048161, 0.79063947,
 0.33488101,  0.78285265, 0.3296602 , 0.77480626, 0.32483117, 0.76651881,
 0.32040504,  0.75800937, 0.31639198, 0.7492975 , 0.31280122, 0.74040325,
 0.30964103,  0.73134706, 0.30691867, 0.72214978, 0.3046404 , 0.71283255,
 0.30281147,  0.7034168 , 0.30143607, 0.6939242 , 0.30051738, 0.68437656,
 0.30005751,  0.67479586, 0.30005751, 0.66520414, 0.30051738, 0.65562344,
 0.30143607,  0.6460758 , 0.30281147, 0.6365832 , 0.3046404 , 0.62716745,
 0.30691867,  0.61785022, 0.30964103, 0.60865294, 0.31280122, 0.59959675,
 0.31639198,  0.5907025 , 0.32040504, 0.58199063, 0.32483117, 0.57348119,
 0.3296602 ,  0.56519374, 0.33488101, 0.55714735, 0.34048161, 0.54936053,
 0.3464491 ,  0.54185118, 0.35276976, 0.53463658, 0.35942906, 0.52773331,
 0.36641167,  0.52115727, 0.37370155, 0.51492357, 0.38128191, 0.50904654,
 0.38913533,  0.50353972, 0.39724374, 0.49841576, 0.40558849, 0.49368645,
 0.4141504 ,  0.48936267, 0.42290975, 0.48545436, 0.43184642, 0.48197051,
 0.44093985,  0.47891913, 0.45016912, 0.47630725, 0.45951299, 0.47414086,
 0.46894999,  0.47242496, 0.47845841, 0.47116348, 0.48801637, 0.47035934,
 0.4976019 ,  0.47001438, 0.50719294, 0.47012939, 0.51676743, 0.47070411,
 0.52630336,  0.47173721, 0.53577879, 0.47322633, 0.54517193, 0.47516803,
 0.55446118,  0.47755785, 0.56362516, 0.4803903 , 0.5726428 , 0.48365885,
 0.58149336,  0.48735599, 0.59015648, 0.49147323, 0.59861224, 0.49600107,
 0.60684119,  0.50092913, 0.6148244 , 0.50624605, 0.62254351, 0.51193961,
 0.62998077,  0.51799671, 0.63711906, 0.52440342, 0.64394198, 0.53114502,
 0.65043383,  0.53820598, 0.65657967, 0.54557007, 0.66236538, 0.55322036,
 0.66777764,  0.56113924, 0.672804  , 0.56930851, 0.67743292, 0.57770937,
 0.68165373,  0.5863225 , 0.68545673, 0.5951281 , 0.68883317, 0.6041059 ,
 0.6917753 ,  0.61323526, 0.69427633, 0.62249518, 0.69633052, 0.63186436,
 0.69793314,  0.64132125, 0.69908052, 0.65084411, 0.69977   , 0.66041103,
 0.7       ,  0.67      , 0.69977   , 0.67958897};
 
 double point1[] ={ 
 0.69908052, 0.68715589, 0.69793314, 0.69667875, 0.69633052, 0.70613564,
 0.69427633, 0.71550482, 0.6917753 , 0.72476474, 0.68883317, 0.7338941 ,
 0.68545673, 0.7428719 , 0.68165373, 0.7516775 , 0.67743292, 0.76029063,
 0.672804  , 0.76869149, 0.66777764, 0.77686076, 0.66236538, 0.78477964,
 0.65657967, 0.79242993, 0.65043383, 0.79979402, 0.64394198, 0.80685498,
 0.63711906, 0.81359658, 0.62998077, 0.82000329, 0.62254351, 0.82606039,
 0.6148244 , 0.83175395, 0.60684119, 0.83707087, 0.59861224, 0.84199893,
 0.59015648, 0.84652677, 0.58149336, 0.85064401, 0.5726428 , 0.85434115,
 0.56362516, 0.8576097 , 0.55446118, 0.86044215, 0.54517193, 0.86283197,
 0.53577879, 0.86477367, 0.52630336, 0.86626279, 0.51676743, 0.86729589,
 0.50719294, 0.86787061, 0.4976019 , 0.86798562, 0.48801637, 0.86764066,
 0.47845841, 0.86683652, 0.46894999, 0.86557504, 0.45951299, 0.86385914,
 0.45016912, 0.86169275, 0.44093985, 0.85908087, 0.43184642, 0.85602949,
 0.42290975, 0.85254564, 0.4141504 , 0.84863733, 0.40558849, 0.84431355,
 0.39724374, 0.83958424, 0.38913533, 0.83446028, 0.38128191, 0.82895346,
 0.37370155, 0.82307643, 0.36641167, 0.81684273, 0.35942906, 0.81026669,
 0.35276976, 0.80336342, 0.3464491 , 0.79614882, 0.34048161, 0.78863947,
 0.33488101, 0.78085265, 0.3296602 , 0.77280626, 0.32483117, 0.76451881,
 0.32040504, 0.75600937, 0.31639198, 0.7472975 , 0.31280122, 0.73840325,
 0.30964103, 0.72934706, 0.30691867, 0.72014978, 0.3046404 , 0.71083255,
 0.30281147, 0.7014168 , 0.30143607, 0.6919242 , 0.30051738, 0.68237656,
 0.30005751, 0.67279586, 0.30005751, 0.66320414, 0.30051738, 0.65362344,
 0.30143607, 0.6440758 , 0.30281147, 0.6345832 , 0.3046404 , 0.62516745,
 0.30691867, 0.61585022, 0.30964103, 0.60665294, 0.31280122, 0.59759675,
 0.31639198, 0.5887025 , 0.32040504, 0.57999063, 0.32483117, 0.57148119,
 0.3296602 , 0.56319374, 0.33488101, 0.55514735, 0.34048161, 0.54736053,
 0.3464491 , 0.53985118, 0.35276976, 0.53263658, 0.35942906, 0.52573331,
 0.36641167, 0.51915727, 0.37370155, 0.51292357, 0.38128191, 0.50704654,
 0.38913533, 0.50153972, 0.39724374, 0.49641576, 0.40558849, 0.49168645,
 0.4141504 , 0.48736267, 0.42290975, 0.48345436, 0.43184642, 0.47997051,
 0.44093985, 0.47691913, 0.45016912, 0.47430725, 0.45951299, 0.47214086,
 0.46894999, 0.47042496, 0.47845841, 0.46916348, 0.48801637, 0.46835934,
 0.4976019 , 0.46801438, 0.50719294, 0.46812939, 0.51676743, 0.46870411,
 0.52630336, 0.46973721, 0.53577879, 0.47122633, 0.54517193, 0.47316803,
 0.55446118, 0.47555785, 0.56362516, 0.4783903 , 0.5726428 , 0.48165885,
 0.58149336, 0.48535599, 0.59015648, 0.48947323, 0.59861224, 0.49400107,
 0.60684119, 0.49892913, 0.6148244 , 0.50424605, 0.62254351, 0.50993961,
 0.62998077, 0.51599671, 0.63711906, 0.52240342, 0.64394198, 0.52914502,
 0.65043383, 0.53620598, 0.65657967, 0.54357007, 0.66236538, 0.55122036,
 0.66777764, 0.55913924, 0.672804  , 0.56730851, 0.67743292, 0.57570937,
 0.68165373, 0.5843225 , 0.68545673, 0.5931281 , 0.68883317, 0.6021059 ,
 0.6917753 , 0.61123526, 0.69427633, 0.62049518, 0.69633052, 0.62986436,
 0.69793314, 0.63932125, 0.69908052, 0.64884411, 0.69977   , 0.65841103,
 0.7       , 0.668     , 0.69977   , 0.67758897};
 double s = std::stof(argv[2]);
 double s1 = std::stof(argv[3]);
  for(uint32_t i = 1; i < 131*2; i+=2)
  {
    point1[i-1] -= s1;
    point1[i] -= s;
    point1[i-1] *= 0.2;
    point1[i] -= 0.2;
    point0[i-1] *= 0.2;
    point0[i] -= 0.2;
  }

  int NSS = 131;
  bool is_fixed_point[NSS] = {1};
  for (int i = 0; i < NSS; i++) 
  {
    is_fixed_point[i] = true;
  }
  int segment[NSS+1] = {0};
  for (uint32_t i = 0; i < 131; i++) 
  {
    segment[i] = i;
  }

  InterfaceParameter i0{point0, is_fixed_point, segment, NSS, NSS+1};
  InterfaceParameter i1{point1, is_fixed_point, segment, NSS, NSS+1};

  const static int N = 1000000;
  double * point_out0 = new double[N*2];
  int * halfedge_out0 = new int[N*6];
  double * point_out1 = new double[N*2];
  int * halfedge_out1 = new int[N*6];
  int * inner_cell = new int[N];
  int * idx0 = new int[N];
  int * idx1 = new int[N];
  int NN[6] = {0};

  OutParameter out{point_out0, halfedge_out0, inner_cell, point_out1, halfedge_out1, idx0, idx1, NN};
  get_cut_mesh2(mp, i0, i1, out);
  clock_t t1 = clock();
  std::cout << (double)(t1-t0)/CLOCKS_PER_SEC << std::endl;
  //for(int i = 0; i < NN[5]; i++)
  //  std::cout << i << " " << idx0[i] << " " << idx1[i] << std::endl;
  return 0;
}

int main(int, char ** argv)
{
  //int NNN = std::stoi(argv[1]);
  //int test_time = std::stoi(argv[2]);
  //test(NNN, test_time);
  ///test111(argv);
  test(1);
  return 0;
}


#endif /* CUT_MESH_APP */ 
