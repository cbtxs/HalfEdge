
#include "geometry.h"
#include "entity.h"
#include "halfedge_mesh.h"
#include "uniform_mesh.h"
#include "cut_mesh_algorithm.h"

using namespace HEM;
using Mesh = CutMesh<UniformMesh>;
using CutMeshAlg = CutMeshAlgorithm<UniformMesh>;
using Interface = typename CutMeshAlg::Interface;

void cut_mesh(double a, double b, double c, double d, uint32_t nx, uint32_t ny,
              double * point, 
              bool * is_fixed_point, 
              uint32_t * segment, 
              uint32_t NP, 
              uint32_t NS,
              double * point_out,
              uint32_t * halfedge_out)
{
  double hx = (c-a)/nx, hy = (d-b)/ny;

  std::shared_ptr<Mesh> meshptr = std::make_shared<Mesh>(0, 0, hx, hy, nx, ny);
  CutMeshAlg cut(meshptr);

  Interface iface;
  std::vector<Point> & points = iface.points;
  for(uint32_t i = 0; i < NP*2; i+=2)
    points.push_back(Point(point[i], point[i+1]));

  auto & is_fixed_points = iface.is_fixed_points;
  for(uint32_t i = 0; i < NP; i++)
    is_fixed_points.push_back(is_fixed_point[i]);

  auto & segments = iface.segments;
  for(uint32_t i = 0; i < NS; i++)
    segments.push_back(segment[i]);

  std::vector<Interface> ifaces;
  ifaces.emplace_back(iface);
  cut.cut_by_interfaces(ifaces);

}
