#ifndef _UNIFORM_MESH_
#define _UNIFORM_MESH_

#include "entity.h"
#include "halfedge_mesh.h"
namespace HEM
{
class UniformMesh : public HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>
{
public:
  using Self = UniformMesh; 
  using Base = HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>;

  struct Parameter 
  {
    Parameter(double _ox, double _oy, double _hx, 
              double _hy, uint32_t _nx, uint32_t _ny): 
              orignx(_ox), origny(_oy), hx(_hx), 
              hy(_hy), nx(_nx), ny(_ny) 
    {}

    double orignx;
    double origny;
    double hx;
    double hy;
    uint32_t nx;
    uint32_t ny;
  };

public:
  UniformMesh(double orign_x, double orign_y, double hx, double hy, 
      uint32_t nx, uint32_t ny):Base(), param_(orign_x, orign_y, hx, hy, nx, ny)
  {
    uint32_t N = 0;
    double * node = new double[(nx+1)*(ny+1)*2]();
    for(uint32_t i = 0; i < nx+1; i++)
    {
      for(uint32_t j = 0; j < ny+1; j++)
      {
        node[N++] = i*hx;
        node[N++] = j*hy;
      }
    }
    N = 0;
    uint32_t * cell = new uint32_t[nx*ny*4]();
    for(uint32_t i = 0; i < nx; i++)
    {
      for(uint32_t j = 0; j < ny; j++)
      {
        cell[N]   = i*(ny+1)+j;
        cell[N+1] = cell[N]+(ny+1);
        cell[N+2] = cell[N+1]+1;
        cell[N+3] = cell[N]+1;
        N += 4;
      }
    }
    reinit(node, cell, (nx+1)*(ny+1), nx*ny, 4);
  }

  uint32_t find_point_in_uniform_mesh(const Point & p)
  {
    uint32_t x = floor((p.x-param_.orignx)/param_.hx);
    uint32_t y = floor((p.y-param_.origny)/param_.hy);
    return x*param_.ny + y;
  }

private:
  Parameter param_;
};

/**
 * @brief 被切割的一致网格 
 */
class CutUniformMesh : public UniformMesh
{
public:
  using Base = UniformMesh;

public:
  CutUniformMesh(double orign_x, double orign_y, double hx, 
      double hy, uint32_t nx, uint32_t ny):
    Base(orign_x, orign_y, hx, hy, nx, ny), cidx_(nx*ny, 0) { update(); }

  void update()
  {
    cidx_.resize(number_of_cells());
    auto & cell = get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point_in_uniform_mesh(c.barycentary());
      cidx_[idx] = &c;
    }
  }

  void cut_by_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface);

  /** 多个界面的情况 */
  void cut_by_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<std::vector<uint32_t>> & interface)
  {
    for(auto & iface : interface)
      cut_by_interface(point, is_fixed_point, iface);
  }

  Cell * find_point(const Point & p)
  {
    HalfEdge * h = cidx_[Base::find_point_in_uniform_mesh(p)]->halfedge();
    HalfEdge * start = h;
    while(true)
    {
      if(h->is_on_the_left(p))
        h = h->next();
      else
      {
        start = h->opposite();
        h = start->next();
      }
      if(h == start)
        break;
    }
    return h->cell();
  }
private:
  /** p0 在 c 中，找到线段 [p0, p1] 的交点 p 和相交的半边 */
  HalfEdge * _out_cell(Cell * c, const Point & p0, const Point & p1, Point & p)
  {
    bool flag = false;
    HalfEdge * h = c->halfedge(); 
    while(flag)
    {
      h = h->next();
      flag = intersection_point_of_two_segments(p0, p1, 
          h->node()->coordinate(), h->previous()->node()->coordinate(), p); 
    }
    return h;
  }

  void _cut_cell(Cell * c, const Point & p0, const Point & p1, HalfEdge * h0, HalfEdge * h1)
  {
    HalfEdge & h2 = add_halfedge();
    HalfEdge & h3 = add_halfedge();
    HalfEdge & h4 = add_halfedge();
    HalfEdge & h5 = add_halfedge();
  }

  /** 线段 [p0, p1] 与网格相交 */
  void _cut_by_segment(const Point & p0, const Point & p1, Cell * c0, Cell * c1)
  {
    Point p(0.0, 0.0);
    HalfEdge * h = _out_cell(c0, p0, p1, p);
    refine_halfedge(h, p);
    while(true)
    {
      if()
    }
  }

private:
  std::vector<Cell * > cidx_;
};

void CutUniformMesh::cut_by_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface)
{
  uint32_t N = interface.size();
  Point p0(point[2*interface[0]], point[2*interface[0]+1]);
  Point p1(point[2*interface[1]], point[2*interface[1]+1]);
  Cell * c0 = find_point(p0);
  Cell * c1 = find_point(p1); 

  while()
    if(c0==c1)
      p1 = Point(point[2*interface[1]], point[2*interface[1]+1]);
    else

  }

}












}
#endif /* _UNIFORM_MESH_ */ 
