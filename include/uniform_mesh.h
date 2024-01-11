#ifndef _UNIFORM_MESH_
#define _UNIFORM_MESH_

#include <algorithm>

#include "geometry.h"
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
    Parameter(double _ox, double _oy, 
              double _hx, double _hy, 
              uint32_t _nx, uint32_t _ny): 
              orignx(_ox), origny(_oy), 
              hx(_hx), hy(_hy), 
              nx(_nx), ny(_ny) 
    {}

    double orignx;
    double origny;
    double hx;
    double hy;
    uint32_t nx;
    uint32_t ny;
  };

public:
  /**
   * @brief 构造函数 
   */
  UniformMesh(double orign_x, 
              double orign_y, 
              double hx, 
              double hy, 
              uint32_t nx, 
              uint32_t ny);

  uint32_t find_point(const Point & p)
  {
    uint32_t x = floor((p.x-param_.orignx)/param_.hx);
    uint32_t y = floor((p.y-param_.origny)/param_.hy);
    return x*param_.ny + y;
  }

  /**
   * @brief 返回网格中单元的尺寸
   */
  double cell_size()
  {
    return std::sqrt(param_.hx*param_.hy);
  }

private:
  Parameter param_;
};

UniformMesh::UniformMesh(double orign_x, 
                         double orign_y, 
                         double hx, 
                         double hy, 
                         uint32_t nx, 
                         uint32_t ny):
  Base(), param_(orign_x, orign_y, hx, hy, nx, ny)
{
  uint32_t N = 0;
  auto & node = *get_node(); 
  auto & halfedge = *get_halfedge(); 
  for(uint32_t i = 0; i < nx+1; i++)
  {
    for(uint32_t j = 0; j < ny+1; j++)
    {
      Node & n = add_node();
      n.reset(Point(i*hx, j*hy), N++, nullptr);
    }
  }

  N = 0;
  for(uint32_t i = 0; i < nx; i++)
  {
    for(uint32_t j = 0; j < ny; j++)
    {
      HalfEdge & h0 = add_halfedge();
      HalfEdge & h1 = add_halfedge();
      HalfEdge & h2 = add_halfedge();
      HalfEdge & h3 = add_halfedge();
      node[i*(ny+1)+j].set_halfedge(&h0);
      node[(i+1)*(ny+1)+j].set_halfedge(&h1);
      node[(i+1)*(ny+1)+j+1].set_halfedge(&h2);
      node[i*(ny+1)+j+1].set_halfedge(&h3);

      Cell & c = add_cell();
      h0.reset(&h1, &h3, &h0, &c, nullptr, &node[i*(ny+1)+j], 4*N);
      h1.reset(&h2, &h0, &h1, &c, nullptr, &node[(i+1)*(ny+1)+j], 4*N+1);
      h2.reset(&h3, &h1, &h2, &c, nullptr, &node[(i+1)*(ny+1)+j+1], 4*N+2);
      h3.reset(&h0, &h2, &h3, &c, nullptr, &node[i*(ny+1)+j+1], 4*N+3);
      c.reset(N++, &h0);
    }
  }

  N = 0;
  for(uint32_t j = 0; j < ny; j++)
  {
    Edge & e = add_edge();
    HalfEdge & h = halfedge[4*j];
    e.reset(N++, &h);
    h.set_edge(&e);
  }
  for(uint32_t i = 1; i < nx; i++)
  {
    for(uint32_t j = 0; j < ny; j++)
    {
      Edge & e = add_edge();
      HalfEdge & h0 = halfedge[4*(i*ny+j)];
      HalfEdge & h1 = halfedge[4*((i-1)*ny+j)+2];
      e.reset(N++, &h0);
      h0.set_edge(&e);
      h1.set_edge(&e);
      h0.set_opposite(&h1);
      h1.set_opposite(&h0);
    }
  }
  for(uint32_t j = 0; j < ny; j++)
  {
    Edge & e = add_edge();
    HalfEdge & h = halfedge[4*((nx-1)*ny+j)+2];
    e.reset(N++, &h);
    h.set_edge(&e);
  }

  /** x edge */
  for(uint32_t i = 0; i < nx; i++)
  {
    Edge & e = add_edge();
    HalfEdge & h = halfedge[4*i*ny+1];
    e.reset(N++, &h);
    h.set_edge(&e);
  }
  for(uint32_t i = 0; i < nx; i++)
  {
    for(uint32_t j = 1; j < ny; j++)
    {
      Edge & e = add_edge();
      HalfEdge & h0 = halfedge[4*(i*ny+j-1)+3];
      HalfEdge & h1 = halfedge[4*(i*ny+j)+1];
      e.reset(N++, &h0);
      h0.set_edge(&e);
      h1.set_edge(&e);
      h0.set_opposite(&h1);
      h1.set_opposite(&h0);
    }
  }
  for(uint32_t i = 0; i < nx; i++)
  {
    Edge & e = add_edge();
    HalfEdge & h = halfedge[4*(i*ny+ny-1)+3];
    e.reset(N++, &h);
    h.set_edge(&e);
  }
}

}
#endif /* _UNIFORM_MESH_ */ 
