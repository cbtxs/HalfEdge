#ifndef _UNIFORM_MESH_
#define _UNIFORM_MESH_

#include "halfedge_mesh.h"
#include "halfedge_mesh_traits.h"

namespace HEM
{

template<int D>
class UniformMesh : public HalfEdgeMeshBase<DefaultHalfEdgeMeshTraits<D>>
{
public:
  using Self = UniformMesh; 
  using Traits = DefaultHalfEdgeMeshTraits<D>;
  using Base = HalfEdgeMeshBase<Traits>;

  using Node = typename Base::Node;
  using Edge = typename Base::Edge;
  using Cell = typename Base::Cell;
  using HalfEdge = typename Base::HalfEdge;
  using Point  = typename Base::Point;
  using Vector = typename Base::Vector;

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

  uint32_t find_point(const Point & p) const
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

  double number_of_blocks()
  {
    return param_.nx*param_.ny;
  }

private:
  Parameter param_;
};

template<int D>
UniformMesh<D>::UniformMesh(double orign_x, 
                         double orign_y, 
                         double hx, 
                         double hy, 
                         uint32_t nx, 
                         uint32_t ny):
  Base((nx+1)*(ny+1), 2*nx*ny + nx + ny, nx*ny, 4*nx*ny, 1e-5 
      ), param_(orign_x, orign_y, hx, hy, nx, ny)
{
  auto & node = *(this->get_node()); 
  auto & edge = *(this->get_edge());
  auto & cell = *(this->get_cell());
  auto & halfedge = *(this->get_halfedge()); 

  uint32_t NN = 0;
  for(uint32_t i = 0; i < nx+1; i++)
  {
    for(uint32_t j = 0; j < ny+1; j++)
    {
      Node & n = node[NN];
      n.reset(Point(i*hx+orign_x, j*hy+orign_y), NN++, nullptr);
    }
  }

  uint32_t NC = 0;
  for(uint32_t i = 0; i < nx; i++)
  {
    for(uint32_t j = 0; j < ny; j++)
    {
      HalfEdge & h0 = halfedge[4*NC]; 
      HalfEdge & h1 = halfedge[4*NC+1]; 
      HalfEdge & h2 = halfedge[4*NC+2]; 
      HalfEdge & h3 = halfedge[4*NC+3]; 

      Node & n0 = node[i*(ny+1)+j];
      Node & n1 = node[(i+1)*(ny+1)+j];
      Node & n2 = node[(i+1)*(ny+1)+j+1];
      Node & n3 = node[i*(ny+1)+j+1];

      n0.set_halfedge(&h0);
      n1.set_halfedge(&h1);
      n2.set_halfedge(&h2);
      n3.set_halfedge(&h3);

      Cell & c = cell[NC];
      h0.reset(&h1, &h3, &h0, &c, nullptr, &n0, 4*NC);
      h1.reset(&h2, &h0, &h1, &c, nullptr, &n1, 4*NC+1);
      h2.reset(&h3, &h1, &h2, &c, nullptr, &n2, 4*NC+2);
      h3.reset(&h0, &h2, &h3, &c, nullptr, &n3, 4*NC+3);

      c.reset(NC++, &h0);
    }
  }

  uint32_t NE = 0;
  for(uint32_t j = 0; j < ny; j++)
  {
    Edge & e = edge[NE]; 
    HalfEdge & h = halfedge[4*j];
    e.reset(NE++, &h);
    h.set_edge(&e);
  }
  for(uint32_t i = 1; i < nx; i++)
  {
    for(uint32_t j = 0; j < ny; j++)
    {
      Edge & e = edge[NE]; 
      HalfEdge & h0 = halfedge[4*(i*ny+j)];
      HalfEdge & h1 = halfedge[4*((i-1)*ny+j)+2];
      e.reset(NE++, &h0);
      h0.set_edge(&e);
      h1.set_edge(&e);
      h0.set_opposite(&h1);
      h1.set_opposite(&h0);
    }
  }
  for(uint32_t j = 0; j < ny; j++)
  {
    Edge & e = edge[NE]; 
    HalfEdge & h = halfedge[4*((nx-1)*ny+j)+2];
    e.reset(NE++, &h);
    h.set_edge(&e);
  }

  /** x edge */
  for(uint32_t i = 0; i < nx; i++)
  {
    Edge & e = edge[NE]; 
    HalfEdge & h = halfedge[4*i*ny+1];
    e.reset(NE++, &h);
    h.set_edge(&e);
  }
  for(uint32_t i = 0; i < nx; i++)
  {
    for(uint32_t j = 1; j < ny; j++)
    {
      Edge & e = edge[NE]; 
      HalfEdge & h0 = halfedge[4*(i*ny+j-1)+3];
      HalfEdge & h1 = halfedge[4*(i*ny+j)+1];
      e.reset(NE++, &h0);
      h0.set_edge(&e);
      h1.set_edge(&e);
      h0.set_opposite(&h1);
      h1.set_opposite(&h0);
    }
  }
  for(uint32_t i = 0; i < nx; i++)
  {
    Edge & e = edge[NE]; 
    HalfEdge & h = halfedge[4*(i*ny+ny-1)+3];
    e.reset(NE++, &h);
    h.set_edge(&e);
  }
  Base::update();
}

using UniformMesh2D = UniformMesh<2>;


}
#endif /* _UNIFORM_MESH_ */ 
