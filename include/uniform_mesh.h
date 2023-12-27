#ifndef _UNIFORM_MESH_
#define _UNIFORM_MESH_

#include <set>

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
    Parameter(double _ox, double _oy, double _hx, double _hy, uint32_t _nx, uint32_t _ny): 
              orignx(_ox), origny(_oy), hx(_hx), hy(_hy), nx(_nx), ny(_ny) {}

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
      double hy, uint32_t nx, uint32_t ny, double eps=1e-15):
    Base(orign_x, orign_y, hx, hy, nx, ny), cidx_(nx*ny, 0), eps_(eps)
  { 
    eps_ = std::min(hx, hy)*1e-2;
    update(); 
  }

  void update()
  {
    cidx_.resize(number_of_cells());
    auto & cell = *get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point_in_uniform_mesh(c.barycenter());
      cidx_[idx] = &c;
    }
  }

  void cut_by_loop_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface);

  void cut_by_non_loop_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface);

  /** 多个界面的情况 */
  void cut_by_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<std::vector<uint32_t>> & interface)
  {
    for(auto & iface : interface)
    {
      if(iface[0]==iface.back())
        cut_by_loop_interface(point, is_fixed_point, iface);
      else
        cut_by_non_loop_interface(point, is_fixed_point, iface);
    }
  }

  Cell * find_point(const Point & p)
  {
    HalfEdge * h = cidx_[Base::find_point_in_uniform_mesh(p)]->halfedge();
    HalfEdge * start = h;
    HalfEdge * endh = h; /**< 增强鲁棒性，处理 p 在单元的顶点上或边上的情况 */
    while(true)
    {
      if(h->is_on_the_left(p))
        h = h->next();
      else
      {
        start = h->opposite();
        h = start->next();
      }
      if(h == start || h==endh)
        break;
    }
    return h->cell();
  }

private:
  /**
   * @brief 找到节点 n 周围的单元中，节点 n 上向量 v 所在的单元
   */
  HalfEdge * _find_cell_by_vector_on_node(Node * n, const Vector & v)
  {
    HalfEdge * h = n->halfedge();
    bool flag0 = h->tangential().cross(v) > 0;
    bool flag1 = h->next()->tangential().cross(v) > 0; 
    while(!flag0 || !flag1)
    {
      flag0 = !flag1;
      h = h->next()->opposite();
      flag1 = h->next()->tangential().cross(v) > 0; 
    }
    return h;
  }

  /**
   * @brief 判断 p 是否在 c 的某条边上
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   * @note 默认 p 在 c 中, 因此 p 不在 c 中就会出 bug。
   */
  HalfEdge * _is_on_the_edge_of_cell(Cell * c, const Point & p)
  {
    HalfEdge * rh = nullptr;
    uint32_t i = 0;
    for(HalfEdge * h = c->halfedge(); h != c->halfedge() || i == 0; h = h->next())
    {
      auto v = h->tangential();
      const auto & p0 = h->previous()->node()->coordinate();
      double t = (p-p0).dot(v)/v.dot(v);
      if((p-p0-v*t).length()<eps_)
      {
        rh = h; break;
      }
      i = 1;
    }
    if(rh)
    {
      const auto & q0 = rh->previous()->node()->coordinate();
      if((q0-p).length() < eps_) /**< 与起点重合 */
        rh = rh->previous();
    }
    return rh;
  }
  
  /**
   * @brief 获取点 p 周围的单元集合 c1s. 
   *   - 如果 p 在某个单元中，那么返回的 hp1 = nullptr。
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   */
  HalfEdge * _get_cell_of_point(Point & p, std::set<Cell * > & c1s)
  {
    Cell * c1 = find_point(p);
    HalfEdge * hp1 = _is_on_the_edge_of_cell(c1, p);
    if(!hp1) /**< 单元内部 */
      c1s.insert(c1);
    else
    {
      auto q1 = hp1->node()->coordinate();
      if((q1-p).length() < eps_) /**< 与终点重合 */
      {
        uint32_t N = hp1->node()->get_top();
        for(uint32_t ii = 0; ii < N; ii++)
          c1s.insert(hp1->node()->node2cell[ii]);
      }
      else /** 在边内部 */
      {
        c1s.insert(hp1->cell());
        c1s.insert(hp1->opposite()->cell());
      }
    }
    return hp1;
  }


  /** 
   * @brief 现在处于 h0 的单元中，要沿着 [p0, p1] 走出去，路径可能会把单元 cut 
   * 掉。步骤如下:
   *    1. 找到单元所有与 [p0, p1] 相交的半边中，交点距离 p0 最近的半边 h1 以及交点 p;
   *    2. 判断是否交点是单元顶点
   *      2.1. 如果是顶点，先判断能不能 cut cell(如果入射点和出射点在同一个半边
   *        上那就不能 cut) 然后判断出射单元(通过遍历顶点周围单元实现)
   *      2.2. 如果不是顶点，那就以 p 为加密点加密半边 h1，并以连接入射点和出射
   *        点来 cut c
   */
  HalfEdge * _out_cell_plus(Cell * c0, HalfEdge* & h0, 
      const Point & p0, const Point & p1, bool can_be_splite=false)
  {
    /** 步骤 1 */ 
    uint32_t ii = 0;
    double t = -1, tempt = 0.0;
    HalfEdge * h1 = nullptr;
    for(HalfEdge * h = c0->halfedge(); h != c0->halfedge() || ii==0; h = h->next())
    {
      auto & q0 = h->previous()->node()->coordinate();
      auto & q1 = h->node()->coordinate();
      bool flag = intersection_point_of_two_segments(p0, p1, q0, q1, tempt, eps_); 
      if(flag && t < tempt && tempt < 1)
      {
        t = tempt; h1 = h;
      }
      ii = 1;
    }
    Point p = p0*(1-t) + p1*t; /**< 获得 p */

    /** 步骤 2 */
    if((p - h1->previous()->node()->coordinate()).length()<eps_)
      h1 = h1->previous();
    if((p - h1->node()->coordinate()).length()<eps_)/**< 交到顶点上 */
    {

      if((h0 && h1 != h0->next() && h0 != h1->next()) || can_be_splite)
        splite_cell(c0, h0, h1);
      h0 = _find_cell_by_vector_on_node(h1->node(), p1-p0);
    }
    else /**< 没有交到顶点上 */
    {
      splite_halfedge(h1, p);
      h1 = h1->previous();
      if(h0)
        splite_cell(c0, h0, h1);
      h0 = h1->opposite()->previous();
    }
    return h1;
  }

  /** 线段 [p0, p1] 与网格相交 */
  HalfEdge * _cut_by_segment(const Point & p0, const Point & p1, HalfEdge * h0, 
      const std::set<Cell *> & c1s)
  {
    Point p = p0;
    Cell * c0 = h0->cell();
    while(c1s.find(c0)==c1s.end())
    { 
      _out_cell_plus(c0, h0, p, p1);
      c0 = h0->cell();
      p = h0->node()->coordinate();
    }
    return h0;
  }

  void _transition(bool is_first, HalfEdge* & h1, HalfEdge* &h1f, 
      std::vector<uint32_t> & interface, std::vector<Point> & fpc, uint32_t i)
  {
    if(is_first)
    {
      h1f = h1;
      if(i>0)
        interface.push_back(i);
    }
    else
    {
      for(auto & p : fpc)
        splite_halfedge(h1->next(), p);
    }
  }

private:
  std::vector<Cell * > cidx_;
  double eps_;
};

void CutUniformMesh::cut_by_loop_interface(double * point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface)
{
  std::cout << "cuting..." << std::endl;
  uint32_t N = interface.size();
  Point p0(point[2*interface[0]], point[2*interface[0]+1]);
  Cell * c0 = find_point(p0);

  HalfEdge * h0 = nullptr, * h1 = nullptr, * h1f = nullptr;
  std::vector<Point> fpc, fpn;
  for(uint32_t i = 1; i < N; i++)
  {
    Point p1(point[2*interface[i]], point[2*interface[i]+1]);
    std::set<Cell * > c1s;
    HalfEdge * hp1 = _get_cell_of_point(p1, c1s);
    if(!hp1) /**< p1 在单元内部 */
    {
      if(c1s.find(c0) != c1s.end()) /**< p1 在单元内部，且和 p0 同一个单元 */
      {
        if(is_fixed_point[interface[i]])
          fpc.push_back(p1);
        p0 = p1; 
      }
      else
      {
        if(is_fixed_point[interface[i]])
          fpn.push_back(p1);
        /** 转折 */
        bool is_first = h0 == nullptr; /**< 判断是不是第一个单元 */
        h1 =  _out_cell_plus(c0, h0, p0, p1, !fpc.empty());
        _transition(is_first, h1, h1f, interface, fpc, i-1);
        /** 连线 */
        h0 = _cut_by_segment(h0->node()->coordinate(), p1, h0, c1s);
      }
    }
    else /**< p1 在边上或者点上 */
    {
      bool is_first = h0 == nullptr;
      bool is_arrived = c1s.find(c0) != c1s.end();
      if(!is_arrived) /**< 没有到达 */
      {
        /** 转折 */
        h1 =  _out_cell_plus(c0, h0, p0, p1, !fpc.empty());
        _transition(is_first, h1, h1f, interface, fpc, i-1);
        fpc.clear();
        /** 连线 */
        h0 = _cut_by_segment(h0->node()->coordinate(), p1, h0, c1s);
        c0 = h0->cell();
      }

      /** 最后一个单元的处理 */
      auto q1 = hp1->node()->coordinate();
      if((q1-p1).length()<eps_)
      {
        while(hp1->cell() != c0)
          hp1 = hp1->next_oppo();
      }
      if (hp1->cell() != c0)
        hp1 = hp1->opposite();
      if (hp1->cell() != c0)
        hp1 = hp1->opposite();
      auto q0 = hp1->previous()->node()->coordinate();
      if((q0-p1).length() > eps_ && (q1-p1).length() > eps_)
      {
        splite_halfedge(hp1, p1);
        hp1 = hp1->previous();
      }
      if(h0 && ((hp1 != h0->next() && h0 != hp1->next() && h0 != hp1) || !fpc.empty()))
        splite_cell(h0->cell(), h0, hp1); 

      if(is_arrived)
        _transition(is_first, hp1, h1f, interface, fpc, i-1);

      Point _p1(point[2*interface[i+1]], point[2*interface[i+1]+1]);
      h0 = _find_cell_by_vector_on_node(hp1->node(), _p1-p1);
      fpn.clear();
    }
    c0 = h0->cell(); p0 = p1; fpc = fpn; fpn.clear();
  }
  if((h1f != h0->next() && h0 != h1f->next()) || !fpc.empty())
    splite_cell(c0, h1f, h0);
  HalfEdge * h1next = h1f->next();
  for(auto & p : fpc)
    splite_halfedge(h1next, p);
}

void CutUniformMesh::cut_by_non_loop_interface(double *, std::vector<bool> & , 
      std::vector<uint32_t> & )
{
}






}
#endif /* _UNIFORM_MESH_ */ 
