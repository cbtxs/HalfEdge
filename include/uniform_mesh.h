#ifndef _UNIFORM_MESH_
#define _UNIFORM_MESH_

#include <algorithm>
#include <queue>

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
    eps_ = std::min(hx, hy)*3e-2;
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

  void cut_by_loop_interface(std::vector<double> & point, 
      std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface);

  void cut_by_non_loop_interface(std::vector<double> & point, std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface);

  /** 多个界面的情况 */
  void cut_by_interface(std::vector<double> & point, std::vector<bool> & is_fixed_point, 
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

  Cell * find_point(Point & p)
  {
    HalfEdge * h = nullptr;
    Cell * c = cidx_[Base::find_point_in_uniform_mesh(p)];
    while(true)
    {
      bool flag = _is_on_the_polygon(p, c, h);
      if(flag==1)
        break;
      c = h->opposite()->cell();
    }
    return c;
  }

private:
  /**
   * @brief 使用射线法判断点 p 是否在多边形 c 中, 如果没有，就返回一个半边 rh,
   *   rh 是距离 p 最近的半边。
   */
  uint8_t _is_on_the_polygon(Point & p, Cell * c, HalfEdge* & rh)
  {
    uint8_t flag = 0;
    double d = 1e10;
    for(HalfEdge * h = c->halfedge(); h != c->halfedge() || d>1e9; h = h->next())
    {
      if(_is_on_the_halfedge(p, h))
        return 1;
      Point p1 = h->node()->coordinate();
      Point p0 = h->previous()->node()->coordinate();
      if((p.x>p0.x) != (p.x>p1.x))
      {
        double ft = (p0.x-p.x)/(p0.x-p1.x);
        if(p.y < p1.y+ft*(p0.y-p1.y))
          flag += 1; 
      }
      double l2 = (p0*0.5+p1*0.5-p).dot(p0*0.5+p1*0.5-p);
      if(d>l2)
      {
        rh = h;
        d = l2;
      }
    }
    return flag%2;
  }

  /**
   * @brief 判断 segment [p0, p1] 与从 start 到 end 之间的哪条半边相交, 交点为 p。
   */
  HalfEdge * _out_cell_0(HalfEdge * start, HalfEdge * end, 
      const Point & p0, const Point & p1, Point & p)
  {
    uint32_t ii = 0; /**< 防止 start == end */
    double t = -1, tempt = 0.0;
    HalfEdge * h1 = nullptr;
    for(HalfEdge * h = start; h != end || ii==0; h = h->next())
    {
      auto & q0 = h->previous()->node()->coordinate();
      auto & q1 = h->node()->coordinate();
      bool flag = _intersection_point_of_two_segments(p0, p1, q0, q1, tempt); 
      if(flag && t < tempt && tempt < 1)
      {
        t = tempt; h1 = h;
      }
      ii = 1;
    }
    p = p0*(1-t) + p1*t; /**< 获得 p */
    return h1;
  }

  /** 
   * @brief 现在处于单元 c0 中，要沿着 [p, v] 走出去，且我们已知路径可能会与
   *   半边 h1 相交，现在我们要做的是
   *   1. 判断是否交点是单元顶点
   *   2. 如果是顶点，先判断能不能 cut cell, 然后判断出射单元(通过遍历顶点周围单元实现)
   *   3. 如果不是顶点，那就以 p 为加密点加密半边 h1，并以连接入射点和出射
   *      点来 cut c
   * @param h0 : 指向的顶点是单元被 [p, v] 入射时加密半边得到的顶点,
   *   但是其也有可能是空的。
   * @param p0 : segment 的起点。
   * @param p1 : segment 的终点。
   * @param can_be_splite : bool 值，表示当前单元是否可以被强行加密，因为如果 
   *   c0 中有固定点，那么及时入射点和出射点在同一个半边，那我们也可以加密这个单元
   */
  void _out_cell_1(Cell * c0, HalfEdge* & h0, HalfEdge* & h1, Point & p, 
      const Vector & v, bool can_be_splite=false)
  {
    if(_is_same_point(p, h1->previous()->node()->coordinate()))
      h1 = h1->previous();
    if(_is_same_point(p, h1->node()->coordinate()))/**< 交到顶点上 */
    {
      /** 
       * 1. 加密产生面积非常小的单元时，不能加密
       * 2. 但是 can_be_splite 为 true 时可无视上一种情况。
       */
      p = h1->node()->coordinate();
      if((h0 && _is_can_be_splite(h0, h1)) || can_be_splite)
        splite_cell(c0, h0, h1);
      h0 = _find_cell_by_vector_on_node(h1->node(), v);
    }
    else /**< 没有交到顶点上 */
    {
      splite_halfedge(h1, p);
      h1 = h1->previous();
      if(h0)
        splite_cell(c0, h0, h1);
      h0 = h1->opposite()->previous();
    }
  }

  /** 
   * @brief 线段 [p0, p1] 与网格相交,  
   */
  HalfEdge * _cut_by_segment(const Point & p0, const Point & p1, HalfEdge * h0, 
      const std::vector<Cell *> & c1s)
  {
    Point p = p0;
    Cell * c0 = h0->cell();
    while(std::find(c1s.begin(), c1s.end(), c0)==c1s.end())
    { 
      HalfEdge * h1 = _out_cell_0(h0->next()->next(), h0, p0, p1, p);
      _out_cell_1(c0, h0, h1, p, p1-p0, false);
      c0 = h0->cell();
      p = h0->node()->coordinate();
    }
    return h0;
  }

  /** 
   * @brief 找到循环界面的第一个点，这个点是一个边上的点或与网格节点重合的点，
   *   如果没有就加一个。
   */
  uint32_t _find_first_point_in_loop_interface(std::vector<double> & point, 
      std::vector<uint32_t> & interface, HalfEdge* & h0)
  {
    Point p0;
    Cell * c0 = nullptr;
    h0 = nullptr;
    uint32_t N = interface.size();
    for(uint32_t i = 0; i < N; i++)
    {
      Point p(point[2*interface[i]], point[2*interface[i]+1]);
      std::vector<Cell * > cs;
      HalfEdge * h = _get_cell_of_point(p, cs);
      if(!h) /**< p 在单元内部 */
      {
        Cell * c = cs[0];
        if(c0==nullptr)
        {
          c0 = c; 
          p0 = p;
          interface.push_back(i);
        }
        else if(c==c0)
        {
          p0 = p;
          interface.push_back(i);
        }
        else
        {
          HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p, p0);
          _out_cell_1(c0, h0, h1, p0, p-p0);
          interface.push_back(point.size()/2);
          point.push_back(h0->node()->coordinate().x);
          point.push_back(h0->node()->coordinate().y);
          return i;
        }
      }
      else /**< p 在边上或者点上 */
      {
        if(c0==nullptr || std::find(cs.begin(), cs.end(), c0)!=cs.end())
        {
          auto q1 = h->node()->coordinate();
          auto q0 = h->previous()->node()->coordinate();
          if(_is_same_point(q0, p))
            h = h->previous();
          else if(!_is_same_point(q1, p))
          {
            splite_halfedge(h, p);
            h = h->previous();
          }

          p0 = Point(point[2*interface[i+1]], point[2*interface[i+1]+1]);
          h0 = _find_cell_by_vector_on_node(h->node(), p0-p);
          interface.push_back(i);
          return i+1;
        }
        else
        {
          HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p, p0);
          _out_cell_1(c0, h0, h1, p0, p-p0);
          interface.push_back(point.size()/2);
          point.push_back(h0->node()->coordinate().x);
          point.push_back(h0->node()->coordinate().y);
          return i;
        }
      }
    }
    return 0;
  }

  /** 下面是涉及的几个几何断言函数 */

  /**
   * @brief 判断 p0, p1 是不是同一个点 
   */
  bool _is_same_point(const Point & p0, const Point & p1)
  {
    return (p0-p1).length()<eps_;
  }

  /**
   * @brief 判断 
   *   "通过连接 h0 和 h1 指向的顶点来加密单元会不会产生面积非常小的单元"
   */
  bool _is_can_be_splite(HalfEdge * h0, HalfEdge * h1)
  {
    Cell * c = h0->cell();
    double ca = c->area();
    Point ce = c->barycenter();
    Vector v0 = h0->node()->coordinate()-ce;
    Vector v1 = h1->node()->coordinate()-ce;
    double nca = v1.cross(v0);
    for(HalfEdge * h = h0->next(); h != h1->next(); h = h->next())
    {
      v1 = v0;
      v0 = h->node()->coordinate()-ce;
      nca += v1.cross(v0);
    }
    nca*=0.5;
    return (ca-nca)>eps_*eps_ && nca > eps_*eps_;
  }

  /**
   * @brief 计算两个线段 [p0, p1], [q0, q1] 的交点，交点为 (1-t)*p0 + t*p1。
   * @note 线段重合不算相交。
   */
  bool _intersection_point_of_two_segments(const Point & p0, const Point & p1, 
      const Point & q0, const Point & q1, double & t)
  {
    Vector v0 = p1-p0;
    Vector v1 = q0-q1;
    Vector v2 = q0-p0;
    double l = v0.length();
    double v = v0.cross(v1);
    if(std::abs(v/l) < eps_/10)
    {
      return false;
    }
    t = (v2.cross(v1))/v;
    if(t*l<l+eps_ && t*l>-eps_)
    {
      l = v1.length();
      double s = (v0.cross(v2))/v;
      if(s*l<l+eps_ && s*l>-eps_)
        return true;
    }
    return false;
  }

  /**
   * @brief 判断一个向量在第几象限 
   */
  uint8_t _quadrant_of_vector(const Vector & v)
  {
    uint8_t a = v.x<0;
    uint8_t b = v.y<0;
    return a+3*b-2*a*b;
  }

  /**
   * @brief 找到节点 n 周围的单元中，节点 n 上向量 v 所在的单元
   */
  HalfEdge * _find_cell_by_vector_on_node(Node * n, const Vector & v)
  {
    HalfEdge * h = n->halfedge();
    Vector rot = v.normalize();
    Vector v0 = h->next()->tangential().normalize();
    Vector v1 = (h->tangential().normalize())*(-1.0);
    uint8_t k0 = _quadrant_of_vector(v0.rotate(rot.y, rot.x));
    uint8_t k1 = _quadrant_of_vector(v1.rotate(rot.y, rot.x));
    while(true) 
    {
      if(k0>k1 || (k0==k1 && (v0.cross(v1)<0.0)))
        return h;
      else
      {
        h = h->next()->opposite();
        v1 = v0;
        k1 = k0; 
        v0 = h->next()->tangential().normalize();
        k0 = _quadrant_of_vector(v0.rotate(rot.y, rot.x));
      }
    }
    return h;
  }

  /**
   * @brief 判断 p 是否在半边 h 上
   * @note 注意！！！该操作会修改 p : 当点 p 距离某一条边很近时，
   *    会把 p 投影到这个边上.
   */
  uint8_t _is_on_the_halfedge(Point & p, HalfEdge * h)
  {
    auto v = h->tangential();
    const auto & p0 = h->previous()->node()->coordinate();
    const auto & p1 = h->node()->coordinate();
    double l = v.length();
    double t = ((p-p0).dot(v))/(l*l);
    if(_is_same_point(p, p0+v*t) && t*l>-eps_ && (1-t)*l > -eps_)
    {
      p = p0+v*t;
      if(_is_same_point(p, p0))
      {
        p=p0;
        return 1; /**< 在起点上 */
      }
      else if(_is_same_point(p, p1))
      {
        p = p1;
        return 2; /**< 在终点上 */
      }
      else
        return 3; /**< 在边上 */
    }
    return 0; /**< 没有在边上 */
  }

  /**
   * @brief 判断 p 是否在 c 的某条边上
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   */
  HalfEdge * _is_on_the_edge_of_cell(Cell * c, Point & p)
  {
    HalfEdge * rh = nullptr;
    uint8_t flag = 4;
    for(HalfEdge * h = c->halfedge(); h != c->halfedge() || flag == 4; h = h->next())
    {
      flag = _is_on_the_halfedge(p, h);
      if(flag!=0)
      {
        if(flag==1)
          rh = h->previous(); 
        else
          rh = h;
        break;
      }
    }
    return rh;
  }
  
  /**
   * @brief 获取点 p 周围的单元集合 c1s. 
   *   - 如果 p 在某个单元中，那么返回的 hp1 = nullptr。
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   */
  HalfEdge * _get_cell_of_point(Point & p, std::vector<Cell * > & c1s)
  {
    Cell * c1 = find_point(p);
    HalfEdge * h = _is_on_the_edge_of_cell(c1, p);
    if(!h) /**< 单元内部 */
      c1s.push_back(c1);
    else
    {
      auto q1 = h->node()->coordinate();
      if(_is_same_point(q1, p)) /**< 与终点重合 */
      {
        uint32_t N = h->node()->get_top();
        for(uint32_t ii = 0; ii < N; ii++)
          c1s.push_back(h->node()->node2cell[ii]);
      }
      else /** 在边内部 */
      {
        c1s.push_back(h->cell());
        c1s.push_back(h->opposite()->cell());
      }
    }
    return h;
  }

private:
  std::vector<Cell * > cidx_;
  double eps_;
};

void CutUniformMesh::cut_by_loop_interface(std::vector<double> & point, 
      std::vector<bool> & is_fixed_point, 
      std::vector<uint32_t> & interface)
{
  std::cout << "cuting..." << std::endl;
  interface.pop_back();

  HalfEdge * h0 = nullptr;
  std::vector<Point> fpc, fpn;

  /** 获取第一个点的信息 */
  uint32_t start = _find_first_point_in_loop_interface(point, interface, h0);
  Cell * c0 = h0->cell();
  Point p0 = h0->node()->coordinate();

  uint32_t N = interface.size();
  for(uint32_t i = start; i < N; i++)
  {
    Point p1(point[2*interface[i]], point[2*interface[i]+1]);
    std::vector<Cell * > c1s;
    HalfEdge * hp1 = _get_cell_of_point(p1, c1s);
    if(!hp1) /**< p1 在单元内部 */
    {
      /** p1 和 p0 同一个单元 */
      if(std::find(c1s.begin(), c1s.end(), c0) != c1s.end()) 
      {
        if(is_fixed_point[interface[i]])
          fpc.push_back(p1);
        p0 = p1; 
        continue;
      }
      /** p1 和 p0 在不同单元 */
      else
      {
        if(is_fixed_point[interface[i]])
          fpn.push_back(p1);
        /** 转折 */
        Point p;
        HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p1, p);
        _out_cell_1(c0, h0, h1, p, p1-p0, !fpc.empty());
        for(auto & p : fpc)
          splite_halfedge(h1->next(), p);
        fpc.clear();
        /** 连线 */
        h0 = _cut_by_segment(h0->node()->coordinate(), p1, h0, c1s);
      }
    }
    else /**< p1 在边上或者点上 */
    {
      bool is_arrived = std::find(c1s.begin(), c1s.end(), c0) != c1s.end();
      if(!is_arrived) /**< 没有到达 */
      {
        /** 转折 */
        Point p;
        HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p1, p);
        _out_cell_1(c0, h0, h1, p, p1-p0, !fpc.empty());
        for(auto & p : fpc)
          splite_halfedge(h1->next(), p);
        fpc.clear();
        /** 连线 */
        h0 = _cut_by_segment(h0->node()->coordinate(), p1, h0, c1s);
        c0 = h0->cell();
      }

      /** 最后一个单元的处理 */
      auto q1 = hp1->node()->coordinate();
      if(_is_same_point(q1, p1))
      {
        while(hp1->cell() != c0)
          hp1 = hp1->next_oppo();
      }
      else if (hp1->cell() != c0)
        hp1 = hp1->opposite();

      q1 = hp1->node()->coordinate();
      auto q0 = hp1->previous()->node()->coordinate();
      if(!_is_same_point(q0, p1) && !_is_same_point(q1, p1))
      {
        splite_halfedge(hp1, p1);
        hp1 = hp1->previous();
      }
      if((_is_can_be_splite(h0, hp1) && h0 != hp1) || !fpc.empty())
        splite_cell(c0, h0, hp1); 
      for(auto & p : fpc)
        splite_halfedge(hp1->next(), p);

      if(i<N-1)
      {
        Point _p1(point[2*interface[i+1]], point[2*interface[i+1]+1]);
        h0 = _find_cell_by_vector_on_node(hp1->node(), _p1-p1);
        fpn.clear();
      }
    }
    c0 = h0->cell(); p0 = p1; fpc = fpn; fpn.clear();
  }
}

void CutUniformMesh::cut_by_non_loop_interface(std::vector<double> &, 
      std::vector<bool> & , 
      std::vector<uint32_t> & )
{
}






}
#endif /* _UNIFORM_MESH_ */ 
