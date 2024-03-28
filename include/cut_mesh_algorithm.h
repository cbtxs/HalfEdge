#ifndef _CUT_MESH_ALGORITHM_
#define _CUT_MESH_ALGORITHM_

#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <stack>
#include <string>

#include "geometry.h"

namespace HEM
{

/**
 * @brief 被切割的网格
 */
template<typename Mesh>
class CutMesh : public Mesh
{
public:
  using Base = Mesh;
  using Cell = typename Mesh::Cell;
  using Edge = typename Mesh::Edge;
  using Node = typename Mesh::Node;
  using HalfEdge = typename Mesh::HalfEdge; 
  template<typename Data>
  using Array = typename Base::template Array<Data>;

  std::shared_ptr<Array<uint8_t>> i3f;

public:
  template<typename... Args>
  CutMesh(Args&&... args) : Base(std::forward<Args>(args)...) 
  {
    update_cidx();
    eps_ = Base::cell_size()*1e-4;
    i3f = this->template add_cell_data<uint8_t>("is_in_the_interface");
  }

  void update_cidx()
  {
    uint32_t NC = Base::number_of_cells();
    subcell_.resize(NC);

    std::fill(cidx_.begin(), cidx_.end(), 0.0);
    cidx_.resize(NC+1, 0);
    auto & cell = *Base::get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      cidx_[idx+1]++;
    }

    for(uint32_t i = 1; i < NC+1; i++)
      cidx_[i] += cidx_[i-1];

    std::vector<uint32_t> I(NC, 0);
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      subcell_[cidx_[idx]+I[idx]] = &c;
      I[idx]++;
    }
  }

  Cell * find_point(Point & p, bool maybe_on_the_edge = true)
  {
    HalfEdge * h = nullptr;
    uint32_t idx = Base::find_point(p);
    Cell * c = nullptr;
    for(uint32_t i = cidx_[idx]; i < cidx_[idx+1]; i++)
    {
      c = subcell_[i]; 
      bool flag = is_on_the_polygon(p, c, h, maybe_on_the_edge);
      if(flag==1)
        break;
    }
    return c;
  }

  double eps()
  {
    return eps_;
  }

  /**
   * @brief 使用射线法判断点 p 是否在多边形 c 中, 如果没有，就返回一个半边 rh,
   *   rh 是距离 p 最近的半边。
   */
  uint8_t is_on_the_polygon(Point & p, Cell * c, HalfEdge* & rh, bool maybe_on_the_edge);

  /**
   * @brief 判断 p0, p1 是不是同一个点 
   */
  bool is_same_point(const Point & p0, const Point & p1)
  {
    return (p0-p1).length()<eps_;
  }

  /**
   * @brief 判断 
   *   "通过连接 h0 和 h1 指向的顶点来加密单元会不会产生面积非常小的单元"
   */
  uint8_t is_can_be_splite(HalfEdge * h0, HalfEdge * h1);

  /**
   * @brief 计算两个线段 [p0, p1], [q0, q1] 的交点，交点为 (1-t)*p0 + t*p1。
   * @note 线段重合不算相交。
   */
  bool intersection_point_of_two_segments(
      const Point & p0, const Point & p1, const Point & q0, const Point & q1, double & t);

  /**
   * @brief 判断一个向量在第几象限 
   */
  uint8_t quadrant_of_vector(const Vector & v);

  /**
   * @brief 找到节点 n 周围的单元中，节点 n 上向量 v 所在的单元
   */
  HalfEdge * find_cell_by_vector_on_node(Node * n, const Vector & v);

  /**
   * @brief 判断 p 是否在半边 h 上
   * @note 注意！！！该操作会修改 p : 当点 p 距离某一条边很近时，
   *    会把 p 投影到这个边上.
   */
  uint8_t is_on_the_halfedge(Point & p, HalfEdge * h);

  /**
   * @brief 判断 p2 在线段 [p0, p1] 上
   */
  uint8_t is_collinear(const Point & p0, const Point & p1, const Point & p2)
  {
    auto v2 = p2 - p0;
    auto v1 = p1 - p0;
    double t = v2.dot(v1)/v1.dot(v1);
    Point pt = p0*(1-t) + p1*t;
    return ((pt-p2).length()<eps_) && (t>0) && (t<1);
  }

  /**
   * @brief 判断 p 是否在 c 的某条边上
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   */
  HalfEdge * is_on_the_edge_of_cell(Cell * c, Point & p);
  
  /**
   * @brief 获取点 p 周围的单元集合 c1s. 
   *   - 如果 p 在某个单元中，那么返回的 hp1 = nullptr。
   *   - 如果 p 在某个边上, 那么返回的 hp1 是 p 所在的半边。
   *   - 如果 p 和某个顶点重合，那么返回的 hp1 是指向与 p 重合的顶点的半边的。
   */
  HalfEdge * get_cell_of_point(Point & p, std::vector<Cell * > & c1s);


private:
  std::vector<Cell * > subcell_;
  std::vector<uint32_t> cidx_;
  double eps_;
};

/**
 * @brief 使用射线法判断点 p 是否在多边形 c 中, 如果没有，就返回一个半边 rh,
 *   rh 是距离 p 最近的半边。
 */
template<typename BaseMesh>
uint8_t CutMesh<BaseMesh>::is_on_the_polygon(Point & p, Cell * c, HalfEdge* & rh, bool maybe_on_the_edge)
{
  uint8_t flag = 0;
  double d = 1e10;
  for(HalfEdge * h = c->halfedge(); h != c->halfedge() || d>1e9; h = h->next())
  {
    if(maybe_on_the_edge && is_on_the_halfedge(p, h))
      return 1;
    Point p1 = h->node()->coordinate();
    Point p0 = h->previous()->node()->coordinate();
    if((p.x>p0.x) != (p.x>p1.x))
    {
      double ft = (p0.x-p.x)/(p0.x-p1.x);
      if(p.y < p0.y-ft*(p0.y-p1.y))
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
 * @brief 判断 
 *   "通过连接 h0 和 h1 指向的顶点来加密单元会不会产生面积非常小的单元"
 */
template<typename BaseMesh>
uint8_t CutMesh<BaseMesh>::is_can_be_splite(HalfEdge * h0, HalfEdge * h1)
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
  if((ca-nca)<eps_*eps_)
    return 0;
  else if(nca < eps_*eps_)
    return 1;
  else
    return 2;
}

/**
 * @brief 计算两个线段 [p0, p1], [q0, q1] 的交点，交点为 (1-t)*p0 + t*p1。
 * @note 线段重合不算相交。
 */
template<typename BaseMesh>
bool CutMesh<BaseMesh>::intersection_point_of_two_segments(
    const Point & p0, const Point & p1, const Point & q0, const Point & q1, double & t)
{
  Vector v0 = p1-p0;
  Vector v1 = q0-q1;
  Vector v2 = q0-p0;
  double l = v0.length();
  double v = v0.cross(v1);
  if(std::abs(v/l) < eps_)
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
template<typename BaseMesh>
uint8_t CutMesh<BaseMesh>::quadrant_of_vector(const Vector & v)
{
  uint8_t a = v.x<0;
  uint8_t b = v.y<0;
  return a+3*b-2*a*b;
}

/**
 * @brief 找到节点 n 周围的单元中，节点 n 上向量 v 所在的单元
 */
template<typename BaseMesh>
typename BaseMesh::HalfEdge * CutMesh<BaseMesh>::find_cell_by_vector_on_node(
    Node * n, const Vector & v)
{
  HalfEdge * h = n->halfedge();
  Vector rot = v.normalize();
  Vector v0 = h->next()->tangential().normalize();
  Vector v1 = (h->tangential().normalize())*(-1.0);
  uint8_t k0 = quadrant_of_vector(v0.rotate(rot.y, rot.x));
  uint8_t k1 = quadrant_of_vector(v1.rotate(rot.y, rot.x));
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
      k0 = quadrant_of_vector(v0.rotate(rot.y, rot.x));
    }
  }
  return h;
}

/**
 * @brief 判断 p 是否在半边 h 上
 * @note 注意！！！该操作会修改 p : 当点 p 距离某一条边很近时，
 *    会把 p 投影到这个边上.
 */
template<typename BaseMesh>
uint8_t CutMesh<BaseMesh>::is_on_the_halfedge(Point & p, HalfEdge * h)
{
  auto v = h->tangential();
  const auto & p0 = h->previous()->node()->coordinate();
  const auto & p1 = h->node()->coordinate();
  double l = v.length();
  double t = ((p-p0).dot(v))/(l*l);
  if(is_same_point(p, p0+v*t) && t*l>-eps_/10.0 && (1-t)*l > -eps_/10.0)
  {
    p = p0+v*t;
    if(is_same_point(p, p0))
    {
      p=p0;
      return 1; /**< 在起点上 */
    }
    else if(is_same_point(p, p1))
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
template<typename BaseMesh>
typename BaseMesh::HalfEdge * CutMesh<BaseMesh>::is_on_the_edge_of_cell(
    Cell * c, Point & p)
{
  HalfEdge * rh = nullptr;
  uint8_t flag = 4;
  for(HalfEdge * h = c->halfedge(); h != c->halfedge() || flag == 4; h = h->next())
  {
    flag = is_on_the_halfedge(p, h);
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
template<typename BaseMesh>
typename BaseMesh::HalfEdge * CutMesh<BaseMesh>::get_cell_of_point(
    Point & p, std::vector<Cell * > & c1s)
{
  Cell * c1 = find_point(p);
  HalfEdge * h = is_on_the_edge_of_cell(c1, p);
  if(!h) /**< 单元内部 */
    c1s.push_back(c1);
  else
  {
    auto q1 = h->node()->coordinate();
    if(is_same_point(q1, p)) /**< 与终点重合 */
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


/**
 * @brief 一个使用界面直接 Cut 网格的算法
 */
template<typename BaseMesh>
class CutMeshAlgorithm
{
public:
  using Mesh = CutMesh<BaseMesh>;

  using Cell = typename BaseMesh::Cell;
  using Edge = typename BaseMesh::Edge;
  using Node = typename BaseMesh::Node;
  using HalfEdge = typename BaseMesh::HalfEdge; 

  template<typename Data>
  using Array = typename BaseMesh::template Array<Data>;

  /**
   * @brief 表示一个二维界面
   */
  struct Interface
  {
    std::vector<bool> is_fixed_points;
    std::vector<Point> points;
    std::vector<uint32_t> segments;

    bool is_loop_interface() 
    {
      return segments[0]==segments.back();
    }
  };

public:

  /**
   * @brief 默认构造函数
   */
  CutMeshAlgorithm() {}

  /**
   * @brief 构造函数
   */
  CutMeshAlgorithm(std::shared_ptr<Mesh> mesh): mesh_(mesh) {}

  void set_mesh(std::shared_ptr<Mesh> mesh)
  {
    mesh_ = mesh;
  }

  void cut_by_loop_interface(Interface & interfaces);

  void cut_by_non_loop_interface(Interface & interfaces);

  /**
   * @brief 多个界面 cut 网格
   */
  void cut_by_interfaces(std::vector<Interface> & interfaces)
  {
    for(auto & iface : interfaces)
    {
      if(iface.is_loop_interface())
        cut_by_loop_interface(iface);
      else
        cut_by_non_loop_interface(iface);
    }
  }

private:
  /**
   * @brief 判断 segment [p0, p1] 与从 start 到 end 之间的哪条半边相交, 交点为 p。
   */
  HalfEdge * _out_cell_0(HalfEdge * start, HalfEdge * end, const Point & p0, const Point & p1, Point & p);

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
      const Point & p0, const Point p1, 
      bool can_be_splite=false, const std::vector<Cell *> * c1s=nullptr);

  /** 
   * @brief 线段 [p0, p1] 与网格相交,  
   */
  HalfEdge * _cut_by_segment(const Point & p0, const Point & p1, HalfEdge * h0, 
      const std::vector<Cell *> & c1s);

  /**
   * @brief 获取内部单元
   */
  void get_inner_cell(Array<uint8_t> & is_in_the_interface);

  /** 
   * @brief 找到循环界面的第一个点，这个点是一个边上的点或与网格节点重合的点，
   *   如果没有就加一个。
   */
  uint32_t _find_first_point_in_loop_interface(Interface & interface, HalfEdge* & h0);


private:
  std::shared_ptr<Mesh> mesh_;
};

/**
 * @brief 判断 segment [p0, p1] 与从 start 到 end 之间的哪条半边相交, 交点为 p。
 */
template<typename BaseMesh>
typename BaseMesh::HalfEdge * CutMeshAlgorithm<BaseMesh>::_out_cell_0(
    HalfEdge * start, HalfEdge * end, const Point & p0, const Point & p1, Point & p)
{
  uint32_t ii = 0; /**< 防止 start == end */
  double t = 2, tempt = 0.0;
  HalfEdge * h1 = nullptr;
  double l = (p1-p0).length();
  for(HalfEdge * h = start; h != end || ii==0; h = h->next())
  {
    auto & q0 = h->previous()->node()->coordinate();
    auto & q1 = h->node()->coordinate();
    bool flag = mesh_->intersection_point_of_two_segments(p0, p1, q0, q1, tempt); 
    if(flag && t > tempt && tempt > mesh_->eps()/l)
    //if(flag && t < tempt && tempt < 1)
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
template<typename BaseMesh>
void CutMeshAlgorithm<BaseMesh>::_out_cell_1(
    Cell * c0, HalfEdge* & h0, HalfEdge* & h1, Point & p, 
    const Point & p0, const Point p1, bool can_be_splite, const std::vector<Cell *> * c1s)
{
  Vector v = p1 - p0;
  if(mesh_->is_same_point(p, h1->previous()->node()->coordinate()))
    h1 = h1->previous();
  if(mesh_->is_same_point(p, h1->node()->coordinate()))/**< 交到顶点上 */
  {
    /** 
     * 1. 加密产生面积非常小的单元时，不能加密
     * 2. 但是 can_be_splite 为 true 时可无视上一种情况。
     */
    p = h1->node()->coordinate();
    if(can_be_splite)
    {
      mesh_->splite_cell(c0, h0, h1);
      (*(mesh_->i3f))[h0->cell()->index()] = 1;
      (*(mesh_->i3f))[h1->cell()->index()] = 2;
    }
    else if(h0)
    {
      uint8_t flag = mesh_->is_can_be_splite(h0, h1);
      if(flag==2)
      {
        mesh_->splite_cell(c0, h0, h1);
        (*(mesh_->i3f))[h0->cell()->index()] = 1;
        (*(mesh_->i3f))[h1->cell()->index()] = 2;
      }
      else if(flag==1)
      {
        (*(mesh_->i3f))[h0->cell()->index()] = 1;
        (*(mesh_->i3f))[h1->opposite()->cell()->index()] = 2;
      }
      else if(flag==0)
      {
        (*(mesh_->i3f))[h0->opposite()->cell()->index()] = 1;
        (*(mesh_->i3f))[h1->cell()->index()] = 2;
      }
    }
    h0 = mesh_->find_cell_by_vector_on_node(h1->node(), v);
    /** 处理特殊情况 1 */
    if(c1s && std::find(c1s->begin(), c1s->end(), h0->cell()) != c1s->end())
    {
      const Point & p2 = h0->previous()->node()->coordinate();
      const Point & p3 = h0->next()->node()->coordinate();
      if(mesh_->is_collinear(p, p1, p2))
        h0 = h0->opposite()->previous();
      else if(mesh_->is_collinear(p, p1, p3))
        h0 = h0->next()->opposite();
    }
  }
  else /**< 没有交到顶点上 */
  {
    mesh_->splite_halfedge(h1, p);
    h1 = h1->previous();
    if(h0)
    {
      mesh_->splite_cell(c0, h0, h1);
      (*(mesh_->i3f))[h0->cell()->index()] = 1;
      (*(mesh_->i3f))[h1->cell()->index()] = 2;
    }
    h0 = h1->opposite()->previous();
  }
}

/** 
 * @brief 线段 [p0, p1] 与网格相交,  
 */
template<typename BaseMesh>
typename BaseMesh::HalfEdge * CutMeshAlgorithm<BaseMesh>::_cut_by_segment(
    const Point & p0, const Point & p1, HalfEdge * h0, const std::vector<Cell *> & c1s)
{
  Point p = p0;
  Cell * c0 = h0->cell();
  while(std::find(c1s.begin(), c1s.end(), c0)==c1s.end())
  { 
    HalfEdge * h1 = _out_cell_0(h0->next()->next(), h0, p, p1, p);
    _out_cell_1(c0, h0, h1, p, p0, p1, false, &c1s);
    //_out_cell_1(c0, h0, h1, p, p0, p1, false);
    c0 = h0->cell();
    p = h0->node()->coordinate();
  }
  return h0;
}

/** 
 * @brief 找到循环界面的第一个点，这个点是一个边上的点或与网格节点重合的点，
 *   如果没有就加一个。
 */
template<typename BaseMesh>
uint32_t CutMeshAlgorithm<BaseMesh>::_find_first_point_in_loop_interface(
    Interface & interface, HalfEdge* & h0)
{
  auto & segments = interface.segments;
  auto & points = interface.points;
  Point p0;
  Cell * c0 = nullptr;
  h0 = nullptr;
  uint32_t N = segments.size();
  for(uint32_t i = 0; i < N; i++)
  {
    Point p = points[segments[i]];
    std::vector<Cell * > cs;
    HalfEdge * h = mesh_->get_cell_of_point(p, cs);
    if(!h) /**< p 在单元内部 */
    {
      Cell * c = cs[0];
      if(c0==nullptr)
      {
        c0 = c; 
        p0 = p;
        segments.push_back(segments[i]);
      }
      else if(c==c0)
      {
        p0 = p;
        segments.push_back(segments[i]);
      }
      else
      {
        HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p, p0);
        _out_cell_1(c0, h0, h1, p0, p0, p);
        segments.push_back(points.size());
        points.push_back(h0->node()->coordinate());
        return i;
      }
    }
    else /**< p 在边上或者点上 */
    {
      if(c0==nullptr || std::find(cs.begin(), cs.end(), c0)!=cs.end())
      {
        auto q1 = h->node()->coordinate();
        auto q0 = h->previous()->node()->coordinate();
        if(mesh_->is_same_point(q0, p))
          h = h->previous();
        else if(!mesh_->is_same_point(q1, p))
        {
          mesh_->splite_halfedge(h, p);
          h = h->previous();
        }

        p0 = points[segments[i+1]];
        h0 = mesh_->find_cell_by_vector_on_node(h->node(), p0-p);
        segments.push_back(segments[i]);
        return i+1;
      }
      else
      {
        HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p, p0);
        _out_cell_1(c0, h0, h1, p0, p0, p);
        segments.push_back(points.size());
        points.push_back(h0->node()->coordinate());
        return i;
      }
    }
  }
  return 0;
}

template<typename BaseMesh>
void CutMeshAlgorithm<BaseMesh>::get_inner_cell(Array<uint8_t> & is_in_the_interface)
{
  /** 处理内部单元标记 */
  std::stack<typename Mesh::Cell*> inner_cell;
  auto & cell = *(mesh_->get_cell());
  for(auto & c : cell)
  {
    if(is_in_the_interface[c.index()]==1)
      inner_cell.push(&c);
  }
  while(!inner_cell.empty())
  {
    typename Mesh::Cell * c = inner_cell.top();
    inner_cell.pop();
    uint32_t N = c->get_top();
    for(uint8_t i = 0; i < N; i++)
    {
      typename Mesh::Cell * ci = c->cell2cell[i];
      if(is_in_the_interface[ci->index()]==0)
      {
        inner_cell.push(ci);
        is_in_the_interface[ci->index()]=1;
      }
    }
  }
  for(auto & c : cell)
    is_in_the_interface[c.index()] = is_in_the_interface[c.index()]==1;
}

template<typename BaseMesh>
void CutMeshAlgorithm<BaseMesh>::cut_by_loop_interface(Interface & interface)
{
  std::cout << "cuting..." << std::endl;

  /** 设置单元状态为在界面外部 */
  auto & is_in_the_interface = *(mesh_->i3f);
  for(auto & t : is_in_the_interface)
    t = 0;

  auto & points = interface.points;
  auto & segments = interface.segments;
  auto & is_fixed_points = interface.is_fixed_points;
  segments.pop_back();

  HalfEdge * h0 = nullptr;
  std::vector<Point> fpc, fpn;

  /** 获取第一个点的信息 */
  uint32_t start = _find_first_point_in_loop_interface(interface, h0);

  Cell * c0 = h0->cell();
  Point p0 = h0->node()->coordinate();

  uint32_t N = segments.size();
  for(uint32_t i = start; i < N; i++)
  {
    Point p1 = points[segments[i]];
    std::vector<Cell * > c1s;
    HalfEdge * hp1 = mesh_->get_cell_of_point(p1, c1s);
    if(!hp1) /**< p1 在单元内部 */
    {
      /** p1 和 p0 同一个单元 */
      if(std::find(c1s.begin(), c1s.end(), c0) != c1s.end()) 
      {
        if(is_fixed_points[segments[i]])
          fpc.push_back(p1);
        p0 = p1; 
        continue;
      }
      /** p1 和 p0 在不同单元 */
      else
      {
        if(is_fixed_points[segments[i]])
          fpn.push_back(p1);
        /** 转折 */
        Point p;
        HalfEdge * h1 = _out_cell_0(c0->halfedge(), c0->halfedge(), p0, p1, p);
        _out_cell_1(c0, h0, h1, p, p0, p1, !fpc.empty(), &c1s);
        for(auto & p : fpc)
          mesh_->splite_halfedge(h1->next(), p);
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
        _out_cell_1(c0, h0, h1, p, p0, p1, !fpc.empty(), &c1s);
        for(auto & p : fpc)
          mesh_->splite_halfedge(h1->next(), p);
        fpc.clear();
        /** 连线 */
        h0 = _cut_by_segment(h0->node()->coordinate(), p1, h0, c1s);
        c0 = h0->cell();
      }

      /** 最后一个单元的处理 */
      auto q1 = hp1->node()->coordinate();
      if(mesh_->is_same_point(q1, p1))
      {
        while(hp1->cell() != c0)
          hp1 = hp1->next_oppo();
      }
      else if (hp1->cell() != c0)
        hp1 = hp1->opposite();

      q1 = hp1->node()->coordinate();
      auto q0 = hp1->previous()->node()->coordinate();
      if(!mesh_->is_same_point(q0, p1) && !mesh_->is_same_point(q1, p1))
      {
        mesh_->splite_halfedge(hp1, p1);
        hp1 = hp1->previous();
      }

      if(!fpc.empty())
      {
        mesh_->splite_cell(c0, h0, hp1);
        (*(mesh_->i3f))[h0->cell()->index()] = 1;
        (*(mesh_->i3f))[hp1->cell()->index()] = 2;
      }
      else if(h0 != hp1)
      {
        uint8_t flag = mesh_->is_can_be_splite(h0, hp1);
        if(flag==2)
        {
          mesh_->splite_cell(c0, h0, hp1);
          (*(mesh_->i3f))[h0->cell()->index()] = 1;
          (*(mesh_->i3f))[hp1->cell()->index()] = 2;
        }
        else if(flag==1)
        {
          (*(mesh_->i3f))[h0->cell()->index()] = 2;
          (*(mesh_->i3f))[hp1->opposite()->cell()->index()] = 1;
        }
        else if(flag==0)
        {
          (*(mesh_->i3f))[h0->opposite()->cell()->index()] = 1;
          (*(mesh_->i3f))[hp1->cell()->index()] = 2;
        }
      }
      for(auto & p : fpc)
        mesh_->splite_halfedge(hp1->next(), p);

      if(i<N-1)
      {
        Point _p1 = points[segments[i+1]];
        h0 = mesh_->find_cell_by_vector_on_node(hp1->node(), _p1-p1);
        fpn.clear();
      }
    }
    c0 = h0->cell(); p0 = p1; fpc = fpn; fpn.clear();
  }
  get_inner_cell(is_in_the_interface);
  mesh_->update_cidx();
}

template<typename BaseMesh>
void CutMeshAlgorithm<BaseMesh>::cut_by_non_loop_interface(Interface & )
{
  //TODO
}

}

#endif /* _CUT_MESH_ALGORITHM_ */ 
