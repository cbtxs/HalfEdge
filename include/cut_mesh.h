#ifndef CUT_MESH_H
#define CUT_MESH_H

#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <stack>
#include <string>

#include "irregular_array2d.h"
#include "geometry_utils.h"

namespace HEM
{

/**
 * @brief 被切割的网格
 * @param Mesh : 网格类型, 该网格必须是 HalfEdgeMesh 的子类
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

  using Point  = typename Mesh::Point;
  using Vector = typename Mesh::Vector;

  template<typename Data>
  using Array = typename Base::template Array<Data>;

  /** 判断点与边的关系 */
  enum PEStatus 
  {
    OUTSIDE_EDGE = 0,
    ON_THE_EDGE = 1, 
    ON_THE_START_POINT = 2,
    ON_THE_END_POINT = 3
  };

  /** 判断点与单元的关系 */
  enum PCStatus
  {
    OUTSIDE_CELL = 0,
    INSIDE_CELL = 1,
    ON_THE_EDGE_OF_CELL = 2,
    ON_THE_VERTEX_OF_CELL = 3
  };

public:
  /** 单元是否在界面内部的标记， 0 表示在外部，1 表示在内部 */
  std::shared_ptr<Array<uint8_t>> cellmarker; 


private:
  /** 
   * @brief 背景网格中单元的子单元
   */
  IrregularArray2D<Cell *> subcell_;

public:
  /**
   * @brief 默认构造函数
   */
  template<typename... Args>
  CutMesh(Args&&... args) : Base(std::forward<Args>(args)...) 
  {
    update_cidx();
    eps_ = Base::cell_size()*1e-4;
    cellmarker = this->template add_cell_data<uint8_t>("is_in_the_interface");
  }

  /**
   * @brief 更新 subcell_
   */
  void update_cidx()
  {
    uint32_t NC = Base::number_of_cells();
    auto & data = subcell_.get_data();
    auto & start = subcell_.get_start_pos();

    data.resize(NC);

    std::fill(start.begin(), start.end(), 0.0);
    start.resize(NC+1, 0);
    auto & cell = *Base::get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      start[idx+1]++;
    }

    for(uint32_t i = 1; i < NC+1; i++)
      start[i] += start[i-1];

    std::vector<uint32_t> I(NC, 0);
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      data[start[idx]+I[idx]] = &c;
      I[idx]++;
    }
  }

  Cell * find_point(Point & p, bool maybe_on_the_edge = true)
  {
    HalfEdge * h = nullptr;
    uint32_t idx = Base::find_point(p);
    auto cellc = subcell_[idx];
    for(auto & c : cellc)
    {
      bool flag = is_on_the_polygon(p, c, h, maybe_on_the_edge);
      if(flag==1)
        return c;
    }
    return nullptr;
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

};

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
      Cell * node2cell[32];
      uint32_t N = h->node()->adj_cell(node2cell);
      for(uint32_t ii = 0; ii < N; ii++)
        c1s.push_back(node2cell[ii]);
    }
    else /** 在边内部 */
    {
      c1s.push_back(h->cell());
      c1s.push_back(h->opposite()->cell());
    }
  }
  return h;
}




}

#endif // CUT_MESH_H
