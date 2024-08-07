#ifndef _CUT_MESH_ALGORITHM_
#define _CUT_MESH_ALGORITHM_

#include <array>
#include <deque>
#include <queue>
#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_set>
#include <assert.h>

#include "interface.h"

namespace HEM
{

/**
 * @brief 一个使用界面直接 Cut 网格的算法
 */
template<typename Mesh>
class CutMeshAlgorithm
{
public:
  using Cell = typename Mesh::Cell;
  using Edge = typename Mesh::Edge;
  using Node = typename Mesh::Node;
  using HalfEdge = typename Mesh::HalfEdge; 

  using Point = typename Mesh::Point;
  using Vector = typename Mesh::Vector;

  using Interface = InterfaceCut<Mesh>;
  using InterfacePoint = typename Interface::InterfacePoint;

  template<typename Data>
  using Array = typename Mesh::template Array<Data>;

  /**
   * @brief a struct to represent the intersection of segment and edge of mesh
   * @param point: the intersection point
   * @param e:     the edge that the intersection point is on
   * @param type:  the type of the intersection 
   *         0 if the two segments intersect at the vertex of the halfedge `h`,
   *         1 if the two segments intersect at a point on the halfedge `h`,
   */
  struct Intersection
  {
    Point point;
    HalfEdge * h;
    uint8_t type = 0;
    HalfEdge * out;
    HalfEdge * in;
    bool operator==(const Intersection & other)
    {
      return point.x == other.point.x && point.y == other.point.y;
    }
  };

public:
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

  /**
   * @brief 计算一个 segment 和网格的交点
   * @param p0: 起点
   * @param p1: 终点
   * @param intersections: 交点列表
   */
  void find_intersections_of_segment(const InterfacePoint & ip0, 
      const InterfacePoint & ip1, std::vector<Intersection> & intersections);

private:

  /**
   * @brief 找到同一个单元的界面角点
   * @param iface: 界面
   * @param corners: 返回角点列表
   */
  std::deque<std::array<uint32_t, 3> > _find_corners_of_same_cell(
      Interface & iface);

  void _find_halfedge(Intersection & a, Intersection & b, Cell * c = nullptr);

  void _find_out_and_in(Intersection & a, Intersection & b);

  /**
   * @brief 连接两个相同的 segment 上的交点
   * @param a: 交点 a
   * @param b: 交点 b
   */
  HalfEdge * _link_two_intersections_in_same_segment(Intersection & a, Intersection & b);

  /**
   * @brief 连接两个交点
   * @param a: 交点 a
   * @param b: 交点 b
   */
  HalfEdge * _link_two_intersections(Intersection & a, Intersection & b, Cell * c=nullptr);

private:
  std::shared_ptr<Mesh> mesh_;
};

/**
 * @brief 计算一个 segment 和网格的交点，算法流程：
 *     1. 从起点和终点所在的单元开始，将所有的边插入到一个队列中
 *     2. 从队列中取出一个边，计算它和 segment 的交点，如果有交点，将交点插入到交点列表中
 *     3. 将边的相邻单元插入到队列中
 *     4. 重复 2-3 步骤，直到队列为空
 *
 * @param p0: 起点
 * @param p1: 终点
 * @param intersections: 交点列表
 */
template<typename Mesh>
void CutMeshAlgorithm<Mesh>::find_intersections_of_segment(const InterfacePoint & ip0, 
        const InterfacePoint & ip1, std::vector<Intersection> & intersections)
{
  auto & p0 = ip0.point;
  auto & p1 = ip1.point;
  auto & geometry_utils = mesh_->geometry_utils();

  std::unordered_set<Edge *> edge_traversed;/** 已经遍历过的边 */
  std::queue<Edge *> edges;       /** 待遍历的边 */

  /** 插入一个单元的所有边 */
  auto insert_cell = [&edges, &edge_traversed](Cell * c)
  {
    auto c2e = c->adj_edges();
    for(Edge & e : c2e)
    {
      if(edge_traversed.find(&e) == edge_traversed.end())
      {
        edges.push(&e);
        edge_traversed.insert(&e);
      }
    }
  };

  for(auto & c : ip0.cells)
    insert_cell(c);

  for(auto & c : ip1.cells)
    insert_cell(c);

  if (ip0.type != 2)
    intersections.push_back(Intersection(p0, ip0.h, ip0.type));

  if (ip1.type != 2)
    intersections.push_back(Intersection(p1, ip1.h, ip1.type));

  Cell  * c[2];
  Point * p[2];
  HalfEdge * h[2];
  while(!edges.empty())
  {
    Edge * e = edges.front();
    edges.pop();
    e->vertices(p);   
    e->adj_cell(c);

    h[0] = e->halfedge();
    h[1] = h[0]->previous();

    Point ip; /**< 交点 */
    uint8_t flag = geometry_utils.relative_position_of_two_segments(p0, p1, *(p[0]), *(p[1]), ip);
    if (flag == 0) /**< 两条线段重合 */
    {
      intersections.push_back(Intersection(*(p[0]), h[1], 0));
      intersections.push_back(Intersection(*(p[1]), h[0], 0));
    }
    else if(flag == 1) /** 交于 e 的起点 */
    {
      intersections.push_back(Intersection(ip, h[1], 0));
    }
    else if(flag == 2) /** 交于 e 的终点 */
    {
      intersections.push_back(Intersection(ip, h[0], 0));
    }
    else if(flag == 3) /** 交于 e 的内部 */
    {
      intersections.push_back(Intersection(ip, h[0], 1));
    }

    edge_traversed.insert(e);

    if(flag != 4)
    {
      insert_cell(c[0]);
      insert_cell(c[1]);
    }
  }

  /** 移除重复的交点 */
  std::sort(intersections.begin(), intersections.end(), 
      [&p0](const Intersection & a, const Intersection & b)
      {
        return (a.point-p0).length() < (b.point-p0).length();
      });
  auto last = std::unique(intersections.begin(), intersections.end());
  intersections.erase(last, intersections.end());
}

/**
 * @brief 连接交点 a 和点 b，要求 a 和 b 都在顶点上, 而且是相邻的
 *        1. 判断两个点是否在同一条边上
 *        2. 如果在同一条边上，可以跳过
 *        3. 如果不在同一条边上，找到两个点的公共单元
 *        4. 在公共单元中连接两个点
 * @param a: 交点 a
 * @param b: 交点 b
 * @return: a 指向 b 的半边
 */
template<typename Mesh>
Mesh::HalfEdge * CutMeshAlgorithm<Mesh>::_link_two_intersections_in_same_segment(
    Intersection & a, Intersection & b)
{  
  assert(a.type == 0 && b.type == 0);

  /** 判断两个点是否在同一条边上 */
  Node * n0 = a.h->node();
  Node * n1 = b.h->node();
  auto n2e  = n0->adj_edges();
  for(auto & e_adj : n2e)
  {
    /** 两个交点在同一条边上 */
    if(e_adj.has_node(n1))
    {
      HalfEdge * h = e_adj.halfedge();
      if(h->node() == n1)
        return h;
      else
        return h->opposite();
    }
  } 
  /** 两个交点不在同一条边上 */
  return _link_two_intersections(a, b);
}

/**
 * @brief 连接两个交点 a 和 b，要求 a 和 b 都在顶点上
 */
template<typename Mesh>
Mesh::HalfEdge * CutMeshAlgorithm<Mesh>::_link_two_intersections(Intersection & a, Intersection & b, Cell * c)
{
  assert(a.type == 0 && b.type == 0);

  if(c==nullptr)
    mesh_->find_point((a.point+b.point)/2, c);

  Node * n0 = a.h->node();
  Node * n1 = b.h->node();
  HalfEdge * out = nullptr;
  HalfEdge * in  = nullptr;
  for(auto & h : c->adj_halfedges())
  {
    if(h.node() == n0)
      out = &h;
    if(h.node() == n1)
      in = &h;
  }
  mesh_->splite_cell(c, out, in);
  return out->next();
}

/**
 * @brief 找到同一个单元的界面角点
 * @param iface: 界面
 * @return: 返回角点列表
 *          每个元素是一个数组，
 *          第一个元素是角点的起始点编号，
 *          第二个元素是角点的终点编号,
 *          第三个元素是角点是否是固定点
 */
template<typename Mesh>
std::deque<std::array<uint32_t, 3> > CutMeshAlgorithm<Mesh>::_find_corners_of_same_cell(
    Interface & iface)
{
  auto & ipoints  = iface.points();

  uint32_t NP = ipoints.size(); /** 界面的点数 */
  std::deque<std::array<uint32_t, 3> > corners; 

  Cell * c = nullptr;
  std::array<uint32_t, 3> corn{0, 0, 0};
  for(uint32_t i = 0; i < NP; i++)
  {
    const auto & ip0 = ipoints[i];
    if(ip0.type == 2)
    {
      if (c != ip0.cells[0])
      {
        if (c != nullptr)
        {
          corners.push_back(corn);
          corn = {0, 0, 0};
        }
        c = ip0.cells[0];
        corn = {i, i, ip0.is_fixed_point};
      }
      else if(c == ip0.cells[0])
      {
        corn[1] = i;
        corn[2] = ip0.is_fixed_point || corn[2];
      }
    }
  }
  corners.push_back(corn);

  if(iface.is_loop_interface() && ipoints[0].type==2 && ipoints.back().type==2
    && ipoints[0].cells[0] == ipoints.back().cells[0])
  {
    corners.front()[0] = corners.back()[0]; 
    corners.front()[2] = corners.back()[2] || corners.front()[2];
    corners.pop_back();
  }
  else if(!iface.is_loop_interface() && ipoints[0].type==2)
  {
    corners.pop_front();
  }
  else if(!iface.is_loop_interface() && ipoints.back().type==2)
  {
    corners.pop_back();
  }
  return corners;
}

template<typename Mesh>
void CutMeshAlgorithm<Mesh>::_find_out_and_in(Intersection & a, Intersection & b)
{
  assert(a.type == 0 && b.type == 0);

  /** 判断两个点是否在同一条边上 */
  Node * n0 = a.h->node();
  Node * n1 = b.h->node();
  auto n2e  = n0->adj_edges();

  for(auto & e_adj : n2e)
  {
    /** 两个交点在同一条边上 */
    if(e_adj.has_node(n1))
    {
      HalfEdge * h = e_adj.halfedge();
      if(h->node() == n1)
        b.in = h;
      else
        b.in = h->opposite();
      a.out = nullptr;
      return;
    }
  } 

  /** 两个交点不在同一条边上 */
  _find_halfedge(a, b);
}

template<typename Mesh>
void CutMeshAlgorithm<Mesh>::_find_halfedge(Intersection & a, Intersection & b, Cell * c)
{
  assert(a.type == 0 && b.type == 0);

  /** 判断两个点是否在同一条边上 */
  Node * n0 = a.h->node();
  Node * n1 = b.h->node();

  if(c==nullptr)
    mesh_->find_point((a.point+b.point)/2, c);
  for(auto & h : c->adj_halfedges())
  {
    if(h.node() == n0)
      a.h = &h;
    if(h.node() == n1)
      b.h = &h;
  }
}


/**
 * @brieg 计算一个界面和网格的所有交点
 * @param iface: 界面
 * @param intersections: 交点列表
 */
template<typename Mesh>
void CutMeshAlgorithm<Mesh>::cut_by_loop_interface(Interface & iface)
{
  const auto & geometry_utils = mesh_->geometry_utils();

  std::vector<std::vector<Intersection> > intersections;
  auto & ipoints  = iface.points();

  uint32_t NP = ipoints.size(); /** 界面的点数 */
  uint32_t NS = NP - 1 + iface.is_loop_interface(); /** 界面的 segment 数 */
  intersections.resize(NS);

  /** 0. 找到同一个单元的界面角点 */
  auto corners = _find_corners_of_same_cell(iface);

  /** 1. 计算每个 segment 和网格的交点 */
  for(uint32_t i = 0; i < NS; i++)
  {
    uint32_t i0 = i; 
    uint32_t i1 = (i+1)%NP; 
    auto & ip0 = ipoints[i0];
    auto & ip1 = ipoints[i1];

    /** 1.1 找到 segment[i] 与网格的交点 */
    find_intersections_of_segment(ip0, ip1, intersections[i]);

    /** 1.2 边上的点加密 */
    for(auto & ins : intersections[i])
    {
      if(ins.type == 1)
      {
        mesh_->splite_halfedge(ins.h, ins.point);
        ins.h = ins.h->previous();
        ins.type = 0;
      }
    }

    /** 1.3 如果 ip0, ip1 在边上，那么它们的状态会因为 1.2 而改变*/
    if(ip0.type == 1)
    {
      ip0.h = intersections[i].front().h;
      ip0.type = 0;
    }
    if(ip1.type == 1)
    {
      ip1.h = intersections[i].back().h;
      ip1.type = 0;
    }
  }

  /** 2. 处理转折点 out in */
  uint32_t NC = corners.size();
  std::vector<bool> need_insert_point;
  need_insert_point.reserve(NC);
  for(auto & corn : corners)
  {
    auto & start = corn[0];
    auto & end   = corn[1];
    auto & fixed = corn[2];
    auto ins0 = intersections[(NP+start-1)%NP].back();
    auto ins1 = intersections[end].front();
    Cell * c = ipoints[start].cells[0];
    bool need_insert = false;
    /** ins0 和 ins1 一定要连接的，对于有固定点的情况，连接了以后加密连接边
     * 对于没有固定点的情况，要判断需不需要插入点*/
    if(fixed) /** 固定转折点的情况 */
    {
      _find_out_and_in_halfedge(ins0, ins1, c);
    }
    else /** 不是固定转折点的情况, 这种情况下要判断 ins0 和 ins1 要不要连接 
     *  1. 如果 ins0 和 ins1 组成的线段与 c 的边相交多于 4 次，那么 ins0 和 ins1 要连接，而且要插入一个点
     *  2. 
     *
     * */
    {
      uint8_t count = 0;
      Point * ep[2];
      auto c2e = c->adj_edges();
      for(const auto & e : c2e)
      {
        e.vertices(ep);
        uint8_t flag = geometry_utils.relative_position_of_two_segments(
            ins0.point, ins1.point, *(ep[0]), *(ep[1]));
        if(flag==0)
        {
          count = 5;
          break;
        }
        else if (flag != 4)
        {
          count++;
        }
      }
      assert (count >=4);
      HalfEdge * h = _find_halfedge(ins0, ins1, c);
      need_insert = count>4;
    }
  }
  for(auto & corn : corners)
  {
    auto & start = corn[0];
    auto & end   = corn[1];
    auto & fixed = corn[2];
    auto ins0 = intersections[(NP+start-1)%NP].back();
    auto ins1 = intersections[end].front();
    Cell * c = ipoints[start].cells[0];
    if(fixed) /** 固定转折点的情况 */
    {
      HalfEdge * h = _link_two_intersections(ins0, ins1, c);
      for(int i = start; i != (end+1)%NP; i = (i+1)%NP) //TODO
      {
        mesh_->splite_halfedge(h, ipoints[i].point);
      }
    }
    else
    {
      uint8_t count = 0;
      Point * ep[2];
      auto c2e = c->adj_edges();
      for(const auto & e : c2e)
      {
        e.vertices(ep);
        uint8_t flag = geometry_utils.relative_position_of_two_segments(
            ins0.point, ins1.point, *(ep[0]), *(ep[1]));
        if(flag==0)
        {
          count = 5;
          break;
        }
        else if (flag != 4)
        {
          count++;
        }
      }
      assert (count >=4);
      HalfEdge * h = _link_two_intersections(ins0, ins1, c);
      if(count>4)
      {
        Point p(0, 0);
        count = 0;
        for(uint32_t i = start; i != (end+1)%NP; i = (i+1)%NP)
        {
          count++;
          p += ipoints[i].point;
        }
        p = p/count;
        mesh_->splite_halfedge(h, p);
      }
    }
  }

  /** 3. 连接每个 segment 的交点*/
  for(uint32_t i = 0; i < NS; i++)
  {
    auto & ips = intersections[i];
    int Ni = ips.size();
    for(int j = 0; j < Ni-1; j++)
      _link_two_intersections_in_same_segment(ips[j], ips[j+1]);
  }
}




}








#endif /* _CUT_MESH_ALGORITHM_ */ 
