#ifndef _CUT_MESH_ALGORITHM_
#define _CUT_MESH_ALGORITHM_

#include <queue>
#include <locale>
#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_set>
#include <stdexcept>
#include <unordered_map>

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
  void find_intersections_of_segment(const InterfacePoint & ip0, const InterfacePoint & ip1, 
      std::vector<Intersection> & intersections);

  /**
   * @brief 计算一个 interface 和网格的所有交点
   * @param iface: 界面
   * @param intersections: 交点列表
   */
  void find_intersections_of_interface(const Interface & iface, 
      std::vector<std::vector<Intersection> > & intersections);

private:

  /**
   * @brief 将两个交点连接起来
   * @param a: 交点 a
   * @param b: 交点 b
   */
  void _link_two_intersections(Intersection & a, Intersection & b);

  /**
   * @brief 计算两个交点连线所在的单元
   * @param a: 交点 a
   * @param b: 交点 b
   * @return: 公共单元
   */
  Cell * _commom_cell_of_two_intersections(Intersection & a, Intersection & b);

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

  if (ip0.type != 3)
    intersections.push_back(Intersection(p0, ip0.h, ip0.type));

  if (ip1.type != 3)
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

template<typename Mesh>
void CutMeshAlgorithm<Mesh>::find_intersections_of_interface(const Interface & iface, 
  std::vector<std::vector<Intersection> > & intersections)
{
  auto & ipoints  = iface.points();
  auto & segments = iface.segments();

  uint32_t N = segments.size()-1;
  intersections.resize(N);

  for(uint32_t i = 0; i < N; i++)
  {
    uint32_t i0 = segments[i];
    uint32_t i1 = segments[i+1];
    auto & p0 = ipoints[i0].point;
    auto & p1 = ipoints[i1].point;

    find_intersections_of_segment(ipoints[i0], ipoints[i1], intersections[i]);
  }
}

template<typename Mesh>
void CutMeshAlgorithm<Mesh>::_link_two_intersections(
    Intersection & a, Intersection & b)
{  

}

template<typename Mesh>
typename Mesh::Cell * CutMeshAlgorithm<Mesh>::_commom_cell_of_two_intersections(
    Intersection & ip0, Intersection & ip1)
{
  HalfEdge * h0 = ip0.h;
  HalfEdge * h1 = ip1.h;

  Edge * e0 = h0->edge();
  Edge * e1 = h1->edge();

  Node * n0 = h0->node();
  Node * n1 = h1->node();

  Point p0 = ip0.point;
  Point p1 = ip1.point;

  uint32_t index = 0;
  Cell * c;
  uint8_t flag = mesh_->find_point(p0, c, index);
  if(flag == 1)/** 在边上 */
  {
  }
  else if(flag == 2)/** 在单元内部 */
  {
  }
  else 
  {
    std::cout << "ERROR in XXX" << std::endl;
  }




  //bool is_in_same_edge; /**< 两个交点是否在同一条边上 */

  //if(ip0.type == 0 && ip1.type == 0) /** 两个点 */
  //{
  //  auto n2n = n0->adj_nodes();
  //  for(auto n_adj : n2n)
  //  {
  //    if(n_adj == n1)
  //    {
  //      is_in_same_edge = true;
  //      break;
  //    }
  //  } 
  //}
  //else if(ip0.type == 0 && ip1.type == 1) /** 一个点一个边 */
  //{
  //  is_in_same_edge = e1->has_node(n0);
  //}
  //else if(ip0.type == 1 && ip1.type == 0) /** 一个边一个点 */
  //{
  //  is_in_same_edge = e0->has_node(n1);
  //}
  //else if(ip0.type == 1 && ip1.type == 1) /** 两个边 */
  //{
  //  is_in_same_edge = (e0 == e1);
  //}

}



}











#endif /* _CUT_MESH_ALGORITHM_ */ 
