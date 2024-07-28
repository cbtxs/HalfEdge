#ifndef INTERFACE_H
#define INTERFACE_H

#include <vector>
#include <memory>
#include <numeric>

/**
 * @brief 表示一个二维界面, 这个界面由一些点组成，这些点可能在边上，
 * 也可能在节点上， 经过判断，会对点进行很小的移动，使得点在边上或者节点上。
 */
template<typename Mesh>
class InterfaceCut
{
public: /** 类型定义 */
  using Cell = typename Mesh::Cell;
  using Edge = typename Mesh::Edge;
  using Node = typename Mesh::Node;
  using HalfEdge = typename Mesh::HalfEdge;

  using Point = typename Mesh::Point;

  /**
   * @brief 表示一个界面点
   * @param point 点的坐标
   * @param is_fixed_point 是否是固定点
   * @param cells 点所在的单元
   * @param edge 点所在的半边
   * @param type 0: 半边的节点;
   *             1: 半边内部; 
   *             2: 单元内的点;
   */
  struct InterfacePoint
  {
    Point point;
    bool is_fixed_point;
    std::vector<Cell * > cells;
    HalfEdge * h;             
    uint8_t type;             
  };

private: /** 属性 */
  std::vector<InterfacePoint> points_;
  std::vector<uint32_t> segments_;
  std::shared_ptr<Mesh> mesh_;
  bool is_loop_;

public: /** 方法 */

  /**
   * @brief 构造函数
   */
  InterfaceCut(std::vector<Point> & points, 
               std::vector<bool>  & is_fixed_points, 
               std::shared_ptr<Mesh> mesh,
               bool _is_loop=true): mesh_(mesh), is_loop_(_is_loop)
  {
    uint32_t N = points.size();
    for(uint32_t i = 0; i < N; i++)
    {
      InterfacePoint ip;
      ip.is_fixed_point = is_fixed_points[i];
      point_to_interface_point(points[i], ip);
      points_.push_back(ip);
    }
    segments_.resize(N);
    std::iota(segments_.begin(), segments_.end(), 0);
    if (is_loop_)
      segments_.push_back(0);
  }

  std::vector<InterfacePoint> & points()
  {
    return points_;
  }

  const std::vector<InterfacePoint> & points() const
  {
    return points_;
  }

  std::vector<uint32_t> & segments()
  {
    return segments_;
  }

  const std::vector<uint32_t> & segments() const
  {
    return segments_;
  }

  /**
   * @brief 获取一个点所在的单元，
   * 如果点在边上，返回两个单元，如果点在节点上，返回相邻的单元, 
   * 同时更新点的坐标。
   */
  void point_to_interface_point(Point & point, InterfacePoint & inerface_point);

  /**
   * @brief 判断界面是否是闭的
   */
  bool is_loop_interface() const
  {
    return is_loop_; 
  }

};

/**
 * @brief 获取一个点所在的单元，
 * 如果点在边上，返回两个单元，如果点在节点上，返回相邻的单元, 
 * 同时更新点的坐标。
 */
template<typename Mesh>
void InterfaceCut<Mesh>::point_to_interface_point(Point & point, InterfacePoint & ip)
{
  auto & cells = ip.cells;
  auto & geometry_utils = mesh_->geometry_utils();
  Cell * c = nullptr;
  uint32_t index = 0;
  uint8_t flag = mesh_->find_point(point, c, index);
  if(flag == 0)/** 在第 index 个节点上 */
  {
    Node * node = c->adj_node(index);
    HalfEdge * h = c->halfedge()->previous()->next(index);

    cells.resize(32);
    uint32_t N = node->adj_cell(cells.data());
    cells.resize(N);

    point   = node->coordinate(); /**< 更新点的坐标 */
    ip.h = h; 
  }
  else if(flag == 1)/** 在第 index 条边上 */
  {
    Point * p[2];
    Edge  * edge = c->adj_edge(index);

    cells.resize(2);
    edge->adj_cell(cells.data());
    edge->vertices(p);
    geometry_utils.project_point_to_line(*(p[0]), *(p[2]), point); /**< 更新点的坐标 */

    ip.h = edge->halfedge(); 
  }
  else if(flag == 2)/** 在单元内部 */
  {
    cells.push_back(c);
  }
  ip.point = point;
  ip.type  = flag;
}

#endif // INTERFACE_H
