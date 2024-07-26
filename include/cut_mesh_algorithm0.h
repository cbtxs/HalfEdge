#ifndef _CUT_MESH_ALGORITHM_
#define _CUT_MESH_ALGORITHM_

#include <queue>
#include <locale>
#include <unordered_map>
#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <stack>
#include <string>
#include <unordered_set>

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

  template<typename Data>
  using Array = typename Mesh::template Array<Data>;

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

  /**
   * @brief 一个交点
   * @param point 交点的坐标
   * @param h 交点所在的半边
   * @param type 交点的类型
   * 0: 交点在边上
   * 1: 交点在半边的顶点上
   */
  struct Intersection
  {
    Point point;
    HalfEdge * h;
    uint8_t type;
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

  void find_intersect_points(const Interface & iface, std::vector<Point> & intersect_points);

  /**
   * @brief 获取一个点所在的单元，
   * 如果点在边上，返回两个单元，如果点在节点上，返回相邻的单元
   */
  uint8_t get_cell_of_point(Point & point, std::queue<Cell *> & cells)
  {
    Cell * c = nullptr;
    uint32_t index = 0;
    uint8_t flag = mesh_->find_point(point, c, index);
    if(flag == 3)/** 在第 index 个节点上 */
    {
      Node * node = c->adj_node(index);
      cells.resize(32);
      uint32_t N = node->adj_cell(cells.data());
      cells.resize(N);
    }
    else if(flag == 2)/** 在第 index 条边上 */
    {

      Edge * edge = c->adj_edge(index);
      cells.resize(2);
      edge->adj_cell(cells.data());
    }
    else if(flag == 1)/** 在单元内部 */
    {
      cells.push_back(c);
    }
    return flag;
  }

  void cut_mesh_with_segment(const Point & p0, const Point & p1, 
      std::vector<Intersection> & intersections);


private:
  std::shared_ptr<Mesh> mesh_;
};

template<typename Mesh>
void CutMeshAlgorithm<Mesh>::cut_mesh_with_segment(const Point & p0, const Point & p1, 
    std::vector<Intersection> & intersections)
{
  std::unordered_set<Cell *, bool> cell_traversed;
  std::queue<Cell *> cells;

  uint8_t flag0 = get_cell_of_point(p0, cells);
  uint8_t flag1 = get_cell_of_point(p1, cells);

  std::vector<Point *> vertices(0);
  while(!cells.empty())
  {
    vertices.resize(32);

    Cell * c = cells.pop();
    uint8_t N = c->vertices(vertices.data());
    vertices.resize(N);

    std::vector<Point> intersection_points;
    std::vector<uint32_t> intersection_index;
    std::vector<uint8_t> intersection_type;
    intersection_point_of_two_segments(vertices, p0, p1, 
        intersection_points, intersection_index, intersection_type);
    
    N = intersection_points.size();
    for(int i = 0; i < N; i++)
    {
      Intersection intersection;
      intersection.point = intersection_points[i];
      intersection.h = c->halfedge()->next(intersection_index[i]);
      intersection.type = intersection_type[i];
      intersections.push_back(intersection);

      if intersection.type == 0
      {
        Cell * adj_cell = intersection.h->opposite()->cell();
        if(cell_traversed.find(adj_cell) == cell_traversed.end())
          cell_traversed.insert(adj_cell);
      }
      else if intersection.type == 1
      {
        Node * node = intersection.h->vertex();
        uint32_t M = node->adj_cell(cells.data());
        for(int j = 0; j < M; j++)
        {
          Cell * adj_cell = cells[j];
          if(cell_traversed.find(adj_cell) == cell_traversed.end())
            cell_traversed.insert(adj_cell);
        }
      }

    }
    cell_traversed.push_back(c);
  }







}

}
#endif /* _CUT_MESH_ALGORITHM_ */ 
