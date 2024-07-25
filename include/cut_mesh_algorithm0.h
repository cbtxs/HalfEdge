#ifndef _CUT_MESH_ALGORITHM_
#define _CUT_MESH_ALGORITHM_

#include <utility>
#include <vector>
#include <algorithm>
#include <memory>
#include <iostream>
#include <stack>
#include <string>

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
   * @brief 获取一个点所在的单元，如果点在边上，返回两个单元，如果点在节点上，返回相邻的单元
   */
  void get_cell_of_point(const Point & point, std::vector<Cell *> & cells)
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
  }


private:
  std::shared_ptr<Mesh> mesh_;
};


}

#endif /* _CUT_MESH_ALGORITHM_ */ 
