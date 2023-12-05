#ifndef _HALFEDGE_
#define _HALFEDGE_

#include <stdint.h>
#include <type_traits>

namespace HEM
{

class Cell;
class Edge;
class Node;

class HalfEdge
{
public:
  HalfEdge(uint32_t index): index_(index) {}

  uint32_t & index() { return index_;}

  void set_index(uint32_t index) {index_=index;}

  const uint32_t & index() const { return index_;}

  HalfEdge * next() {return next_;}

  HalfEdge * previous() {return prev_;}

  HalfEdge * opposite() {return oppo_;}

  HalfEdge * halfedge() {return this;}

  uint32_t & cell() {return cell_;}

  uint32_t & edge() {return edge_;}

  uint32_t & node() {return node_;}

  template<typename Entity>
  uint32_t & entity()
  {
      if constexpr (std::is_same_v<Entity, HalfEdge>)
          return index_;
      else if constexpr (std::is_same_v<Entity, Cell>)
          return cell_;
      else if constexpr (std::is_same_v<Entity, Edge>)
          return edge_;
      else if constexpr (std::is_same_v<Entity, Node>)
          return node_;
  }

  HalfEdge & operator=(const HalfEdge& other)
  {
    if (this != &other) // 避免自我赋值
      index_ = other.index_;
    return *this;
  }

private:
  /** 下一条半边, 上一条半边, 对边 */
  HalfEdge * next_, * prev_, * oppo_;

  /** 半边的存储编号，所属单元的存储编号， 所在边的存储编号，指向顶点的存储编号 */
  uint32_t index_, cell_, edge_, node_;
};

}
#endif /* _HALFEDGE_ */ 
