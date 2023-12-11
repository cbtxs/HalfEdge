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
  HalfEdge(uint32_t index=0) : 
    next_(nullptr), prev_(nullptr), oppo_(this), cell_(nullptr), 
    edge_(nullptr), node_(nullptr), index_(index)  
  {}

  /** 数据接口 */
  uint32_t & index() { return index_;}

  HalfEdge * next() {return next_;}

  HalfEdge * previous() {return prev_;}

  HalfEdge * opposite() {return oppo_;}

  HalfEdge * halfedge() {return this;}

  Cell * cell() {return cell_;}

  Edge * edge() {return edge_;}

  Node * node() {return node_;}

  /** const 数据接口 */
  const uint32_t & index() const { return index_;}

  const HalfEdge * next() const {return next_;}

  const HalfEdge * previous() const {return prev_;}

  const HalfEdge * opposite() const {return oppo_;}

  const HalfEdge * halfedge() const {return this;}

  const Cell * cell() const {return cell_;}

  const Edge * edge() const {return edge_;}

  const Node * node() const {return node_;}

  template<typename Entity>
  Entity & entity()
  {
    if constexpr (std::is_same_v<Entity, Cell>)
      return cell_;
    else if constexpr (std::is_same_v<Entity, Edge>)
      return edge_;
    else if constexpr (std::is_same_v<Entity, Node>)
      return node_;
  }

  void set_next(HalfEdge * next) {next_ = next;}

  void set_previous(HalfEdge * previous) {prev_ = previous;}

  void set_opposite(HalfEdge * opposite) {oppo_ = opposite;}

  void set_node(Node * node) {node_ = node;}

  void set_edge(Edge * edge) {edge_ = edge;}

  void set_cell(Cell * cell) {cell_ = cell;}

  void set_index(uint32_t index) {index_=index;}

  HalfEdge & operator=(const HalfEdge& other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      next_ = other.next_;
      prev_ = other.prev_;
      oppo_ = other.oppo_;
      index_ = other.index_;
      cell_ = other.cell_;
      edge_ = other.edge_;
      node_ = other.node_;
    }
    return *this;
  }

  bool is_boundary() { return this==oppo_;}

private:
  /** 下一条半边, 上一条半边, 对边 */
  HalfEdge * next_, * prev_, * oppo_;

  /** 半边的存储编号，所属单元的存储编号， 所在边的存储编号，指向顶点的存储编号 */
  Cell * cell_;
  Edge * edge_;
  Node * node_;
  uint32_t index_;
};

}
#endif /* _HALFEDGE_ */ 
