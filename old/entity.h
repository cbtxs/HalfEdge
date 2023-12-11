#ifndef _ENTITY_
#define _ENTITY_

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

class Node
{
public:
  Node(uint32_t index, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  Node & operator=(const Node & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
  uint32_t index_;
};

class Edge
{
public:
  static uint32_t edge2node[2];
  static uint32_t edge2cell[4];
public:
  Edge(uint32_t index, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  /** @brief 获取边的邻接关系 */
  void get_top();

  Edge & operator=(const Edge & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
  uint32_t index_;
};

class Cell
{
public:
  static uint32_t N;
  static uint32_t cell2node[16];
  static uint32_t cell2edge[16];
  static uint32_t cell2cell[16];

public:
  Cell(uint32_t index, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  uint8_t halfedge_to_cell_location_number(HalfEdge * h);

  /** @brief 获取单元的邻接关系 */
  void get_top();

  Cell & operator=(const Cell & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
  uint32_t index_;
};

inline uint8_t Cell::halfedge_to_cell_location_number(HalfEdge * h)
{
  uint8_t n = 0;
  HalfEdge * s = h->cell()->halfedge();
  while(s != start_)
  {
    s = s->next();
    n++;
  }
  return n;
}

inline void Edge::get_top()
{
  edge2node[0] = start_->previous()->node()->index();
  edge2node[1] = start_->node()->index();
  edge2cell[0] = start_->cell()->index(); 
  edge2cell[1] = start_->cell()->index(); 
  edge2cell[2] = start_->halfedge_to_cell_location_number();
  edge2cell[3] = start_->opposite()->halfedge_to_cell_location_number();
}

inline void Cell::get_top()
{
  N = 0;
  for(HalfEdge * h = start_; h != start_->previous(); h = h->next())
  {
    cell2node[N] = h->node()->index();
    cell2edge[N] = h->edge()->index(); 
    cell2cell[N++] = h->cell()->index(); 
  }
}


}


#endif /* _ENTITY_ */ 
