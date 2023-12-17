#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
#include <type_traits>
#include "geometry.h"

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
  HalfEdge * next() {return next_;}

  HalfEdge * previous() {return prev_;}

  HalfEdge * opposite() {return oppo_;}

  HalfEdge * halfedge() {return this;}

  Cell * cell() {return cell_;}

  Edge * edge() {return edge_;}

  Node * node() {return node_;}

  uint32_t & index() { return index_;}

  /** const 数据接口 */
  const HalfEdge * next() const {return next_;}

  const HalfEdge * previous() const {return prev_;}

  const HalfEdge * opposite() const {return oppo_;}

  const HalfEdge * halfedge() const {return this;}

  const Cell * cell() const {return cell_;}

  const Edge * edge() const {return edge_;}

  const Node * node() const {return node_;}

  const uint32_t & index() const { return index_;}

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

  /** 设置数据 */
  void set_next(HalfEdge * next) {next_ = next;}

  void set_previous(HalfEdge * prev) {prev_ = prev;}

  void set_opposite(HalfEdge * oppo) {oppo_ = oppo;}

  void set_node(Node * node) {node_ = node;}

  void set_edge(Edge * edge) {edge_ = edge;}

  void set_cell(Cell * cell) {cell_ = cell;}

  void set_index(uint32_t index) {index_=index;}

  void reset(HalfEdge * next, HalfEdge * prev, HalfEdge * oppo, 
      Cell * cell, Edge * edge, Node * node, uint32_t index)
  {
    next_ = next; prev_ = prev; oppo_ = oppo;
    cell_ = cell; edge_ = edge; node_ = node;
    index_ = index;
  }

  HalfEdge & operator=(const HalfEdge& other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      next_ = other.next_; prev_ = other.prev_; oppo_ = other.oppo_;
      cell_ = other.cell_; edge_ = other.edge_; node_ = other.node_;
      index_ = other.index_;
    }
    return *this;
  }

  bool is_boundary() { return this==oppo_;}

  /** h 到 self 需要 next 的次数 */
  uint8_t distance(HalfEdge * h)
  {
    //assert(h->cell() == cell_)
    uint8_t n = 0;
    while(h != this) { h = h->next(); n++; }
    return n;
  }

  bool is_on_the_left(const Point & p);

  double length(); 

  Vector tangential();

  Vector normal();

  Point barycentary();

private:
  /** 下一条半边, 上一条半边, 对边 */
  HalfEdge * next_, * prev_, * oppo_;

  Cell * cell_; /**< 所属单元*/ 
  Edge * edge_; /**< 所在边*/
  Node * node_; /**< 指向顶点 */
  uint32_t index_; /**< 半边的存储编号 */ 
};

class Node
{
public:
  Node(uint32_t index = -1): start_(nullptr), index_(index), coordinate_(0.0, 0.0) {}

  Node(Point & coordinate, uint32_t index, HalfEdge * h = nullptr): 
    start_(h), index_(index), coordinate_(coordinate)
  {}

  /** 获取数据 */
  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  Point & coordinate() { return coordinate_;}

  const HalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  const Point & coordinate() const { return coordinate_;}


  /** 设置数据 */
  void set_coordinate(const Point & p) { coordinate_ = p;}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(const Point & p, uint32_t index, HalfEdge * h) 
  { 
    index_ = index; 
    start_ = h; 
    coordinate_ = p;
  }

  Node & operator=(const Node & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
      coordinate_ = other.coordinate_;
    }
    return *this;
  }

private:
  HalfEdge * start_;
  uint32_t index_;
  Point coordinate_;
};

class Edge
{
public:
  static Node * edge2node[2];
  static Cell * edge2cell[2];
  static uint8_t edge2cellIdx[2];
public:
  Edge(uint32_t index=-1, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, HalfEdge * h) { index_ = index; start_ = h; }

  /** 获取边的邻接关系 */
  void get_top();

  Edge & operator=(const Edge & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  Vector tangential();

  Vector normal();

  Point barycentary();

  double length(); 

private:
  HalfEdge * start_;
  uint32_t index_;
};

class Cell
{
public:
  static uint8_t N;
  static Node * cell2node[32];
  static Edge * cell2edge[32];
  static Cell * cell2cell[32];

public:
  Cell(uint32_t index = -1, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, HalfEdge * h) { index_ = index; start_ = h; }

  /** 获取单元的邻接关系 */
  void get_top();

  Cell & operator=(const Cell & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  Point barycentary() const
  {
    uint8_t n = 1;
    Point p = start_->node()->coordinate();
    for(HalfEdge * h = start_->next(); h != start_; h = h->next(), n++)
      p += h->node()->coordinate(); 
    return p/n; 
  }

  double area();

private:
  HalfEdge * start_;
  uint32_t index_;
};

Node * Edge::edge2node[2] = {nullptr};
Cell * Edge::edge2cell[2] = {nullptr};
uint8_t Edge::edge2cellIdx[2] = {255};

uint8_t Cell::N = 0;
Node * Cell::cell2node[32] = {nullptr};
Edge * Cell::cell2edge[32] = {nullptr};
Cell * Cell::cell2cell[32] = {nullptr};

inline bool HalfEdge::is_on_the_left(const Point & p)
{
  Vector v0 = node_->coordinate()-prev_->node()->coordinate();
  Vector v1 = p-prev_->node()->coordinate();
  return v0.cross(v1)>0;
}

inline void Edge::get_top()
{
  edge2node[0] = start_->previous()->node();
  edge2node[1] = start_->node();
  edge2cell[0] = start_->cell(); 
  edge2cell[1] = start_->cell(); 
  edge2cellIdx[0] = start_->distance(start_->cell()->halfedge());
  edge2cellIdx[0] = start_->opposite()->distance(start_->opposite()->cell()->halfedge());
}

inline void Cell::get_top()
{
  N = 0;
  cell2node[N] = start_->node();
  cell2edge[N] = start_->edge(); 
  cell2cell[N++] = start_->cell(); 
  for(HalfEdge * h = start_->next(); h != start_; h = h->next())
  {
    cell2node[N] = h->node();
    cell2edge[N] = h->edge(); 
    cell2cell[N++] = h->cell(); 
  }
}

inline Point HalfEdge::barycentary()
{
  return (node()->coordinate() + previous()->node()->coordinate())*0.5;
}

inline Vector HalfEdge::tangential() 
{
  return node()->coordinate()-previous()->node()->coordinate();
}

inline Vector HalfEdge::normal()
{
  Vector t = tangential();
  return Vector(-t.y, t.x);
}

inline double HalfEdge::length()
{
  return tangential().length();
}

inline double Edge::length()
{
  return start_->length();
}

inline Vector Edge::tangential()
{
  return start_->tangential();
}

inline Vector Edge::normal()
{
  return start_->normal(); 
}

inline Point Edge::barycentary()
{
  return start_->barycentary(); 
}

inline double Cell::area()
{
  Vector v0 = start_->previous()->tangential();
  Vector v1 = start_->tangential();
  double a = v1.cross(v0);
  for(HalfEdge * h = start_; h != start_->previous(); h = h->next())
  {
    v0 = v1;
    v1 = h->next()->tangential();
    a += v0.cross(v1);
  }
  return a/2;
}


}


#endif /* _ENTITY_ */ 
