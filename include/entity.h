#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
#include <iostream>
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

  HalfEdge * next_oppo() {return next_->opposite();}

  HalfEdge * oppo_prev() {return oppo_->previous();}

  /** const 数据接口 */
  const HalfEdge * next() const {return next_;}

  const HalfEdge * previous() const {return prev_;}

  const HalfEdge * opposite() const {return oppo_;}

  const HalfEdge * halfedge() const {return this;}

  const Cell * cell() const {return cell_;}

  const Edge * edge() const {return edge_;}

  const Node * node() const {return node_;}

  const uint32_t & index() const { return index_;}

  const HalfEdge * next_oppo() const {return next_->opposite();}

  const HalfEdge * oppo_prev() const {return oppo_->previous();}

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

  Vector tangential() const ;

  Vector normal() const;

  Point barycenter();

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
  static Cell * node2cell[32];
  static Edge * node2edge[32];
  static Node * node2node[32];

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

  uint32_t get_top();

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

Cell * Node::node2cell[32] = {nullptr};
Edge * Node::node2edge[32] = {nullptr};
Node * Node::node2node[32] = {nullptr};

/**
 * @brief Edge 类
 */
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

  Vector tangential() const;

  Vector normal();

  Point barycenter();

  double length(); 

private:
  HalfEdge * start_;
  uint32_t index_;
};

Node * Edge::edge2node[2] = {nullptr};
Cell * Edge::edge2cell[2] = {nullptr};
uint8_t Edge::edge2cellIdx[2] = {255};


/**
 * @brief Cell 类
 */
class Cell
{
public:
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
  uint32_t get_top();

  Cell & operator=(const Cell & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  Point barycenter() const
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

Node * Cell::cell2node[32] = {nullptr};
Edge * Cell::cell2edge[32] = {nullptr};
Cell * Cell::cell2cell[32] = {nullptr};



/** HalfEdge 的一些内联函数 */

/** 判断一个点是不是在半边的左边 */
inline bool HalfEdge::is_on_the_left(const Point & p)
{
  Vector v0 = node_->coordinate()-prev_->node()->coordinate();
  Vector v1 = p-prev_->node()->coordinate();
  return v0.cross(v1)>0;
}

inline Point HalfEdge::barycenter()
{
  return (node()->coordinate() + previous()->node()->coordinate())*0.5;
}

inline Vector HalfEdge::tangential() const 
{
  return node()->coordinate()-previous()->node()->coordinate();
}

inline Vector HalfEdge::normal() const
{
  Vector t = tangential();
  return Vector(-t.y, t.x);
}

inline double HalfEdge::length()
{
  return tangential().length();
}


/** Node 的一些内联函数 */
inline uint32_t Node::get_top()
{
  uint32_t N = 1;
  node2node[0] = start_->node();
  node2edge[0] = start_->edge(); 
  node2cell[0] = start_->cell(); 
  for(HalfEdge * h = start_->next_oppo(); h != start_ && !h->is_boundary(); 
      h = h->next_oppo())
  {
    node2node[N] = h->node();
    node2edge[N] = h->edge(); 
    node2cell[N++] = h->cell(); 
  }
  return N;
}

/** Edge 的一些内联函数 */
inline void Edge::get_top()
{
  edge2node[0] = start_->previous()->node();
  edge2node[1] = start_->node();
  edge2cell[0] = start_->cell(); 
  edge2cell[1] = start_->cell(); 
  edge2cellIdx[0] = start_->distance(start_->cell()->halfedge());
  edge2cellIdx[0] = start_->opposite()->distance(start_->opposite()->cell()->halfedge());
}

inline double Edge::length()
{
  return start_->length();
}

inline Vector Edge::tangential() const
{
  return start_->tangential();
}

inline Vector Edge::normal()
{
  return start_->normal(); 
}

inline Point Edge::barycenter()
{
  return start_->barycenter(); 
}


/** Cell 的一些内联函数 */
inline uint32_t Cell::get_top()
{
  uint32_t N = 0;
  cell2node[N] = start_->node();
  cell2edge[N] = start_->edge(); 
  cell2cell[N++] = start_->cell(); 
  for(HalfEdge * h = start_->next(); h != start_; h = h->next())
  {
    cell2node[N] = h->node();
    cell2edge[N] = h->edge(); 
    cell2cell[N++] = h->cell(); 
  }
  return N;
}

inline double Cell::area()
{
  Point p0 = start_->previous()->node()->coordinate();
  Vector v0 = start_->node()->coordinate()-p0;
  Vector v1 = start_->next()->node()->coordinate()-p0;
  double a = v0.cross(v1);
  for(HalfEdge * h = start_->next()->next(); h != start_->previous(); h = h->next())
  {
    v0 = v1;
    v1 = h->node()->coordinate()-p0;
    a += v0.cross(v1);
  }
  return a/2;
}


}


#endif /* _ENTITY_ */ 
