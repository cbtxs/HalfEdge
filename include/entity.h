#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
#include <iostream>
#include "geometry.h"

namespace HEM
{

template<typename HalfEdgeMesh_Traits>
class THalfEdge
{
public:
  THalfEdge(uint32_t index=0) : 
    next_(nullptr), prev_(nullptr), oppo_(this), cell_(nullptr), 
    edge_(nullptr), node_(nullptr), index_(index)  
  {}

  /** 数据接口 */
  THalfEdge * next() {return next_;}

  THalfEdge * previous() {return prev_;}

  THalfEdge * opposite() {return oppo_;}

  THalfEdge * halfedge() {return this;}

  Cell * cell() {return cell_;}

  Edge * edge() {return edge_;}

  Node * node() {return node_;}

  uint32_t & index() { return index_;}

  THalfEdge * next_oppo() {return next_->opposite();}

  THalfEdge * oppo_prev() {return oppo_->previous();}

  /** const 数据接口 */
  const THalfEdge * next() const {return next_;}

  const THalfEdge * previous() const {return prev_;}

  const THalfEdge * opposite() const {return oppo_;}

  const THalfEdge * halfedge() const {return this;}

  const Cell * cell() const {return cell_;}

  const Edge * edge() const {return edge_;}

  const Node * node() const {return node_;}

  const uint32_t & index() const { return index_;}

  const THalfEdge * next_oppo() const {return next_->opposite();}

  const THalfEdge * oppo_prev() const {return oppo_->previous();}

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
  void set_next(THalfEdge * next) {next_ = next;}

  void set_previous(THalfEdge * prev) {prev_ = prev;}

  void set_opposite(THalfEdge * oppo) {oppo_ = oppo;}

  void set_node(Node * node) {node_ = node;}

  void set_edge(Edge * edge) {edge_ = edge;}

  void set_cell(Cell * cell) {cell_ = cell;}

  void set_index(uint32_t index) {index_=index;}

  void reset(THalfEdge * next, THalfEdge * prev, THalfEdge * oppo, 
      Cell * cell, Edge * edge, Node * node, uint32_t index)
  {
    next_ = next; prev_ = prev; oppo_ = oppo;
    cell_ = cell; edge_ = edge; node_ = node;
    index_ = index;
  }

  THalfEdge & operator=(const THalfEdge& other)
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
  uint8_t distance(THalfEdge * h)
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
  THalfEdge * next_, * prev_, * oppo_;

  Cell * cell_; /**< 所属单元*/ 
  Edge * edge_; /**< 所在边*/
  Node * node_; /**< 指向顶点 */
  uint32_t index_; /**< 半边的存储编号 */ 
};

class Node
{
public:
  static Cell * node2cell[128];
  static Edge * node2edge[128];
  static Node * node2node[128];

public:
  Node(uint32_t index = -1): start_(nullptr), index_(index), coordinate_(0.0, 0.0) {}

  Node(Point & coordinate, uint32_t index, THalfEdge * h = nullptr): 
    start_(h), index_(index), coordinate_(coordinate)
  {}

  /** 获取数据 */
  THalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  Point & coordinate() { return coordinate_;}

  const THalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  const Point & coordinate() const { return coordinate_;}

  /** 设置数据 */
  void set_coordinate(const Point & p) { coordinate_ = p;}

  void set_halfedge(THalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(const Point & p, uint32_t index, THalfEdge * h) 
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
  THalfEdge * start_;
  uint32_t index_;
  Point coordinate_;
};

Cell * Node::node2cell[128] = {nullptr};
Edge * Node::node2edge[128] = {nullptr};
Node * Node::node2node[128] = {nullptr};

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
  Edge(uint32_t index=-1, THalfEdge * h = nullptr): start_(h), index_(index) {}

  THalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const THalfEdge * halfedge() const { return start_; }

  void set_halfedge(THalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, THalfEdge * h) { index_ = index; start_ = h; }

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
  THalfEdge * start_;
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
  Cell(uint32_t index = -1, THalfEdge * h = nullptr): start_(h), index_(index) {}

  THalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const THalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  void set_halfedge(THalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, THalfEdge * h) { index_ = index; start_ = h; }

  /** 获取单元的邻接关系 */
  uint32_t get_top();

  /** 获取单元的邻接关系 */
  uint32_t cell_to_edge(Edge ** cell2edge);

  /** 获取单元的邻接关系 */
  uint32_t cell_to_node(Node ** cell2node);

  Cell & operator=(const Cell & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  /** 返回一个内点 */
  Point inner_point() const;

  Point barycenter() const
  {
    uint8_t n = 1;
    Point p = start_->node()->coordinate();
    for(THalfEdge * h = start_->next(); h != start_; h = h->next(), n++)
      p += h->node()->coordinate(); 
    return p/n; 
  }

  double area();

private:
  THalfEdge * start_;
  uint32_t index_;
};

Node * Cell::cell2node[32] = {nullptr};
Edge * Cell::cell2edge[32] = {nullptr};
Cell * Cell::cell2cell[32] = {nullptr};

/** THalfEdge 的一些内联函数 */

/** 判断一个点是不是在半边的左边 */
inline bool THalfEdge::is_on_the_left(const Point & p)
{
  Vector v0 = node_->coordinate()-prev_->node()->coordinate();
  Vector v1 = p-prev_->node()->coordinate();
  return v0.cross(v1)>0;
}

inline Point THalfEdge::barycenter()
{
  return (node()->coordinate() + previous()->node()->coordinate())*0.5;
}

inline Vector THalfEdge::tangential() const 
{
  return node()->coordinate()-previous()->node()->coordinate();
}

inline Vector THalfEdge::normal() const
{
  Vector t = tangential();
  return Vector(-t.y, t.x);
}

inline double THalfEdge::length()
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
  for(THalfEdge * h = start_->next_oppo(); h != start_ && !h->is_boundary(); 
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
  for(THalfEdge * h = start_->next(); h != start_; h = h->next())
  {
    cell2node[N] = h->node();
    cell2edge[N] = h->edge(); 
    cell2cell[N++] = h->opposite()->cell(); 
  }
  return N;
}

/** Cell 的一些内联函数 */
inline uint32_t Cell::cell_to_edge(Edge ** c2e)
{
  uint32_t N = 0;
  c2e[N++] = start_->edge(); 
  for(THalfEdge * h = start_->next(); h != start_; h = h->next())
    c2e[N] = h->edge(); 
  return N;
}

/** Cell 的一些内联函数 */
inline uint32_t Cell::cell_to_node(Node ** c2n)
{
  uint32_t N = 0;
  c2n[N++] = start_->node(); 
  for(THalfEdge * h = start_->next(); h != start_; h = h->next())
    c2n[N] = h->node(); 
  return N;
}

inline double Cell::area()
{
  Point p0 = start_->previous()->node()->coordinate();
  Vector v0 = start_->node()->coordinate()-p0;
  Vector v1 = start_->next()->node()->coordinate()-p0;
  double a = v0.cross(v1);
  for(THalfEdge * h = start_->next()->next(); h != start_->previous(); h = h->next())
  {
    v0 = v1;
    v1 = h->node()->coordinate()-p0;
    a += v0.cross(v1);
  }
  return a/2;
}

/** 返回一个内点 */
inline Point Cell::inner_point() const
{
  Point p(0, 0);
  for(THalfEdge* h = start_->next(); h != start_; h = h->next())
  {
    if(h->tangential().cross(h->next()->tangential())>0)
    {
      p  = h->previous()->node()->coordinate();
      p += h->node()->coordinate();
      p += h->next()->node()->coordinate();
      return p/3.0;
    }
  }
  return p;
}


}


#endif /* _ENTITY_ */ 
