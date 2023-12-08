#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
#include "halfedge.h"

namespace HEM
{

class Node
{
public:
  Node(HalfEdge * h = nullptr): start_(h) {}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  HalfEdge * halfedge() { return start_; }

  const HalfEdge * halfedge() const { return start_; }

  Node & operator=(const Node & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
};


class Edge
{
public:
  static uint32_t edge2node[2];
  static uint32_t edge2cell[2];
public:
  Edge(HalfEdge * h = nullptr): start_(h) {}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  HalfEdge * halfedge() { return start_; }

  const HalfEdge * halfedge() const { return start_; }

  /** @brief 获取边的邻接关系 */
  void get_top()
  {
    edge2node[0] = start_->previous()->node();
    edge2node[1] = start_->node();
    edge2cell[0] = start_->cell(); 
    edge2cell[1] = start_->cell(); 
  }

  Edge & operator=(const Edge & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
};

class Cell
{
public:
  static uint32_t N;
  static uint32_t cell2node[16];
  static uint32_t cell2edge[16];
  static uint32_t cell2cell[16];

public:
  Cell(HalfEdge * h = nullptr): start_(h) {}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  HalfEdge * halfedge() { return start_; }

  const HalfEdge * halfedge() const { return start_; }

  /** @brief 获取单元的邻接关系 */
  void get_top()
  {
    N = 0;
    for(HalfEdge * h = start_; h != start_->previous(); h = h->next())
    {
      cell2node[N] = h->node();
      cell2edge[N] = h->edge(); 
      cell2cell[N++] = h->cell(); 
    }
  }

  Cell & operator=(const Cell & other)
  {
    if (this != &other) // 避免自我赋值
      start_ = other.start_;
    return *this;
  }

private:
  HalfEdge * start_;
};



}


#endif /* _ENTITY_ */ 
