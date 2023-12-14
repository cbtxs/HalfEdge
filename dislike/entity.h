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
  Edge(HalfEdge * h = nullptr): start_(h) {}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  HalfEdge * halfedge() { return start_; }

  const HalfEdge * halfedge() const { return start_; }

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
  Cell(HalfEdge * h = nullptr): start_(h) {}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  HalfEdge * halfedge() { return start_; }

  const HalfEdge * halfedge() const { return start_; }

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
