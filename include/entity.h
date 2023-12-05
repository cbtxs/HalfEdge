#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
namespace HEM
{

class HalfEdge
{
public:
  HalfEdge(uint32_t index): index_(index) {}

  uint32_t & index() { return index_;}

  void set_index(uint32_t index) {index_=index;}

  const uint32_t & index() const { return index_;}

private:
  uint32_t index_;
};

class Cell
{
public:
  Cell(HalfEdge * h = nullptr): start(h) {}

  void set_halfedge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

  const HalfEdge * halfedge() const { return start; }

private:
  HalfEdge * start;
};

class Edge
{
public:
  Edge(HalfEdge * h = nullptr): start(h) {}

  void set_halfedge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

  const HalfEdge * halfedge() const { return start; }

private:
  HalfEdge * start;
};

class Node
{
public:
  Node(HalfEdge * h = nullptr): start(h) {}

  void set_halfedge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

  const HalfEdge * halfedge() const { return start; }

private:
  HalfEdge * start;
};


}


#endif /* _ENTITY_ */ 
