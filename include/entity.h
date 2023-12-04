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

  const uint32_t & index() const { return index_;}

private:
  uint32_t index_;
};

class Cell
{
public:
  Cell(HalfEdge * h = nullptr): start(h) {}

  void SetHalfEdge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

  const HalfEdge * halfedge() const { return start; }

private:
  HalfEdge * start;
};

using Edge = Cell;

using Node = Cell;

}

#endif /* _ENTITY_ */ 
