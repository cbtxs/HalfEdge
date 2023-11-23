#ifndef _EDGE_
#define _EDGE_

#include "Geometry.h"
#include "HalfEdge.h"

namespace HEM 
{

class Edge
{
public:
  Edge(HalfEdge * h = nullptr): start(h) {}

  void SetHalfEdge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

  Vector tangential() {return *(start->vertex())-*(start->opposite()->vertex());}

  Vector normal() { return tangential().rotcw();}


  bool is_boundary() noexcept
  {
    return start->is_boundary() || start->opposite()->is_boundary();
  }

private:
  HalfEdge * start;

};

}

#endif /* _EDGE_ */ 
