#ifndef _VERTEX_
#define _VERTEX_

#include "Geometry.h"

namespace HEM
{
class HalfEdge;

class Vertex : public Vector
{
public:
  Vertex (Point & p, HalfEdge * h = nullptr): Vector(p), start(h) {}

  void SetHalfEdge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

private:
  HalfEdge * start;
};

}

#endif /* _VERTEX_ */ 
