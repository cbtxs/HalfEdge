#ifndef _VERTEX_
#define _VERTEX_

namespace HEM
{
class HalfEdge;

class Node
{
public:
  Node(HalfEdge * h = nullptr): start(h) {}

  void SetHalfEdge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

private:
  HalfEdge * start;
};

}

#endif /* _VERTEX_ */ 
