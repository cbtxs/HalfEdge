#ifndef _CELL_
#define _CELL_

namespace HEM 
{
class HalfEdge;

class Cell
{
public:
  Cell() {}

  Cell(HalfEdge * h): start(h) {}

  void SetHalfEdge(HalfEdge * h) { start = h;}

  HalfEdge * halfedge() { return start; }

private:
  HalfEdge * start;

};

}

#endif /* _CELL_ */ 
