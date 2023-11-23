#ifndef _HALF_EDGE_
#define _HALF_EDGE_

namespace HEM {

class Vertex;
class Edge;
class Cell;

class HalfEdge
{
public:
  HalfEdge() {}

  HalfEdge(Vertex * v, Edge * e, Cell * c, HalfEdge * n, HalfEdge * p, HalfEdge * o):
    ver(v), edg(e), cel(c), nex(n), pre(p), opp(o) {}

  void Init(Vertex * v, Edge * e, Cell * c, HalfEdge * n, HalfEdge * p, HalfEdge * o)
  {
    ver = v; edg = e; cel = c; nex = n; pre = p; opp = o;
  }

  /** 半边的基本数据 */
  Vertex * vertex() {return ver;} 

  Edge * edge() { return edg; }

  Cell * cell() { return cel; }

  HalfEdge * next() { return nex; }

  HalfEdge * previous() { return pre; }

  HalfEdge * opposite() { return opp; }

  /** const 版 */
  const Vertex * vertex() const {return ver;} 

  const Edge * edge() const { return edg; }

  const Cell * cell() const { return cel; }

  const HalfEdge * next() const { return nex; }

  const HalfEdge * previous() const { return pre; }

  const HalfEdge * opposite() const { return opp; }

  /** 设置数据 */
  void setVertex(Vertex * v) { ver = v; }

  void setEdge(Edge * e) { edg = e; }

  void setCell(Cell * c) { cel = c; }

  void setNext(HalfEdge * h) { nex = h; }

  void setPrevious(HalfEdge * h) { pre = h; }

  void setOpposite(HalfEdge * h) { opp = h; }

  /** 基本函数 */
  bool is_boundary() { return !cell(); }

private:
  Vertex   * ver; 
  Edge     * edg;
  Cell     * cel;
  HalfEdge * nex;
  HalfEdge * pre;
  HalfEdge * opp;
};

}
#endif /* _HALF_EDGE_ */ 
