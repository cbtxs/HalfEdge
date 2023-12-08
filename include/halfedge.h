#ifndef _HALFEDGE_
#define _HALFEDGE_

#include <stdint.h>
#include <type_traits>

namespace HEM
{

class Cell;
class Edge;
class Node;

class HalfEdge
{
public:
  HalfEdge(uint32_t index=0) : index_(index) {} 

  uint32_t & index() { return index_;}

  const uint32_t & index() const { return index_;}

  HalfEdge & operator=(const HalfEdge& other)
  {
    if (this != &other) // 避免自我赋值
      index_ = other.index_;
    return *this;
  }

private:

  /** 半边的存储编号，所属单元的存储编号， 所在边的存储编号，指向顶点的存储编号 */
  uint32_t index_;
};

}
#endif /* _HALFEDGE_ */ 
