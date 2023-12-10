#include "halfedge_mesh.h"

namespace HEM
{

class PolygonMesh : public HalfEdgeMeshBase
{
public:
  using Self = PolygonMesh;
  using Base = HalfEdgeMeshBase;

  static uint32_t cell2edge[32];
  static uint32_t cell2node[32];
  static uint32_t edge2node[2];
  static uint32_t edge2cell[4];
  static uint32_t node2cell[32];

public:
  PolygonMesh(): HalfEdgeMeshBase() {} 

  /** 复制构造函数 */
  HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh);

  /** 以单元为中心的网格 */
  HalfEdgeMeshBase(double * node, uint32_t * cell, uint32_t NN, uint32_t NC, uint32_t NV);

  uint8_t node_to_cell(uint32_t n)
  {
    uint32_t N = 0;
    while true  
    return N;
  }

};

}
