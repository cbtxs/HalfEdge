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
  /** 使用基类的构造函数 */
  using HalfEdgeMeshBase::HalfEdgeMeshBase;

  uint8_t node_to_cell(uint32_t n)
  {
    uint8_t N = 0;
    uint32_t h = Base::get_node()[n].halfedge()->index();
    while(true)
    {
      node2cell[N++] = Base::halfedge_to_cell()[h];
      h = Base::next_halfedge()[h];
      h = Base::oppo_halfedge()[h];
      if(h==Base::get_node()[n].halfedge()->index() || h==Base::oppo_halfedge()[h])
        break;
    }
    return N;
  }

  uint8_t cell_to_node(uint32_t c)
  {
    uint8_t N = 0;
    uint32_t h = Base::get_cell()[c].halfedge()->index();
    while(true)
    {
      cell2node[N++] = Base::halfedge_to_node()[h];
      cell2edge[N++] = Base::halfedge_to_edge()[h];
      h = Base::next_halfedge()[h];
      if(h==Base::get_cell()[c].halfedge()->index())
        break;
    }
    return N;
  }

  void edge_to_node(uint32_t e)
  {
    uint32_t h = Base::get_edge()[e].halfedge()->index();
    edge2node[0] = halfedge_to_node()[prev_halfedge()[h]];
    edge2node[1] = halfedge_to_node()[h];
  }

  void edge_to_cell(uint32_t e)
  {
    uint32_t h = Base::get_edge()[e].halfedge()->index();
    edge2cell[0] = halfedge_to_cell()[prev_halfedge()[h]];
    edge2cell[1] = halfedge_to_cell()[h];
  }

};

}
