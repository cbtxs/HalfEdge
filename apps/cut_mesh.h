
#include "entity.h"
#include "halfedge_mesh.h"

namespace HEM
{
using HalfEdgeMesh = HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>;

void cut_mesh(HalfEdgeMesh * mesh, 
              double * point, 
              std::vector<bool> & is_fixed_point, 
              std::vector<uint32_t> & interface)
{

  
}

}
