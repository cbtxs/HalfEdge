#ifndef _HALFEDGE_MESH_BASE_
#define _HALFEDGE_MESH_BASE_

#include "data_set.h"
#include <memory_resource>

namespace HEM
{
template<typename V, typename E, typename H, typename C>
class HalfEdgeMeshBase
{
public:
  HalfEdgeMeshBase(): node(spr), edge(spr), cell(spr), halfedge(spr) {}

  HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh): node(spr), edge(spr), cell(spr), halfedge(spr) 
  {
    int NHE = number_of_halfedges(); 
    halfedge.reserve(NHE);
  }

  template<typename Entity, typename... Args>
  Entity * add_entity(Args &&... args)
  {
    return get_entity<Entity>().add(std::forward<Args>(args)...);
  }

  template<typename Entity>
  data_set<Entity> & get_entity()
  {
    if constexpr (std::is_same_v<Entity, H>)
      return halfedge;
    else if constexpr (std::is_same_v<Entity, V>)
      return node;
    else if constexpr (std::is_same_v<Entity, E>)
      return edge;
    else{
      static_assert(std::is_same_v<Entity, C>);
      return cell;
    }
  }

  int number_of_cells() { return cell.size();}

  int number_of_edges() { return edge.size();}

  int number_of_nodes() { return node.size();}
  
  int number_of_halfedges() { return halfedge.size();}

private:
  std::pmr::synchronized_pool_resource spr;
  data_set<V> node;
  data_set<E> edge;
  data_set<C> cell;
  data_set<H> halfedge;
};

}


#endif /* _HALFEDGE_MESH_BASE_ */ 
