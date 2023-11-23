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
  HalfEdgeMeshBase(): hnode(spr), hedge(spr), hcell(spr), halfedge(spr) {}

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
      return hnode;
    else if constexpr (std::is_same_v<Entity, E>)
      return hedge;
    else{
      static_assert(std::is_same_v<Entity, C>);
      return hcell;
    }
  }

private:
  std::pmr::synchronized_pool_resource spr;
  data_set<V> hnode;
  data_set<E> hedge;
  data_set<C> hcell;
  data_set<H> halfedge;
};

}


#endif /* _HALFEDGE_MESH_BASE_ */ 
