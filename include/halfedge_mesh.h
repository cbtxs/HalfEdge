#ifndef _HALFEDGE_MESH_BASE_
#define _HALFEDGE_MESH_BASE_

#include <memory>

#include "entity.h"
#include "data_container.h"

namespace HEM
{
class HalfEdgeMeshBase
{
public:
  template<typename T>
  using Array = DataContainer::DataArray<T>;

public:
  HalfEdgeMeshBase(): NN(0), NE(0), NC(0), NBE(0)
  {
    node_ = data_[0].add_data<Node>("node");
    edge_ = data_[1].add_data<Edge>("edge");
    cell_ = data_[2].add_data<Cell>("cell");
    halfedge_ = data_[3].add_data<HalfEdge>("halfedge");
  }

  template<typename Entity, typename... Args>
  Entity * add_entity(Args &&... args)
  {
    return get_entity<Entity>().add(std::forward<Args>(args)...);
  }

  template<typename Entity>
  std::shared_ptr<Array<Entity> > & get_entity()
  {
    if constexpr (std::is_same_v<Entity, HalfEdge>)
      return halfedge_;
    else if constexpr (std::is_same_v<Entity, Node>)
      return node_;
    else if constexpr (std::is_same_v<Entity, Edge>)
      return edge_;
    else
    {
      static_assert(std::is_same_v<Entity, Cell>);
      return cell_;
    }
  }

  uint32_t number_of_cells() { return NC;}

  uint32_t number_of_edges() { return NE;}

  uint32_t number_of_nodes() { return NN;}
  
  uint32_t number_of_boundary_edges() { return NBE;}

  uint32_t number_of_halfedges() { return NE*2-NBE;}

private:
  /** NN : 节点个数; NE : 边的个数; NC: 单元个数; NBE: 边界边个数*/
  uint32_t NN, NE, NC, NBE;

  /** 每个实体的数据 */
  DataContainer data_[4];

  /** 实体列表 */
  std::shared_ptr<Array<Node>> node_;
  std::shared_ptr<Array<Edge>> edge_;
  std::shared_ptr<Array<Cell>> cell_;
  std::shared_ptr<Array<HalfEdge>> halfedge_;

  /** 每个实体的实际编号, 其中半边不需要实际编号 */
  //std::shared_ptr<Array<uint32_t>> indices[3];
};

}


#endif /* _HALFEDGE_MESH_BASE_ */ 
