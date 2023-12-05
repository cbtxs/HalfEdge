#ifndef _HALFEDGE_MESH_BASE_
#define _HALFEDGE_MESH_BASE_

#include <memory>
#include <tuple>

#include "entity.h"
#include "data_container.h"

namespace HEM
{
class HalfEdgeMeshBase
{
public:
  template<typename T>
  using Array = DataContainer::DataArray<T>;

  //using Entity = std::tuple<Node, Edge, Cell, HalfEdge>;
public:
  HalfEdgeMeshBase(): NN(0), NE(0), NC(0), NHE(0) 
  {
    node_ = node_data_.add_data<Node>("node");
    edge_ = edge_data_.add_data<Edge>("edge");
    cell_ = cell_data_.add_data<Cell>("cell");
    halfedge_ = halfedge_data_.add_data<HalfEdge>("halfedge");
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

  template<typename Entity>
  void delete_entity(Entity * e)
  {
    if constexpr (std::is_same_v<Entity, HalfEdge>)
      halfedge_data_.delete_entity(e);
    else if constexpr (std::is_same_v<Entity, Node>)
      node_data_.delete_entity(e);
    else if constexpr (std::is_same_v<Entity, Edge>)
      edge_data_.delete_entity(e);
    else
    {
      static_assert(std::is_same_v<Entity, Cell>);
      cell_data_.delete_entity(e);
    }
  }

  template<typename Entity>
  Entity & add_entity()
  {
    if constexpr (std::is_same_v<Entity, HalfEdge>)
      return halfedge_data_.add_entity();
    else if constexpr (std::is_same_v<Entity, Node>)
      return node_data_.add_entity();
    else if constexpr (std::is_same_v<Entity, Edge>)
      return edge_data_.add_entity();
    else
    {
      static_assert(std::is_same_v<Entity, Cell>);
      return cell_data_.add_entity();
    }
  }

  uint32_t number_of_cells() { return NC;}

  uint32_t number_of_edges() { return NE;}

  uint32_t number_of_nodes() { return NN;}
  
  uint32_t number_of_boundary_edges() { return NE*2-NHE;}

  uint32_t number_of_halfedges() { return NHE;}

private:
  /** NN : 节点个数; NE : 边的个数; NC: 单元个数; NBE: 边界边个数*/
  uint32_t NN, NE, NC, NHE;

  /** 每个实体的数据, node, edge, cell, halfedge */
  EntityDataContainer<Node> node_data_;
  EntityDataContainer<Edge> edge_data_;
  EntityDataContainer<Cell> cell_data_;
  EntityDataContainer<HalfEdge> halfedge_data_;

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
