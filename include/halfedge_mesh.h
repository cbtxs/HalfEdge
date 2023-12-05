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

  using Entity = std::tuple<Node, Edge, Cell, HalfEdge>;

public:
  HalfEdgeMeshBase(): NN(0), NE(0), NC(0), NHE(0)
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

  template<typename Entity>
  void delete_entity(Entity * e)
  {
    if constexpr (std::is_same_v<Entity, HalfEdge>)
      data_[3].delete_index(e->index());
    else if constexpr (std::is_same_v<Entity, Node>)
      data_[0].delete_index(h2node_[e->halfedge()]);
    else if constexpr (std::is_same_v<Entity, Edge>)
      data_[1].delete_index(h2edge_[e->halfedge()]);
    else
    {
      static_assert(std::is_same_v<Entity, Cell>);
      data_[2].delete_index(h2cell_[e->halfedge()]);
    }
  }

  HalfEdge & add_halfedge()
  {
    uint32_t index = data_[3].add_index();
    HalfEdge & h = halfedge_->get(index);
    h.set_index(index);
    return h;
  }

  template<typename Entity>
  Entity & add_entity(HalfEdge * h)
  {
    if constexpr (std::is_same_v<Entity, Node>)
    {
      uint32_t index = data_[0].add_index();
      Node & node = node_->get(index); 
      node.set_halfedge(h);
      h2node_->get(h->index()) = &node; 
      return node;
    }
    else if constexpr (std::is_same_v<Entity, Edge>)
    {
      uint32_t index = data_[1].add_index();
      Edge & edge = edge_->get(index); 
      edge.set_halfedge(h);
      h2edge_->get(h->index()) = &edge; 
      return edge;
    }
    else
    {
      static_assert(std::is_same_v<Entity, Cell>);
      uint32_t index = data_[0].add_index();
      Cell & cell = cell_->get(index); 
      cell.set_halfedge(h);
      h2cell_->get(h->index()) = &cell; 
      return cell;
    }
  }

  /** halfedge 的一些操作*/
  HalfEdge * next(HalfEdge * h) { return next_->get(h->index()); }

  HalfEdge * prev(HalfEdge * h) { return prev_->get(h->index()); }

  HalfEdge * oppt(HalfEdge * h) { return oppt_->get(h->index()); }

  Node * node(HalfEdge * h) { return h2node_->get(h->index()); }  

  Edge * edge(HalfEdge * h) { return h2edge_->get(h->index()); }  

  Cell * cell(HalfEdge * h) { return h2cell_->get(h->index()); }  

  uint32_t number_of_cells() { return NC;}

  uint32_t number_of_edges() { return NE;}

  uint32_t number_of_nodes() { return NN;}
  
  uint32_t number_of_boundary_edges() { return NE*2-NHE;}

  uint32_t number_of_halfedges() { return NHE;}

private:
  /** NN : 节点个数; NE : 边的个数; NC: 单元个数; NBE: 边界边个数*/
  uint32_t NN, NE, NC, NHE;

  /** 每个实体的数据, node, edge, cell, halfedge */
  DataContainer data_[4];

  /** 实体列表 */
  std::shared_ptr<Array<Node>> node_;
  std::shared_ptr<Array<Edge>> edge_;
  std::shared_ptr<Array<Cell>> cell_;
  std::shared_ptr<Array<HalfEdge>> halfedge_;

  /** 半边对应的节点, 边，单元，下一条半边，上一条半边，对边*/
  std::shared_ptr<Array<Node*> >   h2node_;
  std::shared_ptr<Array<Edge*> >   h2edge_;
  std::shared_ptr<Array<Cell*> >   h2cell_;
  std::shared_ptr<Array<HalfEdge*> > next_;
  std::shared_ptr<Array<HalfEdge*> > prev_;
  std::shared_ptr<Array<HalfEdge*> > oppt_;

  /** 每个实体的实际编号, 其中半边不需要实际编号 */
  //std::shared_ptr<Array<uint32_t>> indices[3];
};

}


#endif /* _HALFEDGE_MESH_BASE_ */ 
