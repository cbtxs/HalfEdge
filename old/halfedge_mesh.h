#ifndef _HALFEDGE_MESH_BASE_
#define _HALFEDGE_MESH_BASE_

#include <memory>
#include <tuple>

#include "entity.h"
#include "data_container.h"
#include "tuple_type_index.h"
#include "geometry.h"

namespace HEM
{
class HalfEdgeMeshBase
{
public:
  template<typename T>
  using Array = DataContainer::DataArray<T>;
  using ENTITY = std::tuple<Node, Edge, Cell, HalfEdge>;
  using Container = std::tuple<std::shared_ptr<EntityDataContainer<Node>> & , 
                               std::shared_ptr<EntityDataContainer<Edge>> & ,
                               std::shared_ptr<EntityDataContainer<Cell>> & ,
                               std::shared_ptr<EntityDataContainer<HalfEdge>>&  >;
public:
  HalfEdgeMeshBase(): entity_data_(node_data_, edge_data_, cell_data_, halfedge_data_) 
  {
    node_data_ = std::make_shared<EntityDataContainer<Node> >();
    edge_data_ = std::make_shared<EntityDataContainer<Edge> >();
    cell_data_ = std::make_shared<EntityDataContainer<Cell> >();
    halfedge_data_ = std::make_shared<EntityDataContainer<HalfEdge> >();
    node_data_->add_data<Point>("corrdinate");
  }

  HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh);

  template<typename Entity>
  std::shared_ptr<Array<Entity> > & get_entity()
  {
    auto * entity_data = std::get<tuple_type_index<Entity, ENTITY>>(entity_data_);
    return entity_data->get_entiy();
  }

  template<typename Entity>
  void delete_entity(Entity * e)
  {
    auto * entity_data = std::get<tuple_type_index<Entity, ENTITY>>(entity_data_);
    entity_data->delete_entity(e);
  }

  template<typename Entity>
  Entity & add_entity()
  {
    auto * entity_data = std::get<tuple_type_index<Entity, ENTITY>>(entity_data_);
    return entity_data->add_entity();
  }

  template<typename Entity>
  uint32_t number_of_entity()
  {
    auto * entity_data = std::get<tuple_type_index<Entity, ENTITY>>(entity_data_);
    return entity_data->number_of_data();
  }

  uint32_t number_of_boundary_edges() 
  { 
    uint32_t NE = edge_data_->number_of_data();
    uint32_t NHE = halfedge_data_->number_of_data();
    return NE*2-NHE;
  }

  /** @brief 清空网格, 但实际上没有释放内存 */
  void clear()
  {
    node_data_->clear();
    edge_data_->clear();
    cell_data_->clear();
    halfedge_data_->clear();
  }

  /** @brief 清空网格并释放内存 */
  void release()
  {
    node_data_->release();
    edge_data_->release();
    cell_data_->release();
    halfedge_data_->release();
  }

  void swap(HalfEdgeMeshBase & other)
  {
    std::swap(node_data_, other.node_data_);
    std::swap(edge_data_, other.edge_data_);
    std::swap(cell_data_, other.cell_data_);
    std::swap(halfedge_data_, other.halfedge_data_);
  }

  HalfEdgeMeshBase & operator = (const HalfEdgeMeshBase & other);

private:
  std::shared_ptr<EntityDataContainer<Node>> node_data_; 
  std::shared_ptr<EntityDataContainer<Edge>> edge_data_;
  std::shared_ptr<EntityDataContainer<Cell>> cell_data_;
  std::shared_ptr<EntityDataContainer<HalfEdge>> halfedge_data_;

  /** 所有的实体数据存在这里 */
  Container entity_data_;
};

HalfEdgeMeshBase::HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh): HalfEdgeMeshBase() 
{
  (*this) = mesh;
}

HalfEdgeMeshBase & HalfEdgeMeshBase::operator = (const HalfEdgeMeshBase & other)
{
  if(this != &other)
  {
    //clear();
    *node_data_ = *other.node_data_;
    *edge_data_ = *other.edge_data_;
    *cell_data_ = *other.cell_data_;
    *halfedge_data_ = *other.halfedge_data_;

    auto & halfedge_ = halfedge_data_->get_entity();
    for(auto it = halfedge_->begin(); it != halfedge_->end(); ++it)
    {
      (*it).set_next(&(halfedge_->get((*it).next()->index())));
      (*it).set_previous(&(halfedge_->get((*it).previous()->index())));
      (*it).set_opposite(&(halfedge_->get((*it).opposite()->index())));
    }

    auto & node_ = node_data_->get_entity(); 
    for(auto it = node_->begin(); it != node_->end(); ++it)
      (*it).set_halfedge(&(halfedge_->get((*it).halfedge()->index())));

    auto & edge_ = edge_data_->get_entity(); 
    for(auto it = edge_->begin(); it != edge_->end(); ++it)
      (*it).set_halfedge(&(halfedge_->get((*it).halfedge()->index())));

    auto & cell_ = cell_data_->get_entity(); 
    for(auto it = cell_->begin(); it != cell_->end(); ++it)
      (*it).set_halfedge(&(halfedge_->get((*it).halfedge()->index())));
  }
  return *this;
}


}


#endif /* _HALFEDGE_MESH_BASE_ */ 
