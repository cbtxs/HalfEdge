/** 
 * Author : Chunyu Chen 
 * Date   : 2023-12-3
 * Brief  : A base class for halfedge mesh data structure.
 */
#ifndef _HalfEdge_MESH_BASE_
#define _HalfEdge_MESH_BASE_

#include <functional>
#include <map>
#include <memory>

#include "geometry.h"
#include "data_container.h"

namespace HEM
{

/**
 * @brief 半边网格基类
 */
template<typename Traits>
class HalfEdgeMeshBase
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;
  using Point  = typename Traits::Point;
  using Vector = typename Traits::Vector;

  constexpr static const int Dim    = Traits::Dim;

  using Self = HalfEdgeMeshBase<Traits>;

  using NodeDataContainer = EntityDataContainer<Node, 1024u>;
  using EdgeDataContainer = EntityDataContainer<Edge, 1024u>;
  using CellDataContainer = EntityDataContainer<Cell, 1024u>;
  using HalfEdgeDataContainer = EntityDataContainer<HalfEdge, 1024u>;

  template<typename T>
  using Array = typename NodeDataContainer::Base::template DataArray<T>;

public:
  HalfEdgeMeshBase() 
  {
    node_data_ptr_ = std::make_shared<NodeDataContainer>();
    edge_data_ptr_ = std::make_shared<EdgeDataContainer>();
    cell_data_ptr_ = std::make_shared<CellDataContainer>();
    halfedge_data_ptr_ = std::make_shared<HalfEdgeDataContainer>();
  }

  /** 复制构造函数 */
  HalfEdgeMeshBase(const Self & mesh);

  /** 以单元为中心的网格 */
  HalfEdgeMeshBase(double * node, uint32_t * cell, uint32_t NN, uint32_t NC, uint32_t NV):
    HalfEdgeMeshBase()
  {
    reinit(node, cell, NN, NC, NV);
  }

  /** 以单元为中心的网格重新初始化 */
  void reinit(double * node, uint32_t * cell, uint32_t NN, uint32_t NC, uint32_t NV);

  /** 实体接口 */
  std::shared_ptr<Array<Node>> get_node() { return node_data_ptr_->get_entity(); }

  std::shared_ptr<Array<Edge>> get_edge() { return edge_data_ptr_->get_entity(); }

  std::shared_ptr<Array<Cell>> get_cell() { return cell_data_ptr_->get_entity(); }

  std::shared_ptr<Array<HalfEdge>> get_halfedge() { return halfedge_data_ptr_->get_entity(); }

  template<typename Data>
  std::shared_ptr<Array<Data>> get_node_data(std::string dname) 
  { 
    return node_data_ptr_->template get_data<Data>(dname); 
  }

  template<typename Data>
  std::shared_ptr<Array<Data>> get_edge_data(std::string dname) 
  { 
    return edge_data_ptr_->template get_data<Data>(dname); 
  }

  template<typename Data>
  std::shared_ptr<Array<Data>> get_cell_data(std::string dname) 
  { 
    return cell_data_ptr_->template get_data<Data>(dname); 
  }

  template<typename Data>
  std::shared_ptr<Array<Data>> add_node_data(std::string dname) 
  { 
    return node_data_ptr_->template add_data<Data>(dname); 
  }

  template<typename Data>
  std::shared_ptr<Array<Data>> add_edge_data(std::string dname) 
  { 
    return edge_data_ptr_->template add_data<Data>(dname); 
  }

  template<typename Data>
  std::shared_ptr<Array<Data>> add_cell_data(std::string dname) 
  { 
    return cell_data_ptr_->template add_data<Data>(dname); 
  }

  template<typename Entity>
  std::shared_ptr<Array<Entity>> get_entity() 
  { 
    if constexpr (std::is_same_v<Entity, Cell>)
      return cell_data_ptr_->get_entity();
    else if constexpr (std::is_same_v<Entity, Edge>)
      return edge_data_ptr_->get_entity();
    else if constexpr (std::is_same_v<Entity, Node>)
      return node_data_ptr_->get_entity();
    else if constexpr (std::is_same_v<Entity, HalfEdge>)
      return halfedge_data_ptr_->get_entity();
  }

  /** 删除实体 */
  void delete_node(Node & n) { node_data_ptr_->delete_index(n.index());}

  void delete_edge(Edge & e) { edge_data_ptr_->delete_index(e.index());}

  void delete_cell(Cell & c) { cell_data_ptr_->delete_index(c.index());}

  void delete_halfedge(HalfEdge & h) { halfedge_data_ptr_->delete_index(h.index()); }

  /** 添加实体 */
  Node & add_node() { return node_data_ptr_->add_entity(); }

  Edge & add_edge() { return edge_data_ptr_->add_entity(); }

  Cell & add_cell() { return cell_data_ptr_->add_entity(); }

  HalfEdge & add_halfedge() { return halfedge_data_ptr_->add_entity(); }

  /** 获取实体编号 */
  std::shared_ptr<Array<uint32_t>> get_node_indices() { return node_data_ptr_->get_entity_indices(); }

  std::shared_ptr<Array<uint32_t>> get_edge_indices() { return edge_data_ptr_->get_entity_indices(); }

  std::shared_ptr<Array<uint32_t>> get_cell_indices() { return cell_data_ptr_->get_entity_indices(); }

  std::shared_ptr<Array<uint32_t>> get_halfedge_indices() { return halfedge_data_ptr_->get_entity_indices(); }

  /** @brief 实体迭代函数 */
  template<typename Entity>
  void for_each_entity(const std::function<bool(Entity & )> & );

  /** @brief 并行的实体迭代函数 */
  template<typename Entity>
  void parallel_for_each_entity(const std::function<bool(Entity & )> & );

  /** 实体个数 */
  uint32_t number_of_nodes() { return node_data_ptr_->number_of_data(); }

  uint32_t number_of_edges() { return edge_data_ptr_->number_of_data(); }

  uint32_t number_of_cells() { return cell_data_ptr_->number_of_data(); }

  uint32_t number_of_halfedges() { return halfedge_data_ptr_->number_of_data(); }

  uint32_t number_of_boundary_edges() 
  { 
    uint32_t NE = edge_data_ptr_->number_of_data();
    uint32_t NHE = halfedge_data_ptr_->number_of_data();
    return NE*2-NHE;
  }

  /**
   * @brief 获取一个 Box 将整个网格框起来
   */
  std::array<double, 4> get_box()
  {
    std::array<double, 4> box = {-1e-100};
    auto & node = *get_node();
    for(const auto & n : node)
    {
      box[0] = box[0]<n.coordinate().x ? box[0] : n.coordinate().x; 
      box[2] = box[2]>n.coordinate().x ? box[2] : n.coordinate().x; 
      box[1] = box[1]<n.coordinate().y ? box[1] : n.coordinate().y; 
      box[3] = box[3]>n.coordinate().y ? box[3] : n.coordinate().y; 
    }
    box[2] = box[2] - box[0]; box[3] = box[3] - box[1];
    box[0] = box[0]-box[2]*0.05; box[1] = box[1]-box[3]*0.05;
    box[2] *= 1.1; box[3] *= 1.1;
    return box;
  }

  /** 加密半边 */
  void splite_halfedge(HalfEdge * h)
  {
    splite_halfedge(h, (h->node()->coordinate() + 
          h->previous()->node()->coordinate())*0.5);
  }

  /** 加密半边 */
  void splite_halfedge(HalfEdge * h, const Point & p)
  {
    Node & n = add_node();
    Edge & e = add_edge();
    HalfEdge & h0 = add_halfedge();

    e.set_halfedge(&h0);
    n.reset(p, n.index(), &h0);

    h0.reset(h, h->previous(), &h0, h->cell(), &e, &n, h0.index());

    h->previous()->set_next(&h0);
    h->set_previous(&h0);
    h->edge()->set_halfedge(h);
    if(!h->is_boundary())
    {
      HalfEdge & h1 = add_halfedge();
      HalfEdge * o = h->opposite();

      h->set_opposite(&h1);
      h0.set_opposite(o);
      h1.reset(o, o->previous(), h, o->cell(), o->edge(), &n, h1.index());

      o->previous()->set_next(&h1);
      o->set_previous(&h1);
      o->set_opposite(&h0);
      o->set_edge(&e);
    }
  }

  /** 连接 h0 和 h1 的顶点分割单元 c */
  void splite_cell(Cell * c0, HalfEdge * h0, HalfEdge * h1)
  {
    HalfEdge * nh0 = &add_halfedge();
    HalfEdge * nh1 = &add_halfedge();
    Edge * e = &add_edge();
    Cell * c1 = &add_cell();

    c0->set_halfedge(nh0);
    c1->set_halfedge(nh1);
    e->set_halfedge(nh0);

    nh0->reset(h1->next(), h0, nh1, c0, e, h1->node(), nh0->index());
    nh1->reset(h0->next(), h1, nh0, c1, e, h0->node(), nh1->index());

    h0->next()->set_previous(nh1);
    h1->next()->set_previous(nh0);
    h0->set_next(nh0);
    h1->set_next(nh1);
    for(HalfEdge * h = nh1->next(); h != nh1; h = h->next())
      h->set_cell(c1); 
  }

  /** @brief 清空网格, 但实际上没有释放内存 */
  void clear()
  {
    node_data_ptr_->clear();
    edge_data_ptr_->clear();
    cell_data_ptr_->clear();
    halfedge_data_ptr_->clear();
  }

  /** 
   * @brief 清空网格并释放内存 
   */
  void release()
  {
    node_data_ptr_->release();
    edge_data_ptr_->release();
    cell_data_ptr_->release();
    halfedge_data_ptr_->release();
  }

  void swap(HalfEdgeMeshBase & other)
  {
    std::swap(node_data_ptr_, other.node_data_);
    std::swap(edge_data_ptr_, other.edge_data_);
    std::swap(cell_data_ptr_, other.cell_data_);
    std::swap(halfedge_data_ptr_, other.halfedge_data_);
  }

  void update()
  {
    node_data_ptr_->update();
    edge_data_ptr_->update();
    cell_data_ptr_->update();
    halfedge_data_ptr_->update();
  }

  Self & operator = (const Self & other);

private:
  /** 实体数据集合 */
  std::shared_ptr<NodeDataContainer> node_data_ptr_; 
  std::shared_ptr<EdgeDataContainer> edge_data_ptr_;
  std::shared_ptr<CellDataContainer> cell_data_ptr_;
  std::shared_ptr<HalfEdgeDataContainer> halfedge_data_ptr_;
};

}

#include "imp/halfedge_mesh.inl"

#endif /* _HalfEdge_MESH_BASE_ */ 
