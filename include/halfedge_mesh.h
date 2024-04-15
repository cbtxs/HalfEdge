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

#include "data_container.h"
#include "geometry.h"

namespace HEM
{
template<typename Traits>
class HalfEdgeMeshBase
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;
  using Dim  = typename Traits::Dim;
  using Point = typename Node::Point;

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

template<typename N, typename E, typename C, typename H, int D>
class HalfEdgeMesh_Traits
{
  using Node = N; 
  using Edge = N; 
  using Cell = N; 
  using HalfEdge = H;
  constexpr static const int Dim  = D; 
};

/** 复制构造函数 */
template<typename Traits>
HalfEdgeMeshBase<Traits>::HalfEdgeMeshBase(
    const HalfEdgeMeshBase & mesh): HalfEdgeMeshBase() 
{
  //clear();
  *node_data_ptr_     = *mesh.node_data_ptr_;
  *edge_data_ptr_     = *mesh.edge_data_ptr_;
  *cell_data_ptr_     = *mesh.cell_data_ptr_;
  *halfedge_data_ptr_ = *mesh.halfedge_data_ptr_;

  auto & node_ = *get_node(); 
  auto & edge_ = *get_edge(); 
  auto & cell_ = *get_cell(); 
  auto & halfedge_ = *get_halfedge();

  for(auto it = node_.begin(); it != node_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto it = edge_.begin(); it != edge_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto it = cell_.begin(); it != cell_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto & h : halfedge_)
  {
    h.reset(&halfedge_[h.next()->index()], 
            &halfedge_[h.previous()->index()], 
            &halfedge_[h.opposite()->index()], 
            &cell_[h.cell()->index()], 
            &edge_[h.edge()->index()], 
            &node_[h.node()->index()], 
            h.index()); 
  }
  update();
}

/** 
 * @brief 以单元为中心的网格数据为参数的构造函数
 * @note 关键在于生成半边的顶点
 */
template<typename Traits>
void HalfEdgeMeshBase<Traits>::reinit(double * node, uint32_t * cell, 
    uint32_t NN, uint32_t NC, uint32_t NV)
{
  clear();
  auto & node_ = *(get_node());
  for(uint32_t i = 0; i < NN; i++)
  {
    Node * n = &add_node();
    node_[i].set_coordinate(Point(node[2*i], node[2*i+1]));
  }

  std::vector<std::map<uint32_t, HalfEdge*> > n2n(NN);
  std::vector<HalfEdge*> c2h(NV);
  std::vector<uint32_t> previdx(NV);
  std::vector<uint32_t> nextidx(NV);
  for(uint32_t i = 0; i < NV; i++)
  {
    previdx[i] = (NV-1+i)%NV;
    nextidx[i] = (NV+1+i)%NV;
  }
  /** 生成 halfedge_to_cell, halfedge_to_node, next_halfedge **/
  for(uint32_t i = 0; i < NC; i++)
  {
    Cell & c = add_cell();
    for(uint32_t j = 0; j < NV; j++)
      c2h[j] = &add_halfedge();
    for(uint32_t j = 0; j < NV; j++)
    {
      c2h[j]->reset(c2h[nextidx[j]], 
                    c2h[previdx[j]], 
                    c2h[j], 
                    &c, 
                    nullptr, 
                    &node_[cell[i*NV+j]], 
                    c2h[j]->index());
      n2n[cell[i*NV+previdx[j]]].emplace(std::pair<uint32_t, HalfEdge * >(cell[i*NV+j], c2h[j]));
      node_[cell[i*NV+j]].set_halfedge(c2h[j]);
    }
    c.reset(i, c2h[0]);
  }

  /** 生成 edge, opposite_halfedge, halfedge_to_edge, node 的半边 */
  for(uint32_t i = 0; i < NN; i++)
  {
    for(auto & it : n2n[i])
    {
      if(it.second->edge() == nullptr)
      {
        const auto & oppoit = n2n[it.first].find(i);
        if(oppoit != n2n[it.first].end())
        {
          it.second->set_opposite(oppoit->second);
          oppoit->second->set_opposite(it.second);
        }
        Edge & e = add_edge();
        e.reset(number_of_edges()-1, it.second);
        it.second->set_edge(&e);
        it.second->opposite()->set_edge(&e);
      }
    }
  }
  for(auto & n : node_)
  {
    HalfEdge * h = n.halfedge(); 
    if(!h->is_boundary())
    {
      h = h->opposite()->previous();
      while(!h->is_boundary() && h != n.halfedge())
        h = h->opposite()->previous();
    }
    n.set_halfedge(h);
  }
  update();
}

template<typename Traits>
HalfEdgeMeshBase<Traits> & HalfEdgeMeshBase<Traits>::operator = (const Self & other)
{
  if(this != &other)
  {
    HalfEdgeMeshBase m(other);
    swap(m);
  }
  return *this;
}

/**
 * @brief 节点迭代函数
 */
template<typename Traits>
template<typename Entity>
void HalfEdgeMeshBase<Traits>::for_each_entity(
    const std::function<bool(Entity & )> & f)
{
  auto & entitys = *(get_entity<Entity>());
  for(auto & e : entitys)
  {
    f(e);
  }
}

/**
 * @brief 节点迭代函数
 */
template<typename Traits>
template<typename Entity>
void HalfEdgeMeshBase<Traits>::parallel_for_each_entity(
    const std::function<bool(Entity & )> & f)
{
  auto & entitys = *get_entity<Entity>();
  #pragma omp parallel for
  for(auto & e : entitys)
    f(e);
}
}




#endif /* _HalfEdge_MESH_BASE_ */ 
