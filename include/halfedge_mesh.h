#ifndef _HalfEdge_MESH_BASE_
#define _HalfEdge_MESH_BASE_

#include <memory>
#include <map>

#include "data_container.h"
#include "tuple_type_index.h"
#include "geometry.h"

namespace HEM
{
template<typename N, typename E, typename C, typename H>
class HalfEdgeMeshBase
{
public:
  using Node = N;
  using Edge = E;
  using Cell = C;
  using HalfEdge = H;

  using Self = HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>;
  using DataContainer = DataContainer<1024u>;
  template<typename T>
  using Array = DataContainer::DataArray<T>;

public:
  HalfEdgeMeshBase()
  {
    nodes_ptr_ = node_data_ptr_->add_data<Node>("node");
    edges_ptr_ = edge_data_ptr_->add_data<Edge>("edge");
    cells_ptr_ = cell_data_ptr_->add_data<Cell>("cell");
    halfedges_ptr_ = halfedge_data_ptr_->add_data<HalfEdge>("halfedge");

    node_indices_ptr_ = node_data_ptr_->add_data<uint32_t>("indices");
    edge_indices_ptr_ = edge_data_ptr_->add_data<uint32_t>("indices");
    cell_indices_ptr_ = cell_data_ptr_->add_data<uint32_t>("indices");
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
  Array<Node> & get_node() { return *nodes_ptr_; }
  Array<Edge> & get_edge() { return *edges_ptr_; }
  Array<Cell> & get_cell() { return *cells_ptr_; }
  Array<HalfEdge> & get_halfedge() { return *halfedges_ptr_; }

  /** 删除实体 */
  void delete_node(Node & n) { node_data_ptr_->delete_index(n.index());}
  void delete_edge(Edge & e) { edge_data_ptr_->delete_index(e.index());}
  void delete_cell(Cell & c) { cell_data_ptr_->delete_index(c.index());}
  void delete_halfedge(HalfEdge & h) { halfedge_data_ptr_->delete_index(h.index()); }

  /** 添加实体 */
  Node & add_node() { return nodes_ptr_->get(node_data_ptr_->add_index()); }
  Edge & add_edge() { return edges_ptr_->get(edge_data_ptr_->add_index()); }
  Cell & add_cell() { return cells_ptr_->get(cell_data_ptr_->add_index()); }
  HalfEdge & add_halfedge() { return halfedges_ptr_->get(halfedge_data_ptr_->add_index()); }

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

  /** 加密半边 */
  void splite_halfedge(HalfEdge * h)
  {
    splite_halfedge(h, (h->node()->coordinate() + 
          h->previous()->node()->coordinate())*0.5);
  }

  /** 加密半边 */
  void splite_halfedge(HalfEdge * h, const Point & p)
  {
    HalfEdge * o = h->opposite();

    Node & n = add_node();
    Edge & e = add_edge();
    HalfEdge & h0 = add_halfedge();
    HalfEdge & h1 = add_halfedge();

    e.set_halfedge(o);
    n.reset(n.index(), &h0, p);

    h0.reset(h, h->previous(), o, h->cell(), &e, &n, h0.index());
    h1.reset(o, o->previous(), h, o->cell(), o->edge(), &n, h1.index());

    h->previous()->set_next(&h0);
    h->set_previous(&h0);
    h->set_opposite(&h1);
    o->previous()->set_next(&h1);
    o->set_previous(&h1);
    o->set_opposite(&h0);
  }

  /** 连接 h0 和 h1 的顶点分割单元 c */
  void splite_cell(Cell * c0, HalfEdge * h0, HalfEdge * h1)
  {
    HalfEdge * nh0 = &add_halfedge();
    HalfEdge * nh1 = &add_halfedge();
    Edge * e = &add_edge();
    Cell * c1 = &add_cell();

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

  /** @brief 清空网格并释放内存 */
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
    update();
  }

  void update()
  {
    nodes_ptr_ = node_data_ptr_->get_data<Node>("node");
    edges_ptr_ = edge_data_ptr_->get_data<Edge>("edge");
    cells_ptr_ = cell_data_ptr_->get_data<Cell>("cell");
    halfedges_ptr_ = halfedge_data_ptr_->get_data<HalfEdge>("halfedge");
    node_indices_ptr_ = node_data_ptr_->get_data<uint32_t>("indices");
    edge_indices_ptr_ = edge_data_ptr_->get_data<uint32_t>("indices");
    cell_indices_ptr_ = cell_data_ptr_->get_data<uint32_t>("indices");
  }

  Self & operator = (const Self & other);

private:
  /** 实体数据集合 */
  std::shared_ptr<DataContainer> node_data_ptr_; 
  std::shared_ptr<DataContainer> edge_data_ptr_;
  std::shared_ptr<DataContainer> cell_data_ptr_;
  std::shared_ptr<DataContainer> halfedge_data_ptr_;

  /** 方便使用的指针 */
  std::shared_ptr<Array<Node>> nodes_ptr_;
  std::shared_ptr<Array<Edge>>  edges_ptr_;
  std::shared_ptr<Array<Cell>>   cells_ptr_;
  std::shared_ptr<Array<HalfEdge>> halfedges_ptr_;
  std::shared_ptr<Array<uint32_t>> node_indices_ptr_;
  std::shared_ptr<Array<uint32_t>>  edge_indices_ptr_;
  std::shared_ptr<Array<uint32_t>>  cell_indices_ptr_;
};

/** 复制构造函数 */
template<typename Node, typename Edge, typename Cell, typename HalfEdge>
HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>::HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh): 
  HalfEdgeMeshBase() 
{
  //clear();
  *node_data_ptr_ = *mesh.node_data_ptr_;
  *edge_data_ptr_ = *mesh.edge_data_ptr_;
  *cell_data_ptr_ = *mesh.cell_data_ptr_;
  *halfedge_data_ptr_ = *mesh.halfedge_data_ptr_;

  auto & node_ = get_node(); 
  auto & edge_ = get_edge(); 
  auto & cell_ = get_cell(); 
  auto & halfedge_ = get_halfedge();
  for(auto it = node_.begin(); it != node_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto it = edge_.begin(); it != edge_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto it = cell_.begin(); it != cell_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  for(auto & h : halfedge_)
  {
    h.reset(halfedge_[h.next()->index()], 
            halfedge_[h.previous()->index()], 
            halfedge_[h.opposite()->index()], 
            cell_[h.cell()->index()], 
            edge_[h.edge()->index()], 
            node_[h.node()->index()], 
            h.index()); 
  }
  update();
}

/** 
 * @brief 以单元为中心的网格数据为参数的构造函数
 * @note 关键在于生成半边的顶点
 */
template<typename Node, typename Edge, typename Cell, typename HalfEdge>
void HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>::reinit(double * node, uint32_t * cell, 
    uint32_t NN, uint32_t NC, uint32_t NV)
{
  clear();
  auto & node_ = get_node();
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
      c2h[j]->reset(c2h[nextidx[j]], c2h[previdx[j]], c2h[j], &c, nullptr, &node_[i*NV+j], i*NV+j);
      n2n[cell[i*NV+previdx[j]]].emplace(std::pair<uint32_t, HalfEdge * >(j, c2h[j]));
    }
    c.reset(i, c2h[0]);
  }

  /** 生成 edge, opposite_halfedge, halfedge_to_edge, node 的半边 */
  for(uint32_t i = 0; i < NN; i++)
  {
    node_[n2n[i].begin()->first].set_halfedge(n2n[i].begin()->second);
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
}

template<typename Node, typename Edge, typename Cell, typename HalfEdge>
HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge> & HalfEdgeMeshBase<Node, Edge, Cell, HalfEdge>::operator = (const Self & other)
{
  if(this != &other)
  {
    HalfEdgeMeshBase m(other);
    swap(m);
  }
  return *this;
}


}


#endif /* _HalfEdge_MESH_BASE_ */ 
