#ifndef _HALFEDGE_MESH_BASE_
#define _HALFEDGE_MESH_BASE_

#include <memory>
#include <map>
#include <algorithm>

#include "entity.h"
#include "data_container.h"
#include "geometry.h"

namespace HEM
{
class HalfEdgeMeshBase
{
public:
  using DataContainer = DataContainer<1024u>;
  template<typename T>
  using Array = DataContainer::DataArray<T>;

public:
  HalfEdgeMeshBase()
  {
    nodes_ = node_data_->add_data<Node>("node");
    node_indices_ = node_data_->add_data<uint32_t>("indices");
    node_coordinate_ = node_data_->add_data<Point>("coordinate");

    edges_ = edge_data_->add_data<Edge>("edge");
    edge_indices_ = edge_data_->add_data<uint32_t>("indices");

    cells_ = cell_data_->add_data<Cell>("cell");
    cell_indices_ = cell_data_->add_data<uint32_t>("indices");

    halfedges_ = halfedge_data_->add_data<HalfEdge>("halfedge");
    next_halfedge_ =  halfedge_data_->add_data<uint32_t>("next");
    prev_halfedge_ = halfedge_data_->add_data<uint32_t>("previous");
    oppo_halfedge_ =  halfedge_data_->add_data<uint32_t>("opposite");
    halfedge_to_cell_ = halfedge_data_->add_data<uint32_t>("cell");
    halfedge_to_edge_ = halfedge_data_->add_data<uint32_t>("edge");
    halfedge_to_node_ = halfedge_data_->add_data<uint32_t>("node");
  }

  /** 复制构造函数 */
  HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh);

  /** 以单元为中心的网格 */
  HalfEdgeMeshBase(double * node, uint32_t * cell, uint32_t NN, uint32_t NC, uint32_t NV);

  /** 实体接口 */
  Array<Node> & get_node() { return *nodes_; }

  Array<Edge> & get_edge() { return *edges_; }

  Array<Cell> & get_cell() { return *cells_; }

  Array<HalfEdge> & get_halfedge() { return *(halfedges_); }

  /** 删除实体 */
  void delete_node(Node & n)
  {
    auto & h2n = halfedge_to_node(); 
    node_data_->delete_index(h2n[n.halfedge()->index()]);
  }

  void delete_edge(Edge & e)
  {
    auto & h2e = halfedge_to_edge(); 
    edge_data_->delete_index(h2e[e.halfedge()->index()]);
  }

  void delete_cell(Cell & c)
  {
    auto & h2c = halfedge_to_cell(); 
    cell_data_->delete_index(h2c[c.halfedge()->index()]);
  }

  void delete_halfedge(HalfEdge & h) { halfedge_data_->delete_index(h.index()); }

  /** 添加实体 */
  Node & add_node() { return nodes_->get(node_data_->add_index()); }
  Edge & add_edge() { return edges_->get(edge_data_->add_index()); }
  Cell & add_cell() { return cells_->get(cell_data_->add_index()); }
  HalfEdge & add_halfedge() { return halfedges_->get(halfedge_data_->add_index()); }

  /** 实体个数 */
  uint32_t number_of_nodes() { return node_data_->number_of_data(); }
  uint32_t number_of_edges() { return edge_data_->number_of_data(); }
  uint32_t number_of_cells() { return cell_data_->number_of_data(); }
  uint32_t number_of_halfedges() { return halfedge_data_->number_of_data(); }
  uint32_t number_of_boundary_edges() 
  { 
    uint32_t NE = edge_data_->number_of_data();
    uint32_t NHE = halfedge_data_->number_of_data();
    return NE*2-NHE;
  }

  /** 半边数据的接口 */
  Array<uint32_t> & next_halfedge() { return *next_halfedge_; }
  Array<uint32_t> & prev_halfedge() { return *prev_halfedge_; }
  Array<uint32_t> & oppo_halfedge() { return *oppo_halfedge_; }
  Array<uint32_t> & halfedge_to_cell() { return *halfedge_to_cell_;}
  Array<uint32_t> & halfedge_to_edge() { return *halfedge_to_edge_;}
  Array<uint32_t> & halfedge_to_node() { return *halfedge_to_node_;}

  Array<Point> & node_coordinate() { return *node_coordinate_;}

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
    update();
  }

  void update()
  {
    nodes_ = node_data_->get_data<Node>("node");
    node_indices_ = node_data_->get_data<uint32_t>("indices");
    node_coordinate_ = node_data_->get_data<Point>("coordinate");

    edges_ = edge_data_->get_data<Edge>("edge");
    edge_indices_ = edge_data_->get_data<uint32_t>("indices");

    cells_ = cell_data_->get_data<Cell>("cell");
    cell_indices_ = cell_data_->get_data<uint32_t>("indices");

    halfedges_ = halfedge_data_->get_data<HalfEdge>("halfedge");
    next_halfedge_ =  halfedge_data_->get_data<uint32_t>("next");
    prev_halfedge_ = halfedge_data_->get_data<uint32_t>("previous");
    oppo_halfedge_ =  halfedge_data_->get_data<uint32_t>("opposite");
    halfedge_to_cell_ = halfedge_data_->get_data<uint32_t>("cell");
    halfedge_to_edge_ = halfedge_data_->get_data<uint32_t>("edge");
    halfedge_to_node_ = halfedge_data_->get_data<uint32_t>("node");
  }

  HalfEdgeMeshBase & operator = (const HalfEdgeMeshBase & other);

private:

  /** 实体数据集合 */
  std::shared_ptr<DataContainer> node_data_; 
  std::shared_ptr<DataContainer> edge_data_;
  std::shared_ptr<DataContainer> cell_data_;
  std::shared_ptr<DataContainer> halfedge_data_;

  std::shared_ptr<Array<Node>> nodes_;
  std::shared_ptr<Array<uint32_t>> node_indices_;
  std::shared_ptr<Array<Point>>  node_coordinate_;

  std::shared_ptr<Array<Edge>>  edges_;
  std::shared_ptr<Array<uint32_t>>  edge_indices_;

  std::shared_ptr<Array<Cell>>   cells_;
  std::shared_ptr<Array<uint32_t>>  cell_indices_;

  std::shared_ptr<Array<HalfEdge>> halfedges_;
  std::shared_ptr<Array<uint32_t>> next_halfedge_;
  std::shared_ptr<Array<uint32_t>> prev_halfedge_;
  std::shared_ptr<Array<uint32_t>> oppo_halfedge_;
  std::shared_ptr<Array<uint32_t>> halfedge_to_cell_;
  std::shared_ptr<Array<uint32_t>> halfedge_to_edge_;
  std::shared_ptr<Array<uint32_t>> halfedge_to_node_;
};

HalfEdgeMeshBase::HalfEdgeMeshBase(const HalfEdgeMeshBase & mesh): HalfEdgeMeshBase() 
{
  //clear();
  *node_data_ = *mesh.node_data_;
  *edge_data_ = *mesh.edge_data_;
  *cell_data_ = *mesh.cell_data_;
  *halfedge_data_ = *mesh.halfedge_data_;

  auto & halfedge_ = get_halfedge();
  auto & node_ = get_node(); 
  for(auto it = node_.begin(); it != node_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  auto & edge_ = get_edge(); 
  for(auto it = edge_.begin(); it != edge_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));

  auto & cell_ = get_cell(); 
  for(auto it = cell_.begin(); it != cell_.end(); ++it)
    it->set_halfedge(&(halfedge_[it->halfedge()->index()]));
}

/** 
 * @brief 以单元为中心的网格数据为参数的构造函数
 * @note 关键在于生成半边的顶点
 */
HalfEdgeMeshBase::HalfEdgeMeshBase(double * node, uint32_t * cell, 
    uint32_t NN, uint32_t NC, uint32_t NV):HalfEdgeMeshBase()
{
  auto & node_ = get_node();
  auto & coordinate = node_coordinate();
  for(uint32_t i = 0; i < NN; i++)
  {
    add_node();
    coordinate[i] = Point(node[i*2], node[i*2+1]);
  }

  auto & h2n = halfedge_to_node();
  auto & h2e = halfedge_to_edge();
  auto & h2c = halfedge_to_cell();
  auto & next = next_halfedge();
  auto & prev = prev_halfedge();
  auto & oppo = oppo_halfedge();
  
  std::vector<std::map<uint32_t, HalfEdge*> > n2n(NN);
  std::vector<HalfEdge*> c2h(NV);
  std::vector<uint32_t> preidx(NV);
  for(uint32_t i = 0; i < NV; i++)
    preidx[i] = (NV-1+i)%NV;
  /** 生成 halfedge_to_cell, halfedge_to_node, next_halfedge **/
  for(uint32_t i = 0; i < NC; i++)
  {
    Cell & c = add_cell();
    for(uint32_t j = 0; j < NV; j++)
    {
      c2h[j] = &add_halfedge();
      h2c[c2h[j]->index()] = i;
      h2n[c2h[j]->index()] = cell[i*NV+j];
      node_[cell[i*NV+j]].set_halfedge(c2h[j]);
    }
    for(uint32_t j = 0; j < NV; j++)
    {
      next[c2h[j]->index()] = c2h[preidx[j]]->index();
      n2n[cell[i*NV+preidx[j]]].emplace(std::pair<uint32_t, HalfEdge * >(j, c2h[j]));
    }
    c.set_halfedge(c2h[0]);
  }
  /** 生成 previous halfedge */
  auto & halfedge_ = get_halfedge();
  for(HalfEdge & h : halfedge_)
    prev[next[h.index()]] = h.index();

  uint32_t uint32_max = -1;
  h2e.set_value(uint32_max);
  /** 生成 edge, opposite_halfedge, halfedge_to_edge, node 的半边 */
  for(uint32_t i = 0; i < NN; i++)
  {
    node_[n2n[i].begin()->first] = n2n[i].begin()->second;
    for(auto & it : n2n[i])
    {
      if(h2e[it.second->index()] == uint32_max)
      {
        const auto & oppoit = n2n[it.first].find(i);
        if(oppoit != n2n[it.first].end())
        {
          oppo[it.second->index()] = oppoit->second->index();
          oppo[oppoit->second->index()] = it.second->index();
        }
        else
          oppo[it.second->index()] = oppoit->second->index();
        Edge & e = add_edge();
        e.set_halfedge(it.second);
        h2e[it.second->index()] = number_of_edges()-1;
        h2e[oppoit->second->index()] = h2e[it.second->index()];
      }
    }
  }
}

HalfEdgeMeshBase & HalfEdgeMeshBase::operator = (const HalfEdgeMeshBase & other)
{
  if(this != &other)
  {
    HalfEdgeMeshBase m(other);
    swap(m);
  }
  return *this;
}



}


#endif /* _HALFEDGE_MESH_BASE_ */ 
