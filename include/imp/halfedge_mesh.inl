#include<map>

namespace HEM{

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
