#ifndef _ENTITY_
#define _ENTITY_

#include <stdint.h>
#include <iostream>
#include "geometry.h"
#include "view.h"

namespace HEM
{

template<typename Traits>
class THalfEdge
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;

  using Point  = typename Traits::Point;
  using Vector = typename Traits::Vector;

public:
  THalfEdge(uint32_t index=0) : 
    next_(nullptr), prev_(nullptr), oppo_((HalfEdge*)this), cell_(nullptr), 
    edge_(nullptr), node_(nullptr), index_(index)  
  {}

  /** 数据接口 */
  HalfEdge * next() {return next_;}

  HalfEdge * previous() {return prev_;}

  HalfEdge * opposite() {return oppo_;}

  HalfEdge * halfedge() {return this;}

  HalfEdge * next(uint32_t i) 
  {
    HalfEdge * h = (HalfEdge*)this;
    for(uint32_t j = 0; j < i; j++) h = h->next();
    return h;
  }

  HalfEdge * previous(uint32_t i) 
  {
    HalfEdge * h = (HalfEdge*)this;
    for(uint32_t j = 0; j < i; j++) h = h->previous();
    return h;
  }

  Cell * cell() {return cell_;}

  Edge * edge() {return edge_;}

  Node * node() {return node_;}

  uint32_t & index() { return index_;}

  HalfEdge * next_oppo() {return next_->opposite();}

  HalfEdge * oppo_prev() {return oppo_->previous();}

  /** const 数据接口 */
  const HalfEdge * next() const {return next_;}

  const HalfEdge * previous() const {return prev_;}

  const HalfEdge * opposite() const {return oppo_;}

  const HalfEdge * halfedge() const {return this;}

  const Cell * cell() const {return cell_;}

  const Edge * edge() const {return edge_;}

  const Node * node() const {return node_;}

  const uint32_t & index() const { return index_;}

  const HalfEdge * next_oppo() const {return next_->opposite();}

  const HalfEdge * oppo_prev() const {return oppo_->previous();}

  template<typename Entity>
  Entity & entity()
  {
    if constexpr (std::is_same_v<Entity, Cell>)
      return cell_;
    else if constexpr (std::is_same_v<Entity, Edge>)
      return edge_;
    else if constexpr (std::is_same_v<Entity, Node>)
      return node_;
  }

  /** 设置数据 */
  void set_next(HalfEdge * next) {next_ = next;}

  void set_previous(HalfEdge * prev) {prev_ = prev;}

  void set_opposite(HalfEdge * oppo) {oppo_ = oppo;}

  void set_node(Node * node) {node_ = node;}

  void set_edge(Edge * edge) {edge_ = edge;}

  void set_cell(Cell * cell) {cell_ = cell;}

  void set_index(uint32_t index) {index_=index;}

  void reset(HalfEdge * next, HalfEdge * prev, HalfEdge * oppo, 
      Cell * cell, Edge * edge, Node * node, uint32_t index)
  {
    next_ = next; prev_ = prev; oppo_ = oppo;
    cell_ = cell; edge_ = edge; node_ = node;
    index_ = index;
  }

  HalfEdge & operator=(const HalfEdge& other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      next_ = other.next_; prev_ = other.prev_; oppo_ = other.oppo_;
      cell_ = other.cell_; edge_ = other.edge_; node_ = other.node_;
      index_ = other.index_;
    }
    return *this;
  }

  bool is_boundary() { return this==oppo_;}

  /** h 到 self 需要 next 的次数 */
  uint8_t distance(HalfEdge * h)
  {
    //assert(h->cell() == cell_)
    uint8_t n = 0;
    while(h != this) { h = h->next(); n++; }
    return n;
  }

  bool is_on_the_left(const Point & p);

  double length(); 

  Vector tangential() const ;

  Vector normal() const;

  Point barycenter();

private:
  /** 下一条半边, 上一条半边, 对边 */
  HalfEdge * next_, * prev_, * oppo_;

  Cell * cell_; /**< 所属单元*/ 
  Edge * edge_; /**< 所在边*/
  Node * node_; /**< 指向顶点 */
  uint32_t index_; /**< 半边的存储编号 */ 
};

template<typename Traits>
class TNode
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;

  using Point  = typename Traits::Point;
  using Vector = typename Traits::Vector;

  /** 定义邻接实体的迭代子 */
  template<typename H, typename N>
  class AdjNodeIterator : public AdjEntityIteratorBase<AdjNodeIterator<H, N>, H, N>
  {
  public:
    N & entity_imp(H * h) const { return *(h->previous()->node()); }
    H * next_imp(H * h) 
    { 
      if(h->node() != this->start_->node())
        return nullptr;
      h = h->next_oppo();
      if (h == this->start_) 
        return nullptr;
      else if (h->is_boundary())
        return h->next();
      else
        return h;
    }
  };

  template<typename H, typename E>
  class AdjEdgeIterator : public AdjEntityIteratorBase<AdjEdgeIterator<H, E>, H, E>
  {
  public:
    E & entity_imp(H * h) const { return *(h->edge()); }
    H * next_imp(H * h) 
    { 
      if(h->node() != this->start_->node())
        return nullptr;
      h = h->next_oppo();
      if (h == this->start_) 
        return nullptr;
      else
        return h;
    }
  };

  template<typename H, typename C>
  class AdjCellIterator : public AdjEntityIteratorBase<AdjCellIterator<H, C>, H, C>
  {
  public:
    C & entity_imp(H * h) const { return *(h->cell()); }
    H * next_imp(H * h) 
    { 
      h = h->next_oppo(); 
      if (h == this->start_ || h->is_boundary()) 
        return nullptr;
      else
        return h;
    }
  };

  /** 定义邻接实体的视图 */
  using AdjNodeView = AdjEntityViewBase<AdjNodeIterator<HalfEdge, Node>>;
  using ConstAdjNodeView = AdjEntityViewBase<AdjNodeIterator<const HalfEdge, const Node>>;

  using AdjEdgeView = AdjEntityViewBase<AdjEdgeIterator<HalfEdge, Edge>>;
  using ConstAdjEdgeView = AdjEntityViewBase<AdjEdgeIterator<const HalfEdge, const Edge>>;

  using AdjCellView = AdjEntityViewBase<AdjCellIterator<HalfEdge, Cell>>;
  using ConstAdjCellView = AdjEntityViewBase<AdjCellIterator<const HalfEdge, const Cell>>;

public:
  TNode(uint32_t index = -1): start_(nullptr), index_(index), coordinate_(0.0, 0.0) {}

  TNode(Point & coordinate, uint32_t index, HalfEdge * h = nullptr): 
    start_(h), index_(index), coordinate_(coordinate)
  {}

  /** 获取数据 */
  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  Point & coordinate() { return coordinate_;}

  const HalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  const Point & coordinate() const { return coordinate_;}

  /** 设置数据 */
  void set_coordinate(const Point & p) { coordinate_ = p;}

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(const Point & p, uint32_t index, HalfEdge * h) 
  { 
    index_ = index; 
    start_ = h; 
    coordinate_ = p;
  }

  uint32_t adj_cell(Cell ** n2c);

  AdjNodeView adj_nodes() const { return AdjNodeView(start_); }

  AdjEdgeView adj_edges() const { return AdjEdgeView(start_); }

  AdjCellView adj_cells() const { return AdjCellView(start_); }

  Node & operator=(const Node & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
      coordinate_ = other.coordinate_;
    }
    return *this;
  }

private:
  HalfEdge * start_;
  uint32_t index_;
  Point coordinate_;
};

/**
 * @brief Edge 类
 */
template<typename Traits>
class TEdge
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;

  using Point  = typename Traits::Point;
  using Vector = typename Traits::Vector;

public:
  TEdge(uint32_t index=-1, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, HalfEdge * h) { index_ = index; start_ = h; }

  /** 获取边的邻接关系 */
  void get_top(Node ** e2n, Cell ** e2c, uint8_t * e2cidx);

  uint32_t adj_cell(Cell ** n2c);

  uint32_t adj_node(Node ** n2n);

  void vertices(Point ** vertices);

  bool has_node(Node * n) const
  {
    return start_->node() == n || start_->previous()->node() == n;
  }

  Edge & operator=(const Edge & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  Vector tangential() const;

  Vector normal() const;

  Point barycenter() const;

  double length(); 

private:
  HalfEdge * start_;
  uint32_t index_;
};

/**
 * @brief Cell 类
 */
template<typename Traits>
class TCell
{
public:
  using Node = typename Traits::Node;
  using Edge = typename Traits::Edge;
  using Cell = typename Traits::Cell;
  using HalfEdge = typename Traits::HalfEdge;

  using Point  = typename Traits::Point;
  using Vector = typename Traits::Vector;

  /** 定义邻接实体的迭代子 */
  template<typename H, typename N>
  class AdjNodeIterator : public AdjEntityIteratorBase<AdjNodeIterator<H, N>, H, N>
  {
  public:
    N & entity_imp(H * h) const { return *(h->node()); }
    H * next_imp(H * h) 
    { 
      h = h->next();
      if (h == this->start_) 
        return nullptr;
      else
        return h;
    }
  };

  template<typename H, typename E>
  class AdjEdgeIterator : public AdjEntityIteratorBase<AdjEdgeIterator<H, E>, H, E>
  {
  public:
    E & entity_imp(H * h) const { return *(h->edge()); }
    H * next_imp(H * h) 
    { 
      h = h->next();
      if (h == this->start_) 
        return nullptr;
      else
        return h;
    }
  };

  template<typename H, typename C>
  class AdjCellIterator : public AdjEntityIteratorBase<AdjCellIterator<H, C>, H, C>
  {
  public:
    C & entity_imp(H * h) const { return *(h->opposite()->cell()); }
    H * next_imp(H * h) 
    { 
      h = h->next();
      if (h == this->start_) 
        return nullptr;
      else
        return h;
    }
  };

  /** 定义邻接实体的视图 */
  using AdjNodeView = AdjEntityViewBase<AdjNodeIterator<HalfEdge, Node>>;
  using ConstAdjNodeView = AdjEntityViewBase<AdjNodeIterator<const HalfEdge, const Node>>;

  using AdjEdgeView = AdjEntityViewBase<AdjEdgeIterator<HalfEdge, Edge>>;
  using ConstAdjEdgeView = AdjEntityViewBase<AdjEdgeIterator<const HalfEdge, const Edge>>;

  using AdjCellView = AdjEntityViewBase<AdjCellIterator<HalfEdge, Cell>>;
  using ConstAdjCellView = AdjEntityViewBase<AdjCellIterator<const HalfEdge, const Cell>>;

public:
  TCell(uint32_t index = -1, HalfEdge * h = nullptr): start_(h), index_(index) {}

  HalfEdge * halfedge() { return start_; }

  uint32_t & index() { return index_; }

  const HalfEdge * halfedge() const { return start_; }

  const uint32_t & index() const { return index_; }

  void set_halfedge(HalfEdge * h) { start_ = h;}

  void set_index(uint32_t index) { index_ = index;}

  void reset(uint32_t index, HalfEdge * h) { index_ = index; start_ = h; }

  /** 获取单元的邻接关系 */
  uint32_t get_top(Node ** c2n, Edge ** c2e, Cell ** c2c);

  /** 获取单元的邻接关系 */
  uint32_t adj_edge(Edge ** cell2edge);

  /** 获取单元的邻接关系 */
  uint32_t adj_node(Node ** cell2node);

  /** 获取单元的邻接关系 */
  uint32_t adj_cell(Cell ** cell2cell);

  /** 获取单元的顶点 */
  uint32_t vertices(Point ** vertices);

  /** 单元第 i 个邻接边 */
  Edge * adj_edge(uint32_t i) const 
  { 
    return start_->next(i)->edge(); 
  }

  /** 单元第 i 个邻接点 */
  Node * adj_node(uint32_t i) const 
  { 
    return start_->previous()->next(i)->node(); 
  }

  /** 单元第 i 个邻接单元 */
  Cell * adj_cell(uint32_t i) const
  { 
    return start_->next(i)->opposite()->cell(); 
  }

  /** 单元第 i 个顶点的坐标 */
  Point * vertex(uint32_t i) const 
  { 
    return &(adj_node(i)->coordinate());
  }

  AdjNodeView adj_nodes() { return AdjNodeView(start_->previous()); }

  AdjEdgeView adj_edges() { return AdjEdgeView(start_); }

  AdjCellView adj_cells() { return AdjCellView(start_); }

  ConstAdjNodeView adj_nodes() const { return ConstAdjNodeView(start_->previous()); }

  ConstAdjEdgeView adj_edges() const { return ConstAdjEdgeView(start_); }

  ConstAdjCellView adj_cells() const { return ConstAdjCellView(start_); }

  Cell & operator=(const Cell & other)
  {
    if (this != &other) /**< 避免自我赋值 */
    {
      start_ = other.start_;
      index_ = other.index_;
    }
    return *this;
  }

  /** 返回一个内点 */
  Point inner_point() const;

  Point barycenter() const
  {
    uint8_t n = 1;
    Point p = start_->node()->coordinate();
    for(HalfEdge * h = start_->next(); h != start_; h = h->next(), n++)
      p += h->node()->coordinate(); 
    return p/n; 
  }

  double area();

private:
  HalfEdge * start_;
  uint32_t index_;
};

/** HalfEdge 的一些内联函数 */

/** 判断一个点是不是在半边的左边 */
template<typename Traits>
inline bool THalfEdge<Traits>::is_on_the_left(const Point & p)
{
  Vector v0 = node_->coordinate()-prev_->node()->coordinate();
  Vector v1 = p-prev_->node()->coordinate();
  return v0.cross(v1)>0;
}

template<typename Traits>
inline typename THalfEdge<Traits>::Point THalfEdge<Traits>::barycenter()
{
  return (node()->coordinate() + previous()->node()->coordinate())*0.5;
}

template<typename Traits>
inline typename THalfEdge<Traits>::Vector THalfEdge<Traits>::tangential() const 
{
  return node()->coordinate()-previous()->node()->coordinate();
}

template<typename Traits>
inline typename THalfEdge<Traits>::Vector THalfEdge<Traits>::normal() const 
{
  Vector t = tangential();
  return Vector(-t.y, t.x);
}

template<typename Traits>
inline double THalfEdge<Traits>::length()
{
  return tangential().length();
}


template<typename Traits>
inline uint32_t TNode<Traits>::adj_cell(Cell ** n2c)
{
  uint32_t N = 1;
  n2c[0] = start_->cell(); 
  for(HalfEdge * h = start_->next_oppo(); h != start_ && !h->is_boundary(); 
      h = h->next_oppo())
  {
    n2c[N++] = h->cell(); 
  }
  return N;
}

/** Edge 的一些内联函数 */
template<typename Traits>
inline void TEdge<Traits>::get_top(Node ** e2n, Cell ** e2c, uint8_t * e2cidx)
{
  e2n[0] = start_->previous()->node();
  e2n[1] = start_->node();
  e2c[0] = start_->cell(); 
  e2c[1] = start_->cell(); 
  e2cidx[0] = start_->distance(start_->cell()->halfedge());
  e2cidx[1] = start_->opposite()->distance(start_->opposite()->cell()->halfedge());
}

template<typename Traits>
inline uint32_t TEdge<Traits>::adj_cell(Cell ** e2c)
{
  e2c[0] = start_->cell(); 
  e2c[1] = start_->opposite()->cell(); 
  return 2;
}

template<typename Traits>
inline uint32_t TEdge<Traits>::adj_node(Node ** e2n)
{
  e2n[0] = start_->previous()->node();
  e2n[1] = start_->node();
  return 2;
}

template<typename Traits>
inline void TEdge<Traits>::vertices(Point ** vertices)
{
  vertices[0] = &start_->previous()->node()->coordinate();
  vertices[1] = &start_->node()->coordinate();
}

template<typename Traits>
inline double TEdge<Traits>::length()
{
  return start_->length();
}

template<typename Traits>
inline typename TEdge<Traits>::Vector TEdge<Traits>::tangential() const 
{
  return start_->tangential();
}

template<typename Traits>
inline typename TEdge<Traits>::Vector TEdge<Traits>::normal() const 
{
  return start_->normal(); 
}

template<typename Traits>
inline typename TEdge<Traits>::Point TEdge<Traits>::barycenter() const
{
  return start_->barycenter(); 
}


/** Cell 的一些内联函数 */
template<typename Traits>
uint32_t TCell<Traits>::get_top(Node ** c2n, Edge ** c2e, Cell ** c2c)
{
  uint32_t N = 0;
  c2n[N] = start_->node();
  c2e[N] = start_->edge(); 
  c2c[N++] = start_->cell(); 
  for(HalfEdge * h = start_->next(); h != start_; h = h->next())
  {
    c2n[N] = h->node();
    c2e[N] = h->edge(); 
    c2c[N++] = h->opposite()->cell(); 
  }
  return N;
}

template<typename Traits>
uint32_t TCell<Traits>::adj_edge(Edge ** c2e)
{
  uint32_t N = 0;
  c2e[N++] = start_->edge(); 
  for(HalfEdge * h = start_->next(); h != start_; h = h->next())
    c2e[N++] = h->edge(); 
  return N;
}

template<typename Traits>
uint32_t TCell<Traits>::adj_node(Node ** c2n)
{
  uint32_t N = 0;
  c2n[N++] = start_->previous()->node(); 
  for(HalfEdge * h = start_; h != start_->previous(); h = h->next())
    c2n[N++] = h->node(); 
  return N;
}

template<typename Traits>
uint32_t TCell<Traits>::adj_cell(Cell ** c2c)
{
  uint32_t N = 0;
  c2c[N++] = start_->cell(); 
  for(HalfEdge * h = start_->next(); h != start_; h = h->next())
    c2c[N++] = h->cell(); 
  return N;
}

template <typename Traits>
uint32_t TCell<Traits>::vertices(Point ** vertices)
{
  uint32_t N = 0;
  vertices[N++] = &start_->previous()->node()->coordinate();
  for(HalfEdge * h = start_; h != start_->previous(); h = h->next())
    vertices[N++] = &h->node()->coordinate();
  return N;
}

template<typename Traits>
double TCell<Traits>::area()
{
  Point p0 = start_->previous()->node()->coordinate();
  Vector v0 = start_->node()->coordinate()-p0;
  Vector v1 = start_->next()->node()->coordinate()-p0;
  double a = v0.cross(v1);
  for(HalfEdge * h = start_->next()->next(); h != start_->previous(); h = h->next())
  {
    v0 = v1;
    v1 = h->node()->coordinate()-p0;
    a += v0.cross(v1);
  }
  return a/2;
}

/** 返回一个内点 */
template<typename Traits>
typename TCell<Traits>::Point TCell<Traits>::inner_point() const
{
  Point p(0, 0);
  for(HalfEdge* h = start_->next(); h != start_; h = h->next())
  {
    if(h->tangential().cross(h->next()->tangential())>0)
    {
      p  = h->previous()->node()->coordinate();
      p += h->node()->coordinate();
      p += h->next()->node()->coordinate();
      return p/3.0;
    }
  }
  return p;
}


}


#endif /* _ENTITY_ */ 
