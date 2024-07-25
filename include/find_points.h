#ifndef FIND_POINT_H
#define FIND_POINT_H

#include <memory>
#include <vector>

#include "irregular_array2d.h"

namespace HEM
{

template <typename Derived, typename Mesh>
class FindPointAlgorithmBase
{
public:
  using Point = typename Mesh::Point;

public:
  FindPointAlgorithmBase(std::shared_ptr<Mesh> mesh) : mesh_(mesh) {}
public:

  void find_point(const Point & query_point)
  {
    static_cast<Derived*>(this)->find_point_imp(query_point);
  }

  void update()
  {
    static_cast<Derived*>(this)->update_imp();
  }

  std::shared_ptr<Mesh> get_mesh() { return mesh_; }

protected:
  std::shared_ptr<Mesh> mesh_;
};

/**
 * @brief The KDTreeFindPointAlg class
 * This class is used to find the nearest point in a mesh using a KDTree
 * @todo TODO
 */
template <typename Mesh>
class KDTreeFindPointAlg : public FindPointAlgorithmBase<KDTreeFindPointAlg<Mesh>, Mesh>
{
public:
  using Base = FindPointAlgorithmBase<KDTreeFindPointAlg<Mesh>, Mesh>;

public:
  KDTreeFindPointAlg(const Mesh& mesh) : Base(mesh) {}
};

/**
 * @brief The UniformMesh2dFindPointAlg class
 * This class is used to find the nearest point in a 2D uniform mesh
 */
template<typename Mesh>
class UniformMesh2dFindPointAlg : public HEM::FindPointAlgorithmBase<UniformMesh2dFindPointAlg<Mesh>, Mesh> 
{
public:
  using Base = FindPointAlgorithmBase<UniformMesh2dFindPointAlg<Mesh>, Mesh>; 
  using Self = UniformMesh2dFindPointAlg<Mesh>;

  using Cell = typename Mesh::Cell;
  using Edge = typename Mesh::Edge;
  using Node = typename Mesh::Node;
  using HalfEdge = typename Mesh::HalfEdge; 

  using Point  = typename Mesh::Point;
  using Vector = typename Mesh::Vector;

private:
  /** 
   * @brief 背景网格中单元的子单元
   */
  IrregularArray2D<Cell *> subcell_;

public:

  /**
   * @brief 默认构造函数
   */
  UniformMesh2dFindPointAlg(std::shared_ptr<Mesh> mesh) : Base(mesh), 
    subcell_() 
  {
    update_imp();
  }

  
  uint32_t find_point_in_uniform_mesh(const Point & p)
  {
    return mesh_->find_point(p);
  }

  /**
   * @brief 更新 subcell_
   */
  void update_imp()
  {
    uint32_t NC = mesh_->number_of_cells();
    auto & data = subcell_.get_data();
    auto & start = subcell_.get_start_pos();

    data.resize(NC);

    std::fill(start.begin(), start.end(), 0.0);
    start.resize(NC+1, 0);
    auto & cell = *(mesh_->get_cell());
    for(auto & c : cell)
    {
      uint32_t idx = meshfind_point(c.barycenter());
      start[idx+1]++;
    }

    for(uint32_t i = 1; i < NC+1; i++)
      start[i] += start[i-1];

    std::vector<uint32_t> I(NC, 0);
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      data[start[idx]+I[idx]] = &c;
      I[idx]++;
    }
  }

  Cell * find_point(Point & p, bool maybe_on_the_edge = true)
  {
    HalfEdge * h = nullptr;
    uint32_t idx = Base::find_point(p);
    auto cellc = subcell_[idx];
    for(auto & c : cellc)
    {
      bool flag = is_on_the_polygon(p, c, h, maybe_on_the_edge);
      if(flag==1)
        return c;
    }
    return nullptr;
  }


}

} // namespace HEM

#endif // FIND_POINT_H
