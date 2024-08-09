#ifndef UNIFORM_MESH_CUT_H
#define UNIFORM_MESH_CUT_H

#include "uniform_mesh.h"
#include "geometry_utils.h"
#include <vector>

namespace HEM 
{

template<int D>
class UniformMeshCut : public UniformMesh<D>
{
public:
  using Base = UniformMesh<D>;
  using Self = UniformMeshCut<D>;

  using Cell = typename Base::Cell;
  using Edge = typename Base::Edge;
  using Node = typename Base::Node;
  using HalfEdge = typename Base::HalfEdge; 

  using Point  = typename Base::Point;
  using Vector = typename Base::Vector;

  using SubCellArray = std::vector<std::vector<Cell *> >;
        
public:
  /**
   * @brief 默认构造函数
   */
  template<typename... Args>
  UniformMeshCut(Args&&... args) : Base(std::forward<Args>(args)...), geo_(), 
  subcell_(Base::number_of_blocks()) 
  {
    update_subcell();
  } 

  /**
   * @brief 更新 subcell_
   */
  void update_subcell()
  {
    uint32_t NB = Base::number_of_blocks();

    subcell_.clear();
    subcell_.resize(NB);

    auto & cell = *Base::get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      subcell_[idx].push_back(&c);
    }
  }

  /**
   * @brief 查找点所在的单元
   * @param p 点
   * @param out 输出单元
   * @param index 若点在单元边上，返回边的索引，若点在单元顶点上，返回顶点的索引
   * @return 点相对于单元的位置
   * @retval 0 点在顶点上
   * @retval 1 点在单元边上
   * @retval 2 点在单元内部
   * @retval 3 点在网格外面
   */
  uint32_t find_point(const Point & p, Cell* & out, uint32_t & index) const 
  {
    uint32_t idx = Base::find_point(p);
    const auto & cellc = subcell_[idx];

    for(auto & c : cellc)
    {
      std::vector<Point *> points(32, nullptr);
      int N = c->vertices(points.data());
      points.resize(N);
      uint32_t flag = geo_.relative_position_of_point_and_polygon(points, p, index);
      if(flag!=3)
      {
        out = c;
        return flag;
      }
    }
    return 3;
  }

  /**
   * @brief 查找点所在的单元, 不需要返回 index
   */
  uint32_t find_point(const Point & p, Cell* & out) const
  {
      uint32_t dummy_index; /**< 定义一个临时的 index 变量 */
      return find_point(p, out, dummy_index);
  }

  Cell * find_point(const Point & p) const
  {
    Cell * out;
    find_point(p, out);
    return out;
  }

private:
  /** 
   * @brief 背景网格中单元的子单元
   */
  SubCellArray subcell_;

  GeometryUtils2D geo_;

};

} // namespace HEM


#endif // UNIFORM_MESH_CUT_H
