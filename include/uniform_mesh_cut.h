#ifndef UNIFORM_MESH_CUT_H
#define UNIFORM_MESH_CUT_H

#include "uniform_mesh.h"
#include "irregular_array2d.h"
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
        
public:
  /**
   * @brief 默认构造函数
   */
  template<typename... Args>
  UniformMeshCut(Args&&... args) : Base(std::forward<Args>(args)...), geo_() 
  {
    update();
  } 

  /**
   * @brief 更新 subcell_
   */
  void update()
  {
    uint32_t NC = Base::number_of_cells();
    uint32_t NB = Base::number_of_blocks();
    auto & data = subcell_.get_data();
    auto & start = subcell_.get_start_pos();

    data.resize(NC);

    std::fill(start.begin(), start.end(), 0.0);
    start.resize(NB+1, 0);
    auto & cell = *Base::get_cell();
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      start[idx+1]++;
    }

    for(uint32_t i = 1; i < NB+1; i++)
      start[i] += start[i-1];

    std::vector<uint32_t> I(NC, 0);
    for(auto & c : cell)
    {
      uint32_t idx = Base::find_point(c.barycenter());
      data[start[idx]+I[idx]] = &c;
      I[idx]++;
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
  uint32_t find_point(Point & p, Cell* & out, uint32_t & index) const 
  {
    uint32_t idx = Base::find_point(p);
    auto cellc = subcell_.get_row(idx);
    uint32_t L = subcell_.row_size(idx);

    std::vector<Point *> points(32, nullptr);
    for(uint32_t i = 0; i < L; i++)
    {
      Cell * c = cellc[i];
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

private:
  /** 
   * @brief 背景网格中单元的子单元
   */
  IrregularArray2D<Cell *> subcell_;

  GeometryUtils2D geo_;

};

} // namespace HEM


#endif // UNIFORM_MESH_CUT_H
