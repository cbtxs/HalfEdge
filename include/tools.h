#ifndef TOOLS_H
#define TOOLS_H

#include <ios>
#include <vector>
#include <span>
#include <iostream>

namespace HEM
{

/*
 * @brief 一个不规则二维数组的类
 */
template <typename T>
class IrregularArray2D 
{
public:
  /** 
   * @brief 构造函数
   */
  IrregularArray2D() = default;

  /** 
   * @brief 获取行数
   */
  size_t len() const 
  {
    return start_pos.size() - 1;
  }

  /**
   * @brief 获取指定行的大小
   */
  size_t row_size(size_t row) const
  {
    return start_pos[row + 1] - start_pos[row];
  }

  /** 
   * @brief 获取指定行的视图 
   */
  std::span<T> get_row(size_t row) 
  {
    auto loc = start_pos[row];
    auto rowsize = row_size(row); 
    return std::span<T>(data.data() + loc, rowsize);
  }

  /** 
   * @brief 获取指定行的视图 
   */
  std::span<T> operator[](size_t rowIndex) 
  {
    return get_row(rowIndex);
  }

  /** 
   * @brief 获取所有数据
   */
  std::vector<T> & get_data()
  {
    return data;
  }

  /** 
   * @brief 获取每行的起始位置
   */
  std::vector<size_t> & get_start_pos()
  {
    return start_pos;
  }

private:
    std::vector<T> data;                 /**< 存储所有数据 */
    std::vector<size_t> start_pos;    /**< 存储每行的起始位置 */
};
}
#endif // TOOLS_H
