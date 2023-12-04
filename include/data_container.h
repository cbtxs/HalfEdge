#ifndef _DATA_CONTAINER_
#define _DATA_CONTAINER_

#include <vector>
#include <memory>
#include <algorithm>

#include "chunk_array.h"

class DataContainer
{
public:
  template<typename T>
  using DataArray = ChunkArray<T>; 

  using MarkArray = ChunkArrayBool<1024u>; 

public:
  /**
   * @brief 构造函数
   */
  DataContainer(unsigned int size=0): 
    is_free_("is_free", size), free_index_(), datas_() 
  {
    is_free_.set_true();
  }

  uint32_t size() {return is_free_.size();}

  /**
   * @brief 添加一个数据
   */
  template<typename T>
  std::shared_ptr<DataArray<T> > & add_data(std::string name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto & it = std::find_if(datas_.begin(), datas_.end(), ff);
    if(it==datas_.end())
    {
      std::shared_ptr<DataArray<T> > data = std::make_shared<DataArray<T> >(name, is_free_);
      data->resize(size());
      datas_.push_back(data);
      return data;
    }
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

  /**
   * @brief 删除一个数据
   */
  template<typename T>
  void remove_data(std::string & name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto & it = std::find_if(datas_.begin(), datas_.end(), ff);
    assert(it!=datas_.end() && "Container don't have data named \" " + name + "\" ");
    *it = datas_.back();
    datas_.pop_back();
  }

  /**
   * @brief 获取某个数据
   */
  template<typename T>
  std::shared_ptr<DataArray<T> > & get_data(std::string & name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto & it = std::find_if(datas_.begin(), datas_.end(), ff);
    assert(it!=datas_.end() && "Container don't have data named \" " + name + "\" ");
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

private:
  /** is_used_[i] = false 代表第 i 个位置已经被释放了*/
  MarkArray is_free_;

  /** 存储已经被释放的位置的编号 */
  std::vector<uint32_t> free_index_;

  /** 数据使用 shared_ptr 管理 */
  std::vector<std::shared_ptr<ArrayBase>> datas_;
};

#endif /* _DATA_CONTAINER_ */ 
