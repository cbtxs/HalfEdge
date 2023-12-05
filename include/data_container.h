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
    is_free_("is_free", size), free_index_(), data_() 
  {
    is_free_.set_false();
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
    auto & it = std::find_if(data_.begin(), data_.end(), ff);
    if(it==data_.end())
    {
      std::shared_ptr<DataArray<T> > data = std::make_shared<DataArray<T> >(name, is_free_);
      data->resize(size());
      data_.push_back(data);
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
    auto & it = std::find_if(data_.begin(), data_.end(), ff);
    assert(it!=data_.end() && "Container don't have data named \" " + name + "\" ");
    *it = data_.back();
    data_.pop_back();
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
    auto & it = std::find_if(data_.begin(), data_.end(), ff);
    assert(it!=data_.end() && "Container don't have data named \" " + name + "\" ");
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

  void delete_index(uint32_t idx)
  {
    is_free_.set_true(idx); 
    free_index_.emplace_back(idx);
  }

  uint32_t add_index()
  {
    /** 当 free_index_ 为空时，没有可用指标所以要重新分配内存 */
    if(free_index_.empty())
    {
      for(auto & ptr : data_)
        ptr->resize(size()+1);
      is_free_.push_back_false();
      return size();
    }
    /** 当 free_index_ 非空时，启用最后一个可用指标 */
    else
    {
      uint32_t index = free_index_.back();
      is_free_.set_false(index);
      free_index_.pop_back();
      return index; 
    }
  }

private:
  /** is_free_[i] = true 代表第 i 个位置已经被释放了*/
  MarkArray is_free_;

  /** 存储已经被释放的位置的编号 */
  std::vector<uint32_t> free_index_;

  /** 数据使用 shared_ptr 管理 */
  std::vector<std::shared_ptr<ArrayBase>> data_;
};

template<typename Entity>
class EntityDataContainer : public DataContainer
{
public:
  EntityDataContainer(uint32_t size) : DataContainer(size) 
  {
    entity = std::make_shared<DataArray<Entity> >(this->is_free_, "entity", size);
  }

  void add_entity(Entity & e)
  {
    entity.push_back(e);
    this->add_index();
  }

private:
  std::shared_ptr<DataArray<Entity> > entity;
};

#endif /* _DATA_CONTAINER_ */ 
