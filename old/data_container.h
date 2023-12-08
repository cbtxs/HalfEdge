#ifndef _DATA_CONTAINER_
#define _DATA_CONTAINER_

#include <vector>
#include <memory>
#include <algorithm>

#include "chunk_array.h"
#include "entity.h"

namespace HEM {

class DataContainer
{
public:
  template<typename T>
  using DataArray = ChunkArrayWithMark<T, 1024>; 

  using MarkArray = ArrayBase::MarkArray; 
public:
  /**
   * @brief 构造函数
   */
  DataContainer(unsigned int size=0): 
    is_free_("is_free", size), free_index_(), data_() 
  {
    is_free_.set_value(0);
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
      std::shared_ptr<DataArray<T> > data = std::make_shared<DataArray<T> >(name, is_free_, size());
      data_.push_back(data);
      return data;
    }
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

  /**
   * @brief 删除一个数据
   */
  template<typename T>
  void delete_data(std::string & name)
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

  void clear()
  {
    for(auto & data : data_)
      data->clear();
    free_index_.clear();
    data_number_ = 0;
  }

  void release()
  {
    for(auto & data : data_)
      data.reset();
    free_index_.clear();
    data_number_ = 0;
  }

  DataContainer & operator = (DataContainer & other)
  {
    if (this != &other)
    {
      clear();
      data_number_ = other.data_number_;
      is_free_ = other.is_free_;
      free_index_ = other.free_index_;
      for(auto & otherData : other.data_)
      {
        auto data = otherData->copy_self(is_free_); 
        data_.emplace_back(data);
      }
    }
    return *this;
  }

  /** @brief 交换数据所有权 TODO 需要 Test*/ 
  void swap(DataContainer & other)
  {
    std::swap(data_number_, other.data_number_);
    is_free_.swap(other.is_free_);
    std::swap(free_index_, other.free_index_);
    std::swap(data_, other.data_);
  }

  MarkArray & is_free() {return is_free_;}

  /** 获取 data 的接口 */
  std::vector<std::shared_ptr<ArrayBase>> & data() {return data_;}

  void delete_index(uint32_t idx)
  {
    is_free_[idx] = 1; 
    free_index_.emplace_back(idx);
    data_number_--;
  }

  uint32_t add_index()
  {
    data_number_++;
    /** 当 free_index_ 为空时，没有可用指标所以要重新分配内存 */
    if(free_index_.empty())
    {
      for(auto & ptr : data_)
        ptr->resize(size()+1);
      is_free_.push_back(0);
      return size();
    }
    /** 当 free_index_ 非空时，启用最后一个可用指标 */
    else
    {
      uint32_t index = free_index_.back();
      is_free_[index] = 0;
      free_index_.pop_back();
      return index; 
    }
  }

  uint32_t number_of_data() { return data_number_;}

private:
  /** 实际上 data 的大小 */ 
  uint32_t data_number_;

  /** is_free_[i] = 1 代表第 i 个位置已经被释放了*/
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
  using Base = DataContainer;
  using Self = EntityDataContainer<Entity>;
public:
  EntityDataContainer(uint32_t size=0) : DataContainer(size) 
  {
    entity = Base::add_data<DataArray<Entity> >("entity");
    indices = Base::add_data<DataArray<Entity> >("indices");
  }

  Entity & add_entity()
  {
    uint32_t index = Base::add_index();
    return entity->get(index);
  }

  void add_entity(Entity & e)
  {
    uint32_t index = Base::add_index();
    entity->get(index) = e;
  }

  void delete_entity(Entity * e)
  {
    Base::delete_index(e->halfedge()->template entity<Entity>());
  }

  void update_indices()
  {
    uint32_t k = 0;
    for(auto it = indices->begin(); it != indices->end(); ++it)
      *it = k++;
  }

  Self & operator=(const Self & other)
  {
    Base::operator=(other);
    indices = Base::get_data<uint32_t>("indices");
    entity = Base::get_data<Entity>("entity");
  }

  std::shared_ptr<DataArray<Entity> > & get_entity() {return entity;}

private:
  std::shared_ptr<DataArray<uint32_t> > indices;
  std::shared_ptr<DataArray<Entity> > entity;
};

}

#endif /* _DATA_CONTAINER_ */ 
