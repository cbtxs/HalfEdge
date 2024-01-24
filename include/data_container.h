#ifndef _DATA_CONTAINER_
#define _DATA_CONTAINER_

#include <vector>
#include <memory>
#include <algorithm>

#include "chunk_array.h"
#include "entity.h"

namespace HEM {

template<uint32_t CHUNK_SIZE = 1024u>
class DataContainer
{
public:
  template<typename T>
  using DataArray = ChunkArrayWithMark<T, CHUNK_SIZE>; 
  using MarkArray = ArrayBase::MarkArray; 

public:
  /**
   * @brief 构造函数
   */
  DataContainer(unsigned int size=0): data_number_(size), 
    is_free_(std::make_shared<MarkArray>(size)), 
    free_index_(), data_() 
  {
    //is_free_->set_false();
    is_free_->set_value(0);
  }

  uint32_t size() {return is_free_->size();}

  /**
   * @brief 添加一个数据
   */
  template<typename T>
  std::shared_ptr<DataArray<T> > add_data(std::string name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto it = std::find_if(data_.begin(), data_.end(), ff);
    if(it==data_.end())
    {
      std::shared_ptr<DataArray<T> > data = std::make_shared<DataArray<T> >(name, is_free_);
      data_.push_back(data);
      return data;
    }
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

  /**
   * @brief 删除一个数据
   */
  template<typename T>
  void delete_data(std::string name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto it = std::find_if(data_.begin(), data_.end(), ff);
    assert(it!=data_.end());
    *it = data_.back();
    data_.pop_back();
  }

  /**
   * @brief 获取某个数据
   */
  template<typename T>
  std::shared_ptr<DataArray<T> > get_data(std::string name)
  {
    auto ff = [&](std::shared_ptr<ArrayBase> & data) -> bool
    {
      return data->get_name().compare(name)==0;
    };
    auto it = std::find_if(data_.begin(), data_.end(), ff);
    assert(it!=data_.end());
    return std::dynamic_pointer_cast<DataArray<T> >(*it);
  }

  void clear()
  {
    data_.clear();
    free_index_.clear();
    data_number_ = 0;
  }

  void release()
  {
    data_.clear();
    free_index_.clear();
    data_number_ = 0;
  }

  DataContainer & operator = (DataContainer & other)
  {
    if (this != &other)
    {
      clear();
      data_number_ = other.data_number_;
      is_free_ = std::make_shared<MarkArray>(*(other.is_free_));
      free_index_ = other.free_index_;
      for(auto & otherData : other.data_)
      {
        std::shared_ptr<ArrayBase> data = std::make_shared<DataArray<uint8_t>>();
        otherData->copy_self(is_free_, data);
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

  std::shared_ptr<MarkArray> & is_free() {return is_free_;}

  /** 获取 data 的接口 */
  std::vector<std::shared_ptr<ArrayBase>> & data() {return data_;}

  void delete_index(uint32_t idx)
  {
    //is_free_->set_true(idx); 
    is_free_->get(idx) = 1; 
    free_index_.emplace_back(idx);
    data_number_--;
  }

  uint32_t add_index()
  {
    data_number_++;
    /** 当 free_index_ 为空时，没有可用指标所以要重新分配内存 */
    if(!free_index_.empty())
    {
      uint32_t index = free_index_.back();
      //is_free_->set_false(index);
      is_free_->get(index) = 0;
      free_index_.pop_back();
      return index; 
    }
    /** 当 free_index_ 非空时，启用最后一个可用指标 */
    else
    {
      for(auto & ptr : data_)
        ptr->resize(size()+1);
      //is_free_->push_back_false();
      is_free_->push_back(0);
      return size()-1;
    }
  }

  uint32_t number_of_data() { return data_number_;}

private:
  /** 实际上 data 的大小 */ 
  uint32_t data_number_;

  /** is_free_[i] = 1 代表第 i 个位置已经被释放了*/
  std::shared_ptr<MarkArray> is_free_;

  /** 存储已经被释放的位置的编号 */
  std::vector<uint32_t> free_index_;

  /** 数据使用 shared_ptr 管理 */
  std::vector<std::shared_ptr<ArrayBase>> data_;
};

template<typename Entity, uint32_t CHUNK_SIZE>
class EntityDataContainer : public DataContainer<CHUNK_SIZE>
{
public:
  using Base = DataContainer<CHUNK_SIZE>;
  using Self = EntityDataContainer<Entity, CHUNK_SIZE>;
  template<typename T>
  using DataArray = typename Base::template DataArray<T>;

public:
  EntityDataContainer(uint32_t size=0): Base(size)
  {
    entity_ = Base::template add_data<Entity>("entity");
    indices_ = Base::template add_data<uint32_t>("indices");
  }

  std::shared_ptr<DataArray<Entity> > get_entity() { return entity_;}

  const std::shared_ptr<DataArray<Entity> > get_entity() const  { return entity_;}

  std::shared_ptr<DataArray<uint32_t> > get_entity_indices() { return indices_;}

  const std::shared_ptr<DataArray<uint32_t> > get_entity_indices() const  { return indices_;}

  Entity & add_entity()
  {
    uint32_t idx = Base::add_index();
    entity_->get(idx).set_index(idx);
    return entity_->get(idx);
  }

  void delete_entity(Entity & e)
  {
    delete_index(e->index());
  }

  void update()
  {
    entity_ = Base::template get_data<Entity>("entity");
    indices_ = Base::template get_data<uint32_t>("indices");
    uint32_t N = 0;
    for(auto & idx : *indices_)
      idx = N++;
  }

  void clear()
  {
    Base::clear();
    entity_ = Base::template add_data<Entity>("entity");
    indices_ = Base::template add_data<uint32_t>("indices");
  }

  Self & operator = (Self & other)
  {
    Base::operator = (other);
    update();
    return *this;
  }

private:
  std::shared_ptr<DataArray<Entity> > entity_;
  std::shared_ptr<DataArray<uint32_t> > indices_;
};

}

#endif /* _DATA_CONTAINER_ */ 
