#ifndef Vector_hpp
#define Vector_hpp


#include <memory>
#include <string>
#include <vector>

#include "mesh_data_container.hpp"

namespace OF
{

template <typename T>
class Vector : public DataSetBase_
{

private:

  std::vector<T> data_; //数据容器

  /**
   * @brief 检查新指标是否超过数据容器的大小，如果大于等于则增加内存
   */
  inline void manage_index(uint32_t index) override
  {
    while (index >= uint32_t(data_.size()))
      data_.push_back(T());
  }

  uint32_t number_of_datas()
  {

  }

public:

  Vector(DataSetContainerBase* container, const std::string& name) :
    DataSetBase_(container, name)
  {
    data_.reserve(512u);
  }

  ~Vector() override
  {
  }

  inline T& operator[](uint32_t index)
  {
    //TODO: assert
    return data_[index];
  }

  inline const T& operator[](uint32_t index) const
  {
    //TODO: assert
    return data_[index];
  }

  inline void fill(const T& value)
  {
    std::fill(data_.begin(), data_.end(), value);
  }

  inline void swap(Vector<T>* vector)
  {
    if (vector->container_ == this->container_)//only swap from same container
      data_.swap(vector->data_);
  }

  inline void copy(Vector<T>* vector)
  {
    if (vector->container_ == this->container_)//only swap from same container
      data_ = vector->data_;
  }

  inline void clear() override
  {
    data_.clear();
    data_.reserve(512u);
    data_.shrink_to_fit();
  }

  inline std::shared_ptr<DataSetBase_> create_in(DataSetContainerBase& dst) const override
  {
    using DataSetContainer = DataSetContainerT<Vector>;
    DataSetContainer* pcontainer = dynamic_cast<DataSetContainer*>(&dst);
    if (pcontainer)
    { // 如果是有效指针
      auto mdata = pcontainer->get_mesh_data<T>(name_);
      if (!mdata) // 如果不存在，则增加网格数据
        mdata = pcontainer->add_mesh_data<T>(name_);
      return mdata; 
    }
    return nullptr;
  }

  inline void copy(const DataSetBase_& src) override
  {
    const Vector<T>* src_vector = dynamic_cast<const Vector<T>*>(&src);
    if (src_vector)
    {
      data_ = src_vector->data_;
    }
  }

  inline const void* data_pointer() const
  {
    return &data_[0];
  }

  class const_iterator
  {
    const Vector<T>* vector_;
    uint32_t index_;
  public:
    inline const_iterator(const Vector<T>* vector, uint32_t index) :
      vector_(vector), index_(index)
    {
    }

    inline const_iterator(const const_iterator& it) : 
      vector_(it.vector_), index_(it.index_)
    {
    }

    inline const_iterator& operator=(const const_iterator& it)
    {
      vector_ = it.vector_;
      index_ = it.index_;
      return *this;
    }

    inline bool operator!=(const_iterator it) const
    { //TODO: 断言 vector 相同
      return index_ != it.index_; 
    }

    inline const_iterator& operator++()
    {
      index_ = vector_->container_->next_index(index_);
      return *this;
    }

    inline const T& operator*() const
    {
      return vector_->operator[](index_);
    }

    inline uint32_t index()
    {
      return index_;
  };

  inline const_iterator begin() const
  {
    return const_iterator(this, this->container_->begin_index());
  }

  inline const_iterator end() const
  {
    return const_iterator(this, this->container_->end_index());
  }

 	class iterator
	{
		Vector<T>* vector_;
		uint32_t index_;
	public:
		inline iterator(Vector<T>* vector, uint32_t index) 
      : vector_(vector), index_(index)
		{
		}

		inline iterator(const iterator& it) 
      : vector_(it.vector_), index_(it.index_)
		{
		}

		inline iterator& operator=(const iterator& it)
		{
			vector_ = it.vector_;
			index_ = it.index_;
			return *this;
		}

		inline bool operator!=(iterator it) const
		{
			// assert vector_ == it.vector_
			return index_ != it.index_;
		}

		inline iterator& operator++()
		{
			index_ = vector_->container_->next_index(index_);
			return *this;
		}

		inline T& operator*() const
		{
			return vector_->operator[](index_);
		}

		inline uint32_t index()
		{
			return index_;
		}
	};

	inline iterator begin()
	{
		return iterator(this, this->container_->begin_index());
	}

	inline iterator end()
	{
		return iterator(this, this->container_->end_index());
	}
};

} // end of namespace OF

#endif // end of Vector_hpp
