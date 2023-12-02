#ifndef BLOCK_VECTOR_HPP_
#define BLOCK_VECTOR_HPP_

#include <memory>
#include <string>
#include <vector>

#include "data_set_container.hpp"

namespace OF 
{

/**
 * @class BlockVector 
 *
 * @brief 分块序列容器
 *
 * @note 采用一个 std::vector 存储一组内存块的指针，通过索引的求商、求余来索引
 * 数据。
 *
 */
template <typename T>
class BlockVector : public DataSetBase_
{
public:
	static const uint32_t BLOCK_SIZE = 1024u; // 这里也许可以定制

private:
	std::vector<T*> blocks_;
	uint32_t capacity_;

  /**
   * @brief 检查新指标是否超过数据容器的容量，如果大于等于则增加内存
   */
	inline void manage_index(uint32_t index) override
	{
		while (index >= capacity_)
		{
			blocks_.push_back(new T[BLOCK_SIZE]());
			capacity_ = uint32_t(blocks_.size()) * BLOCK_SIZE;
		}
	}

public:
	BlockVector(DataSetContainerBase* container, const std::string& name)
		: DataSetBase_(container, name)
	{
		blocks_.reserve(512u);
		capacity_ = 0u;
	}

	~BlockVector() override
	{
		for (auto block: blocks_)
			delete[] block; // 释放每个 block 内存
	}

	inline T& operator[](uint32_t index)
	{
		// assert index < capacity_, "index out of bounds");
		return blocks_[index / BLOCK_SIZE][index % BLOCK_SIZE];
	}

	inline const T& operator[](uint32_t index) const
	{
		// assert (index < capacity_, "index out of bounds");
		return blocks_[index / BLOCK_SIZE][index % BLOCK_SIZE];
	}

	inline void fill(const T& value)
	{
		for (auto block : blocks_)
			std::fill(block, block + BLOCK_SIZE, value);
	}

	inline void swap(BlockVector<T>* v)
	{
		if (v->container_ == this->container_) // only swap from same container
			blocks_.swap(v->blocks_);
	}

	inline void copy(BlockVector<T>* v)
	{
		if (v->container_ == this->container_) // only copy from same container
			for (uint32_t i = 0; i < uint32_t(blocks_.size()); ++i)
				std::copy(v->blocks_[i], v->blocks_[i] + BLOCK_SIZE, blocks_[i]);
	}

	inline void clear() override
	{
		for (auto block : blocks_)
			delete[] block;
		blocks_.clear();
		capacity_ = 0;
	}

	inline std::shared_ptr<DataSetBase_> create_in(DataSetContainerBase& dst) const override
	{
		using DataSetContainer = DataSetContainerT<BlockVector>;
		DataSetContainer* container = dynamic_cast<DataSetContainer*>(&dst);
		if (container)
		{
			auto mdata = container->get_data_set<T>(name_);
			if (!mdata)
				mdata = container->add_data_set<T>(name_);
			return mdata;
		}
		return nullptr;
	}

	inline void copy(const DataSetBase_& v) override
	{
		const BlockVector<T>* pv = dynamic_cast<const BlockVector<T>*>(&v);
		if (pv)
		{
			// assert(src_ca->capacity_ == capacity_, "Copy from src with different capacity");
			for (uint32_t i = 0; i < uint32_t(pv->blocks_.size()); ++i)
				std::copy(pv->blocks_[i], pv->blocks_[i] + BLOCK_SIZE, blocks_[i]);
		}
	}

	inline uint32_t number_of_blocks() const
	{
		return uint32_t(blocks_.size());
	}

	inline std::vector<const void*> block_pointers() const
	{
		std::vector<const void*> pointers;
		pointers.reserve(uint32_t(blocks_.size()));
		for (auto block : blocks_)
			pointers.push_back(block);
		return pointers;
	}

  /**
   * @class const_iterator
   * @brief
   * @note 它是一个 BlockVector  指针和一个编码的封装
   */
	class const_iterator
	{
		const BlockVector<T>* pv_;
		uint32_t index_;

	public:

		inline const_iterator(const BlockVector<T>* pv, uint32_t index) 
      : pv_(pv), index_(index)
		{
		}

		inline const_iterator(const const_iterator& it) 
      : pv_(it.pv_), index_(it.index_)
		{
		}

		inline const_iterator& operator=(const const_iterator& it)
		{
			pv_ = it.pv_;
			index_ = it.index_;
			return *this;
		}

		inline bool operator!=(const_iterator it) const
		{
			// 断言是否为同一个数据，这里可能要考虑性能，少一次比较
			return index_ != it.index_;
		}

		inline const_iterator& operator++()
		{
			index_ = pv_->container_->next_index(index_);
			return *this;
		}

		inline const T& operator*() const
		{
			return pv_->operator[](index_);
		}
		inline uint32_t index()
		{
			return index_;
		}
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
		BlockVector<T>* pv_;
		uint32_t index_;

	public:
		inline iterator(BlockVector<T>* pv, uint32_t index) : pv_(pv), index_(index)
		{
		}

		inline iterator(const iterator& it) : pv_(it.pv_), index_(it.index_)
		{
		}

		inline iterator& operator=(const iterator& it)
		{
			pv_ = it.pv_;
			index_ = it.index_;
			return *this;
		}

		inline bool operator!=(iterator it) const
		{
      // 一般情况下不会比较不同 BlockVector 的迭代子，增加断言以便在
      // 调试模式下发现问题
			return index_ != it.index_;
		}

		inline iterator& operator++()
		{
			index_ = pv_->container_->next_index(index_);
			return *this;
		}

		inline T& operator*() const
		{
			return pv_->operator[](index_);
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

#endif // end of BLOCK_VECTOR_HPP_
