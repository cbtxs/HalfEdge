#ifndef _CHUNK_ARRAY_
#define _CHUNK_ARRAY_

#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>

class ChunkArrayBase {};

/**
 * @brief 分块存储的数组
 * @param ChunkSize : 每个块的元素个数
 * @param size_     : 当前数组的长度
 * @param chunks_   : 每个块的指针
 */
template <typename T, size_t ChunkSize = 1024u>
class ChunkArray : public ChunkArrayBase 
{
public:
  // 构造函数
  ChunkArray() : size_(0), chunks_(0) 
  { 
    chunks_.reserve(1024); 
  }

  // 析构函数，释放分配的内存
  ~ChunkArray() 
  {
    for (T* chunk : chunks_) 
      delete[] chunk;
  }

  // 获取元素
  T& operator[](size_t index) 
  {
    assert(index < size_ && "Index out of range");
    size_t chunkIndex = index / ChunkSize;
    size_t offset = index % ChunkSize;
    return chunks_[chunkIndex][offset];
  }

  // 添加元素
  void push_back(const T& value) 
  {
    if (size_ == chunks_.size() * ChunkSize) 
      chunks_.push_back(new T[ChunkSize]);

    size_t chunkIndex = size_ / ChunkSize;
    size_t offset = size_ % ChunkSize;

    chunks_[chunkIndex][offset] = value;
    size_++;
  }

  // emplace_back 类似于 std::vector
  template <typename... Args>
  void emplace_back(Args&&... args)
  {
    if (size_ == capacity())
      chunks_.push_back(new T[ChunkSize]);

    size_t chunkIndex = size_ / ChunkSize;
    size_t offset = size_ % ChunkSize;

    chunks_[chunkIndex][offset] = T(std::forward<Args>(args)...);
    size_++;
  }

  // 获取当前元素个数
  size_t size() const { return size_;}

  // 获取当前分配的内存容量
  size_t capacity() const { return chunks_.size() * ChunkSize;}

  // 交换两个 ChunkedVector 的内容
  void swap(ChunkArray& other) 
  {
    std::swap(chunks_, other.chunks_);
    std::swap(size_, other.size_);
  }

  void clear() 
  {
    size_ = 0;
  }

  // 复制另一个 ChunkArray 的内容
  void copy(const ChunkArray& other) 
  {
    clear();
    reserve(other.size_);
    int N_chunk = other.chunks_.size();
    for (size_t i = 0; i < N_chunk; ++i) 
      std::copy(other.chunks_[i], other.chunks_[i]+ChunkSize, chunks_[i]);
  }

  // 预留空间，使得至少可以容纳指定数量的元素
  void reserve(size_t newCapacity) 
  {
    if (newCapacity > capacity()) 
    {
      size_t requiredChunks = (newCapacity + ChunkSize - 1) / ChunkSize;
      chunks_.resize(requiredChunks, nullptr);
      for (size_t i = capacity()/ChunkSize; i < requiredChunks; ++i) 
        chunks_[i] = new T[ChunkSize];
    }
  }

  void resize(size_t newSize) 
  {
    if (newSize <= capacity()) 
    {
      size_ = newSize;  // 如果新大小小于等于当前大小，直接更新 size_
    } 
    else 
    {
      reserve(newSize);  // 如果新大小大于当前大小，调用 reserve 函数
      size_ = newSize;
    }
  }


public:
  // 迭代子
  class iterator {
  public:
    iterator(ChunkArray& array, size_t index)
        : array_(array), index_(index) {}

    bool operator!=(const iterator& other) const { return index_ != other.index_; }

    iterator& operator++() 
    {
      index_++;
      return *this;
    }

    T& operator*() { return array_[index_];}

  private:
    ChunkArray& array_;
    size_t index_;
  };

  // 迭代子的起始位置
  iterator begin() { return iterator(*this, 0); }

  // 迭代子的结束位置
  iterator end() { return iterator(*this, size_);}

private:
    size_t size_;  // 元素个数
    std::vector<T*> chunks_;  // 存储块的指针
};

#endif /* _CHUNK_ARRAY_ */ 
