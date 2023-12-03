#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>

#include "chunk_array.h"

template <size_t ChunkSize = 1024u>
class ChunkArrayBool : public ChunkArrayBase 
{
public:
  // 构造函数
  ChunkArrayBool() : size_(0), chunks_(0) {}

  // 析构函数，释放分配的内存
  ~ChunkArrayBool() {
      for (unsigned int* chunk : chunks_)
          delete[] chunk;
  }

  // 设置指定位置的 bool 值
  void set_true(size_t index) 
  {
    assert(index < size_ && "Index out of range");
    size_t chunkIndex = index / ChunkSize;  // 计算是第几个 chunk
    size_t offset = index % ChunkSize;  // 计算在这个 chunk 中的偏移位置
    size_t intIndex = offset / 32;
    offset %= 32;
    chunks_[chunkIndex][intIndex] |= (1u << offset);  // 设置为 1
  }

  // 设置指定位置的 bool 值
  void set_fasle(size_t index) 
  {
    assert(index < size_ && "Index out of range");
    size_t chunkIndex = index / ChunkSize;  // 计算是第几个 chunk
    size_t offset = index % ChunkSize;  // 计算在这个 chunk 中的偏移位置
    size_t intIndex = offset / 32;
    offset %= 32;
    chunks_[chunkIndex][intIndex] &= ~(1u << offset);  // 设置为 0
  }

  // 获取指定位置的 bool 值
  bool get(size_t index) const 
  {
    assert(index < size_ && "Index out of range");
    size_t chunkIndex = index / ChunkSize;  // 计算是第几个 chunk
    size_t offset = index % ChunkSize;  // 计算在这个 chunk 中的偏移位置
    size_t intIndex = offset / 32;
    offset %= 32;
    return (chunks_[chunkIndex][intIndex] & (1u << offset)) != 0;  // 检查指定位置是否为 1
  }

  void push_back() { push_back_true(); }

  void push_back_true()
  {
    if (size_ == capacity())
      chunks_.emplace_back(new unsigned int[ChunkSize/32]);
    set_true(size_);
    size_++;
  }

  void push_back_false()
  {
    if (size_ == capacity())
      chunks_.emplace_back(new unsigned int[ChunkSize/32]);
    set_false(size_);
    size_++;
  }

  // 交换两个 ChunkArray 的内容
  void swap(ChunkArrayBool& other) 
  {
    std::swap(chunks_, other.chunks_);
    std::swap(size_, other.size_);
  }

  // 获取元素
  bool operator[](size_t index) { return get(index); }

  // 获取当前元素个数
  size_t size() const { return size_;}

  // 获取当前分配的内存容量
  size_t capacity() const { return chunks_.size() * ChunkSize;}

  void clear() 
  {
    size_ = 0;
  }

  // 复制另一个 ChunkArray 的内容
  void copy(const ChunkArrayBool& other) 
  {
    clear();
    reserve(other.size_);
    int N_chunk = other.chunks_.size();
    for (size_t i = 0; i < N_chunk; ++i) 
      std::copy(other.chunks_[i], other.chunks_[i]+ChunkSize/32, chunks_[i]);
  }

  // 预留空间，使得至少可以容纳指定数量的元素
  void reserve(size_t newCapacity) 
  {
    if (newCapacity > capacity()) 
    {
      size_t requiredChunks = (newCapacity + ChunkSize - 1) / ChunkSize;
      chunks_.resize(requiredChunks, nullptr);
      for (size_t i = capacity()/ChunkSize; i < requiredChunks; ++i) 
        chunks_[i] = new unsigned int[ChunkSize/32];
    }
  }

  void resize(size_t newSize) 
  {
    reserve(newSize);  // 如果新大小大于当前大小，调用 reserve 函数
    size_ = newSize;
  }

private:
    size_t size_;  // 元素个数
    std::vector<unsigned int*> chunks_;  // 存储块的指针
};

