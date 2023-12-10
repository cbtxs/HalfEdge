#ifndef _MARK_ARRAY_
#define _MARK_ARRAY_

#include <cassert>
#include <vector>
#include <stdexcept>
#include <memory>

namespace HEM
{
/**
 * @brief bool 类型的 chunk array, 使用 uint32_t 表示 32 位的 bool 值. 
 */
template <size_t ChunkSize = 1024u>
class ChunkArrayBool
{
public:
  // 构造函数
  ChunkArrayBool(): size_(0), chunks_(0) {}

  // 构造函数
  ChunkArrayBool(int size) : size_(0), chunks_(0) 
  {
    resize(size);
  }

  // 构造函数
  ChunkArrayBool(const ChunkArrayBool & other)
  {
    this->copy(other);
  }

  // 析构函数，释放分配的内存
  ~ChunkArrayBool() 
  {
    for (uint32_t* chunk : chunks_)
      delete[] chunk;
  }

  void release()
  {
    for (uint32_t* chunk : chunks_)
      delete[] chunk;
  }

  void set_true()
  {
    for(auto & chunk : chunks_)
      std::fill(chunk, chunk+(ChunkSize/32), 4294967295);
  }

  void set_false()
  {
    for(auto & chunk : chunks_)
      std::fill(chunk, chunk+(ChunkSize/32), 0);
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
  void set_false(size_t index) 
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
      chunks_.emplace_back(new  uint32_t[ChunkSize/32]);
    set_true(size_);
    size_++;
  }

  void push_back_false()
  {
    if (size_ == capacity())
      chunks_.emplace_back(new  uint32_t[ChunkSize/32]);
    set_false(size_);
    size_++;
  }

  // 交换两个 ChunkArray 的内容
  void swap(ChunkArrayBool& other) 
  {
    std::swap(chunks_, other.chunks_);
    std::swap(size_, other.size_);
  }

  // 获取当前元素个数
  size_t size() const { return size_;}

  // 获取当前分配的内存容量
  size_t capacity() const { return chunks_.size() * ChunkSize;}

  void clear() { size_ = 0;}

  // 复制另一个 ChunkArrayBool 的内容
  void copy(const ChunkArrayBool& other) 
  {
    if(this != &other)
    {
      resize(other.size_);
      uint32_t N_chunk = other.chunks_.size();
      for (size_t i = 0; i < N_chunk; ++i) 
        std::copy(other.chunks_[i], other.chunks_[i]+ChunkSize/32, chunks_[i]);
    }
  }

  // 获取元素
  bool operator[](uint32_t index) const { return get(index); }

  /** @brief operator =  */
  ChunkArrayBool & operator = (const ChunkArrayBool & other)
  {
    this->copy(other);
    return *this;
  }

  // 预留空间，使得至少可以容纳指定数量的元素
  void reserve(size_t newCapacity) 
  {
    uint32_t cap = capacity();
    if (newCapacity > capacity()) 
    {
      size_t requiredChunks = (newCapacity + ChunkSize - 1) / ChunkSize;
      chunks_.resize(requiredChunks, nullptr);
      for (size_t i = cap/ChunkSize; i < requiredChunks; ++i) 
        chunks_[i] = new uint32_t[ChunkSize/32]();
    }
  }

  void resize(size_t newSize)
  {
    reserve(newSize);  // 如果新大小大于当前大小，调用 reserve 函数
    size_ = newSize;
  }

private:
    size_t size_;  // 元素个数
    std::vector<uint32_t*> chunks_;  // 存储块的指针
};

}

#endif /* _MARK_ARRAY_ */ 
