#ifndef _DATA_SET_
#define _DATA_SET_

#include <vector>
#include <unordered_map>
#include <memory_resource>

namespace HEM
{

template<typename T>
class data_set
{
public:
  data_set(std::pmr::synchronized_pool_resource & spr): alloc(&spr) {}

  /**
   * @brief 添加 1 个元素
   */
  template <typename... Args>
  T * add(Args&&... args)
  {
    T * data = alloc.allocate(1);
    new (data) T(std::forward<Args>(args)...);
    insert(data);
    return data;
  }

  /**
   * @brief 在尾部插入 1 个元素 data
   */
  void insert(T * data)
  {
    ptr2idx.emplace(data, size());
    ptrs.push_back(data);
  }

  /**
   * @brief 删除一个元素 data
   */
  void remove(T * data)
  {
    auto it = ptr2idx.find(data); /**< 找到元素的 idx */
    assert(it != ptr2idx.end());
    if(it->second != size()-1) /**< 若不是最后一个，就用最后一个替换 */
    {
      ptr2idx[ptrs.back()] = it->second;
      ptrs[it->second] = ptrs.back();
    }
    alloc.destroy(data); /**< 析构对象 */
    alloc.deallocate(data, 1);
    ptrs.pop_back();
    ptr2idx.erase(it);
  }

  /**
   * @brief 删除最后一个元素 
   */
  void pop_back()
  {
    T * data = ptrs.back();
    alloc.deallocate(data, 1);
    ptrs.pop_back();
    ptr2idx.erase(data);
  }

  /**
   * @brief 预开辟内存
   */
  void reserve(std::size_t n) 
  {
    ptrs.reserve(n);
    ptr2idx.reserve(n);
  }

  /**
   * @brief 清空内存
   */
  void clear() noexcept
  {
    ptrs.clear();
    ptr2idx.clear();
  }

  std::size_t size() const { return ptrs.size(); }

  T * operator[](std::size_t index) const noexcept { return ptrs[index]; }

  T * at(std::size_t index) { return ptrs[index]; }

  auto begin() noexcept { return ptrs.begin(); }

  auto begin() const noexcept { return ptrs.begin(); }

  auto end() noexcept { return ptrs.end(); }

  auto end() const noexcept { return ptrs.end(); }

  std::size_t idx(T* e) const { return ptr2idx.at(e); }

  bool contains(T* e) const { return ptr2idx.find(e) != ptr2idx.end(); }

  bool empty() const noexcept { return ptrs.empty(); }

private:
  std::vector<T*> ptrs;
  std::unordered_map<T*, int> ptr2idx;
  std::pmr::polymorphic_allocator<T> alloc;
};

}

#endif // end of _DATA_SET_
