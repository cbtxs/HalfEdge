#ifndef DATA_SET_CONTAINER_HPP_
#define DATA_SET_CONTAINER_HPP_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>


namespace OF {

class DataSetContainerBase;
template <template <typename> class DataSet>
class DataSetContainerT;

/**
 * @class DataSetBase
 * @brief 数据集的基类，存储数据集所在容器的指针，规定了子类的接口
 *  和功能，如指标的管理，以及在容器中如何创建自己的
 */
class DataSetBase_
{
public:
  DataSetBase_(DataSetContainerBase * container, const std::string& name);
  virtual ~DataSetBase_();

  inline const std::string& name() const
  {
    return name_;
  }

  /**
   * @
   * @brief 返回数据中用到的索引上界指标
   * @note 实际由所在容器来管理
   */
  uint32_t end_index() const;

  /**
   * @brief 返回数据集中实际存储的数据个数
   * @note 实际由所在容器来管理
   */
  uint32_t number_of_datas() const;

  /**
   * @brief 返回数据集中实际存储的数据个数
   */
  uint32_t size() const;

protected:
    DataSetContainerBase* container_; // 数据所在的容器
    std::string name_; // 数据名字
private:
    friend DataSetContainerBase;
    template <template <typename> class DataSet>
    friend class DataSetContainerT;

    /**
     * @brief 根据指标做内存的管理
     */
    virtual void manage_index(uint32_t index) = 0;
    virtual void clear() = 0;
    virtual std::shared_ptr<DataSetBase_> create_in(DataSetContainerBase& container) const = 0;
    virtual void copy(const DataSetBase_& src) = 0;
};

/**
 * @class DataSetContainerBase
 * @brief 指标管理 
 */

class DataSetContainerBase
{
public:
  using const_iterator = std::vector<std::shared_ptr<DataSetBase_>>::const_iterator;
public:

  DataSetContainerBase();

  /**
   * @brief 初始化一个有 nd 数据的容器
   */
  DataSetContainerBase(std::uint32_t nd);

  /**
   * @brief 
   */
  void init(std::uint32_t nd);
  virtual ~DataSetContainerBase();

  inline std::uint32_t number_of_datas() const
  {
    return nb_datas_;
  }

  /**
   * @brief 获取一个可用的索引指标
   * @note 
   */
  std::uint32_t new_index();
  void release_index(std::uint32_t index);

  void remove_data_set(const std::shared_ptr<DataSetBase_> & ds);
  void remove_data_set(DataSetBase_* ds);

  /**
   * @brief 清理容器中所有数据集的数据
   */
  void clear_data_sets();

  /**
   * @brief 移除数据集
   * @note 如果存在共享指针，一些数据集可能不会被移除
   */
  void remove_data_sets();

  inline const_iterator begin() const
  {
    return sp_data_sets_.begin();
  }

  inline const_iterator end() const
  {
    return sp_data_sets_.end();
  }

  inline uint32_t begin_index() const
  {
    uint32_t index = 0u;
    while (index < end_index_ && number_of_refs(index) == 0)
      ++index;
    return index;
  }

  inline uint32_t end_index() const
  {
    return end_index_;
  }

  inline uint32_t next_index(uint32_t index) const
  {
    do
    {
      ++index;
    } while (index < end_index_ && number_of_refs(index) == 0);
    return index;
  }

protected:

  std::vector<DataSetBase_*> data_sets_; // 这里为什么不用 Map 呢？
  /**
   * @note 这里用于数据内存资源的自动释放
   */
  std::vector<std::shared_ptr<DataSetBase_>> sp_data_sets_; // 同样问上面的问题

  std::vector<uint32_t> available_indices_;
  uint32_t end_index_; // 容器中每个数据集所用的索引指标上界
  uint32_t nb_datas_; // 容器中每个数据集中存储的数据个数

  friend DataSetBase_;

  /**
   * @brief 从 data_sets_ 中删除数据集的指针
   * @note 只有当数据集析构时，才会调用这个函数。当最后一个共享指针被销毁时，
   *  就会调用这个函数。而且这个函数是被保护的，用户不能调用
   *  
   */
  void delete_data_set(DataSetBase_* ds);

  virtual void init_ref_counter(uint32_t index) = 0;
  virtual void reset_ref_counter(uint32_t index) = 0;
  virtual uint32_t number_of_refs(uint32_t index) const = 0;
};


/**
 * @class DataSetContainerT
 * @brief 主要用来管理指标引用计数 
 * @note
 */
template <template <typename> class DATASET>
class DataSetContainerT : public DataSetContainerBase
{
public:
  template <typename T>
  using DataSet = DATASET<T>;
  using DataSetBase = DataSetBase_;

protected:
  std::unique_ptr<DataSet<uint32_t>> ref_counter_; // 存储每个数据位置的引用计数

  /**
   * @brief 初始化 index 被引用的次数
   */
  inline void init_ref_counter(uint32_t index) override
  {
    // 管理内存
    static_cast<DataSetBase*>(ref_counter_.get())->manage_index(index); 
    // 初始化
    (*ref_counter_)[index] = 1u;
  }

  /**
   * @brief 重置一个指标的引用次数为 0 
   */
  inline void reset_ref_counter(uint32_t index) override
  {
    // assert 
    (*ref_counter_)[index] = 0u;
  }

  /**
   * @brief 当前指标被引用的次数
   */
  inline uint32_t number_of_refs(uint32_t index) const override
  {
    //assert
    return (*ref_counter_)[index];
  }


public:

  /**
   * @brief 构造函数
   */
  DataSetContainerT() : DataSetContainerBase()
  {
    // 这个数据唯一存在，且不属于任何容器
    ref_counter_ = std::make_unique<DataSet<uint32_t>>(nullptr, "__refs");
  }

  /**
   * @brief 析构函数
   */
  ~DataSetContainerT()
  {
  }

  /**
   * @brief 增加一个新的数据
   * @note 
   */
  template <typename T>
  std::shared_ptr<DataSet<T>> add_data_set(const std::string& name)
  {
    // 查找是否已经存在数据
    auto it = std::find_if(
       data_sets_.begin(), 
       data_sets_.end(),
       [&](DataSetBase* ds){
         return ds->name().compare(name) == 0;
       });
    if (it == data_sets_.end())
    { // 如果不存在
      std::shared_ptr<DataSet<T>> spds = std::make_shared<DataSet<T>>(this, name);
      DataSet<T>* pds = spds.get();
      static_cast<DataSetBase*>(pds)->manage_index(end_index_);//初始化内存
      data_sets_.push_back(pds);
      sp_data_sets_.push_back(spds);
      return  spds;
    }
    return std::shared_ptr<DataSet<T>>(); //为什么要返回这个数据？
 }

 /**
  * @
  * @brief 获取一个数据的共享指针
  * @note 
  */
 template <typename T>
 std::shared_ptr<DataSet<T>> get_data_set(const std::string& name) const
 {
   auto it = std::find_if(
       sp_data_sets_.begin(), 
       sp_data_sets_.end(),
       [&](const auto& ds){
         return ds->name().compare(name) == 0;
       });
   if ( it != sp_data_sets_.end())
     return std::dynamic_pointer_cast<DataSet<T>>(*it);

   return std::shared_ptr<DataSet<T>>();//如果没有这个数据，就返回空指针
 }

 /**
  * @
  * @brief 索引指标引用增加 1
  * @note
  */
 inline void ref_index(uint32_t index)
 {
   (*ref_counter_)[index]++;
 }

 /**
  * @
  * @brief 指标索引自减 1
  * @note
  */
 inline bool unref_index(uint32_t index)
 {
   (*ref_counter_)[index]--;
   if ((*ref_counter_)[index] == 1u) // 这里是一个 Bug 吗？应该为 0 吧？
   {
     release_index(index);
     return true;
   }
   return false;
 }

};

} // end of namespace OF

#endif // end of DATA_SET_CONTAINER_HPP_ 
