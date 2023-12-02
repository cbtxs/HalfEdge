#include "data_set_container.hpp"

namespace OF {


DataSetBase_::DataSetBase_(
    DataSetContainerBase* container, 
    const std::string& name): container_(container), name_(name)
{
}

DataSetBase_::~DataSetBase_()
{
  if (container_)
    container_->delete_data_set(this); // 数据所在容器负责释放数据，
}

uint32_t DataSetBase_::end_index() const
{
  if (container_)
    return container_->end_index(); // 数据所在容器负责管理索引指标
  return 0;
}

uint32_t DataSetBase_::number_of_datas() const
{
  if (container_)
    return container_->number_of_datas(); // 数据中数据的个数
  return 0;
}

uint32_t DataSetBase_::size() const
{
  if (container_)
    return container_->number_of_datas(); // 数据中数据的个数
  return 0;
}

DataSetContainerBase::DataSetContainerBase() : 
  end_index_(0), nb_datas_(0)
{
  data_sets_.reserve(32);
  sp_data_sets_.reserve(32);
  available_indices_.reserve(1024);
}

DataSetContainerBase::~DataSetContainerBase()
{
  // data_sets_ 和 sp_data_sets_ 会自动销毁
  //
  // sp_data_sets_ 中存储的是共享指针，
  // 当最后一个共享指针被销毁时，会自动调用数据本身的析构函数，实现分配内存
  // 的回收
}

/*
 * @note 增加一个新的指标
 */
uint32_t DataSetContainerBase::new_index()
{
  std::uint32_t index;
  if (std::uint32_t(available_indices_.size()) > 0)
  {
    index = available_indices_.back();
    available_indices_.pop_back(); // 记录可用的指标，如果有可用指标，取最后一个可用指标
  }
  else
  {
    index = end_index_++; // 如果没有可用的指标，就取指标上界，然后指标上界自加
  }

  // 每个数据管理这个指标，根据指标调整每一个数据容器的大小 
  for (DataSetBase_* ds : data_sets_)
  {
    ds->manage_index(index); // 根据指标调整每一个数据容器的实际大小
  }

  init_ref_counter(index); // 初始这个指标的引用次数为 1
  ++nb_datas_; // 数据中的总数增加 1

  return index;
}

void DataSetContainerBase::release_index(uint32_t index)
{
  // assert(nb_refs(index) > 0)
  // Trying to release an unused index
  available_indices_.push_back(index); // 放入可用指标
  reset_ref_counter(index);
  --nb_datas_;
}

void DataSetContainerBase::remove_data_set(
    const std::shared_ptr<DataSetBase_>& ds)
{
  auto it = std::find(sp_data_sets_.begin(), sp_data_sets_.end(), ds);
  if (it != sp_data_sets_.end()) 
  {
    *it = sp_data_sets_.back(); // 把最后一个交换到这个位置
    sp_data_sets_.pop_back(); // 并把最后一个数据移除
  }
}

void DataSetContainerBase::remove_data_set(DataSetBase_* ds)
{
  auto it = std::find_if(
      sp_data_sets_.begin(), 
      sp_data_sets_.end(),
      [&](const auto& md) {return md.get() == ds;}
      );
  if (it != sp_data_sets_.end())
  {
    *it = sp_data_sets_.back(); // 把最后一个交换到这个位置
    sp_data_sets_.pop_back(); // 并把最后一个数据移除
  }
}

void DataSetContainerBase::clear_data_sets()
{
  for (DataSetBase_* md : data_sets_)
    md->clear();

  available_indices_.clear();
  nb_datas_ = 0;
  end_index_ = 0;
}

void DataSetContainerBase::remove_data_sets()
{
  sp_data_sets_.clear(); // 共享指针清除，所有数据就都被移除了
}


void DataSetContainerBase::delete_data_set(DataSetBase_* ds)
{
  auto it = std::find(data_sets_.begin(), data_sets_.end(), ds);
  if (it != data_sets_.end())
  { 
    *it = data_sets_.back(); // 把最后一个交换到这个位置
    data_sets_.pop_back(); //并把最后一个数据移除
  }
}

} // end of namespace OF
