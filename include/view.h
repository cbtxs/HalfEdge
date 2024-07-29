#ifndef HEM_VIEW_H
#define HEM_VIEW_H

#include<ranges>

/**
 * @brief 邻接实体迭代器基类，用于遍历邻接实体。
 * 其子类需要实现以下两个方法： 
 * 1. entity_imp: 返回当前半边对应的实体
 * 2. next_imp : 返回下一个半边, 如果返回nullptr则表示遍历结束
 */
template<typename Imp, typename HalfEdge, typename Entity>
class AdjEntityIteratorBase
{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Entity;
    using difference_type = std::ptrdiff_t;

    AdjEntityIteratorBase(HalfEdge * start__): current_(start__), start_(start__) {}

    Entity & entity() 
    { 
      return static_cast<Imp*>(this)->entity_imp(current_); 
    }

    HalfEdge * next()
    {
      return static_cast<Imp*>(this)->next_imp(current_);
    }

    Entity & operator*() { return entity();}

    Entity & operator->() { return entity();}

    Imp & operator++() 
    {
      //auto This = static_cast<Imp*>(this);
      //current_ = This->next_half_edge(current_);
      //if (current_ == start_ || This->check_end_condition(current_)) 
      //{
      //  current_ = nullptr;  // End of iteration
      //}
      current_ = next();
      return *this;
    }

    bool operator!=(const AdjEntityIteratorBase & other) const 
    {
      return current_ != other.current_;
    }

private:
    HalfEdge* current_;
    HalfEdge* start_;
};


/**
 * @brief 邻接实体视图基类，用于遍历邻接实体。
 */
template<typename Iterator>
class AdjEntityViewBase : public std::ranges::view_base
{
public:
  using HalfEdge = typename Iterator::HalfEdge;

public:
  AdjEntityViewBase(HalfEdge * start) : start_(start) {}

  Iterator begin() const { return Iterator(start_); }

  Iterator end() const { return Iterator(nullptr); }

private:
    HalfEdge * start_;
};



















#endif // HEM_VIEW_H
