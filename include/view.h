#ifndef HEM_VIEW_H
#define HEM_VIEW_H

#include<ranges>

template<typename Imp, typename HalfEdge, typename Entity>
class AdjEntityIteratorBase
{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Entity;
    using difference_type = std::ptrdiff_t;

    AdjEntityIteratorBase(HalfEdge * start_): current(start_), start(start_) {}

    Entity & entity() 
    { 
      auto This = static_cast<Imp*>(this);
      return This->entity_imp(current); 
    }

    Entity & next()
    {
      return static_cast<Imp*>(this)->next_imp(current);
    }

    Entity & operator*() { return entity();}

    Entity & operator->() { return entity();}

    AdjEntityIteratorBase & operator++() 
    {
      //auto This = static_cast<Imp*>(this);
      //current = This->next_half_edge(current);
      //if (current == start || This->check_end_condition(current)) 
      //{
      //  current = nullptr;  // End of iteration
      //}
      current = next();
      return *this;
    }

    bool operator!=(const AdjEntityIteratorBase & other) const 
    {
      return current != other.current;
    }

private:
    HalfEdge* current;
    HalfEdge* start;
};


template<typename Imp, typename Iterator>
class AdjEntityViewBase: public std::ranges::view_base
{
public:


};




















#endif // HEM_VIEW_H
