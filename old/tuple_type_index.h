#ifndef _TUPLE_TYPE_INDEX_
#define _TUPLE_TYPE_INDEX_

#include <cstddef>
#include <tuple>

template <typename T, typename Tuple>
struct tuple_type_index;

// 基本情况：如果元组的第一个元素的类型匹配，返回 0
template <typename T, typename... Ts>
struct tuple_type_index<T, std::tuple<T, Ts...>> : std::integral_constant<std::size_t, 0> {};

// 递归情况：如果元组的第一个元素的类型不匹配，继续在剩余的元组中查找
template <typename T, typename U, typename... Ts>
struct tuple_type_index<T, std::tuple<U, Ts...>> : std::integral_constant<std::size_t, 1 + tuple_type_index<T, std::tuple<Ts...>>::value> {};


#endif /* _TUPLE_TYPE_INDEX_ */ 
