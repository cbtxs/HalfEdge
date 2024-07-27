#ifndef IRREGULAR_ARRAY2D_H
#define IRREGULAR_ARRAY2D_H

#include <vector>
#include <iostream>

namespace HEM
{

template <typename T>
class IrregularArray2D 
{
public:
 IrregularArray2D() = default;

 uint32_t len() const 
 {
   return start_pos.size() - 1;
 }

 uint32_t row_size(uint32_t row) const
 {
   return start_pos[row + 1] - start_pos[row];
 }

 T* get_row(uint32_t row) 
 {
   auto loc = start_pos[row];
   return data.data() + loc;
 }

 const T* get_row(uint32_t row) const 
 {
   auto loc = start_pos[row];
   return data.data() + loc;
 }

 T * operator[](uint32_t rowIndex)
 {
   return get_row(rowIndex);
 }

 const T * operator[](uint32_t rowIndex) const
 {
     return get_row(rowIndex);
 }

 std::vector<T>& get_data()
 {
     return data;
 }

 std::vector<uint32_t>& get_start_pos()
 {
     return start_pos;
 }

private:
   std::vector<T> data;                 
   std::vector<uint32_t> start_pos;   
};

} // namespace HEM

#endif // IRREGULAR_ARRAY2D_H
