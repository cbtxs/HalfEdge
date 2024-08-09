#include <iostream>
#include <chrono>

#include "chunk_array.h"

using namespace std::chrono;
using namespace HEM;


void test_chunkarray_iterator()
{
  ChunkArray<int> l("l");
  l.resize(1000);
  l.set_value(100);
  for(auto & val : l)
    std::cout << "val : " << val << std::endl;
}

void test_chunkarray_operaotr()
{
  ChunkArray<int, 8> l("l");
  l.resize(17);
  for(uint32_t i = 0; i < l.size(); i++)
    l[i] = i;
  for(auto & val : l)
    std::cout << "val : " << val << std::endl;
  long start = (long)&(l[0]);
  for(auto & val : l)
    std::cout << "val : " << (long)(&val)-start << " " << ((long)(&val)-start)/4 << std::endl;
  std::cout << "size of l : " << l.size() << std::endl;
  std::cout << "capacity of l : " << l.capacity() << std::endl;

  ChunkArray<int, 8> ll("l");
  ll = l;
  for(auto & val : ll)
    std::cout << "val : " << val << " " << (long)(&val)-start << " " << ((long)(&val)-start)/4 << std::endl;
  std::cout << "size of l : " << ll.size() << std::endl;
  std::cout << "capacity of l : " << ll.capacity() << std::endl;
  l.resize(25);
  for(auto & val : l)
    std::cout << "val : " << val << " " << (long)(&val)-start << " " << ((long)(&val)-start)/4 << std::endl;
}

int main()
{
  test_chunkarray_operaotr();
  return 0;
}

