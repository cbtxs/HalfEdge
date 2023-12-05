#include <vector>
#include <iostream>
#include <memory>
#include <memory_resource>
#include "chunk_array.h"
#include <type_traits>

template<typename T>
class vector_test
{
public:
  vector_test(): alloc(&pool) {}

public:
  std::pmr::unsynchronized_pool_resource pool;
  std::pmr::polymorphic_allocator<T> alloc;
};

using int_vector_test = vector_test<int>;
using double_vector_test = vector_test<double>;

void test_vector(int N)
{
  constexpr size_t BUF_SIZE = 1000000000;
  std::pmr::pool_options options;
  options.max_blocks_per_chunk = 40;
  options.largest_required_pool_block = 640;

  char* buffer = new char[BUF_SIZE];
  std::pmr::monotonic_buffer_resource pool{buffer, BUF_SIZE};
  //std::pmr::unsynchronized_pool_resource mem (options,&pool);
  //std::pmr::unsynchronized_pool_resource pool;
  std::pmr::vector<int> v0{&pool};
  std::vector<int> v1;
  std::vector<int> v2;

  auto s1 = std::chrono::high_resolution_clock::now();
  v0.emplace_back(2);
  std::cout << (long long)&v0[0] << std::endl;
  for(int i=0; i<N; i++)
    v1.emplace_back(i);
  std::cout << (long long)&v0[0] << std::endl;
  auto s2 = std::chrono::high_resolution_clock::now();
  auto d0 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
  std::cout << "0 : " << d0.count()/1e6 << std::endl;

  s1 = std::chrono::high_resolution_clock::now();
  for(int i=0; i<N; i++)
    v0.emplace_back(i);
  s2 = std::chrono::high_resolution_clock::now();
  d0 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
  std::cout << "1 : " << d0.count()/1e6 << std::endl;

  s1 = std::chrono::high_resolution_clock::now();
  v2.reserve(N);
  for(int i=0; i<N; i++)
    v2[i] = i;
  s2 = std::chrono::high_resolution_clock::now();
  d0 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
  std::cout << "1 : " << d0.count()/1e6 << std::endl;

  /**
  ChunkArray<int, 2048> vv("vv");
  s1 = std::chrono::high_resolution_clock::now();
  for(int i=0; i<N; i++)
    vv.push_back(i);
  s2 = std::chrono::high_resolution_clock::now();
  d0 = std::chrono::duration_cast<std::chrono::microseconds>(s2-s1);
  std::cout << "1 : " << d0.count()/1e6 << std::endl;

  ChunkArray<int> vint("vint");
  for(int i = 0; i < 10000; i++)
  {
    vint.emplace_back(i);
  }
  std::cout <<  vint.capacity()/1024 << std::endl;
  vint.clear();
  for(int i = 0; i < 10000; i++)
  {
    vint.emplace_back(i);
  }
  std::cout <<  vint.capacity()/1024 << std::endl;
  for (int i = 0; i < 10000; i++) 
  {
    if(vint[i] != i)
      std::cout << (int)(vint[i]==i) << std::endl;
  }
  */

}

void test_bety_operator()
{
  int * s = new int[10];
  s[1] = 1;
  *s |= 1 << 31;
  std::cout << s[0] << std::endl;
  delete [] s;
}

class ABase
{
public:
  ABase(int aa):a(aa) {}
  int & get_val() {return a;}

  virtual int add(int a) = 0;
private:
  int a;
};

template<typename T>
class B: public ABase
{
public:
  B(int a_): ABase(a_) {}

  int add(int b) override
  {
    return this->get_val()+b;
  }
private:
  std::vector<T> s;
};

using C = ABase;

template<typename Entity>
auto entity()
{
  return std::conditional_t<std::is_same_v<Entity, int>,  ,
         std::conditional_t<std::is_same_v<Entity, double>, double,
         std::conditional_t<std::is_same_v<Entity, uint32_t>, uint32_t, int>>>
}

void test_inherit()
{
  B<int> b(0);
  ABase * c = &b;
  auto ss = c->add(3);
  std::cout << ss << std::endl;

  using TE = std::tuple<C, B<int>, ABase>;

  std::cout << entity<double>() << std::endl;
}

int main(int , char**argv)
{
  //int N = std::stoi(argv[1]);
  //test_vector(N);
  test_inherit();
}
