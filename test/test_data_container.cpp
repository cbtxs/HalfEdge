#include <iostream>

#include "geometry.h"
#include "data_container.h"

using namespace HEM;

void test_data_container_free_index()
{
  DataContainer<8> ddd(10);
  {
    auto p = ddd.add_data<Point>("p");
    std::cout << p->size()  << std::endl;
    for(auto pp : (*p))
    {
      std::cout << pp.x << " " << pp.y << std::endl;
    }
  }
  auto p = ddd.get_data<Point>("p");

}

void test_data_container_copy()
{
  DataContainer<8> ddd(10);
  auto p = ddd.add_data<Point>("p");
  std::cout << p->size()  << std::endl;
  p->get(1) = Point(10, 123);
  for(auto pp : (*p))
    std::cout << pp.x << " " << pp.y << std::endl;

  DataContainer<8> fff(20);
  auto fffp = fff.add_data<Point>("p");

  std::cout << (int)(fffp==p) << std::endl;

  std::cout << "size  :  " << fffp->size() << std::endl;
  fff = ddd;
  std::cout << "size  :  " << fffp->size() << std::endl;
  std::cout << "size  :  " << p->size() << std::endl;

  for(auto pp : (*fffp))
    std::cout << pp.x << " " << pp.y << std::endl;

  auto fffp1 = fff.get_data<Point>("p");
  std::cout << "size  :  " << fffp1->size() << std::endl;

  for(auto pp : (*fffp1))
    std::cout << pp.x << " " << pp.y << std::endl;
}

void test_data_container_delete_data()
{
  DataContainer<8> ddd(10);
  auto p = ddd.add_data<Point>("p");
  std::cout << p->size()  << std::endl;
  p->get(1) = Point(10, 123);
  for(auto pp : (*p))
    std::cout << pp.x << " " << pp.y << std::endl;
  ddd.delete_data<Point>("p");
  p = ddd.get_data<Point>("p");
}

void test_data_container_delete_entity()
{
  DataContainer<8> ddd(10);
  auto p = ddd.add_data<int>("p");
  std::cout << p->size()  << std::endl;
  int i = 0;
  for(auto & pp : (*p))
    pp = i++;

  ddd.delete_index(5);
  for(auto & pp : *p)
    std::cout << pp << std::endl;
  uint32_t idx = ddd.add_index();
  std::cout << idx << std::endl;
  uint32_t a = -1;
  std::cout << (int)(a==(uint32_t)-1) << std::endl;
}

int main()
{
  //test_data_container_free_index();
  //test_data_container_copy();
  //test_data_container_delete_data();
  test_data_container_delete_entity();
  return 0;
}





