#include "tools.h"

int main()
{
  HEM::IrregularArray2D<int> a;
  auto & data = a.get_data();
  auto & start = a.get_start_pos();

  data.resize(6);
  start.resize(3);
  start = {0, 2, 6};

  data = {1, 2, 3, 4, 5, 6};
  auto row0 = a[1];
  for (auto & x : row0)
    std::cout << x << " "<< std::endl;
  return 0;
}
