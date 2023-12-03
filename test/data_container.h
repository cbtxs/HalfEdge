#ifndef _DATA_CONTAINER_
#define _DATA_CONTAINER_

#include <vector>
#include <memory>
#include <map>

#include "chunk_array.h"
#include "chunk_array_bool.h"

template<typename T>
class DataContainer
{
public:
  DataContainer(): data_() {}
private:
  std::map<std::string, std::shared_ptr<ChunkArrayBase>> data_;
  ChunkArrayBool<1024u> is_used_;
  std::vector<unsigned int> free_index;
};

#endif /* _DATA_CONTAINER_ */ 
