#ifndef INDEX
#define INDEX

#include <fstream>

class Index{
private:
  uint32_t title_length;
  char* title;
  uint32_t number_of_sequences;
  uint64_t residue_count; // STORED IN LITTLE ENDIAN
  uint32_t max_sequence_length;
  uint32_t* header_offset_table;
  uint32_t* sequence_offset_table;
public:
  // still need to make constructor and destructor
  Index();
  ~Index();

  void read_data(std::ifstream database_index);
};

#endif
