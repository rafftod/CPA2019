#ifndef HEADER
#define HEADER
#include <fstream>

class Header{
private:
  char* title = NULL;
  uint8_t length;
  int msb_find(uint8_t in_byte);
  int length_find(std::ifstream& databatse_header, uint8_t nb_bytes);


public:
  Header();
  ~Header();
  
  void read_data(std::ifstream& database_header,int offset);
  char* get_title() const;
  uint8_t get_length() const;
};

#endif
