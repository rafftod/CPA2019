#include "Index.h"
#include <iostream>
#include <fstream>

using namespace std;

Index::Index(){

}

Index::~Index(){

}

void Index::read_data(ifstream database_index){
  // reads needed data in database_index binary file
  database_index.seekg(2*sizeof(uint32_t)); // skip version and db type
  database_index.read((char*)&title_length,sizeof(uint32_t)); // read int32 in db
  title_length = __bswap_32(title_length); // conversion from little endian to big endian
  title = new char[title_length];
  database_index.read((char*)title,sizeof(char)*title_length); // read char[title_length] in db
  uint32_t timestamp_length; // we need to read it to seekg far enough in file
  database_index.read((char*)&timestamp_length,sizeof(uint32_t));
  timestamp_length = __bswap_32(timestamp_length);
  // tellg() provides current position in database
  // as seekg needs an absolute position, we sum tellg() and
  // the part we need to skip
  database_index.seekg((int)database_index.tellg()+sizeof(char)*timestamp_length);
  database_index.read((char*)&number_of_sequences, sizeof(uint32_t));
  number_of_sequences = __bswap_32(number_of_sequences);
  database_index.read((char*)&residue_count, sizeof(uint64_t));
  database_index.read((char*)&max_sequence_length, sizeof(uint32_t));
  max_sequence_length = __bswap_32(max_sequence_length);
  header_offset_table = new uint32_t[number_of_sequences+1];
  sequence_offset_table = new uint32_t[number_of_sequences+1];
  database_index.read((char*)header_offset_table, sizeof(uint32_t)*(number_of_sequences+1));
  database_index.read((char*)sequence_offset_table, sizeof(uint32_t)*(number_of_sequences+1));
  for(int i = 0; i < number_of_sequences + 1; i++){
    header_offset_table[i] = __bswap_32(header_offset_table[i]);
    sequence_offset_table[i] = __bswap_32(sequence_offset_table[i]);
  }
}
