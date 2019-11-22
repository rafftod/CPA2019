#include "Index.h"
#include <fstream>

Index::Index(){
  title = NULL;
  header_offset_table = NULL;
  sequence_offset_table = NULL;
}

Index::~Index(){
  delete[] title;
  delete[] header_offset_table;
  delete[] sequence_offset_table;
}

void Index::read_data(std::ifstream& database_index){
/**
 * Reads data from database index file opened stream.
 * 
 * For every int32, we need to convert from big endian to little endian to get correct values.
 *
 * @param database_index Database file opened ifstream.
 */
    // reads needed data in database_index binary file
    database_index.seekg(2*sizeof(uint32_t)); // skip version and db type
    database_index.read((char*)&title_length,sizeof(uint32_t)); // read int32 in database
    title_length = __bswap_32(title_length); // conversion from big endian to little endian
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
    database_index.read((char*)&residue_count, sizeof(uint64_t)); // already in little endian; no byteswap needed
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

uint32_t Index::get_max_sequence_length() const{
  return max_sequence_length;
}

uint64_t Index::get_residue_count() const{
  return residue_count;
}

char* Index::get_title() const{
  return title;
}

uint32_t Index::get_title_length() const{
  return title_length;
}

uint32_t Index::get_number_of_sequences() const{
  return number_of_sequences;
}

uint32_t* Index::get_header_offset_table() const{
  return header_offset_table;
}

uint32_t* Index::get_sequence_offset_table() const{
  return sequence_offset_table;
}