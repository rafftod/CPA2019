#include "Header.h"
#include <iostream>
#include <fstream>

Header::Header(){ title = NULL; } 

Header::~Header(){ delete[] title; }

void Header::read_data(std::ifstream& database_header, int offset){
  //reads the title in database_header binary file
  //database_header.seekg(offset*sizeof(uint8_t));//skip to the wanted entry in .phr file
  //database_header.seekg(7*sizeof(uint8_t)); //skip directly to the title
  database_header.seekg(offset+7);
  database_header.read((char*)&length,sizeof(uint8_t));//read the length of the title
  if (! msb_find(length)){//most significant bit is off
    //the byte is the length of the title string
    title = new char[length];
    database_header.read(title,length*sizeof(char));
  } else{//most significant bit is on
      //the byte is the number of bytes containing length of the title
      int title_length = length_find(database_header,length);
      if(title_length != -1)
        { title = new char[title_length]; database_header.read((char*)&title, title_length*sizeof(char)); }
      else
        title = (char *)"Empty header";
      //std::cout << title << std::endl;
    }
  }

char* Header::get_title() const {
  return title;
}

uint8_t Header::get_length() const{
  return length;
}

int Header::msb_find(uint8_t in_byte){
  //finds the most significant bit of a byte
  unsigned char msb = in_byte >> 7;//most significant bit is the last one
  return (int) msb;
}

int Header::length_find(std::ifstream& database_header, uint8_t nb_bytes){
  nb_bytes -= 128;
  if (nb_bytes == 1){//the length is encoded on a byte
    uint8_t length_8;
    database_header.read((char*)&length_8, sizeof(uint8_t));
    return length_8;
  }
  else if (nb_bytes == 2){
    uint16_t length_16;
    database_header.read((char*)&length_16, sizeof(uint16_t));
    length_16 = __bswap_16(length_16);//conversion from litlle endian to big endian
    return length_16;
  }
  else if (nb_bytes == 4){
    uint32_t length_32;
    database_header.read((char*)&length_32, sizeof(uint32_t));
    length_32 = __bswap_32(length_32);//conversion from litlle endian to big endian
    return length_32;
  }
  else if (nb_bytes == 8){
    uint64_t length_64;
    database_header.read((char*)&length_64, sizeof(uint64_t));
    length_64 = __bswap_64(length_64);//conversion from litlle endian to big endian
    return length_64;
  }
  else{//Length is unlikely to be coded on another number of bytes
    std::cout << "Incorrect number of bytes"<<std::endl;
    return -1;
  }
}
