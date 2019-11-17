#include <iostream>
#include <fstream>
#include <cstring>

#include "Header.h"


using namespace std;

int msb_find(uint8_t in_byte);

int main(int argc, char const *argv[]) {
  string fname = "P00533.fasta.phr";
  ifstream in (fname,ios::binary);
  if (in.is_open()){
    Header* header = new Header();
    //header->read_data(move(in));
    header->read_data(in,0);
    cout << header->get_title() << endl;
    in.close();
  } else {
    cout << "Could not open " << fname << ": " << strerror(errno) << endl;
  }

  return 0;
}
