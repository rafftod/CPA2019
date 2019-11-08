#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char const *argv[]) {
  // For now, we suppose we don't have optional arguments
  // and only database and query protein
  // First command-line argument is database; second is query protein
  if (argc < 2) {
    cout << "Too few arguments" << '\n';
    exit(EXIT_FAILURE);
  }
  ifstream database_index; ifstream database_header; ifstream database_sequence;
  ifstream protein;
  // we need to open the 3 binary files of the BLAST format
  database_index.open((string) argv[1] + ".pin", ios::binary | ios::in);
  database_header.open((string) argv[1] + ".phr", ios::binary | ios::in);
  database_sequence.open((string) argv[1] + ".psq", ios::binary | ios::in);
  // FASTA can be read as normal text
  protein.open(argv[2]);
  if (database_index.is_open() && database_header.is_open() && database_sequence.is_open()) {
    if (protein.is_open()) {
      // TODO : read database and look for match
    } else {
      cout << "Protein file couldn't be read" << '\n';
      exit(EXIT_FAILURE);
    }
  } else {
    cout << "Database files couldn't be read" << '\n';
    exit(EXIT_FAILURE);
  }
  database_header.close();
  database_index.close();
  database_sequence.close();
  protein.close();
  return 0;
}