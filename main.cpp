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
  ifstream database(argv[1]); ifstream protein(argv[2]); //
  if (database.is_open()) {
    if(protein.is_open()){
      // TODO
    } else {
      cout << "Protein file couldn't be read" << '\n';
      exit(EXIT_FAILURE);
    }
    string line;
    //read database & search for match
    database.close();
    protein.close();
  } else {
    cout << "Database file couldn't be read" << '\n';
    exit(EXIT_FAILURE);
  }
  return 0;
}
