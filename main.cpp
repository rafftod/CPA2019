#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Index.h"
#include "SequenceReader.h"
#include "Header.h"

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
    database_index.open((string) argv[1] + ".pin", ios::binary);
    database_header.open((string) argv[1] + ".phr", ios::binary);
    database_sequence.open((string) argv[1] + ".psq", ios::binary);
    // FASTA can be read as normal text
    protein.open(argv[2]);
    if (database_index.is_open() && database_header.is_open() && database_sequence.is_open()) {
        if (protein.is_open()) {
            Index* index = new Index();
            index->read_data(database_index); // read index
            SequenceReader* seq_reader = new SequenceReader(index);
            int position = seq_reader->exact_match(database_sequence, protein);
            if (position != -1){
                cout << "Exact match at position " << position << ":" << endl; 
                Header* header = new Header();
                header->read_data(database_header, index->get_header_offset_table()[position]);
                for(int i = 0; i < header->get_length(); i++)
                  cout << header->get_title()[i];
                cout << endl;
                delete header;
            } else {
                cout << "There is no matching protein in the database." << endl;
            }
            delete seq_reader;
            delete index;
        } else {
            throw "Protein file couldn't be read.";
            exit(EXIT_FAILURE);
        }
    } else {
        throw "Database files couldn't be read.";
        exit(EXIT_FAILURE);
    }
    database_header.close();
    database_index.close();
    database_sequence.close();
    protein.close();
    return 0;
}
