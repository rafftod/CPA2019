#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <unistd.h>
#include "Index.h"
#include "SequenceReader.h"
#include "Header.h"
#include "Smith_Waterman.h"

using namespace std;

struct Sequence {
    int id;
    int score;
};

int main(int argc, char const *argv[]) {

    /* Determination of optional command-line arguments given */

    int gap_open_penalty = 11;
    int gap_expansion_penalty = 1;
    string blosum_path = "BLOSUM62";

    if (argc < 3) {
        cout << "Too few arguments given." << '\n';
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i < argc; i++){
        string current_arg = (string) argv[i];
        if(current_arg == "-o") {
            // Optional gap open penalty is set
            string next_arg = (string) argv[i+1];
            if(argv[i+1] == NULL || !any_of(next_arg.begin(), next_arg.end(), ::isdigit)) {
                cout << "Invalid gap open penalty argument." << endl;
                exit(EXIT_FAILURE);
            } else {
                gap_open_penalty = atoi(next_arg.c_str());
                continue; // we can skip next argument as it is the gap open penalty value
            } 
        } else if(current_arg == "-e") {
            // Optional gap expansion penalty is set
            string next_arg = (string) argv[i+1];
            if(argv[i+1] == NULL || !any_of(next_arg.begin(), next_arg.end(), ::isdigit)) {
                cout << "Invalid gap expansion penalty argument." << endl;
                exit(EXIT_FAILURE);
            } else {
                gap_expansion_penalty = atoi(next_arg.c_str());
                continue; // we can skip next argument as it is the gap expansion penalty value
            }
        } else if(current_arg == "-b") {
            // Optional BLOSUM matrix is set
            string next_arg = (string) argv[i+1];
            if(argv[i+1] == NULL || access(next_arg.c_str(), F_OK) != -1) {
                cout << "Invalid BLOSUM matrix path." << endl;
                exit(EXIT_FAILURE);
            } else {
                blosum_path = next_arg;
                continue; // we can skip next argument as it is the BLOSUM file path
            }
        }
    }

    /* Open files */

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

            /* Database treatment */

            Index* index = new Index();
            index->read_data(database_index); // read index
            SequenceReader* seq_reader = new SequenceReader(index, database_sequence);
            Smith_Waterman* sw = new Smith_Waterman(gap_open_penalty, gap_expansion_penalty, blosum_path);
            const int n_seq = index->get_number_of_sequences();
            //const int n_seq = 500;
            struct Sequence sequences[n_seq]; // array of sequences, with id and score
            std::vector<int> query_protein = seq_reader->convert_query_sequence(protein);
            const uint8_t *db_seq; int db_seq_length;
            for(int i = 0; i < n_seq; i++)
            {
                db_seq = seq_reader->get_sequence(i);
                db_seq_length = seq_reader->get_sequence_length(i);
                sequences[i].score = sw->compare(db_seq, query_protein, db_seq_length+1, query_protein.size()+1);
                std::cout << "Sequence " << i << " score : " << sequences[i].score << std::endl;
                sequences[i].id = i;
            }
            // sorting sequences by score
            std::sort(sequences, sequences+n_seq-1,
                        [](struct Sequence const & s1, struct Sequence const & s2) -> bool
                        { return s1.score > s2.score; }); // sort by score value with lambda 
            
            Header* header = new Header();
            for(int i = 0; i < 20; i++)
            {
                int sq_offset = index->get_header_offset_table()[sequences[i].id];
                header->read_data(database_header, sq_offset);
                std::cout << header->get_title()+'\0' << "\n Score : " << sequences[i].score << std::endl;
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
