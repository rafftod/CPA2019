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
#include <boost/thread/thread.hpp>

using namespace std;

struct Sequence {
    int id;
    int score;
};

int check_cores()
{
    int numCPU = 0;
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    numCPU = sysinfo.dwNumberOfProcessors
    #endif

    #ifdef __linux__
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    #endif

    #ifdef __APPLE__
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    #endif
}

struct arguments {
    //arguments needed by a thread
    const int n_seq;
    struct Sequence* sequences; // array of sequences, with id and score
    const int* query_protein;
    const uint8_t *db_seq; int db_seq_length;
    const int query_size;
    int offset;
    int offset2;
    SequenceReader* seq_reader;
    Smith_Waterman* sw;

};

void* routine(void* args)
//Process carried out by thread
{
    struct arguments* arguments = (struct arguments*)args;
    for(int i = arguments->offset; i < arguments->n_seq+arguments->offset2; ++i)
            {
                arguments->db_seq = arguments->seq_reader->get_sequence(i);
                arguments->db_seq_length = arguments->seq_reader->get_sequence_length(i);
                arguments->sequences[i-arguments->offset].score = arguments->sw->compare2(arguments->db_seq, arguments->query_protein, arguments->db_seq_length+1, arguments->query_size+1);
                std::cout << "Sequence " << i << " score : " << arguments->sequences[i-arguments->offset].score << std::endl;
                arguments->sequences[i-arguments->offset].id = i;
            }
}

void manage_seq(struct Sequence* & sequences,const int* & query_protein, const uint8_t* & db_seq, int db_seq_length,
 const int query_size, const int offset, const int offset2, Smith_Waterman& sw, const int n_seq, SequenceReader& seq_reader)
{
  for(int i = offset+offset2; i < n_seq+offset+offset2; ++i)
    {
        db_seq = seq_reader.get_sequence(i);
        db_seq_length = seq_reader.get_sequence_length(i);
        sequences[i-offset].score = sw.compare(db_seq, query_protein, db_seq_length+1, query_size+1);
        std::cout << "Sequence " << i << " score : " << sequences[i-offset].score << std::endl;
        sequences[i-offset].id = i;
    }
}

int main(int argc, char const *argv[]) {

    /* Determination of optional command-line arguments given */

    int gap_open_penalty = 11;
    int gap_expansion_penalty = 1;
    string blosum_path = "BLOSUM62";

    if (argc < 3) {
        cout << "Too few arguments given." << '\n';
        exit(EXIT_FAILURE);
    }

    for(int i = 1; i < argc; ++i){
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
            //const int n_seq = index->get_number_of_sequences();
            const int n_seq = 100;
            struct Sequence sequences[n_seq]; // array of sequences, with id and score
            const std::vector<int> query_protein_vec = seq_reader->convert_query_sequence(protein);
            const int* query_protein = &query_protein_vec[0];
            const uint8_t *db_seq; int db_seq_length;
            const int query_size = seq_reader->get_query_size();
            const int offset = 116000;

            //creation of as many threads as there are cores on the machine
            int n = check_cores();
            boost::thread_group threads;
            for (int i = 0; i < n; i++)
            {
                if (i!= n-1)
                {
                    struct arguments args = {n_seq, sequences,query_protein,db_seq,db_seq_length,query_size,offset+i*n_seq/n,(i+1)*n_seq/n -1 , seq_reader,sw};
                    boost::thread th = boost::thread(routine,boost::ref(args));
                    threads.add_thread(&th);
                }
                else
                {
                    struct arguments args = {n_seq, sequences,query_protein,db_seq,db_seq_length,query_size,offset+i*n_seq/n,(i+1)*n_seq/n + n_seq%n, seq_reader,sw};
                    boost::thread th = boost::thread(routine,boost::ref(args));
                    threads.add_thread(&th);                    
                }
            }
            threads.join_all();
            
            /*for(int i = offset; i < n_seq+offset; ++i)
            {
                db_seq = seq_reader->get_sequence(i);
                db_seq_length = seq_reader->get_sequence_length(i);
                sequences[i-offset].score = sw->compare(db_seq, query_protein, db_seq_length+1, query_size+1);
                std::cout << "Sequence " << i << " score : " << sequences[i-offset].score << std::endl;
                sequences[i-offset].id = i;
            }*/

            // sorting sequences by score
            std::sort(sequences, sequences+n_seq-1,
                        [](struct Sequence const & s1, struct Sequence const & s2) -> bool
                        { return s1.score > s2.score; }); // sort by score value with lambda 
            
            Header* header = new Header();
            for(int i = 0; i < 20; ++i)
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
