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
#include <pthread.h>


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

    return numCPU;
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

void* routine(struct arguments& args)
//Process carried out by thread
{
    for(int i = args.offset; i < args.n_seq+args.offset2; ++i)
            {
                args.db_seq = args.seq_reader->get_sequence(i);
                args.db_seq_length = args.seq_reader->get_sequence_length(i);
                args.sequences[i-args.offset].score = args.sw->compare2(args.db_seq, args.query_protein, args.db_seq_length+1, args.query_size+1);
                std::cout << "Sequence " << i << " score : " << args.sequences[i-args.offset].score << std::endl;
                args.sequences[i-args.offset].id = i;
            }
}

void* routine2(void* argz)
//Process carried out by thread
{   arguments* args = (arguments*) argz;
    std::cout << args->n_seq << std::endl;
    for(int i = args->offset; i < args->n_seq+args->offset; ++i)
    {
        args->db_seq = args->seq_reader->get_sequence(i);
        args->db_seq_length = args->seq_reader->get_sequence_length(i);
        args->sequences[i-args->offset2].score = args->sw->compare(args->db_seq, args->query_protein, args->db_seq_length+1, args->query_size+1);
        std::cout << "Sequence " << i << " score : " << args->sequences[i-args->offset].score << std::endl;
        args->sequences[i-args->offset2].id = i;
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
        cout << "Too few arguments given. Minimum : main database protein." << '\n';
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
            ifstream blosum(next_arg);
            if(argv[i+1] == NULL || !blosum.good()) {
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
            const int n_seq = 1001;
            struct Sequence sequences[n_seq]; // array of sequences, with id and score
            const std::vector<int> query_protein_vec = seq_reader->convert_query_sequence(protein);
            const int* query_protein = &query_protein_vec[0];
            const uint8_t *db_seq; int db_seq_length;
            const int query_size = seq_reader->get_query_size();
            const int offset = 116000;

            //creation of as many threads as there are cores on the machine
            int n = check_cores();
            pthread_t threads[n];
            int offsets[n]; int thread_n_seq[n];
            for (int i = 0; i < n; i++)
            {
                if (i != 0)
                {
                    offsets[i] = offsets[i-1] + thread_n_seq[i-1];
                    thread_n_seq[i] = n_seq/n;
                    struct arguments args = {thread_n_seq[i], sequences,query_protein,db_seq,db_seq_length,query_size,offsets[i],offset, seq_reader,sw};
                    pthread_create(&threads[i],NULL,routine2,(void*)&args);
                }
                else
                {
                    // first thread will do 1/n sequences + rest of division (n_seq%n)
                    offsets[i] = offset;
                    thread_n_seq[i] = n_seq/n + n_seq%n;
                    struct arguments args = {thread_n_seq[i], sequences,query_protein,db_seq,db_seq_length,query_size,offsets[i],offset, seq_reader,sw};
                    pthread_create(&threads[i],NULL,routine2,(void*)&args);                  
                }
                pthread_join(threads[i],NULL);
            }

            for (int i = 0; i < n; i++)
            {
                pthread_join(threads[i],NULL);
            }
            
            

            /*
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
            */
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
