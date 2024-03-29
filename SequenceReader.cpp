#include "SequenceReader.h"
#include <iostream>
#include <vector>

SequenceReader::SequenceReader(Index* index, std::ifstream& database_sequence){
    sq_index = index;
    read_data(database_sequence);
}

void SequenceReader::read_data(std::ifstream& database_sequence){
    /* there are residue_count + number_of_sequences + 1 residues in the file
    (residue_count is all non null residues; we have then N+1 separators (1 for each sequence + beginning/end)) */
    int total_residue_number = sq_index->get_residue_count() + sq_index->get_number_of_sequences() + 1;
    sq_buffer = new uint8_t[total_residue_number];
    database_sequence.read((char*)sq_buffer, sizeof(uint8_t)*total_residue_number);
}

SequenceReader::SequenceReader(Index* index){
/**
 * Constructor 
 *
 * Attributes initialization.
 */
    sq_index = index;
    sequences = new std::string[(int)index->get_number_of_sequences()];
}

SequenceReader::~SequenceReader(){
/**
 * Destructor 
 *
 * Frees allocated memory.
 */
    delete[] sequences;
}

int SequenceReader::exact_match(std::ifstream& database_sequence,std::ifstream& query_protein){
/**
 * Compares each sequence in the database and returns match if there is one.
 *
 * @param query_protein Query protein file.
 * @return Index of matching sequence in the database, or -1 if no match.
 */
    std::vector<int> query_sequence;
    std::string line;
    getline(query_protein,line); //skip the first line which contains the header
    char residue;
    while(query_protein.get(residue)){ //read query_protein sequence
        std::map<char,int>::const_iterator it = residue_char_conversion_map.find(residue);
        if(it != residue_char_conversion_map.end()){
            query_sequence.push_back(residue_char_conversion_map.at(residue));
        }
    }
    read_data(database_sequence);
    //int total_residue_number = sq_index->get_residue_count() + sq_index->get_number_of_sequences() + 1;
    bool match = true;
    for(unsigned i = 0; i < sq_index->get_number_of_sequences(); i++){ // skip first 0 byte
        int length = sq_index->get_sequence_offset_table()[i+1] - sq_index->get_sequence_offset_table()[i];
        int offset = sq_index->get_sequence_offset_table()[i];
        for(unsigned j = 0; j < length-1; j++){ // length-1 to avoid last 0
            if(j > query_sequence.size() || (int)sq_buffer[offset+j] != query_sequence[j]){
                match = false;
            }
        }
        if(match){
            return i;
        } else {
            match = true;
        }
    }
    return -1;
}

std::vector<int> SequenceReader::convert_query_sequence(std::ifstream& query_protein)
/**
 * Converts opened stream to vector of integers representing query sequence.
 *
 * @param query_protein Query protein file stream.
 * @return Vector of integers.
 */
{
    std::vector<int> query_sequence;
    std::string line;
    getline(query_protein,line); //skip the first line which contains the header
    char residue;
    while(query_protein.get(residue))
    {   
        //read query_protein sequence
        std::map<char,int>::const_iterator it = residue_char_conversion_map.find(residue);
        if(it != residue_char_conversion_map.end())
            query_sequence.push_back(residue_char_conversion_map.at(residue));
    }
    query_size = query_sequence.size();
    return query_sequence; // return pointer to first element as they are stored contiguously
}

int SequenceReader::get_query_size() const { return query_size; }