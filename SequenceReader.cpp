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

uint8_t* SequenceReader::get_sequence(int i) const
{
/**
 * Finds sequence at given index in database, and throws exception if index is out of bounds.
 *
 * @param i Number of the sequence in the database.
 * @return Sequence found.
 */
    if(i < sq_index->get_number_of_sequences())
    {
        //int length = sq_index->get_sequence_offset_table()[i+1] - sq_index->get_sequence_offset_table()[i];
        int offset = sq_index->get_sequence_offset_table()[i];
        return &sq_buffer[offset];
    }
    throw "Index out of bounds for sequences.";
}

int SequenceReader::get_sequence_length(int i) const
{
    if(i < sq_index->get_number_of_sequences())
    {
        return sq_index->get_sequence_offset_table()[i+1] - sq_index->get_sequence_offset_table()[i]-1; // -1 to remove last 0
    }
    return -1;
}

int SequenceReader::search_sequences(std::ifstream& query_protein){
/**
 * Compares each sequence in the database and returns match if there is one.
 *
 * @param query_protein Query protein file.
 * @return Index of matching sequence in the database, or -1 if no match.
 */

    std::string sequence = "";
    std::string line;
    getline(query_protein,line);//skip the first line which contains the header
    while(getline(query_protein,line)){//read query_protein sequence
        sequence = sequence + line;
    }
    for(int i=0; i < sq_index->get_number_of_sequences();i++){ // search database for exact match
        if(sequences[i]==sequence){
            return i;
        }
    }
    return -1;
}

int SequenceReader::exact_match(std::ifstream& database_sequence,std::ifstream& query_protein){
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
    int total_residue_number = sq_index->get_residue_count() + sq_index->get_number_of_sequences() + 1;
    bool match = true;
    for(int i = 0; i < sq_index->get_number_of_sequences(); i++){ // skip first 0 byte
        int length = sq_index->get_sequence_offset_table()[i+1] - sq_index->get_sequence_offset_table()[i];
        int offset = sq_index->get_sequence_offset_table()[i];
        for(int j = 0; j < length-1; j++){ // length-1 to avoid last 0
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