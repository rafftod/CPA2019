#ifndef SEQUENCE_H
#define SEQUENCE_H
#include "SequenceReader.h"
#include <iostream>
#include <vector>

/* void SequenceReader::read_data(std::ifstream& database_sequence){
/**
 * Reads data from database sequence file opened stream.
 *
 * Separates sequences in file with NUL byte, and stores each in string pointer.
 * 
 * To make the program faster for intermediate deadline, we could only read sequences
 * compatible with query protein (by filtering with size, comparing residues progressively...),
 * but, as we need to evaluate all database sequences score for the final deadline, we already store all
 * sequences even if we don't need all of them.
 *
 * @param database_sequence Database file opened ifstream.
    uint8_t current_residue; // each residue is stored on 8 bits
    std::string current_sequence = "";
    int i = 0;
    while(database_sequence.read((char*)&current_residue, sizeof(uint8_t))){
        if((int)current_residue == 0) { // NUL byte considered as separator
            if(current_sequence != ""){ // to make sure there is no empty sequence @ beginning or end
                sequences[i++]=current_sequence;
                current_sequence = "";
            }
        } else {
            current_sequence.push_back(residue_int_conversion_map.at((int)current_residue)); // at instead of [] to ensure limit check
        }
    }
} */

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

std::string SequenceReader::get_sequence(int i) const{
/**
 * Finds sequence at given index in database, and throws exception if index is out of bounds.
 *
 * @param i Number of the sequence in the database.
 * @return Sequence found.
 */
    if(i < sq_index->get_number_of_sequences()){
        return sequences[i];
    }
    throw "Index out of bounds for sequences vector";
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

/* int SequenceReader::exact_match(std::ifstream& database_sequence, std::ifstream& query_protein){
    // Read query protein 
    std::string query_sequence = "";
    std::string line;
    getline(query_protein,line);//skip the first line which contains the header
    while(getline(query_protein,line)){//read query_protein sequence
        query_sequence += line;
    }
    // Read database and return when exact match 
    uint8_t current_residue; // each residue is stored on 8 bits
    std::string current_sequence = "";
    int i = 0;
    while(database_sequence.read((char*)&current_residue, sizeof(uint8_t))){
        if((int)current_residue == 0) { // NUL byte considered as separator
            if(current_sequence != ""){ // to make sure there is no empty sequence @ beginning or end
                if(current_sequence==query_sequence){
                    return i;
                }
                i++;
                current_sequence = "";
            }
        } else {
            current_sequence.push_back(residue_int_conversion_map.at((int)current_residue)); // at instead of [] to ensure limit check
        }
    }
    return -1;

} */

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

#endif 