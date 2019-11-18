#include "SequenceReader.h"
#include <iostream>

using namespace std;

void SequenceReader::read_data(ifstream database_sequence){
    uint8_t current_residue;
    string current_sequence = "";
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
}

SequenceReader::SequenceReader(Index* index){
    sq_index = index;
    sequences = new string[(int)index->get_number_of_sequences()];
}

SequenceReader::~SequenceReader(){
    
}

string SequenceReader::get_sequence(int i) const{
    if(i < sq_index->get_number_of_sequences()){
        return sequences[i];
    }
    throw "Index out of bounds for sequences vector";
}

int SequenceReader::search_sequences(std::ifstream& query_protein){

    string sequence = "";
    string line;
    getline(query_protein,line);//skip the first line which contains the header
    while(getline(query_protein,line)){//read query_protein sequence
        sequence = sequence + line;
    }

    int i;
    for(i=0; (i < sq_index->get_number_of_sequences()) && (sequences[i]!=sequence);i++){}//search all sequences for exact match
    return i;
}