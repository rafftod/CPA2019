#include "SequenceReader.h"
#include <iostream>

using namespace std;

vector<string> SequenceReader::read_data(ifstream database_sequence){
    int8_t current_residue;
    string current_sequence = "";
    while(database_sequence.read((char*)&current_residue, sizeof(int8_t))){
        if((int)current_residue == 0) { // NUL byte considered as separator
            if(current_sequence != ""){ // to make sure there is no empty sequence @ beginning or end
                sequences.push_back(current_sequence);
                current_sequence = "";
            }
        } else {
            current_sequence.push_back(residue_int_conversion_map.at((int)current_residue)); // at instead of [] to ensure limit check
        }
    }
}

SequenceReader::SequenceReader(){

}

SequenceReader::~SequenceReader(){
    
}

string SequenceReader::get_sequence(int index) const{
    if(index < sequences.size()){
        return sequences.at(index);
    }
    throw "Index out of bounds for sequences vector";
}