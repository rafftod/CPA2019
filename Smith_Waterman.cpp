#include "Smith_Waterman.h"
#include <string>
#include <iostream>

Smith_Waterman::Smith_Waterman(int gap1, int gap2, int m, int nm){
    gap_penalty_open = gap1;
    gap_penalty_exp = gap2;
    match = m;
    no_match = nm;
}

Smith_Waterman::~Smith_Waterman(){}

void Smith_Waterman::compare(std::string& sequence1, std::string& sequence2){
    //compare 2 sequences with the Smith-Waterman algorithm

    //initialize score matrix
    
    matrix =  new int*[sequence1.length()];//sequence1.length() = number of rows
    for (int i = 0; i < sequence1.length(); ++i){
        matrix[i] = new int[sequence2.length()];//sequence.length() = number of columns
    }

    //first row and column are set at 0

    for (int i = 0; i < sequence1.length(); i++){
        matrix[i][0]=0;
    }

    for (int j = 1; j < sequence2.length(); j++){
        matrix[0][j]=0;
    }

    //compare each residu and change score accordingly

    for(int i = 0; i < sequence1.length(); i++){
        for (int j = 0; j < sequence2.length(); j++){

        }
    }
}