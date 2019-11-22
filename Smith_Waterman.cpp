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
    matrix =  new int[sequence1.length()][sequence2.length()];
}