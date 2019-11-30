#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN
#include <iostream>
#include <string>

class Smith_Waterman{
private:
    int** matrix;
    int gap_penalty_open;
    int gap_penalty_exp;
    int match;
    int no_match;
    int blosom_matrix[24][24]; // Al BLOSUM matrix are 24 x 24

public:
    Smith_Waterman(int gap1,int gap2,int m, int nm);
    ~Smith_Waterman();

    void compare(std::string& sequence1,std::string& sequence2);

};
#endif