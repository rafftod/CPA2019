#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN
#include <iostream>
#include <string>
#include <vector>

class Smith_Waterman{
private:
    int** matrix;
    int gap_penalty_open = 11;
    int gap_penalty_exp = 1;
    int match;
    int no_match;
    int blosum_matrix[24][24]; // Al BLOSUM matrix are 24 x 24

public:
    Smith_Waterman(int gap1,int gap2,int m, int nm);
    ~Smith_Waterman();

    void compare(uint8_t* sequence1,std::vector<int>& sequence2, int length1, int length2);
    int find_max(int a, int b, int c);
    int max_column(int i, int j, int** matrix);
    int max_row(int i, int j, int** matrix);

};
#endif