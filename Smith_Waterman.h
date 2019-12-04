#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN
#include <iostream>
#include <string>
#include <vector>
#include <map>
#define BLOSUM_SIZE 24 // All BLOSUM matrix are 24 x 24

class Smith_Waterman{
private:
    int** matrix;
    int gap_penalty_open = 11;
    int gap_penalty_exp = 1;
    int match;
    int no_match;
    int blosum_matrix[BLOSUM_SIZE][BLOSUM_SIZE];
    std::map<int, int> residue_int_to_blosum_pos_map;

public:
    Smith_Waterman(int gap1,int gap2,int m, int nm);
    Smith_Waterman(int gap_open_penalty, int gap_expansion_penalty, std::string blosum_path);
    ~Smith_Waterman();

    void compare(uint8_t* sequence1,std::vector<int>& sequence2, int length1, int length2);
    int find_max(int a, int b, int c);
    int max_column(int i, int j, int** matrix);
    int max_row(int i, int j, int** matrix);
    int build_BLOSUM(std::string path);
};
#endif