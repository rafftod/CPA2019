#ifndef SMITH_WATERMAN
#define SMITH_WATERMAN
#include <iostream>
#include <string>
#include <vector>
#include <map>
#define BLOSUM_SIZE 26 // All BLOSUM matrix are 24 x 24, but ours are 26 x 26 to match with BLAST number of residues

class Smith_Waterman{
private:
    int gap_penalty_open = 11;
    int gap_penalty_exp = 1;
    //std::map<int, int> residue_int_to_blosum_pos_map;
    int blosum_column_to_blast_number(int j) const ;

public:
    Smith_Waterman(int gap1,int gap2,int m, int nm);
    Smith_Waterman(int gap_open_penalty, int gap_expansion_penalty, std::string blosum_path);
    ~Smith_Waterman();

    int compare(const uint8_t* & sequence1, const int*& sequence2, const int length1, const int length2) const;
    int compare2(const uint8_t* & sequence1, const int* & sequence2, int length1, int length2);
    int find_max(int a, int b, int c);
    void max_column(int i, int j, int** matrix, int** max_col_matrix);
    void max_row(int i, int j, int** matrix, int** max_row_matrix);
    int build_BLOSUM(std::string path);
    int get_blosum_element(const int index) const;
};
#endif