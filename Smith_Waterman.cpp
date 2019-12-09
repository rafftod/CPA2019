#include "Smith_Waterman.h"
#include "SequenceReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>

Smith_Waterman::Smith_Waterman(int gap1, int gap2, int m, int nm){
    gap_penalty_open = gap1;
    gap_penalty_exp = gap2;
    match = m;
    no_match = nm;
}

Smith_Waterman::Smith_Waterman(int gap_open_penalty, int gap_expansion_penalty, std::string blosum_path){
    gap_penalty_open = gap_open_penalty;
    gap_penalty_exp = gap_expansion_penalty;
    residue_int_to_blosum_pos_map = {};
    build_BLOSUM(blosum_path);
}

Smith_Waterman::~Smith_Waterman(){}

int Smith_Waterman::compare(uint8_t* sequence1, std::vector<int>& sequence2, int length1, int length2)
{
    //compare 2 sequences with the Smith-Waterman algorithm and returns score
    //sequence1 : sequence from the database
    //sequence2 : query sequence
    //length1 : number of residues in sequence1 + 1
    //length2 : number of residues in sequence2 + 1

    //initialize score matrix
    
    matrix =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        matrix[i] = new int[length2];//length2 = number of columns
    }

    //initialize memoisation matrices
    int** max_col_matrix =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        max_col_matrix[i] = new int[length2];//length2 = number of columns
    }

    int** max_row_matrix =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        max_row_matrix[i] = new int[length2];
    }

    //first row and column are set at 0

    for (int i = 0; i < length1; i++)
    {
        matrix[i][0]=0;
        max_row_matrix[i][0]=0;
    }

    for (int j = 1; j < length2; j++)
    {
        matrix[0][j]=0;
        max_col_matrix[0][j] = 0;
    }

    //scoring matrix

    int score = 0;

    for(int i = 1; i < length1; i++)
    {
        for (int j = 1; j < length2; j++)
        {
            int blosum_i, blosum_j;
            try {
                blosum_i = residue_int_to_blosum_pos_map.at((int)sequence1[i-1]);
            } catch (std::string const & e) { // in case our character is not in the blosum matrix
                blosum_i = 23;
            }
            try {
                blosum_j = residue_int_to_blosum_pos_map.at(sequence2[j-1]);
            } catch (std::string const & e) {
                blosum_j = 23;
            }
            int a = matrix[i-1][j-1] + blosum_matrix[blosum_i][blosum_j];
            /* b is maximum on row i */
            int column_max = 0;
            int b = 0;
            for(int k = 0; k < i; k++)
            {
                if(matrix[k][j] > column_max)
                    b = matrix[k][j] - gap_penalty_open - gap_penalty_exp*abs(i-k);
            }
            /* c is maximum on column j */
            int row_max = 0;
            int c = 0;
            for(int k = 0; k < j; k++)
            {
                if(matrix[i][k] > row_max)
                    c = matrix[i][k] - gap_penalty_open - gap_penalty_exp*abs(j-k);
            }
            matrix[i][j] = std::max({a,b,c,0});
            //max_column(i,j,matrix, max_col_matrix);
            //max_row(i,j,matrix, max_row_matrix);
            //matrix[i][j] = find_max(a,max_row_matrix[i][j],max_col_matrix[i][j]);
            if(matrix[i][j] > score)
            {
                score = matrix[i][j]; 
            }
        }
    }
    return score;
}

int Smith_Waterman::find_max(int a, int b, int c)
{
    /* returns the maximum of {a,b,c}, or 0 if all are negatives */
    int max = 0;
    if (a > max)
    {
        max = a;
    }

    if (b > max)
    {
        max = b;
    }
    if (c > max)
    {
        max = c;
    }
    return max;
}

void Smith_Waterman::max_row(int i, int j, int** matrix, int** max_row_matrix)
//returns the maximum score on a column
{
    max_row_matrix[i][j] = std::max(matrix[i][j-1] - gap_penalty_open - gap_penalty_exp, max_row_matrix[i][j-1] - gap_penalty_exp);
}

void Smith_Waterman::max_column(int i, int j, int** matrix,int** max_col_matrix)
//returns the maximum score on a line
{
    max_col_matrix[i][j] = std::max(matrix[i-1][j] - gap_penalty_open - gap_penalty_exp, max_col_matrix[i-1][j] - gap_penalty_exp);
}

int Smith_Waterman::build_BLOSUM(const std::string path)
{
    /* Builds BLOSUM matrix from given path */
    std::ifstream blosum_file(path);
    std::string s;
    std::vector<std::string> line_vector;
    for(int i = 0; i < 7; i++) {std::getline(blosum_file,s);} // ignore first 6 lines and use 7th to build conversion map
    int i = 0; 
    for(char c : s)
    {
        if(isalpha(c) || c == '*')
        {
            residue_int_to_blosum_pos_map.insert(std::pair<int,int> (residue_char_conversion_map.at(c), i++));
        }
    }
    int sign = 1; char current_char;
    i = 0; int j = 0;
    while(getline(blosum_file,s))
    {
        boost::split(line_vector, s, boost::is_any_of(" "));
        for(std::string sub: line_vector)
        {
            if(!isdigit(sub[sub.size()-1])) { continue; }
            int digit = strtol(sub.c_str(), NULL, 10);
                if(j > BLOSUM_SIZE-1) 
                {
                    i++; 
                    j = 0;
                }
                blosum_matrix[i][j++] = digit;
        }
    }
    return 0;
}