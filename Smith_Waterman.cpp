#include "Smith_Waterman.h"
#include "SequenceReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <boost/algorithm/string.hpp>
#define sw_max2(a,b) ((a)>(b)?(a):(b))
#define sw_max3(a,b,c) sw_max2(sw_max2(a,b),c)
#define sw_max(a,b,c,d) sw_max2(sw_max3(a,b,c),d)

Smith_Waterman::Smith_Waterman(int gap1, int gap2, int m, int nm){
    gap_penalty_open = gap1;
    gap_penalty_exp = gap2;
}

Smith_Waterman::Smith_Waterman(int gap_open_penalty, int gap_expansion_penalty, std::string blosum_path){
    gap_penalty_open = gap_open_penalty;
    gap_penalty_exp = gap_expansion_penalty;
    residue_int_to_blosum_pos_map = {};
    build_BLOSUM(blosum_path);
}

Smith_Waterman::~Smith_Waterman(){}

int Smith_Waterman::compare(const uint8_t* & sequence1, const std::vector<int> & sequence2, const int length1, const int length2) const
{
    //compare 2 sequences with the Smith-Waterman algorithm and returns score
    //sequence1 : sequence from the database
    //sequence2 : query sequence
    //length1 : number of residues in sequence1 + 1
    //length2 : number of residues in sequence2 + 1

    // initialize score matrix on the heap, because stack length limit is 1000
    /*
    int** matrix = new int*[length1];
    int* matrix_data = new int[length1*length2];
    for(int i = 0; i < length1; i++)
        matrix[i] = matrix_data + i*length2;
        */
    int** matrix = new int*[length1];
    for(int i = 0; i < length1; ++i)
        matrix[i] = new int[length2];
    
    // first row and column are set at 0

    for (int i = 0; i < length1; ++i)
    {
        matrix[i][0] = 0;
    }

    for (int j = 1; j < length2; j++)
    {
        matrix[0][j]=0;
    }

    //scoring matrix

    // memoization arrays : we need to memoize each row and each column maximum, and its position
    int maximum_on_rows[length1] = { -1 };
    int maximum_on_rows_index[length1] = { -1 };
    int maximum_on_columns[length2] = { -1 }; // by default, we set values at -1 to know if we memoized
    int maximum_on_columns_index[length2] = { -1 };

    int score = 0;

    for(int i = 1; i < length1; i++)
    {
        for (int j = 1; j < length2; j++)
        {
            int a = matrix[i-1][j-1] + blosum_matrix[get_blosum_element((int)sequence1[i-1])][get_blosum_element(sequence2[j-1])];
            //int a = matrix[index(i-1,j-1)] + blosum_matrix[get_blosum_element((int)sequence1[i-1])][get_blosum_element(sequence2[j-1])];
            /* b is maximum on row i */
            int b = 0;
            if(maximum_on_rows[i] != -1)
            {
                // we memoized the maximum of this row
                int max = maximum_on_rows[i];
                int max_index = maximum_on_rows_index[i];
                b = max - gap_penalty_open - gap_penalty_exp*abs(i-max_index);
            }
            else
            {
                int column_max = 0;
                for(int k = 0; k < i; k++)
                {
                    if(matrix[k][j] > column_max)
                    //if(matrix[index(k,j)] > column_max)
                    {
                        b = matrix[k][j] - gap_penalty_open - gap_penalty_exp*abs(i-k);
                        //b = matrix[index(k,j)] - gap_penalty_open - gap_penalty_exp*abs(i-k);
                        column_max = matrix[k][j];
                        //column_max = matrix[index(k,j)];
                        maximum_on_rows[j] = matrix[k][j];
                        //maximum_on_rows[j] = matrix[index(k,j)];
                        maximum_on_rows_index[j] = k;
                    }
                }
            }
            /* c is maximum on column j */
            int c = 0;
            if(maximum_on_columns[j] != -1)
            {
                // we memoized the maximum of this column
                int max = maximum_on_columns[j];
                int max_index = maximum_on_columns_index[j];
                c = max - gap_penalty_open - gap_penalty_exp*abs(j-max_index);
            }
            else
            {
                int row_max = 0;
                for(int k = 0; k < j; k++)
                {
                    if(matrix[i][k] > row_max)
                    //if(matrix[index(i,k)] > row_max)
                    {
                        c = matrix[i][k] - gap_penalty_open - gap_penalty_exp*abs(j-k);
                        //c = matrix[index(i,k)] - gap_penalty_open - gap_penalty_exp*abs(j-k);
                        row_max = matrix[i][k];
                        //row_max = matrix[index(i,k)];
                        maximum_on_columns[j] = matrix[i][k];
                        //maximum_on_columns[j] = matrix[index(i,k)];
                        maximum_on_columns_index[j] = k;
                    }
                }
            }
            //matrix[i][j] = std::max({a,b,c,0});
            matrix[i][j] = sw_max(a,b,c,0);
            //matrix[index(i,j)] = std::max({a,b,c,0});
            if(matrix[i][j] > score)
            //if(matrix[index(i,j)] > score)
            {
                score = matrix[i][j];
                //score = matrix[index(i,j)]; 
            }
        }
    }
    for(int i = 0; i < length1; i++)
        delete matrix[i];
    delete matrix;
    return score;
}

int Smith_Waterman::compare2(const uint8_t* & sequence1, const std::vector<int>& sequence2, int length1, int length2)
{
    //compare 2 sequences with the Smith-Waterman algorithm and returns score
    //sequence1 : sequence from the database
    //sequence2 : query sequence
    //length1 : number of residues in sequence1 + 1
    //length2 : number of residues in sequence2 + 1

    //initialize score matrix
    
    int** matrix =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        matrix[i] = new int[length2];//length2 = number of columns
    }

    
    //initialize memoisation matrices
    int** E =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        E[i] = new int[length2];//length2 = number of columns
    }

    int** F =  new int*[length1];//length1 = number of rows
    for (int i = 0; i < length1; ++i)
    {
        F[i] = new int[length2];
    }
    

    //first row and column are set at 0

    for (int i = 0; i < length1; i++)
    {
        matrix[i][0]=0;
        E[i][0]=0;
    }

    for (int j = 1; j < length2; j++)
    {
        matrix[0][j]=0;
        F[0][j]=0;
    }

    //scoring matrix
    int score = 0;

    for (int j = 1; j < length2; j++)
    {
        for (int i = 1; i < length1; i++)
        {
            /* int blosum_i, blosum_j;
            try {
                blosum_i = residue_int_to_blosum_pos_map.at((int)sequence1[i-1]);
            } catch (std::string const & e) { // in case our character is not in the blosum matrix
                blosum_i = 23;
            }
            try {
                blosum_j = residue_int_to_blosum_pos_map.at(sequence2[j-1]);
            } catch (std::string const & e) {
                blosum_j = 23;
            } */
            //int a = matrix[i-1][j-1] + blosum_matrix[blosum_i][blosum_j];
            int a = matrix[i-1][j-1] + blosum_matrix[get_blosum_element((int)sequence1[i-1])][get_blosum_element(sequence2[j-1])];
            E[i][j] = std::max(matrix[i][j-1]-gap_penalty_open - gap_penalty_exp,E[i][j-1] - gap_penalty_exp);
            int b = E[i][j];

            F[i][j]= std::max(matrix[i-1][j]- gap_penalty_open -gap_penalty_exp, F[i-1][j] -gap_penalty_exp);
            int c = F[i][j];

            if (matrix[i][j] = std::max({a,b,c,0}) > score)
            {
                score = matrix[i][j];
            }    
        } 
    }
    delete matrix;
    return score;
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
    /*for(char c : s)
    {
        if(isalpha(c) || c == '*')
        {
            residue_int_to_blosum_pos_map.insert(std::pair<int,int> (residue_char_conversion_map.at(c), i++));
        }
    }
    for(std::map<int,int>::const_iterator it = residue_int_to_blosum_pos_map.begin(); it != residue_int_to_blosum_pos_map.end(); ++it)
    {
        std::cout << "Map key : " << it->first << " | Map value : " << it->second << std::endl;
    }*/
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

int Smith_Waterman::get_blosum_element(int index) const
{
    /* this function acts as a map, but as we know the map content
    in advance, this is O(1) and not a O(log n) to find element with this function */
    switch (index)
    {
    case 1:
        return 0;
    case 2:
        return 20;
    case 3:
        return 4;
    case 4:
        return 3;
    case 5:
        return 6;
    case 6:
        return 13;
    case 7:
        return 7;
    case 8:
        return 8;
    case 9:
        return 9;
    case 10:
        return 11;
    case 11:
        return 10;
    case 12:
        return 12;
    case 13:
        return 2;
    case 14:
        return 14;
    case 15:
        return 5;
    case 16:
        return 1;
    case 17:
        return 15;
    case 18:
        return 16;
    case 19:
        return 19;
    case 20:
        return 17;
    case 21:
        return 22;
    case 22:
        return 18;
    case 23:
        return 21;
    case 25:
        return 23;
    default:
        return 23; // consider other case as * character
    }
}