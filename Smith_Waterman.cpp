#include "Smith_Waterman.h"
#include "SequenceReader.h"
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <boost/algorithm/string.hpp>
#define sw_max2(a,b) ((a)>(b)?(a):(b))
#define sw_max3(a,b,c) sw_max2(sw_max2(a,b),c)
#define sw_max(a,b,c,d) sw_max2(sw_max3(a,b,c),d)

int blosum_matrix[BLOSUM_SIZE][BLOSUM_SIZE];

Smith_Waterman::Smith_Waterman(int gap1, int gap2, int m, int nm){
    gap_penalty_open = gap1;
    gap_penalty_exp = gap2;
}

Smith_Waterman::Smith_Waterman(int gap_open_penalty, int gap_expansion_penalty, std::string blosum_path){
/**
 * Constructor.
 *
 * @param gap_open_penalty Gap open penalty. Default : 11.
 * @param gap_expansion_penalty Gap expansion penalty. Default : 1.
 * @param blosum_path Path to BLOSUM matrix. Default : BLOSUM62
 * @return Index of matching sequence in the database, or -1 if no match.
 */
    gap_penalty_open = gap_open_penalty;
    gap_penalty_exp = gap_expansion_penalty;
    build_BLOSUM(blosum_path);
}

Smith_Waterman::~Smith_Waterman(){}

int Smith_Waterman::compare(const uint8_t* & sequence1, const int* & sequence2, const int length1, const int length2) const
/**
 * Compares 2 sequences with Swith-Waterman algorithm and returns score.
 *
 * @param sequence1 Sequence to compare from the database.
 * @param sequence2 Query sequence to compare.
 * @param length1 Number of residues in sequence 1 + 1.
 * @param length2 Number of residues in sequence 2 + 1.
 * @return Alignment bitscore between the 2 sequences.
 */
{
    // initialize score matrix on the heap, because stack length limit is 1000
    
    int** matrix;
    matrix = new int*[length1];
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
    int maximum_on_rows[length1] = { 0 }; 
    int maximum_on_rows_index[length1] = { 0 }; 
    int maximum_on_columns[length2] = { 0 }; 
    int maximum_on_columns_index[length2] = { 0 };

    int score = 0;

    for(int i = 1; i < length1; i++)
    {
        for (int j = 1; j < length2; j++)
        {
            int a = matrix[i-1][j-1] + blosum_matrix[(int)sequence1[i-1]][sequence2[j-1]];

            int b = maximum_on_rows[i] - gap_penalty_open - gap_penalty_exp*abs(j-maximum_on_rows_index[i]);

            /* c is gap penalty value on column j */

            int c = maximum_on_columns[j] - gap_penalty_open - gap_penalty_exp*abs(i-maximum_on_columns_index[j]);
            /* calculation of matrix element */
            matrix[i][j] = sw_max(a,b,c,0);
            if(matrix[i][j] >= maximum_on_rows[i] - gap_penalty_exp*abs(j-maximum_on_rows_index[i]))
            {
                // we made a new maximum we can memoize
                maximum_on_rows[i] = matrix[i][j];
                maximum_on_rows_index[i] = j;
            }
            if(matrix[i][j] > maximum_on_columns[j] - gap_penalty_exp*abs(i-maximum_on_columns_index[j]))
            {
                // we made a new maximum we can memoize
                maximum_on_columns[j] = matrix[i][j];
                maximum_on_columns_index[j] = i;
            }
            if(matrix[i][j] > score)
            {
                score = matrix[i][j];
            }
        }
    }
    /* heap cleanup */
    for(int i = 0; i < length1; i++)
        delete [] matrix[i];
    delete [] matrix;
    int bitscore = (0.267*score + 3.34)/log(2);
    return bitscore;
}

int Smith_Waterman::compare2(const uint8_t* & sequence1, const int* & sequence2, int length1, int length2)
{
    //compare 2 sequences with the Smith-Waterman algorithm and returns score
    //sequence1 : sequence from the database
    //sequence2 : query sequence
    //length1 : number of residues in sequence1 + 1
    //length2 : number of residues in sequence2 + 1

    //initialize score matrix

    int** matrix; int** E; int** F;
    if(length1 >= 1000 || length2 >= 1000) // we can not allocate the matrix on the stack
    {
        matrix = new int*[length1];
        for(int i = 0; i < length1; ++i)
            matrix[i] = new int[length2];
        E =  new int*[length1];//length1 = number of rows
        for (int i = 0; i < length1; ++i)
        {
            E[i] = new int[length2];//length2 = number of columns
        }

        F =  new int*[length1];//length1 = number of rows
        for (int i = 0; i < length1; ++i)
        {
            F[i] = new int[length2];
        }
        
    }
    else
    {
        matrix[length1][length2];
        E[length1][length2];
        F[length1][length2];
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

    for (int j = 1; j < length2; ++j)
    {
        for (int i = 1; i < length1; ++i)
        {
            int a = matrix[i-1][j-1] + blosum_matrix[(int)sequence1[i-1]][sequence2[j-1]];
            //std::cout << "a : " << a << std::endl;
            E[i][j] = sw_max2(matrix[i][j-1]-gap_penalty_open - gap_penalty_exp,E[i][j-1] - gap_penalty_exp);
            int b = E[i][j];
            //std::cout << "b : " << b << std::endl;
            F[i][j]= sw_max2(matrix[i-1][j]- gap_penalty_open -gap_penalty_exp, F[i-1][j] -gap_penalty_exp);
            int c = F[i][j];
            //std::cout << "c : " << c << std::endl;
            matrix[i][j] = sw_max(a,b,c,0);
            //std::cout << "max : " << matrix[i][j] << std::endl;
            if (matrix[i][j] > score)
            {
                score = matrix[i][j];
            }    
        } 
    }
    for(int i = 0; i < length1; ++i)
        delete [] matrix[i];
    delete [] matrix;
    for(int i = 0; i < length1; ++i)
        delete [] E[i];
    delete E;
    for(int i = 0; i < length1; ++i)
        delete [] F[i];
    delete F;
    int bitscore = (0.267*score + 3.34)/log(2);
    return bitscore;
}

int Smith_Waterman::build_BLOSUM(const std::string path)
/**
 * Reads BLOSUM file from given path and fills BLOSUM matrix.
 * 
 * This function supposes that all BLOSUM matrices
 * rows and columns are in the same order as BLOSUM62
 * and have the same format.
 *
 * @param path Path to BLOSUM matrix file.
 * @return Alignment bitscore between the 2 sequences.
 */
{
    /* Builds BLOSUM matrix from given path */
    std::ifstream blosum_file(path);
    std::string s;
    std::vector<std::string> line_vector;
    for(int i = 0; i < 7; i++) {std::getline(blosum_file,s);} // ignore first 6 lines and use 7th to build conversion map
    int sign = 1; char current_char;
    int i = 0; int j = 0;
    while(getline(blosum_file,s))
    {
        boost::split(line_vector, s, boost::is_any_of(" "));
        for(std::string sub: line_vector)
        {
            if(!isdigit(sub[sub.size()-1])) { continue; }
            int digit = strtol(sub.c_str(), NULL, 10);
                if(j > BLOSUM_SIZE-1-2) // we are at end of line
                {
                    i++; // go to next line
                    j = 0; // reset column counter
                }
                blosum_matrix[blosum_column_to_blast_number(i)][blosum_column_to_blast_number(j++)] = digit;
        }
    }
    return 0;
}

int Smith_Waterman::blosum_column_to_blast_number(int j) const
{
    switch (j)
    {
    case 0:
        return 1;
    case 20:
        return 2;
    case 4:
        return 3;
    case 3:
        return 4;
    case 6:
        return 5;
    case 13:
        return 6;
    case 7:
        return 7;
    case 8:
        return 8;
    case 9:
        return 9;
    case 11:
        return 10;
    case 10:
        return 11;
    case 12:
        return 12;
    case 2:
        return 13;
    case 14:
        return 14;
    case 5:
        return 15;
    case 1:
        return 16;
    case 15:
        return 17;
    case 16:
        return 18;
    case 19:
        return 19;
    case 17:
        return 20;
    case 22:
        return 21;
    case 18:
        return 22;
    case 21:
        return 23;
    case 23:
        return 25;
    default:
        return 25; // consider other case as * character
    }
}