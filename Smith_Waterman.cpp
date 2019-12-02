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

void Smith_Waterman::compare(uint8_t* sequence1, std::vector<int>& sequence2, int length1, int length2)
{
    //compare 2 sequences with the Smith-Waterman algorithm
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

    //first row and column are set at 0

    for (int i = 0; i < length1; i++)
    {
        matrix[i][0]=0;
    }

    for (int j = 1; j < length2; j++)
    {
        matrix[0][j]=0;
    }

    //scoring matrix

    for(int i = 1; i < length1; i++)
    {
        for (int j = 1; j < length2; j++)
        {
            int a = matrix[i-1][j-1] + blosum_matrix[sequence1[i]][sequence2[j]];
            int b = max_column(i,j,matrix);
            int c = max_row(i,j,matrix);
            matrix[i][j] = find_max(a,b,c);
        }
    }
}

int Smith_Waterman::find_max(int a, int b, int c)
{
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

int Smith_Waterman::max_column(int i, int j, int** matrix)
//returns the maximum score on a column
{
    int max = matrix[i][j-1] - gap_penalty_exp - gap_penalty_open;
    for (int k = 2; k <= i; k++)
    {
        int score = matrix[i][j - k] - k*gap_penalty_exp - gap_penalty_open;
        if (score > max)
        {
            max = score;
        }  
    }
    return max;
}

int Smith_Waterman::max_row(int i, int j, int** matrix)
//returns the maximum score on a line
{
    int max = matrix[i-1][j] - gap_penalty_exp - gap_penalty_open;
    for (int k = 2; k <= i; k++)
    {
        int score = matrix[i - k][j] - k*gap_penalty_exp - gap_penalty_open;
        if (score > max)
        {
            max = score;
        }  
    }
    return max;
}