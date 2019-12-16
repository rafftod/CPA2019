#ifndef SEQ_READER
#define SEQ_READER
#include <fstream>
#include <string>
#include "Index.h"
#include <map>
#include <vector>

/* all (ordered) maps are binary search trees, which makes searching into it a O(log n) */

const std::map<int, char> residue_int_conversion_map = { // corresponding int and char, following conversion chart of BLAST format
    {1,'A'},{2,'B'},{3,'C'},{4,'D'},{5,'E'},{6,'F'},
    {7,'G'},{8,'H'},{9,'I'},{27,'J'},{10,'K'},{11,'L'},{12,'M'},
    {13,'N'},{26,'O'},{14,'P'},{15,'Q'},{16,'R'},{17,'S'},{18,'T'},
    {24,'U'},{19,'V'},{20,'W'},{21,'X'},{22,'Y'},{23,'Z'},{25,'*'},
};

const std::map<char, int> residue_char_conversion_map = {
    {'A',1},{'B',2},{'C',3},{'D',4},{'E',5},{'F',6},
    {'G',7},{'H',8},{'I',9},{'J',27},{'K',10},{'L',11},{'M',12},
    {'N',13},{'O',26},{'P',14},{'Q',15},{'R',16},{'S',17},{'T',18},
    {'U',24},{'V',19},{'W',20},{'X',21},{'Y',22},{'Z',23},{'*',25},
};

class SequenceReader
{
private:
    std::string* sequences = NULL; // stores all the database sequences
    Index* sq_index = NULL; // database index
    uint8_t* sq_buffer = NULL;
    int query_size;
public:
    SequenceReader(Index* index);
    SequenceReader(Index* index, std::ifstream& database_sequence);
    ~SequenceReader();

    void read_data(std::ifstream& database_sequence);
    uint8_t* get_sequence(int index) const;
    int get_sequence_length(int index) const;
    int exact_match(std::ifstream& database_sequence, std::ifstream& query_protein);
    std::vector<int> convert_query_sequence(std::ifstream& query_protein);
    int get_query_size() const ;
};

#endif