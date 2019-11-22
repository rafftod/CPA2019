#include <fstream>
#include <string>
#include "Index.h"
#include <map>

const std::map<int, char> residue_int_conversion_map = { // corresponding int and char, following conversion chart of BLAST format
    {1,'A'},{2,'B'},{3,'C'},{4,'D'},{5,'E'},{6,'F'},
    {7,'G'},{8,'H'},{9,'I'},{27,'J'},{10,'K'},{11,'L'},{12,'M'},
    {13,'N'},{26,'O'},{14,'P'},{15,'Q'},{16,'R'},{17,'S'},{18,'T'},
    {24,'U'},{19,'V'},{20,'W'},{21,'X'},{22,'Y'},{23,'Z'},{25,'*'},
};

class SequenceReader
{
private:
    std::string* sequences; // stores all the database sequences
    Index* sq_index; // database index
public:
    SequenceReader(Index* index);
    ~SequenceReader();

    void read_data(std::ifstream& database_sequence);
    std::string get_sequence(int index) const;
    int search_sequences(std::ifstream& query_protein);
    int exact_match(std::ifstream& database_sequence, std::ifstream& query_protein);
};
