#include <fstream>
#include <string>
#include <vector>
#include <map>

const std::map<int, char> residue_int_conversion_map = {
    {1,'A'},{2,'B'},{3,'C'},{4,'D'},{5,'E'},{6,'F'},
    {7,'G'},{8,'H'},{9,'I'},{27,'J'},{10,'K'},{11,'L'},{12,'M'},
    {13,'N'},{26,'O'},{14,'P'},{15,'Q'},{16,'R'},{17,'S'},{18,'T'},
    {24,'U'},{19,'V'},{20,'W'},{21,'X'},{22,'Y'},{23,'Z'},{25,'*'},
};

class SequenceReader
{
private:
    vector<string> sequences;
public:
    SequenceReader();
    ~SequenceReader();

    std::vector<std::string> read_data(std::ifstream database_sequence);
    std::string get_sequence(int index) const;
};
