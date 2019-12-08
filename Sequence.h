#ifndef SEQUENCE
#define SEQUENCE

class Sequence
{
private:
    int id;
    int score;
public:
    Sequence();
    ~Sequence();

    int get_id() const;
    int get_score() const;
    void set_id(int id);
    void set_score(int score);
};

#endif