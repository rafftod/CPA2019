#include "Sequence.h"

Sequence::Sequence()
{
    score = 0;
    id = -1;
}

Sequence::~Sequence() {}

int Sequence::get_id() const
    { return id; }

int Sequence::get_score() const
    { return score; }

void Sequence::set_id(int id)
    { this->id = id; }

void Sequence::set_score(int score)
    { this->score = score;}