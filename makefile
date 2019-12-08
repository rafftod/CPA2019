all: main
main: main.cpp Index.cpp Header.cpp SequenceReader.cpp Smith_Waterman.cpp Sequence.cpp
	g++ main.cpp Index.cpp Header.cpp SequenceReader.cpp Smith_Waterman.cpp Sequence.cpp -o main
