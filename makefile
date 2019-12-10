all: main
main: main.cpp Index.cpp Header.cpp SequenceReader.cpp Smith_Waterman.cpp Sequence.cpp
	g++ -o3 -std=c++11 main.cpp Index.cpp Header.cpp SequenceReader.cpp Smith_Waterman.cpp Sequence.cpp -o main
