all: main
main: main.cpp Index.cpp HeaderReader.cpp SequenceReader.cpp
	g++ main.cpp Index.cpp HeaderReader.cpp SequenceReader.cpp -o main
