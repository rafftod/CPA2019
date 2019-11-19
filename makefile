all: main
main: main.cpp Index.cpp Header.cpp SequenceReader.cpp
	g++ main.cpp Index.cpp Header.cpp SequenceReader.cpp -o main
