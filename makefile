all: main
main: main.cpp Index.cpp
	g++ main.cpp Index.cpp -o main
