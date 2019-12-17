CPP = g++
CPPFLAGS = -Ofast -std=c++11
CLSFLAGS = -c -o $@
LDFLAGS=

all: main clean
main: Index.o Header.o SequenceReader.o Smith_Waterman.o main.o
	$(CPP) $(CPPFLAGS) -lpthread -o $@ $^

Index.o: Index.cpp
	$(CPP) $(CPPFLAGS) $(CLSFLAGS) Index.cpp

Header.o: Header.cpp
	$(CPP) $(CPPFLAGS) $(CLSFLAGS) Header.cpp

SequenceReader.o: SequenceReader.cpp
	$(CPP) $(CPPFLAGS) $(CLSFLAGS) SequenceReader.cpp

Smith_Waterman.o: Smith_Waterman.cpp
	$(CPP) $(CPPFLAGS) $(CLSFLAGS) Smith_Waterman.cpp

main.o: main.cpp
	$(CPP) $(CPPFLAGS) $(CLSFLAGS) main.cpp

clean:
	rm *.o
