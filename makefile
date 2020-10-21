objects = main.o FCI.o CI.o File_string.o
CC = g++
CXXFLAGS = -std=c++11

app : $(objects)
	g++ -o app $(objects)

main.o : FCI.h File_string.h
FCI.o : FCI.h CI.h
CI.o : CI.h
File_string.o : File_string.h

.PHONY : clean
clean :
	-rm main $(objects)