EXECUTABLE=DCT
FLAGS= -O3 -Wall -std=c++0x 
POSTFLAGS=-lrt -lm
CC=g++


parallel: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) -fopenmp $(EXECUTABLE).cpp $(POSTFLAGS)
all: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) $(EXECUTABLE).cpp $(POSTFLAGS)
clean:
	@rm -f $(EXECUTABLE)
