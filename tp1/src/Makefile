EXECUTABLE=MN
FLAGS=-O3 -Wall
POSTFLAGS=-lrt -lm
CC=g++
all: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) $(EXECUTABLE).cpp $(POSTFLAGS)

clean:
	@rm -f $(EXECUTABLE)
