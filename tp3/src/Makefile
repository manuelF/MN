EXECUTABLE=OCR
FLAGS=-O3 -std=c++0x
POSTFLAGS=-lrt -lm
CC=g++

#Este build hace para ir haciendo una iteracion a la vez de la matriz de entrenamiento
pasos: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) -fopenmp $(EXECUTABLE).cpp $(POSTFLAGS)
#Este build genera la matriz de entrenamiento de una y luego calcula
entero: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) -fopenmp -DPRECALC $(EXECUTABLE).cpp $(POSTFLAGS)
#Este build hace para ir haciendo una iteracion a la vez de la matriz de entrenamiento (serial)
pasosserial: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) $(EXECUTABLE).cpp $(POSTFLAGS)
#Este build genera la matriz de entrenamiento de una y luego calcula (serial)
enteroserial: clean
	$(CC) -o $(EXECUTABLE) $(FLAGS) -DPRECALC $(EXECUTABLE).cpp $(POSTFLAGS)
clean:
	@rm -f $(EXECUTABLE)
