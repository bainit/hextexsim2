CC := g++
MPICC := mpicxx
CFLAGS := --std=c++11 -g

TARGET := simulations
TARGET_MPI := simulations_mpi

TARGETFLAGS := -lgsl -lm -lgslcblas -fopenmp

SRC := $(wildcard src/*.cpp)
OBJ := $(patsubst %.cpp,%.o,$(SRC))

all : $(TARGET)
mpi : $(TARGET_MPI)

$(TARGET) : $(OBJ)
	$(CC) $(CFLAGS) -c simulations.cpp -o simulations.o $(TARGETFLAGS)
	$(CC) -o $@ $^ simulations.o $(TARGETFLAGS)

$(TARGET_MPI) : $(OBJ)
	$(MPICC) $(CFLAGS) -c simulations.cpp -o simulations.o -DMPI_ENABLED
	$(MPICC) -o $@ $^ simulations.o  $(TARGETFLAGS) 

$(SCR)/%.o : %.cpp
	$(CC) $(CFLAGS) -c $< 

clean:
	rm -rf *.o
	rm -rf src/*.o
	rm simulations

.PHONY: all clean mpi

