FAST= -fast
OMP= -fopenmp
CCINTEL=icpc
CCMPI=mpicxx
CCMPIINTEL=mpiicpc
MPIFAST=-0fast

FLAGS= -O3 
DEBUG= -g -p
FAST = -fast
CC = g++
SRCS = main.cpp mersenne.cpp  
TEST_SRCS = test.cpp mersenne.cpp

OBJS = $(SRCS)
TEST_OBJS = $(TEST_SRCS)

main: $(OBJS) 
	$(CC)  $(OBJS) $(FLAGS) -o nesi -w -std=c++11 

omp: $(OBJS) 
	$(CC)  $(OBJS) $(FLAGS) $(OMP) -o nesi_omp -w -std=c++11 

test: $(TEST_OBJS) 
	$(CC)  $(TEST_OBJS) $(FLAGS) -o test -w -std=c++11 

test_omp: $(TEST_OBJS) 
	$(CC)  $(TEST_OBJS) $(FLAGS) $(OMP) -o test -w -std=c++11 

mpi: $(OBJS)1
	$(CCMPIINTEL) $(OBJS) $(MPIFAST) -lm -w -o dsa

mpi_debug: $(OBJS)
	$(CCMPI)	$(OBJS)	$(DEBUG) -o dsa

intelfast: $(OBJS)
	$(CCINTEL)  $(OBJS) $(FAST) -o dsa

intelstatic: $(OBJS)
	$(CCINTEL)  $(OBJS) $(FAST) -ipo -o dsa

fast: $(OBJS)
	$(CC)  $(OBJS) $(FAST) -o dsa

debug: $(OBJS)
	$(CC)  $(OBJS) $(DEBUG) -w -o debug

%.cpp.o: %.cpp
	$(CC) $(FLAGS) $(DEBUG) -c -o $@ $< 
clean:
	rm *.cpp.o

