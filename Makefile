OUTPUT_DIR = ./
CFLAGS= -g -Wall -O3 -lm -L.

all: ${OUTPUT_DIR}random_walk

RandomWalkPthread.o: RandomWalkPthread.cpp
	g++ $(CFLAGS) -c RandomWalkPthread.cpp

${OUTPUT_DIR}random_walk: RandomWalkPthread.o
	g++ $(CFLAGS) RandomWalkPthread.o -o ${OUTPUT_DIR}"random_walk" -lrt -lpthread

clean:
	rm *.o
	rm ${OUTPUT_DIR}random_walk