NAME=mmh
CC=g++
OBJS=mmh.o loadtxt.o rng.o ezh5.o
#OFLAGS=-O0 -g
OFLAGS=-O3
CFLAGS=-c $(OFLAGS) -std=c++14

%.o : %.cpp
	$(CC) $(CFLAGS) $<

all: $(OBJS) main.cpp
	$(CC) -o exec $(OFLAGS) $(OBJS) main.cpp -lhdf5

clean:
	rm *.o exec
