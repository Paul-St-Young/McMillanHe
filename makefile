NAME=mmh
CC=g++
OBJS=mmh.o loadtxt.o rng.o
#OFLAGS=-O0 -g
OFLAGS=-O3
CFLAGS=-c $(OFLAGS) -std=c++14

%.o : %.cpp
	$(CC) $(CFLAGS) $<

all: $(OBJS) main.cpp
	$(CC) -o exec $(OFLAGS) $(OBJS) main.cpp

clean:
	rm *.o exec
