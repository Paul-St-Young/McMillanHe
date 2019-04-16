NAME=mmh
CC=g++
OBJS=main.o mmh.o loadtxt.o ezh5.o
#OFLAGS=-O0 -g
OFLAGS=-O3
CFLAGS=-c $(OFLAGS) -I. -std=c++14

%.o : %.cpp
	$(CC) $(CFLAGS) $<

all: $(OBJS)
	$(CC) -o exec $(OFLAGS) $(OBJS) -lhdf5
clean:
	rm *.o exec
