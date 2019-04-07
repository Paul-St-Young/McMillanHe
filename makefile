NAME=mmh
CC=g++
OBJS=$(NAME).o
CFLAGS=-c -O3 -std=c++14

%.o : %.cpp
	$(CC) $(CFLAGS) $<

all: $(OBJS) main.cpp
	$(CC) -o exec $(OBJS) main.cpp

clean:
	rm $(OBJS) exec
