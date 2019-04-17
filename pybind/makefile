NAME=mmh
CXX=g++
OFLAGS=-O3
CFLAGS=$(OFLAGS) -shared -fPIC
OBJS=mmh.o pybind_mmh.o
LIB=$(NAME).so

%.o : %.cpp
	$(CXX) $(CFLAGS) -c $<

all: $(OBJS)
	$(CXX) $(CFLAGS) -o $(LIB) $(OBJS)

clean:
	rm $(OBJS) $(LIB)
