CXX=g++
CXXFLAGS=-std=c++11 -O2 -Wall -W

all: shared_ptr_test.o
	$(CXX) $(CXXFLAGS) -o shared_ptr_test shared_ptr_test.o
