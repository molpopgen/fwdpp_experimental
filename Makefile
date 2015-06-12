CXX=g++
CXXFLAGS=-std=c++11 -O2 -Wall -W -DNDEBUG

all: shared_ptr_test.o fwdpp_lite.o
	$(CXX) $(CXXFLAGS) -o shared_ptr_test shared_ptr_test.o
	$(CXX) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o -lgsl -lgslcblas -ltcmalloc

clean:
	rm -f *.o
