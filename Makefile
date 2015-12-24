CXX=g++
DEBUG=-DNDEBUG
#DEBUG=
CXXFLAGS=-std=c++11 -O2 -Wall -W $(DEBUG)

all: shared_ptr_test.o fwdpp_lite.o fwdpp_indexes.o
	$(CXX) $(CXXFLAGS) -o shared_ptr_test shared_ptr_test.o
	$(CXX) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o -lgsl -lgslcblas -ltcmalloc
	$(CXX) $(CXXFLAGS) -o fwdpp_indexes fwdpp_indexes.o -lgsl -lgslcblas
clean:
	rm -f *.o
