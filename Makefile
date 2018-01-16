CXX=g++
DEBUG=-DNDEBUG
DEBUG=
CXXFLAGS=-std=c++11 -O2 -Wall -W $(DEBUG)

all: fwdpp_lite.o fwdpp_lite_intrusive.o fwdpp_simplify.o split_breakpoints.o test_simplify.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite_intrusive fwdpp_lite_intrusive.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_simplify fwdpp_simplify.o  split_breakpoints.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o test_simplify test_simplify.o -lgsl -lgslcblas $(LDFLAGS)

clean:
	rm -f *.o
