CXX=g++
DEBUG=-DNDEBUG
DEBUG=
CXXFLAGS=-std=c++11 -O2 -Wall -W $(DEBUG)

all: fwdpp_lite.o fwdpp_lite_intrusive.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o  -lgsl -lgslcblas -ltcmalloc $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite_intrusive fwdpp_lite_intrusive.o  -lgsl -lgslcblas -ltcmalloc $(LDFLAGS)

clean:
	rm -f *.o
