DEBUG=-DNDEBUG
DEBUG=
OPT=-O2
#OPT=
CXXFLAGS=-std=c++11 $(OPT) -Wall -W $(DEBUG) -g

all: fwdpp_lite.o fwdpp_lite_intrusive.o fwdpp_simplify.o split_breakpoints.o test_simplify.o fwdpp_dump_nodes_edges.o mutation_counting.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite_intrusive fwdpp_lite_intrusive.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_simplify fwdpp_simplify.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_dump_nodes_edges fwdpp_dump_nodes_edges.o  split_breakpoints.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o test_simplify test_simplify.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o mutation_counting mutation_counting.o

fwdpp_dump_nodes_edges.o: node.hpp edge.hpp msprime_algo.hpp
test_simplify.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp
fwdpp_simplify.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp



clean:
	rm -f *.o
