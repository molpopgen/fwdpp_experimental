DEBUG=-DNDEBUG
DEBUG=
OPT=-O2
#OPT=
CXXFLAGS=-std=c++11 $(OPT) -Wall -W $(DEBUG) -g

all: fwdpp_lite.o fwdpp_lite_intrusive.o fwdpp_simplify.o split_breakpoints.o test_simplify.o fwdpp_dump_nodes_edges.o test_mut_simplify.o \
	process_decap.o fwdpp_simplify_gc_interval.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite fwdpp_lite.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_lite_intrusive fwdpp_lite_intrusive.o  -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_simplify fwdpp_simplify.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_simplify_gc_interval fwdpp_simplify_gc_interval.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o fwdpp_dump_nodes_edges fwdpp_dump_nodes_edges.o  split_breakpoints.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o test_simplify test_simplify.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o test_mut_simplify test_mut_simplify.o -lgsl -lgslcblas $(LDFLAGS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o process_decap process_decap.o -lgsl -lgslcblas $(LDFLAGS)


fwdpp_dump_nodes_edges.o: node.hpp edge.hpp msprime_algo.hpp
test_simplify.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp
fwdpp_simplify.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp variant_filler.hpp data_matrix_generator.hpp indexed_edge.hpp
fwdpp_simplify_gc_interval.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp
test_mut_simplify.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp
process_decap.o: node.hpp edge.hpp table_simplifier.hpp table_collection.hpp msprime_algo.hpp



clean:
	rm -f *.o
