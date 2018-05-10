// This file is helping me get an implementation
// of the simplify algorithm from Kelleher et al.
// implemented in C++.  This implemenation
// is based on the set theoretic algorithm
// presented by Jerome Kelleher in the repo
// of that paper.
//
// Almost no intellectual input is from me
// (KRT). All I did was put the Python
// implementation in terms of c++11 code.
#include <cassert>
#include <tuple>
#include <vector>
#include <functional>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <chrono>
#include <fwdpp/sugar/popgenmut.hpp>
#include "edge.hpp"
#include "node.hpp"
#include "table_simplifier.hpp"
using namespace fwdpp::ancestry;

// The following needs to be dealt with before we'd be able to move
// the code into the context of a class to be tied to a simulation.

// All benchmarking is donw with output from dump_big.sh, which
// generates 80004000 nodes and 80997299 edges.

// TODO list based on simplifying large numbers of edges.
// This list all made measurable speedups, and should be re-introduced
// separately for the purposes of git history:
// 1. Replace priority_queue with vector that we sort as needed. DONE
// 2. Replace "alpha" variable with raw data types that we emplace_back. DONE
// 3. Replace all uses of push_back of alpha with emplace_back of raw types. DONE

// TODO explore
// 1. Give node and edge constructors so that we may emplace_back them, too.
//   As with the above list, this will reduce temporaries a lot. Experimenting
//   shows that for these types, emplacement is slower than push_back of a constructor call.
//   The constructor call is more idiomatic than a make_foo function, so we
//   stick with that.  This means we must copy to get data to msprime, but that
//   is probably what we want if we're going all-in on the fwdpp side.
// 2. Replace segment with a plain tuple.  It is a hidden/internal data type,
//   so we really don't need to see it.
// 3. Move X higher in scope and clear() it each iteration.  This will
//    re-use its RAM each time instead of reallocating. Done--this was worth
//    a 25% reduction in simplify time by itself!!!!
// 4. Try to see if reserve on No,Eo helps

// TODO issues
// 1. There's a trick to ancient samples.  One has to add them into output node
// list,
//   but take care not to re-add them.  I may have to steal a look at msprime
//   here...

// TODO API
// 1. separate sorting from simplifying. DONE

// TODO Homework
// 2. Profile with big input on dev server

// Some notes:
// 1. Sorting is much faster here than in msprime.  Simplifying small numbers
// of edges
//   is also faster, but huge numbers favors msprime as Jerome predicted.
//   Some, but
//   not all, of the difference is due to temporary object creation, and
//   refactoring
//   to allow emplace_back wherever possible will help.

// Any idea of reserving memory will be moot in practice.  When simplify
// becomes a member of a C++ class, all of the temporary containers
// can be private variables whose memory gets re-used during each bout
// of simplification.

int
main(int argc, char** argv)
{
    std::string nodefilename, edgefilename, nodeoutfile, edgeoutfile;
    nodefilename = std::string(argv[1]);
    edgefilename = std::string(argv[2]);
    nodeoutfile = std::string(argv[3]);
    edgeoutfile = std::string(argv[4]);
    std::int32_t N = std::atoi(argv[5]);

    std::ifstream in(nodefilename.c_str());
    std::int32_t a, b;
    double x, y;
	table_collection tables(1.0);
    auto start = std::chrono::steady_clock::now();
    while (!in.eof())
        {
            in.read(reinterpret_cast<char*>(&a), sizeof(decltype(a)));
            if (a == -1)
                break;
            in.read(reinterpret_cast<char*>(&x), sizeof(decltype(x)));
            // in >> a >> x >> std::ws;
			tables.emplace_back_node(a,0,x);
        }
    // for(auto & n : nodes)
    //{
    //    std::cout << n.id << ' ' << n.generation << '\n';
    //}
    in.close();
    in.open(edgefilename.c_str());
    while (!in.eof())
        {
            // in >> a >> b >> x >> y >> std::ws;
            in.read(reinterpret_cast<char*>(&a), sizeof(decltype(a)));
            if (a == -1)
                break;
            in.read(reinterpret_cast<char*>(&b), sizeof(decltype(b)));
            in.read(reinterpret_cast<char*>(&x), sizeof(decltype(x)));
            in.read(reinterpret_cast<char*>(&y), sizeof(decltype(y)));
			tables.emplace_back_edge(x,y,a,b);
        }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cerr
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;
    std::cerr << tables.node_table.size() << " nodes, " << tables.edge_table.size() << " edges\n";
    std::vector<std::int32_t> samples;
    for (unsigned i = tables.node_table.size() - 2 * N; i < tables.node_table.size(); ++i)
        {
            samples.push_back(i);
        }
    start = std::chrono::steady_clock::now();
	tables.sort_edges();
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cerr
        << "sort time = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;
	table_simplifier ancestry(1.0);
    start = std::chrono::steady_clock::now();
	std::vector<fwdpp::popgenmut> dummy;
    ancestry.simplify(tables,samples,dummy);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cerr
        << "simplify time = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;

    std::ofstream out(nodeoutfile.c_str());
    for (auto& n : tables.node_table)
        {
            out << n.id << ' ' << n.generation << '\n';
        }
    out.close();
    out.open(edgeoutfile.c_str());
    for (auto& n : tables.edge_table)
        {
            out << n.left << ' ' << n.right << ' ' << n.parent << ' '
                << n.child << '\n';
        }
    out.close();
}
