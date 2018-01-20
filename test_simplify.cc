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
#include "edge.hpp"
#include "node.hpp"

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
//   As with the above list, this will reduce temporaries a lot.
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

struct segment
{
    double left, right;
    std::int32_t node;
    segment() : left{}, right{}, node{} {}
    segment(double l, double r, std::int32_t n)
        : left{ l }, right{ r }, node{ n }
    {
    }
};

void
sort_tables(std::vector<edge>& edge_table,
            const std::vector<node>& node_table) noexcept
{
    // Sort the edge table.  On PARENT birth times.
    // The sorting differs from msprime here. We
    // assume that birth times are recorded forward in
    // time rather than backwards.
    std::sort(edge_table.begin(), edge_table.end(),
              [&node_table](const edge& a, const edge& b) {
                  auto ga = node_table[a.parent].generation;
                  auto gb = node_table[b.parent].generation;
                  return ga > gb
                         || (ga == gb
                             && std::tie(a.parent, a.child, a.left)
                                    < std::tie(b.parent, b.child, b.left));
              });
    assert(std::is_sorted(
        edge_table.begin(), edge_table.end(),
        [&node_table](const edge& a, const edge& b) {
            auto ga = node_table[a.parent].generation;
            auto gb = node_table[b.parent].generation;
            return ga > gb || (ga == gb
                               && std::tie(a.parent, a.child, a.left)
                                      < std::tie(b.parent, b.child, b.left));
        }));
}

bool
sort_queue(bool added2Q, std::vector<segment>& Q) noexcept
{
    if (added2Q)
        {
            std::sort(Q.begin(), Q.end(),
                      [](const segment& a, const segment& b) {
                          return a.left > b.left;
                      });
        }
    return false;
}

std::vector<std::int32_t>
simplify(const std::vector<std::int32_t>& samples,
         std::vector<edge>& edge_table, std::vector<node>& node_table)
{
    std::vector<edge> Eo;
    std::vector<node> No;
    std::vector<std::vector<segment>> Ancestry(node_table.size());

	// Relates input node ids to output node ids
    std::vector<std::int32_t> idmap(node_table.size(), -1);

    // This plays the role of a min queue on segments, meaning
    // that it is always sorted such that Q.back() is the
    // smallest value, according to the lambda in sort_queue
    std::vector<segment> Q;
    // TODO: document a gotcha re: samples not being sorted w.r.to
    // index

	// We take our samples and add them to both the output
	// node list and initialize their ancestry with
	// a segment on [0,L).
    for (auto& s : samples)
        {
            No.push_back(make_node(No.size(), node_table[s].generation, 0));
            Ancestry[s].emplace_back(0, 1,
                                     static_cast<std::int32_t>(No.size() - 1));
            idmap[s] = static_cast<std::int32_t>(No.size() - 1);
        }

    auto edge_ptr = edge_table.begin();
    std::int32_t anode;
    double aleft, aright;
    std::vector<segment> X;
    bool added2Q = false;
    X.reserve(1000); //Arbitrary

	// At this point, our edges are sorted by birth
	// order of parents, from present to past.
	// We can now work our way up the pedigree.
	// This outer loop differs from how we describe it in the 
	// paper, but the strict sorting of edges means that this 
	// equivalent.
    while (edge_ptr < edge_table.end())
        {
            auto u = edge_ptr->parent;
            for (; edge_ptr < edge_table.end() && edge_ptr->parent == u;
                 ++edge_ptr)
                {
					//For each edge corresponding to this parent,
					//we look at all segments from the child.
					//If the two segments overlap, we add the minimal
					//overlap to our queue.
					//This is Step S3.
                    for (auto& seg : Ancestry[edge_ptr->child])
                        {
                            if (seg.right > edge_ptr->left
                                && edge_ptr->right > seg.left)
                                {
                                    Q.emplace_back(
                                        std::max(seg.left, edge_ptr->left),
                                        std::min(seg.right, edge_ptr->right),
                                        seg.node);
                                    added2Q = true;
                                }
                        }
                }
            added2Q = sort_queue(added2Q, Q);
            std::int32_t v = -1;
            while (!Q.empty())
				// Steps S4 through S8 of the algorithm.
                {
                    X.clear();
                    auto l = Q.back().left;
                    double r = 1.0;
                    while (!Q.empty() && Q.back().left == l)
					// This while loop is Step S4. This step 
					// adds to X all segments with left == l
					// and also finds the minimum right for 
					// all segs with left == l.
					// TODO: this can be done with reverse iteration,
					// but testing on 0.5e9 edges didn't seem to 
					// make it worthwhile.
                        {
                            X.emplace_back(Q.back().left, Q.back().right,
                                           Q.back().node);
                            r = std::min(r, Q.back().right);
                            Q.pop_back();
                        }
                    if (!Q.empty())
                        {
                            r = std::min(r, Q.back().left);
                        }
                    if (X.size() == 1)
                        {
                            aleft = X[0].left;
                            aright = X[0].right;
                            anode = X[0].node;
                            if (!Q.empty() && Q.back().left < X[0].right)
                                {
                                    aleft = X[0].left;
                                    aright = Q.back().left;
                                    anode = X[0].node;
                                    Q.emplace_back(Q.back().left, X[0].right,
                                                   X[0].node);
                                    added2Q = true;
                                }
                        }
                    else
                        {
                            if (v == -1)
                                {
									// Overlap/coalescence, and thus
									// a new node. Step S6.
                                    No.push_back(make_node(
                                        static_cast<std::int32_t>(No.size()),
                                        node_table[u].generation, 0));
                                    v = No.size() - 1;
                                    // update sample map
                                    idmap[u] = v;
                                }
                            aleft = l;
                            aright = r;
                            anode = v;
                            for (auto& x : X)
                                {
                                    Eo.push_back(make_edge(l, r, v, x.node));
                                    if (x.right > r)
                                        {
                                            x.left = r;
                                            Q.emplace_back(x.left, x.right,
                                                           x.node);
                                            added2Q = true;
                                        }
                                }
                        }
                    added2Q = sort_queue(added2Q, Q);
                    Ancestry[u].emplace_back(aleft, aright, anode);
                }
        }

    assert(std::count_if(idmap.begin(), idmap.end(),
                         [](const std::int32_t i) { return i != -1; })
           == No.size());
    std::size_t start = 0;

    // Now, we compact the edges,
    // which means removing redundant
    // info due to different edges
    // representing the same ancestry.
    std::vector<edge> E;
    E.swap(Eo);
    assert(Eo.empty());

    std::sort(E.begin(), E.end(), [](const edge& a, const edge& b) {
        return std::tie(a.parent, a.child, a.left, a.right)
               < std::tie(b.parent, b.child, b.left, b.right);
    });

    for (std::size_t j = 1; j < E.size(); ++j)
        {
            bool condition = E[j - 1].right != E[j].left
                             || E[j - 1].parent != E[j].parent
                             || E[j - 1].child != E[j].child;
            if (condition)
                {
                    Eo.push_back(make_edge(E[start].left, E[j - 1].right,
                                           E[j - 1].parent, E[j - 1].child));
                    start = j;
                }
        }
    auto j = E.size();
    Eo.push_back(make_edge(E[start].left, E[j - 1].right, E[j - 1].parent,
                           E[j - 1].child));
    edge_table.swap(Eo);
    node_table.swap(No);
    assert(std::is_sorted(
        edge_table.begin(), edge_table.end(),
        [&node_table](const edge& a, const edge& b) {
            auto ga = node_table[a.parent].generation;
            auto gb = node_table[b.parent].generation;
            return ga > gb || (ga == gb
                               && std::tie(a.parent, a.child, a.left)
                                      < std::tie(b.parent, b.child, b.left));
        }));
    return idmap;
}

int
main(int argc, char** argv)
{
    std::string nodefilename, edgefilename, nodeoutfile, edgeoutfile;
    nodefilename = std::string(argv[1]);
    edgefilename = std::string(argv[2]);
    nodeoutfile = std::string(argv[3]);
    edgeoutfile = std::string(argv[4]);
    std::int32_t N = std::atoi(argv[5]);
    std::vector<node> nodes;
    std::vector<edge> edges;

    std::ifstream in(nodefilename.c_str());
    std::int32_t a, b;
    double x, y;
    auto start = std::chrono::steady_clock::now();
    while (!in.eof())
        {
            in.read(reinterpret_cast<char*>(&a), sizeof(decltype(a)));
            if (a == -1)
                break;
            in.read(reinterpret_cast<char*>(&x), sizeof(decltype(x)));
            // in >> a >> x >> std::ws;
            nodes.push_back(make_node(a, x, 0));
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
            edges.push_back(make_edge(x, y, a, b));
        }
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cerr
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;
    std::cerr << nodes.size() << " nodes, " << edges.size() << " edges\n";
    std::vector<std::int32_t> samples;
    for (unsigned i = nodes.size() - 2 * N; i < nodes.size(); ++i)
        {
            samples.push_back(i);
        }
    start = std::chrono::steady_clock::now();
    sort_tables(edges, nodes);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cerr
        << "sort time = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;
    start = std::chrono::steady_clock::now();
    simplify(samples, edges, nodes);
    end = std::chrono::steady_clock::now();
    diff = end - start;
    std::cerr
        << "simplify time = "
        << std::chrono::duration_cast<std::chrono::milliseconds>(diff).count()
        << " ms" << std::endl;

    std::ofstream out(nodeoutfile.c_str());
    for (auto& n : nodes)
        {
            out << n.id << ' ' << n.generation << '\n';
        }
    out.close();
    out.open(edgeoutfile.c_str());
    for (auto& n : edges)
        {
            out << n.left << ' ' << n.right << ' ' << n.parent << ' '
                << n.child << '\n';
        }
    out.close();
}
