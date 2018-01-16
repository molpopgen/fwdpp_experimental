// This file is helping me get an implementation
// of the simplify algorithm from Ralph et al.
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
#include <queue>
#include <algorithm>
#include <fstream>
#include <iostream>
#include "edge.hpp"
#include "node.hpp"

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
simplify(const std::vector<std::int32_t>& samples,
         std::vector<edge>& edge_table, std::vector<node>& node_table)
{
    // Don't do this, as we're testing w/data input from coalescent
    // reverse time
    // auto maxtime = node_table.back().generation;
    // for (auto& n : node_table)
    //    {
    //        n.generation -= maxtime;
    //        // Note: leads to 0 being -0.  SHOULDFIX
    //        n.generation *= -1.0;
    //    }

    // Sort the edge table.  On PARENT birth times.
    std::sort(edge_table.begin(), edge_table.end(),
              [&node_table](const edge& a, const edge& b) {
                  return std::tie(node_table[a.parent].generation, a.parent,
                                  a.child, a.left)
                         < std::tie(node_table[b.parent].generation, b.parent,
                                    b.child, b.left);
              });
    // std::cerr << "sorted: "
    //          << std::is_sorted(
    //                 edge_table.begin(), edge_table.end(),
    //                 [&node_table](const edge& a, const edge& b) {
    //                     return std::tie(node_table[a.parent].generation,
    //                                     a.parent)
    //                            < std::tie(node_table[b.parent].generation,
    //                                       b.parent);
    //                 }) << '\n';

    std::vector<edge> Eo;
    std::vector<node> No;
    std::vector<std::vector<segment>> Ancestry(node_table.size());
    // The algorithm using a min queue.  The default C++ queue
    // is a max queue.  Thus, we must use > rather than <
    // to generate a min queue;
    const auto segment_sorter_q
        = [](const segment& a, const segment& b) { return a.left > b.left; };
    std::priority_queue<segment, std::vector<segment>,
                        decltype(segment_sorter_q)>
        Q(segment_sorter_q);

    for (auto& s : samples)
        {
            No.push_back(make_node(No.size(), node_table[s].generation, 0));
            Ancestry[s].push_back(
                segment(0, 1, static_cast<std::int32_t>(No.size() - 1)));
        }
    // for(auto &s:samples)
    //{
    //    std::cout << Ancestry[s][0].left << ' ' << Ancestry[s][0].right << '
    //    ' << Ancestry[s][0].node << '\n';
    //}

    // auto last_edge = edge_table.begin();
    segment alpha;
    for (std::int32_t u = 0; u < static_cast<std::int32_t>(node_table.size());
         ++u)
        {
            // auto f
            //    = std::find_if(edge_table.begin(), edge_table.end(),
            //                   [u](const edge& e) { return e.parent == u; });

            // Edges are sorted according to (parent time, parent),
            // meaning was can find the first element containing
            // parent u in log time.
            auto f = std::lower_bound(
                edge_table.begin(), edge_table.end(), u,
                [&node_table](const edge& e, const std::int32_t value) {
                    return std::tie(e.parent, node_table[e.parent].generation)
                           < std::tie(value, node_table[value].generation);
                });
            auto l = std::find_if(f, edge_table.end(), [u](const edge& e) {
                return e.parent != u;
            });
            auto d1 = std::distance(edge_table.begin(), f),
                 d2 = std::distance(edge_table.begin(), l);
            for (; f < l; ++f)
                {
                    assert(f->parent == u);
                    for (auto& seg : Ancestry[f->child])
                        {
                            // if (f->child < 3)
                            //    {
                            //        std::cout
                            //            << "Ancestry of child: " << f->child
                            //            << ' ' << f->left << ' ' << f->right
                            //            << ' ' << seg.left << ' ' <<
                            //            seg.right
                            //            << '\n';
                            //    }
                            if (seg.right > f->left && f->right > seg.left)
                                {
                                    Q.emplace(std::max(seg.left, f->left),
                                              std::min(seg.right, f->right),
                                              seg.node);
                                    // std::cout << "here: " << Q.size() <<
                                    // '\n';
                                }
                        }
                }
            auto d3 = std::distance(edge_table.begin(), f);
            auto check = std::find_if(f, edge_table.end(), [u](const edge& e) {
                return e.parent == u;
            });
            if (check < edge_table.end())
                {
                    auto d4 = std::distance(edge_table.begin(), check);
                    std::cout << u << ' ' << d1 << ' ' << d2 << ' ' << d3
                              << ' ' << d4 << '\n';
                    std::cout << (edge_table.begin() + d1)->parent << ' '
                              << (edge_table.begin() + d2)->parent << ' '
                              << (edge_table.begin() + d3)->parent << ' '
                              << (edge_table.begin() + d4)->parent << '\n';
                    assert(false);
                }
            // if (!Q.empty())
            //    std::cout << "Qsize: " << u << ' ' << Q.size() << '\n';

            std::int32_t v = -1;
            while (!Q.empty())
                {
                    auto l = Q.top().left;
                    double r = 1.0;
                    std::vector<segment> X;
                    while (!Q.empty() && Q.top().left == l)
                        {
                            // Can be done w/fewer lines of code.
                            auto seg = Q.top();
                            Q.pop();
                            r = std::min(r, seg.right);
                            X.push_back(std::move(seg));
                        }
                    if (!Q.empty())
                        {
                            r = std::min(r, Q.top().left);
                        }
                    if (X.size() == 1)
                        {
                            alpha = X[0];
                            auto x = X[0];
                            if (!Q.empty() && Q.top().left < x.right)
                                {
                                    alpha = segment(x.left, Q.top().left,
                                                    x.node);
                                    x.left = Q.top().left;
                                    Q.push(x);
                                }
                        }
                    else
                        {
                            if (v == -1)
                                {
                                    No.push_back(make_node(
                                        static_cast<std::int32_t>(No.size()),
                                        node_table[u].generation, 0));
                                    v = No.size() - 1;
                                }
                            alpha = segment(l, r, v);
                            for (auto& x : X)
                                {
                                    Eo.push_back(make_edge(l, r, v, x.node));
                                    if (x.right > r)
                                        {
                                            x.left = r;
                                            Q.emplace(x);
                                        }
                                }
                        }
                    Ancestry[u].push_back(alpha);
                }
        }

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
}

int
main(int argc, char** argv)
{
    std::vector<node> nodes;
    std::vector<edge> edges;

    std::ifstream in("test_nodes.txt");
    std::int32_t a, b;
    double x, y;
    while (!in.eof())
        {
            in >> a >> x >> std::ws;
            nodes.push_back(make_node(a, x, 0));
        }
    // for(auto & n : nodes)
    //{
    //    std::cout << n.id << ' ' << n.generation << '\n';
    //}
    in.close();
    in.open("test_edges.txt");
    while (!in.eof())
        {
            in >> a >> b >> x >> y >> std::ws;
            edges.push_back(make_edge(x, y, a, b));
        }
    simplify({ 0, 1, 2, 19, 33, 11, 12 }, edges, nodes);
    std::cout << "nodes:\n";
    for (auto& n : nodes)
        {
            std::cout << n.id << ' ' << n.generation << '\n';
        }
    std::cout << "edges:\n";
    for (auto& n : edges)
        {
            std::cout << n.left << ' ' << n.right << ' ' << n.parent << ' '
                      << n.child << '\n';
        }
}
