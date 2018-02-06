#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "table_collection.hpp"

using namespace std;
using namespace fwdpp::ancestry;

vector<edge>::const_reverse_iterator
get_parent_range(const int32_t parent, const vector<edge>& edges)
{
    auto start = lower_bound(
        edges.rbegin(), edges.rend(), parent,
        [](const edge& e, int32_t value) { return e.parent > value; });
    return start;
}

int
count_tips(int32_t parent, const double pos, const vector<edge>& edges)
{
    int c = 0;
    auto range = get_parent_range(parent, edges);
    if (range == edges.rend())
        {
            ++c;
        }
    else
        {
            for (; range < edges.rend() && range->parent == parent; ++range)
                {
                    if (pos >= range->left && pos < range->right)
                        c += count_tips(range->child, pos, edges);
                }
        }
    return c;
}

vector<pair<double, int>>
count_mutations(const table_collection& tables,
                const vector<double>& positions)
//Count the number of times each mutation
//occurs in an edge table.
//
//Works by relying on the fact that output edges
//are sorted by parent ID and the fact that we
//store edges as structs in a contiguous vector.
//
//The algorithm works by doing a binary search
//from the back for the node associated with a mutation.
//
//We rely on the following:
//1. If we cannot find a parent label with a mutation node label,
//   then the mutation is a singleton
//
//TODO: we can end the inner loop earlier.  The simplify
//   algorithm ensures that edges for a parent cannot overlap,
//   meaning that once we find what edge a mutation is in, we can
//   break.
//TODO: more code re-use.  The inner loop below is copied in count_tips
//   above
//TODO: make this work for back mutation
{
    vector<pair<double, int>> rv;

    for (auto& m : tables.mutation_table)
        {
            auto range = get_parent_range(m.first, tables.edge_table);
            if (range == tables.edge_table.rend())
                {
                    rv.emplace_back(positions[m.second], 1);
                }
            else
                {
                    int c = 0;
                    for (; range < tables.edge_table.rend()
                           && range->parent == m.first;
                         ++range)
                        {
                            if (positions[m.second] >= range->left
                                && positions[m.second] < range->right)
                                {
                                    c += count_tips(range->child,
                                                    positions[m.second],
                                                    tables.edge_table);
                                }
                        }
                    rv.emplace_back(positions[m.second], c);
                }
        }
    sort(rv.begin(), rv.end(),
         [](const pair<double, int>& a, const pair<double, int>& b) {
             return a.first < b.first;
         });

    return rv;
}

int
main(int argc, char** argv)
{
    table_collection tables;
    ifstream in("test_nodes.bin");
    int32_t n;
    in.read(reinterpret_cast<char*>(&n), sizeof(int32_t));
    for (int32_t i = 0; i < n; ++i)
        {
            double t;
            in.read(reinterpret_cast<char*>(&t), sizeof(double));
            tables.emplace_back_node(i, t, 0);
        }
    in.close();
    in.open("test_edges.bin");

    in.read(reinterpret_cast<char*>(&n), sizeof(int32_t));

    for (int32_t i = 0; i < n; ++i)
        {
            double l, r;
            int32_t p, c;
            in.read(reinterpret_cast<char*>(&l), sizeof(double));
            in.read(reinterpret_cast<char*>(&r), sizeof(double));
            in.read(reinterpret_cast<char*>(&p), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&c), sizeof(int32_t));
            tables.emplace_back_edge(l, r, p, c);
        }
    in.close();

    in.open("test_sites.bin");
    in.read(reinterpret_cast<char*>(&n), sizeof(int32_t));

    vector<double> positions;
    for (int32_t i = 0; i < n; ++i)
        {
            double pos;
            int32_t id;
            in.read(reinterpret_cast<char*>(&id), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&pos), sizeof(double));
            positions.emplace_back(pos);
            tables.mutation_table.emplace_back(id, i);
        }

    auto mcounts = count_mutations(tables, positions);
    for (auto& i : mcounts)
        {
            cout << i.first << ' ' << i.second << '\n';
        }
}
