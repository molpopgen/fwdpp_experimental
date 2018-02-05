#include <cassert>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <vector>
#include <algorithm>
#include "table_collection.hpp"

using namespace std;
using namespace fwdpp::ancestry;

int
count_tips(int32_t parent, const vector<edge>& edges)
{
    int c = 0;
    auto start = lower_bound(
        edges.rbegin(), edges.rend(), parent,
        [](const edge& e, const int32_t& value) { return e.parent > value; });
    auto end = upper_bound(
        start, edges.rend(), parent,
        [](const int32_t& value, const edge& e) { return value > e.parent; });
	auto d = distance(start,end);
	if(d==0)
	{
		++c;
	}
	else
	{
		assert(d>0);
		for(;start<end;++start)
		{
			c += count_tips(start->child,edges);
		}
	}
	return c;
}

vector<pair<double, int>>
count_mutations(const table_collection& tables,
                const vector<double>& positions)
{
    vector<pair<double, int>> rv;

    for (auto& m : tables.mutation_table)
        {
            auto start = lower_bound(tables.edge_table.rbegin(),
                                     tables.edge_table.rend(), m.first,
                                     [](const edge& e, const int32_t& value) {
                                         return e.parent > value;
                                     });
            auto end = upper_bound(start, tables.edge_table.rend(), m.first,
                                   [](const int32_t& value, const edge& e) {
                                       return value > e.parent;
                                   });
            auto d = distance(start, end);
            if (d == 0)
                {
                    rv.emplace_back(positions[m.second], 1);
                }
            else
                {
                    int c = 0;
                    for (; start < end; ++start)
                        {
							c+=count_tips(start->child,tables.edge_table);
                        }
					rv.emplace_back(positions[m.second],c);
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
