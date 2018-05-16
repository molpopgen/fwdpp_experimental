#include <iostream>
#include <fstream>
#include <cstdint>
#include <limits>
#include "table_collection.hpp"
#include "table_simplifier.hpp"

using namespace std;
using namespace fwdpp;
using namespace fwdpp::ancestry;

struct fake_mut
{
    double pos;
};

int
main(int argc, char** argv)
{
    table_collection tables(1.0);
    std::int32_t nrecords;
    double time, mtime = 0.0, mintime = numeric_limits<double>::max();
    ifstream in("decap_nodes.bin", std::ios::binary);
    in.read(reinterpret_cast<char*>(&nrecords), sizeof(std::int32_t));
    for (int i = 0; i < nrecords; ++i)
        {
            in.read(reinterpret_cast<char*>(&time), sizeof(double));
            tables.push_back_node(i, time, 0);
            mtime = max(mtime, time);
        }
    for (auto& n : tables.node_table)
        {
            n.generation -= mtime;
            n.generation *= -1.0;
            mintime = min(mintime, n.generation);
        }

    vector<int32_t> samples;
    for (int i = 0; i < tables.node_table.size(); ++i)
        {
            if (tables.node_table[i].generation == mtime)
                {
                    samples.push_back(i);
                }
        }
    std::cout << samples.size() << " samples found\n";

    in.close();

    in.open("decap_edges.bin");
    in.read(reinterpret_cast<char*>(&nrecords), sizeof(std::int32_t));
    int32_t p, c;
    double l, r;
    for (int i = 0; i < nrecords; ++i)
        {
            in.read(reinterpret_cast<char*>(&p), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&c), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&l), sizeof(double));
            in.read(reinterpret_cast<char*>(&r), sizeof(double));
            tables.push_back_edge(l, r, p, c);
        }
    in.close();

    vector<fake_mut> mutations;
    in.open("decap_mutations.bin", ios::binary);
    in.read(reinterpret_cast<char*>(&nrecords), sizeof(std::int32_t));
    std::cout << nrecords << "mutations\n";
    for (int i = 0; i < nrecords; ++i)
        {
            in.read(reinterpret_cast<char*>(&p), sizeof(int32_t));
            in.read(reinterpret_cast<char*>(&l), sizeof(double));
            mutations.emplace_back(fake_mut{ l });
            tables.mutation_table.emplace_back(mutation_record{ p, i });
        }

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0);
    auto res = simplifier.simplify(tables, samples, mutations);
    tables.mutation_table.erase(
        remove_if(tables.mutation_table.begin(), tables.mutation_table.end(),
                  [&res, &samples](const mutation_record& mr) {
                      return res.second[mr.key] == samples.size();
                  }),
        tables.mutation_table.end());
    std::cout << "mutation table size = " << tables.mutation_table.size()
              << '\n';
    ofstream out("cpp_counts.txt");
    for (auto mr : tables.mutation_table)
        {
            out << res.second[mr.key] << '\n';
        }
    out.close();
}
