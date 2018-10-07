// Unit tests for mut simplification
#include <iostream>
#include "table_collection.hpp"
#include "table_simplifier.hpp"

using namespace std;
using namespace fwdpp;
using namespace fwdpp::ts;

struct fake_mut
// We need to mock minimal API needed for a fwdpp
// mutation
{
    double pos;
};

void
test1()
//one mutation being passed onto a "complete"
//descendant
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 1., 0, 1);
    tables.emplace_back_edge(0., 1., 0, 2);

    // We measure time forwards, not backwards
    tables.push_back_node(
        0.0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(1.0, 0);
    tables.push_back_node(1.0, 0);

    vector<fake_mut> mutations(1, fake_mut{ 0.5 }); // one mutation at pos 0.5
    tables.mutation_table.emplace_back(
        mutation_record{ 0, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for (auto &s : samples)
        {
            s = res[s];
        }
    tables.build_indexes();
    tables.count_mutations(mutations, samples, mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 1);

    //the mutation must remap to xx.first[2]
    assert(tables.mutation_table[0].node == res[2]);
}

void
test2()
// Add an extra node that will be simplified out.
// The mutation is not in the edge that gets
// transmitted from this node
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 1., 0, 1);
    tables.emplace_back_edge(0., 1., 0, 2);
    tables.emplace_back_edge(0., 0.4, 3, 0); //this one will have the mutations
    // We measure time forwards, not backwards
    tables.push_back_node(
        0.0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(1.0, 0);
    tables.push_back_node(1.0, 0);
    tables.push_back_node(0., 0); // Node that has the mut

    vector<fake_mut> mutations(1, fake_mut{ 0.5 }); // one mutation at pos 0.5
    tables.mutation_table.emplace_back(
        mutation_record{ 3, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for (auto &s : samples)
        {
            s = res[s];
        }
    tables.build_indexes();
    tables.count_mutations(mutations, samples, mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 0);
    assert(tables.mutation_table.empty());
}
void
test3()
// This is test1, but with mutation at 0.4,
// putting it on edge of 1/2-open interval
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 1., 0, 1);
    tables.emplace_back_edge(0., 1., 0, 2);

    // We measure time forwards, not backwards
    tables.push_back_node(
        0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(0, 0);
    tables.push_back_node(0, 0);

    vector<fake_mut> mutations(1, fake_mut{ 0.4 }); // one mutation at pos 0.4
    tables.mutation_table.emplace_back(
        mutation_record{ 0, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for (auto &s : samples)
        {
            s = res[s];
        }
    tables.build_indexes();
    tables.count_mutations(mutations, samples, mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 1);
    //the mutation must remap to xx.first[2]
    assert(tables.mutation_table[0].node == res[2]);
}

void
test4()
//This is test2, but with the mutation
//position == rightmost position of
//most ancient edge
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 1., 0, 1);
    tables.emplace_back_edge(0., 1., 0, 2);
    tables.emplace_back_edge(0., 0.4, 3, 0); //this one will have the mutations
    // We measure time forwards, not backwards
    tables.push_back_node(
        0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(0, 0);
    tables.push_back_node(0, 0);
    tables.push_back_node(0., 0); // Node that has the mut

    vector<fake_mut> mutations(1, fake_mut{ 0.4 }); // one mutation at pos 0.4
    tables.mutation_table.emplace_back(
        mutation_record{ 3, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for (auto &s : samples)
        {
            s = res[s];
        }
    tables.build_indexes();
    tables.count_mutations(mutations, samples, mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 0);
    assert(tables.mutation_table.empty());
}

void
test5()
//Test 4, but we put a mutation at the very right
//of the most ancient node.
//To replicate this one in msprime, you gotta
//pass sequence_length = 1.0 to simplify
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 1.0, 0, 1);
    tables.emplace_back_edge(0., 1.0, 0, 2);
    tables.emplace_back_edge(0., 0.4, 3, 0); //this one will have the mutations
    tables.emplace_back_edge(0.6, 0.8, 3,
                             0); //this one will have the mutations
    // We measure time forwards, not backwards
    tables.push_back_node(
        0.0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(1.0, 0);
    tables.push_back_node(1.0, 0);
    tables.push_back_node(0., 0); // Node that has the mut

    vector<fake_mut> mutations(1, fake_mut{ 0.6 }); // one mutation at pos 0.6
    tables.mutation_table.emplace_back(
        mutation_record{ 3, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for(auto&s:samples){s=res[s];}
    tables.build_indexes();
    tables.count_mutations(mutations,samples,mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 2);
    assert(tables.mutation_table[0].node == res[0]);
}

// The above tests all have gaps, and seem to cover relevant
// cases.  So, let's test cases where an ancestor passes on
// [a,b)[b,c) to descendants and the mutation is at position b.
// Ideally, we'd also integrate this with the other machinery
// for adding edges based on breakpoints, etc., but that's not easy
// right now

void
test6()
//Test 4, but we put a mutation at the very right
//of the most ancient node.
//To replicate this one in msprime, you gotta
//pass sequence_length = 1.0 to simplify
{
    table_collection tables(1.); //max length = 1.0
    tables.emplace_back_edge(0., 0.4, 0, 1);
    tables.emplace_back_edge(0.6, 0.8, 0, 1);
    tables.emplace_back_edge(0., 1.0, 0, 2);
    tables.emplace_back_edge(0., 0.4, 3, 0); //this one will have the mutations
    tables.emplace_back_edge(0.4, 0.8, 3,
                             0); //this one will have the mutations
    // We measure time forwards, not backwards
    tables.push_back_node(
        0.0, 0); // node id, time, pop.  Yes, the ID is redundant...
    tables.push_back_node(1.0, 0);
    tables.push_back_node(1.0, 0);
    tables.push_back_node(0., 0); // Node that has the mut

    vector<fake_mut> mutations(1, fake_mut{ 0.4 }); // one mutation at pos 0.4
    tables.mutation_table.emplace_back(
        mutation_record{ 3, 0 }); //node, index of mutation in mutations

    tables.sort_tables(mutations);
    table_simplifier simplifier(1.0); // max length
    std::vector<std::int32_t> samples{ 1, 2 };
    std::vector<std::uint32_t> mcounts;
    auto res = simplifier.simplify(tables, samples, mutations);
    for(auto&s:samples){s=res[s];}
    tables.build_indexes();
    tables.count_mutations(mutations,samples,mcounts);

    // res is the node remapping,
    // and mcounts contains mutation counts
    assert(mcounts.size() == 1);
    assert(mcounts[0] == 1);
    assert(tables.mutation_table[0].node == res[2]);
}

int
main(int argc, char **argv)
{
    test1();
    test2();
    test3();
    test4();
    test5();
    test6();
}
