#ifndef FWDPP_MSPRIME_ALGO_HPP__
#define FWDPP_MSPRIME_ALGO_HPP__

/// Implements Algorithms T and L from
/// Kelleher et al., 2016, aka "the msprime
/// paper", DOI: 10.1371/journal/pcbi.1004842
/// This implementation is experimental, and
/// is intended for integration into a future
/// release of fwdpp.
/// The most likely fate is that this code
/// end up as member functions of table_collection
/// and/or ancestry_tracker, so that we can
/// re-use memory allocated for I, O, and pi,
/// as we will be calling these functions a LOT.
/// LICENSE: GPL3
/// Author: Kevin Thornton
/// Thanks: Jerome Kelleher

#include <algorithm>
#include <vector>
#include <cstdint>
#include <tuple>
#include <limits>
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ancestry
    {

        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
        fill_I_O(const table_collection& tables)
        // Fill I and O.  This is not described in Algorithm T,
        // but can be sleuthed from msprime/lib/tables.c
        {
            std::vector<std::size_t> I(tables.edge_table.size(), 0),
                O(tables.edge_table.size(), 0);
            using edge_iterator = std::vector<edge>::const_iterator;
            std::vector<edge_iterator> edge_pointers;
            edge_pointers.reserve(tables.edge_table.size());
            for (auto i = std::begin(tables.edge_table);
                 i < std::end(tables.edge_table); ++i)
                {
                    edge_pointers.push_back(i);
                }
            // To fill I, sort the pointers by increasing left
            // and increasing parent time (moving from present back into past).
            // Because, unlike msprime, we are not recording node times
            // backwards into the past, our sort is simple.
            std::sort(
                std::begin(edge_pointers), std::end(edge_pointers),
                [&tables](const edge_iterator i, const edge_iterator j) {
                    return std::tie(i->left,
                                    tables.node_table[i->parent].generation)
                           < std::tie(j->left,
                                      tables.node_table[j->parent].generation);
                });

            for (std::size_t i = 0; i < I.size(); ++i)
                {
                    I[i] = static_cast<std::size_t>(std::distance(
                        std::begin(tables.edge_table), edge_pointers[i]));
                }

            // To fill O, sort by right and decreasing parent time, again
            // moving into the past.  Our one trick is to sort on the
            // reversed node times.
            std::sort(std::begin(edge_pointers), std::end(edge_pointers),
                      [&tables](const edge_iterator i, const edge_iterator j) {
                          auto ig = -tables.node_table[i->parent].generation;
                          auto jg = -tables.node_table[j->parent].generation;

                          return std::tie(i->right, ig)
                                 < std::tie(j->right, jg);
                      });
            for (std::size_t i = 0; i < O.size(); ++i)
                {
                    O[i] = static_cast<std::size_t>(std::distance(
                        std::begin(tables.edge_table), edge_pointers[i]));
                }
            return std::make_pair(std::move(I), std::move(O));
        }

        void
        algorithmT(const table_collection& tables)
        // Assumes tables are not empty.  Probably unsafe.
        {
            std::vector<std::size_t> pi(
                tables.edge_table.back().parent,
                std::numeric_limits<std::size_t>::max());
            auto p = fill_I_O(tables);

            // Move data for variable name convenience
            auto I = std::move(p.first);
            auto O = std::move(p.second);

            std::size_t j = 0, k = 0, M = I.size();
            double x = 0.0;
            // TODO: replace .at with []
            while (j < M)
                {
                    while (j < M
                           && tables.edge_table[I[j]].left == x) // Step T2
                        {
                            auto h = I[j];
                            pi.at(tables.edge_table.at(h).child) = h;
                            ++j;
                        }
                    // At this point, pi refers to the marginal tree
                    // beginning at x for all pi[i] != 
                    // std::numeric_limits<std::size_t>::max()
                    // The useful thing to do here would 
                    // be to define a "visitor function",
                    // taking pi and x as arguments.
                    //if (x != 0.0)
                    //    {
                    //        for (auto& e : tables.edge_table)
                    //            {
                    //                std::cerr << e.parent << ' ' << e.child
                    //                          << ' ' << e.left << ' '
                    //                          << e.right << '\n';
                    //            }
                    //        for (auto& pp : pi)
                    //            {
                    //                if (pp
                    //                    != std::numeric_limits<
                    //                           std::size_t>::max())
                    //                    {
                    //                        std::cout
                    //                            << tables.edge_table.at(pp)
                    //                                   .parent
                    //                            << ' '
                    //                            << tables.edge_table.at(pp)
                    //                                   .child
                    //                            << ' '
                    //                            << tables.edge_table.at(pp)
                    //                                   .left
                    //                            << '\n';
                    //                    }
                    //            }
                    //        std::exit(0);
                    //    }
                    if (j >= M)
                        break;
                    x = tables.edge_table.at(I.at(j)).left; // T3
                    while (k < M
                           && tables.edge_table.at(O[k]).right == x) // T4
                        {
                            auto h = O.at(k);
                            pi.at(tables.edge_table.at(h).child)
                                = std::numeric_limits<std::size_t>::max();
                            ++k;
                        }
                }
        }
    }
}
#endif
