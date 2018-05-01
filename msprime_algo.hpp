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
        // but can be sleuthed from msprime/tests/tsutil.py
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
                                    tables.node_table[i->parent].generation,
                                    i->parent, i->child)
                           < std::tie(j->left,
                                      tables.node_table[j->parent].generation,
                                      j->parent, j->child);
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
                          //std::tie works via references.  So, to sort on -x,
                          //we need to make a copy, else we are trying to tie
                          //temporaries.
                          auto ig = -tables.node_table[i->parent].generation;
                          auto ip = -i->parent;
                          auto ic = -i->child;
                          auto jg = -tables.node_table[j->parent].generation;
                          auto jp = -j->parent;
                          auto jc = -j->child;

                          return std::tie(i->right, ig, ip, ic)
                                 < std::tie(j->right, jg, jp, jc);
                      });
            for (std::size_t i = 0; i < O.size(); ++i)
                {
                    O[i] = static_cast<std::size_t>(std::distance(
                        std::begin(tables.edge_table), edge_pointers[i]));
                }
            return std::make_pair(std::move(I), std::move(O));
        }

        void
        algorithmT(const table_collection& tables, const double maxpos)
        // Assumes tables are not empty.  Probably unsafe.
        {
            std::vector<std::size_t> pi(
                tables.node_table.size(),
                std::numeric_limits<std::size_t>::max());
            auto p = fill_I_O(tables);

            // Move data for variable name convenience
            auto I = std::move(p.first);
            auto O = std::move(p.second);

            std::size_t j = 0, k = 0, M = tables.edge_table.size();
            double x = 0.0;
            // TODO: replace .at with []
            while (j < M || x < maxpos)
                {
                    while (k < M
                           && tables.edge_table.at(O[k]).right == x) // T4
                        {
                            auto& edge_ = tables.edge_table.at(O[k]);
                            pi.at(edge_.child)
                                = std::numeric_limits<std::size_t>::max();
                            ++k;
                        }
                    while (j < M
                           && tables.edge_table[I[j]].left == x) // Step T2
                        {
                            auto& edge_ = tables.edge_table.at(I[j]);
							// The entry for the child refers to
							// the parent's location in the node table.
                            pi.at(edge_.child) = edge_.parent;
                            ++j;
                        }
                    double right = maxpos;
                    if (j < M)
                        {
                            right = std::min(right,
                                             tables.edge_table.at(I[j]).left);
                        }
                    if (k < M)
                        {
                            right = std::min(right,
                                             tables.edge_table.at(O[k]).right);
                        }
                    // At this point, pi refers to the marginal tree
                    // beginning at x for all pi[i] !=
                    // std::numeric_limits<std::size_t>::max()
                    // The useful thing to do here would
                    // be to define a "visitor function",
                    // taking pi, x, and right as arguments.
					// The i-th element in pi indexes its 
					// parent in the node table.

                    //if (x != 0.0)
                    //    {
                    //        for (auto& e : tables.edge_table)
                    //            {
                    //                std::cerr << e.parent << ' ' << e.child
                    //                          << ' ' << e.left << ' '
                    //                          << e.right << '\n';
                    //            }
					//		for(std::size_t x=0;x<pi.size();++x)
					//		{
					//			if(pi[x]!=std::numeric_limits<std::size_t>::max())
					//			{
					//				std::cout<<x <<' '<<pi[x]<<' '<<tables.node_table[x].generation<<' '<<tables.node_table[pi[x]].generation<<'\n';
					//			}
					//		}
                    //       std::exit(0);
                    //    }
                    //if (j >= M)
                    //    break;
                    x = right;
                }
        }
    }
}
#endif
