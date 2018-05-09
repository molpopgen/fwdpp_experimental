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
        struct index_key
        {
            double pos, time;
            std::int32_t parent, child;
            explicit index_key(double pos_, double t, std::int32_t p,
                               std::int32_t c)
                : pos{ pos_ }, time{ t }, parent{ p }, child{ c }
            {
            }
            inline bool
            operator<(const index_key& rhs) const
            {
                return std::tie(pos, time) < std::tie(rhs.pos, rhs.time);
            }
        };

        using index_vector = std::vector<index_key>;
        std::pair<index_vector, index_vector>
        fill_I_O(const table_collection& tables)
        // TODO: better function name
        // Fill I and O.  This is not described in Algorithm T,
        // but can be sleuthed from msprime/tests/tsutil.py
        {
            index_vector I, O;
            I.reserve(tables.edge_table.size());
            O.reserve(tables.edge_table.size());
            for (auto& e : tables.edge_table)
                {
                    I.emplace_back(e.left,
                                   -tables.node_table[e.parent].generation,
                                   e.parent, e.child);
                    O.emplace_back(e.right,
                                   tables.node_table[e.parent].generation,
                                   e.parent, e.child);
                }
            std::sort(I.begin(), I.end());
            std::sort(O.begin(), O.end());
            return std::make_pair(std::move(I), std::move(O));
        }

        void
        algorithmT(const index_vector& input_left,
                   const index_vector& output_right, const std::size_t nnodes,
                   const double maxpos)
        // Assumes tables are not empty.  Probably unsafe.
        {
            std::vector<std::size_t> pi(
                nnodes, std::numeric_limits<std::size_t>::max());

            auto j = input_left.begin(), jM = input_left.end(),
                 k = output_right.begin(), kM = output_right.end();
            double x = 0.0;
            // TODO: replace .at with []
            while (j != jM || x < maxpos)
                {
                    while (k != kM && k->pos == x) // T4
                        {
                            pi[k->child]
                                = std::numeric_limits<std::size_t>::max();
                            ++k;
                        }
                    while (j != jM && j->pos == x) // Step T2
                        {
                            // The entry for the child refers to
                            // the parent's location in the node table.
                            pi[j->child] = j->parent;
                            ++j;
                        }
                    double right = maxpos;
                    if (j != jM)
                        {
                            right = std::min(right, j->pos);
                        }
                    if (k != kM)
                        {
                            right = std::min(right, k->pos);
                        }
                    //for(auto & p : pi){std::cout<<p <<' ';}
                    //std::cout<<'\n';
                    // At this point, pi refers to the marginal tree
                    // beginning at x for all pi[i] !=
                    // std::numeric_limits<std::size_t>::max()
                    // The useful thing to do here would
                    // be to define a "visitor function",
                    // taking pi, x, and right as arguments.
                    // The i-th element in pi indexes its
                    // parent in the node table.

                    //if (x != 0.0)
					//{
					//	for(auto p:pi){std::cout<<p<<' ';}
					//	std::exit(0);
					//}
                    //if (j >= M)
                    //    break;
                    x = right;
                }
        }
    }
}
#endif
