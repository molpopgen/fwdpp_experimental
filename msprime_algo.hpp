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
/// and/or table_simplifier, so that we can
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
                if (pos == rhs.pos)
                    return time < rhs.time;
                return pos < rhs.pos;
            }
        };

        using index_vector = std::vector<index_key>;

        //std::pair<index_vector, index_vector>
        //fill_I_O(const table_collection& tables)
        //// TODO: better function name
        //// Fill I and O.  This is not described in Algorithm T,
        //// but can be sleuthed from msprime/tests/tsutil.py
        //{
        //    index_vector I, O;
        //    I.reserve(tables.edge_table.size());
        //    O.reserve(tables.edge_table.size());
        //    for (auto& e : tables.edge_table)
        //        {
        //            I.emplace_back(e.left,
        //                           -tables.node_table[e.parent].generation,
        //                           e.parent, e.child);
        //            O.emplace_back(e.right,
        //                           tables.node_table[e.parent].generation,
        //                           e.parent, e.child);
        //        }
        //    std::sort(I.begin(), I.end());
        //    std::sort(O.begin(), O.end());
        //    return std::make_pair(std::move(I), std::move(O));
        //}

        void
        algorithmT(const index_vector& input_left,
                   const index_vector& output_right, const std::size_t nnodes,
                   const double maxpos)
        // TODO: needs to take a visitor as argument
        // Assumes tables are not empty.  Probably unsafe.
        {
            std::vector<std::int32_t> pi(nnodes, -1);

            auto j = input_left.begin(), jM = input_left.end(),
                 k = output_right.begin(), kM = output_right.end();
            double x = 0.0;
            // TODO: replace .at with []
            while (j != jM || x < maxpos)
                {
                    while (k != kM && k->pos == x) // T4
                        {
                            pi[k->child] = -1;
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
                    // beginning at x for all pi[i] != -1
                    // The useful thing to do here would
                    // be to define a "visitor function",
                    // taking pi, x, and right as arguments.
                    // The i-th element in pi indexes its
                    // parent in the node table.

                    //if (x != 0.0)
                    //    {
                    //        for (auto p : pi)
                    //            {
                    //                std::cout << p << ' ';
                    //            }
                    //        std::exit(0);
                    //    }
                    //if (j >= M)
                    //    break;
                    x = right;
                }
        }

        struct marginal_tree
        {
            std::vector<std::int32_t> parents, leaf_counts;
            double left, right;
            marginal_tree(std::int32_t nnodes,
                          const std::vector<std::int32_t>& samples)
                : parents(nnodes, -1), leaf_counts(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() },
                  right{ std::numeric_limits<double>::quiet_NaN() }
            {
                for (auto s : samples)
                    {
                        if (s >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        leaf_counts[s] = 1;
                    }
            }
        };

        template <typename visitor>
        void
        algorithmL(const std::vector<index_key>& input_left,
                   const std::vector<index_key>& output_right,
                   const std::vector<std::int32_t>& sample_indexes,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            auto j = input_left.begin(), jM = input_left.end(),
                 k = output_right.begin(), kM = output_right.end();
            double x = 0.0;
            marginal_tree marginal(nnodes, sample_indexes);
            while (j != jM || x < maxpos)
                {
                    // TODO: this asserts may be incorrect
					assert(j<jM);
					assert(k<kM);
                    while (k != kM && k->pos == x) // T4
                        {
                            marginal.parents[k->child] = -1;
                            // Decrement leaf counts for outgoing nodes
                            auto p = k->parent;
                            auto lc = marginal.leaf_counts[k->child];
                            while (p != -1)
                                {
                                    marginal.leaf_counts[p] -= lc;
                                    assert(marginal.leaf_counts[p] >= 0);
                                    p = marginal.parents[p];
                                }
                            ++k;
                        }
                    while (j != jM && j->pos == x) // Step T2
                        {
                            // The entry for the child refers to
                            // the parent's location in the node table.
                            marginal.parents[j->child] = j->parent;
                            // Increment leaf counts for incoming nodes
                            auto p = j->parent;
                            auto lc = marginal.leaf_counts[j->child];
                            while (p != -1)
                                {
                                    marginal.leaf_counts[p] += lc;
                                    p = marginal.parents[p];
                                }
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
                    marginal.left = x;
                    marginal.right = right;
                    //This "yields"
                    //the data for this tree
                    //to the visitor
                    v(marginal);
                    x = right;
                }
        }
    }
}
#endif
