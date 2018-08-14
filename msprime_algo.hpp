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
    namespace ts
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
                    {
                        return time < rhs.time;
                    }
                return pos < rhs.pos;
            }
        };

        struct do_nothing
        {
        };

        struct count_leaves
        {
        };

        using index_vector = std::vector<index_key>;

        struct marginal_tree
        {
            //TODO separate leaf_counts from this type,
            //and require the visitor to take a const &
            //as an arg.
            std::vector<std::int32_t> parents, leaf_counts;
            double left, right;
            marginal_tree(std::int32_t nnodes,
                          const std::vector<std::int32_t>& samples)
                : parents(nnodes, -1), leaf_counts(nnodes, 0),
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  }
            {
                for (auto s : samples)
                    {
                        if (static_cast<std::size_t>(s) >= leaf_counts.size())
                            {
                                throw std::invalid_argument(
                                    "sample index out of range");
                            }
                        leaf_counts[s] = 1;
                    }
            }
            marginal_tree(std::int32_t nnodes)
                : parents(nnodes, -1), leaf_counts{},
                  left{ std::numeric_limits<double>::quiet_NaN() }, right{
                      std::numeric_limits<double>::quiet_NaN()
                  }
            {
            }
        };

        void
        outgoing_leaf_counts(marginal_tree&, const std::int32_t,
                             const std::int32_t, const do_nothing)
        // TODO: internal namespace
        {
        }

        inline void
        incoming_leaf_counts(marginal_tree&, const std::int32_t,
                             const std::int32_t, const do_nothing)
        // TODO: internal namespace
        {
        }

        inline void
        outgoing_leaf_counts(marginal_tree& marginal,
                             const std::int32_t parent,
                             const std::int32_t child, const count_leaves)
        // TODO: internal namespace
        {
            auto p = parent;
            auto lc = marginal.leaf_counts[child];
            while (p != -1)
                {
                    marginal.leaf_counts[p] -= lc;
                    assert(marginal.leaf_counts[p] >= 0);
                    p = marginal.parents[p];
                }
        }

        inline void
        incoming_leaf_counts(marginal_tree& marginal,
                             const std::int32_t parent,
                             const std::int32_t child, const count_leaves)
        // TODO: internal namespace
        {
            auto p = parent;
            auto lc = marginal.leaf_counts[child];
            while (p != -1)
                {
                    marginal.leaf_counts[p] += lc;
                    p = marginal.parents[p];
                }
        }

        template <typename visitor, typename leaf_policy>
        inline void
        iterate_marginal_trees(const leaf_policy lp,
                               const index_vector& input_left,
                               const index_vector& output_right,
                               const double maxpos, marginal_tree& marginal,
                               visitor v)
        // TODO: internal namespace
        {
            auto j = input_left.begin(), jM = input_left.end(),
                 k = output_right.begin(), kM = output_right.end();
            double x = 0.0;
            // TODO: replace .at with []
            while (j < jM || x < maxpos)
                {
                    while (k < kM && k->pos == x) // T4
                        {
                            marginal.parents[k->child] = -1;
                            outgoing_leaf_counts(marginal, k->parent, k->child,
                                                 lp);
                            ++k;
                        }
                    while (j < jM && j->pos == x) // Step T2
                        {
                            // The entry for the child refers to
                            // the parent's location in the node table.
                            marginal.parents[j->child] = j->parent;
                            incoming_leaf_counts(marginal, j->parent, j->child,
                                                 lp);
                            ++j;
                        }
                    double right = maxpos;
                    if (j < jM)
                        {
                            right = std::min(right, j->pos);
                        }
                    if (k < kM)
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
        template <typename visitor>
        void
        algorithmT(const std::vector<index_key>& input_left,
                   const std::vector<index_key>& output_right,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes);
            iterate_marginal_trees(do_nothing(), input_left, output_right,
                                   maxpos, marginal, v);
        }

        template <typename visitor>
        void
        algorithmL(const std::vector<index_key>& input_left,
                   const std::vector<index_key>& output_right,
                   const std::vector<std::int32_t>& sample_indexes,
                   const std::int32_t nnodes, const double maxpos, visitor v)
        {
            marginal_tree marginal(nnodes, sample_indexes);
            iterate_marginal_trees(count_leaves(), input_left, output_right,
                                   maxpos, marginal, v);
        }
    } // namespace ts
} // namespace fwdpp
#endif
