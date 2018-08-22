#ifndef FWDPP_TS_ITERATE_MARGINAL_TREES_HPP
#define FWDPP_TS_ITERATE_MARGINAL_TREES_HPP

#include "marginal_tree.hpp"
#include "indexed_edge.hpp"
#include "iterate_marginal_trees_details.hpp"

namespace fwdpp
{
    namespace ts
    {
        template <typename visitor, typename leaf_policy,
                  typename sample_list_policy>
        inline void
        iterate_marginal_trees(marginal_tree& marginal,
                               const indexed_edge_container& input_left,
                               const indexed_edge_container& output_right,
                               const double maxpos, visitor v,
                               const leaf_policy lp,
                               const sample_list_policy slp)
        {
            auto j = input_left.begin(), jM = input_left.end(),
                 k = output_right.begin(), kM = output_right.end();
            double x = 0.0;
            while (j < jM || x < maxpos)
                {
                    while (k < kM && k->pos == x) // T4
                        {
                            const auto p = k->parent;
                            const auto c = k->child;
                            const auto lsib = marginal.left_sib[c];
                            const auto rsib = marginal.right_sib[c];
                            if (lsib == -1)
                                {
                                    marginal.left_child[p] = rsib;
                                }
                            else
                                {
                                    marginal.right_sib[lsib] = rsib;
                                }
                            if (rsib == -1)
                                {
                                    marginal.right_child[p] = lsib;
                                }
                            else
                                {
                                    marginal.left_sib[rsib] = lsib;
                                }
                            marginal.parents[c] = -1;
                            marginal.left_sib[c] = -1;
                            marginal.right_sib[c] = -1;
							detail::outgoing_leaf_counts(marginal, k->parent, k->child,
                                                 lp);
							detail::update_sample_list(marginal, k->parent, slp);
                            ++k;
                        }
                    while (j < jM && j->pos == x) // Step T2
                        {
                            const auto p = j->parent;
                            const auto c = j->child;
                            const auto rchild = marginal.right_child[p];
                            if (rchild == -1)
                                {
                                    marginal.left_child[p] = c;
                                    marginal.left_sib[c] = -1;
                                    marginal.right_sib[c] = -1;
                                }
                            else
                                {
                                    marginal.right_sib[rchild] = c;
                                    marginal.left_sib[c] = rchild;
                                    marginal.right_sib[c] = -1;
                                }
                            // The entry for the child refers to
                            // the parent's location in the node table.
                            marginal.parents[c] = j->parent;
                            marginal.right_child[p] = c;
							detail::incoming_leaf_counts(marginal, j->parent, j->child,
                                                 lp);
							detail::update_sample_list(marginal, j->parent, slp);
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
    } // namespace ts
} // namespace fwdpp
#endif
