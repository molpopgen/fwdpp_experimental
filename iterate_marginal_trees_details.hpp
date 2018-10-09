#ifndef FWDPP_TS_ITERATE_MARGINAL_TREES_DETAIL_HPP
#define FWDPP_TS_ITERATE_MARGINAL_TREES_DETAIL_HPP

#include <cstdint>
#include <type_traits>
#include <cassert>
#include "marginal_tree.hpp"

namespace fwdpp
{
    namespace ts
    {
        namespace detail
        {
            void
            outgoing_leaf_counts(marginal_tree&, const std::int32_t,
                                 const std::int32_t, const std::false_type)
            {
            }

            inline void
            incoming_leaf_counts(marginal_tree&, const std::int32_t,
                                 const std::int32_t, const std::false_type)
            {
            }

            inline void
            outgoing_leaf_counts(marginal_tree& marginal,
                                 const std::int32_t parent,
                                 const std::int32_t child,
                                 const std::true_type)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                while (p != -1)
                    {
                        marginal.leaf_counts[p] -= lc;
                        marginal.preserved_leaf_counts[p] -= plc;
                        assert(marginal.leaf_counts[p] >= 0);
                        assert(marginal.preserved_leaf_counts[p] >= 0);
                        p = marginal.parents[p];
                    }
            }

            inline void
            incoming_leaf_counts(marginal_tree& marginal,
                                 const std::int32_t parent,
                                 const std::int32_t child,
                                 const std::true_type)
            {
                auto p = parent;
                auto lc = marginal.leaf_counts[child];
                auto plc = marginal.preserved_leaf_counts[child];
                while (p != -1)
                    {
                        marginal.leaf_counts[p] += lc;
                        marginal.preserved_leaf_counts[p] += plc;
                        p = marginal.parents[p];
                    }
            }

            inline void
            update_sample_list(marginal_tree& marginal,
                               const std::int32_t parent, const std::true_type)
            {
                for (auto n = parent; n != -1; n = marginal.parents[n])
                    {
                        if (marginal.is_sample[n] == 1)
                            {
                                marginal.left_sample[n]
                                    = marginal.right_sample[n];
                            }
                        else
                            {
                                marginal.left_sample[n]
                                    = marginal.right_sample[n] = -1;
                            }
                        for (auto v = marginal.left_child[n]; v != -1;
                             v = marginal.right_sib[v])
                            {
                                if (marginal.left_sample[v] != -1)
                                    {
                                        assert(marginal.right_sample[v] != -1);
                                        if (marginal.left_sample[n] == -1)
                                            {
                                                marginal.left_sample[n]
                                                    = marginal.left_sample[v];
                                                marginal.right_sample[n]
                                                    = marginal.right_sample[v];
                                            }
                                        else
                                            {
                                                marginal.next_sample
                                                    [marginal.right_sample[n]]
                                                    = marginal.left_sample[v];
                                                marginal.right_sample[n]
                                                    = marginal.right_sample[v];
                                            }
                                    }
                            }
                    }
            }

            inline void
            update_sample_list(marginal_tree&, const std::int32_t,
                               const std::false_type)
            {
            }
        } // namespace detail
    }     // namespace ts
} // namespace fwdpp

#endif
