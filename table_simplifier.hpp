#ifndef FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP__
#define FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP__

#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstddef>
#include <stdexcept>
#include <unordered_map>
#include "node.hpp"
#include "edge.hpp"
#include "table_collection.hpp"
#include "msprime_algo.hpp"

namespace fwdpp
{
    namespace ts
    {
        class table_simplifier
        {
          private:
            struct segment
            {
                double left, right;
                std::int32_t node;
                segment(double l, double r, std::int32_t n)
                    : left{ l }, right{ r }, node{ n }
                {
                    if (right <= left)
                        {
                            throw std::invalid_argument(
                                "right must be > left");
                        }
                }
            };

            // These are temp tables/buffer
            // for simplification.  We keep
            // their allocated memory persistent.
            edge_vector new_edge_table;
            node_vector new_node_table;
            // segment_queue mimics a min queue of segments w.r.to
            // segment::left. X is a temporary container
            // for storing segments during ancestry merging.
            std::vector<segment> segment_queue, X;
            std::vector<std::vector<segment>> Ancestry;
            /// Temp container used for compacting edges
            edge_vector E;
            // region length
            const double L;

            bool
            add_to_queue(const double left, const double right,
                         const std::int32_t node, const bool added)
            {
                bool rv = added;
                if (!segment_queue.empty() && left > segment_queue.back().left)
                    {
                        return false;
                    }
                segment_queue.emplace_back(left, right, node);
                return rv;
            }

            void
            sort_queue(std::size_t beg) noexcept
            {
                std::sort(segment_queue.begin() + beg, segment_queue.end(),
                          [](const segment& a, const segment& b) {
                              return a.left > b.left;
                          });
                assert(std::is_sorted(segment_queue.begin(),
                                      segment_queue.end(),
                                      [](const segment& a, const segment& b) {
                                          return a.left > b.left;
                                      }));
            }

            void
            cleanup() noexcept
            // Clears out data from
            // temp containers after simplify.
            // Retains container capacity.
            {
                new_edge_table.clear();
                new_node_table.clear();
                E.clear();
                for (auto& ai : Ancestry)
                    {
                        ai.clear();
                    }
            }

            edge_vector::const_iterator
            step_S3(edge_vector::const_iterator edge_ptr,
                    const edge_vector::const_iterator edge_end, std::int32_t u)
            {
                for (; edge_ptr < edge_end && edge_ptr->parent == u;
                     ++edge_ptr)
                    {
                        // For each edge corresponding to this parent,
                        // we look at all segments from the child.
                        // If the two segments overlap, we add the
                        // minimal
                        // overlap to our queue.
                        // This is Step S3.
                        // TODO: the data here are sorted in ascending
                        // order by left, meaning we can process these
                        // data using binary searches if we had an
                        // interval
                        // tree data structure instead of a straight
                        // vector.
                        for (auto& seg : Ancestry[edge_ptr->child])
                            {
                                if (seg.right > edge_ptr->left
                                    && edge_ptr->right > seg.left)
                                    {
                                        segment_queue.emplace_back(
                                            std::max(seg.left, edge_ptr->left),
                                            std::min(seg.right,
                                                     edge_ptr->right),
                                            seg.node);
                                    }
                            }
                    }
                if (!segment_queue.empty())
                    {
                        sort_queue(0);
                    }
                return edge_ptr;
            }

            void
            merge_ancestors(const node_vector& input_node_table,
                            const std::int32_t parent_input_id,
                            std::vector<std::int32_t>& idmap)
            // TODO: will have to be made aware of sample labels.
            {
                std::int32_t anode, znode = -1;
                double aleft, aright, zright = 0;
                std::int32_t v = -1;
                bool defrag_required = false;
                std::size_t current_queue_size = segment_queue.size();
                bool added_to_queue = false;
                E.clear();
                while (!segment_queue.empty())
                    // Steps S4 through S8 of the algorithm.
                    {
                        X.clear();
                        auto l = segment_queue.back().left;
                        double r = L;
                        while (!segment_queue.empty()
                               && segment_queue.back().left == l)
                            // This while loop is Step S4. This step
                            // adds to X all segments with left == l
                            // and also finds the minimum right for
                            // all segs with left == l.
                            // TODO: this can be done with reverse
                            // iteration,
                            // but testing on 0.5e9 edges didn't seem
                            // to
                            // make it worthwhile.
                            {
                                r = std::min(r, segment_queue.back().right);
                                X.emplace_back(
                                    std::move(segment_queue.back()));
                                assert(current_queue_size > 0);
                                segment_queue.pop_back();
                                --current_queue_size;
                            }
                        double next_left = 0.0;
                        if (!segment_queue.empty())
                            {
                                next_left = segment_queue.back().left;
                                r = std::min(r, next_left);
                            }
                        if (X.size() == 1)
                            {
                                if (!segment_queue.empty()
                                    && next_left < X[0].right)
                                    {
                                        aleft = X[0].left;
                                        aright = next_left;
                                        anode = X[0].node;
                                        X[0].left = next_left;
                                        added_to_queue = add_to_queue(
                                            next_left, X[0].right, X[0].node,
                                            added_to_queue);
                                    }
                                else
                                    {
                                        aleft = X[0].left;
                                        aright = X[0].right;
                                        anode = X[0].node;
                                    }
                            }
                        else
                            {
                                if (v == -1)
                                    {
                                        // Overlap/coalescence, and
                                        // thus
                                        // a new node. Step S6.
                                        new_node_table.emplace_back(node{
                                            static_cast<std::int32_t>(
                                                new_node_table.size()),
                                            input_node_table[parent_input_id]
                                                .population,
                                            input_node_table[parent_input_id]
                                                .generation });
                                        v = new_node_table.size() - 1;
                                        // update sample map
                                        idmap[parent_input_id] = v;
                                    }
                                aleft = l;
                                aright = r;
                                anode = v;
                                for (auto& x : X)
                                    {
                                        E.emplace_back(
                                            edge{ l, r, v, x.node });
                                        if (x.right > r)
                                            {
                                                x.left = r;
                                                added_to_queue = add_to_queue(
                                                    x.left, x.right, x.node,
                                                    added_to_queue);
                                            }
                                    }
                            }
                        if (added_to_queue)
                            {
                                sort_queue(current_queue_size);
                                added_to_queue = false;
                            }
                        // need to make sure this variable is up to date
                        current_queue_size = segment_queue.size();
                        Ancestry[parent_input_id].emplace_back(aleft, aright,
                                                               anode);
                        defrag_required |= zright == aleft && znode == anode;
                        zright = aright;
                        znode = anode;
                    }
                assert(current_queue_size == 0);
                if (defrag_required)
                    {
                        defragment(Ancestry[parent_input_id]);
                    }
                // remove redundant info from
                // edge data for this parent
                squash_and_flush_edges();
            }

            void
            squash_and_flush_edges()
            // Implementation copied from msprime.
            // Squashes identical edges on a per-parent
            // basis and adds them to the output list of edges.
            // Based on implementation found in msprime/lib/msprime.c
            // by Jerome Kelleher, version 0.5.0 of that code.
            // http://github.com/jeromekelleher/msprime.
            {
                if (!E.empty())
                    {
                        std::sort(E.begin(), E.end(),
                                  [](const edge& a, const edge& b) {
                                      return std::tie(a.child, a.left)
                                             < std::tie(b.child, b.left);
                                  });
                        std::size_t j = 0, k, l = 0;
                        for (k = 1; k < E.size(); ++k)
                            {
                                assert(E[k - 1].parent == E[k].parent);
                                if (E[j].child != E[k].child
                                    || E[k - 1].right != E[k].left)
                                    {
                                        auto e = E[j];
                                        e.right = E[k - 1].right;
                                        E[l] = e;
                                        j = k;
                                        ++l;
                                    }
                            }
                        auto e = E[j];
                        e.right = E[k - 1].right;
                        E[l] = e;
                        new_edge_table.insert(
                            new_edge_table.end(),
                            std::make_move_iterator(E.begin()),
                            std::make_move_iterator(E.begin() + l + 1));
                        E.clear();
                    }
            }

            void
            defragment(std::vector<segment>& ancestry_segment) noexcept
            /// This implementation is based on Jerome Kelleher's
            /// implementation in msprime/lib/table.c, which is found at
            /// http://github.com/jeromekelleher/msprime.
            {
                assert(std::is_sorted(ancestry_segment.begin(),
                                      ancestry_segment.end(),
                                      [](const segment& a, const segment& b) {
                                          return std::tie(a.left, a.right)
                                                 < std::tie(b.left, b.right);
                                      }));
                auto prev_seg = ancestry_segment.begin();
                auto next_seg = prev_seg + 1;
                while (next_seg < ancestry_segment.end())
                    {
                        assert(prev_seg->node != -1);
                        if (prev_seg->node == next_seg->node
                            && prev_seg->right == next_seg->left)
                            {
                                prev_seg->right = next_seg->right;
                                next_seg->node = -1;
                                ++next_seg;
                            }
                        else
                            {
                                prev_seg = next_seg++;
                            }
                    }
                ancestry_segment.erase(
                    std::remove_if(
                        ancestry_segment.begin(), ancestry_segment.end(),
                        [](const segment& seg) { return seg.node == -1; }),
                    ancestry_segment.end());
            }

            using mutation_map_t = std::unordered_map<
                std::int32_t,
                std::vector<std::pair<std::size_t, std::size_t>>>;
            using mutation_node_map_t = std::vector<std::int32_t>;

            template <typename mcont_t>
            mutation_map_t
            prep_mutation_simplification(
                const mcont_t& mutations,
                const mutation_key_vector& mutation_table) const
            {
                mutation_map_t mutation_map;
                for (std::size_t i = 0; i < mutation_table.size(); ++i)
                    {
                        mutation_map[mutation_table[i].node].emplace_back(
                            mutation_table[i].key, i);
                    }

                for (auto& mm : mutation_map)
                    {
                        std::sort(
                            mm.second.begin(), mm.second.end(),
                            [&mutations](
                                const std::pair<std::size_t, std::size_t>& a,
                                const std::pair<std::size_t, std::size_t>& b) {
                                return mutations[a.first].pos
                                       < mutations[b.first].pos;
                            });
                    }

                return mutation_map;
            }

            template <typename mcont_t>
            void
            simplify_mutations(const mcont_t& mutations,
                               std::vector<std::uint32_t>& mcounts,
                               const std::vector<std::int32_t>& samples,
                               const std::vector<std::int32_t>& idmap,
                               table_collection& tables,
                               mutation_map_t& mutation_map) const
            {
                // Simplify mutations

                // 0. Remap mutation input node ids.  To do this, we use
                // the data stored in Ancestry, which allows us to "push"
                // mutation nodes down the tree.

                for (auto& mr : tables.mutation_table)
                    {
                        mr.node = -1;
                    }
                for (auto& mm : mutation_map)
                    {
                        auto seg = Ancestry[mm.first].cbegin();
                        const auto seg_e = Ancestry[mm.first].cend();
                        auto mut = mm.second.cbegin();
                        const auto mute = mm.second.cend();
                        while (seg < seg_e && mut < mute)
                            {
                                auto pos = mutations[mut->first].pos;
                                if (seg->left <= pos && pos < seg->right)
                                    {
                                        assert(mut < mute);
                                        assert(mut->first < mutations.size());
                                        //The following assert fails when
                                        //processing "decaptitated" trees
                                        //assert(
                                        //    static_cast<std::size_t>(seg->node)
                                        //    < new_edge_table.size());
                                        tables.mutation_table[mut->second].node
                                            = seg->node;
                                        ++mut;
                                    }
                                else if (pos >= seg->right)
                                    {
                                        ++seg;
                                    }
                                else
                                    {
                                        ++mut;
                                    }
                            }
                    }

                // 1. Remove all mutations whose output nodes are simply gone.
                // Note this does not remove mutations where the node still exists
                // somewhere in the pedigree, but the mutation is on a marginal tree
                // where the node is not an ancestor.
                // This is a fast O(n).
                tables.mutation_table.erase(
                    std::remove_if(tables.mutation_table.begin(),
                                   tables.mutation_table.end(),
                                   [&idmap](const mutation_record& mr) {
                                       return mr.node == -1;
                                   }),
                    tables.mutation_table.end());
                assert(std::is_sorted(
                    tables.mutation_table.begin(), tables.mutation_table.end(),
                    [&mutations](const mutation_record& a,
                                 const mutation_record& b) {
                        return mutations[a.key].pos < mutations[b.key].pos;
                    }));

                // 2. Now, we use Kelleher et al. (2016)'s Algorithm L
                // to march through each marginal tree and its leaf
                // counts. At the same time, we march through our mutation
                // table, which is sorted by position.
                std::fill(mcounts.begin(), mcounts.end(), 0);
                mcounts.resize(mutations.size(), 0);

                auto mtable_itr = tables.mutation_table.begin();
                auto mtable_end = tables.mutation_table.end();
                auto mutation_counter = [&mutations, &mtable_itr, mtable_end,
                                         &mcounts](
                                            const marginal_tree& marginal) {
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.right)
                        {
                            assert(mutations[mtable_itr->key].pos
                                   >= marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < marginal.right);
                            mcounts[mtable_itr->key]
                                = marginal.leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                };

                tables.build_indexes();
                std::vector<std::int32_t> remapped_samples(samples.size());
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        remapped_samples[i] = idmap[samples[i]];
                    }
                algorithmL(tables.input_left, tables.output_right,
                           remapped_samples, tables.node_table.size(),
                           tables.L, mutation_counter);

                // 3. Remove any elements from table with count == 0
                // TODO: this has to be updated to handle fixations
                // and losses in one pass!
                tables.mutation_table.erase(
                    std::remove_if(tables.mutation_table.begin(),
                                   tables.mutation_table.end(),
                                   [&mcounts](const mutation_record& mr) {
                                       return mcounts[mr.key] == 0;
                                   }),
                    tables.mutation_table.end());
            }

          public:
            table_simplifier(const double maxpos)
                : new_edge_table{}, new_node_table{},
                  segment_queue{}, X{}, Ancestry{}, E{}, L{ maxpos }
            {
                if(maxpos < 0 || !std::isfinite(maxpos))
                {
                    throw std::invalid_argument("maxpos must be > 0 and finite");
                }
            }

            template <typename mutation_container>
            std::vector<std::int32_t>
            simplify(table_collection& tables,
                     const std::vector<std::int32_t>& samples,
                     const mutation_container& mutations,
                     std::vector<std::uint32_t>& mcounts)
            /// Set theoretic simplify.
            /// TODO: compare against current msprime, 
            /// which has streamlined some steps
            {
                Ancestry.resize(tables.node_table.size());

                // Set some things up for later mutation simplification
                auto mutation_map = prep_mutation_simplification(
                    mutations, tables.mutation_table);

                // Relates input node ids to output node ids
                std::vector<std::int32_t> idmap(tables.node_table.size(), -1);

                // We take our samples and add them to both the output
                // node list and initialize their ancestry with
                // a segment on [0,L).
                for (const auto& s : samples)
                    {
                        new_node_table.emplace_back(node{
                            static_cast<std::int32_t>(new_node_table.size()),
                            tables.node_table[s].population,
                            tables.node_table[s].generation });
                        Ancestry[s].emplace_back(
                            0, L,
                            static_cast<std::int32_t>(new_node_table.size()
                                                      - 1));
                        idmap[s] = static_cast<std::int32_t>(
                            new_node_table.size() - 1);
                    }

                // At this point, our edges are sorted by birth
                // order of parents, from present to past.
                // We can now work our way up the pedigree.
                // This outer loop differs from how we describe it in the
                // paper, but the strict sorting of edges means that this
                // equivalent.
                auto edge_ptr = tables.edge_table.cbegin();
                const auto edge_end = tables.edge_table.cend();
                while (edge_ptr < edge_end)
                    {
                        auto u = edge_ptr->parent;
                        edge_ptr = step_S3(edge_ptr, edge_end, u);
                        merge_ancestors(tables.node_table, u, idmap);
                    }

                assert(static_cast<std::size_t>(std::count_if(
                           idmap.begin(), idmap.end(),
                           [](const std::int32_t i) { return i != -1; }))
                       == new_node_table.size());
                // After swapping, new_node_table
                // contains the input nodes
                tables.edge_table.swap(new_edge_table);
                tables.node_table.swap(new_node_table);
                // TODO: allow for exception instead of assert
                assert(tables.edges_are_sorted());
                tables.update_offset();

                simplify_mutations(mutations, mcounts, samples, idmap, tables,
                                   mutation_map);

                cleanup();
                return idmap;
            }
        };
    } // namespace ancestry
} // namespace fwdpp

#endif
