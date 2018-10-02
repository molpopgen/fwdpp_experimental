#ifndef FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP
#define FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP

#include <iostream>
#include <cmath>
#include <map>
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

            class segment_overlapper
            {
              private:
                std::vector<segment>::const_iterator sbeg, send;

                inline double
                set_partition()
                {
                    double tright = std::numeric_limits<double>::max();
                    auto b = overlapping.begin();
                    for (auto i = overlapping.begin(); i < overlapping_end;
                         ++i)
                        {
                            if (i->right > left)
                                {
                                    *b = *i;
                                    tright = std::min(tright, b->right);
                                    ++b;
                                }
                        }
                    overlapping_end = b;
                    return tright;
                    //overlapping_end = std::stable_partition(
                    //    overlapping.begin(), overlapping_end,
                    //    [this](const segment& seg) {
                    //        return seg.right > left;
                    //    });
                }

                //inline double
                //min_right_overlap()
                //{
                //    return std::min_element(
                //               overlapping.begin(), overlapping_end,
                //               [](const segment& a, const segment& b) {
                //                   return a.right < b.right;
                //               })
                //        ->right;
                //}

              public:
                std::vector<segment> overlapping;
                std::vector<segment>::iterator overlapping_end;
                double left, right;
                segment_overlapper(const std::vector<segment>& segs)
                    // The - 1 for send assumes a "cap"/sentinel value.
                    : sbeg(segs.begin()), send(segs.end() - 1), overlapping{},
                      overlapping_end(overlapping.end()), left(0),
                      right(std::numeric_limits<double>::max())
                {
                }
                bool
                operator()()
                {
                    bool rv = 0;
                    if (sbeg < send)
                        {
                            left = right;
                            auto tright = set_partition();
                            if (num_overlaps() == 0)
                                {
                                    left = sbeg->left;
                                }
                            while (sbeg < send && sbeg->left == left)
                                {
                                    tright = std::min(tright, sbeg->right);
                                    overlapping_end
                                        = overlapping.insert(overlapping_end,
                                                             *sbeg)
                                          + 1;
                                    ++sbeg;
                                }
                            right = std::min(sbeg->left, tright);
                            rv = true;
                        }
                    else
                        {
                            left = right;
                            right = std::numeric_limits<double>::max();
                            auto tright = set_partition();
                            if (num_overlaps() > 0)
                                {
                                    right = tright;
                                    rv = true;
                                }
                        }
                    return rv;
                }
                std::int64_t
                num_overlaps()
                {
                    return std::distance(overlapping.begin(), overlapping_end);
                }
            };
            // These are temp tables/buffer
            // for simplification.  We keep
            // their allocated memory persistent.
            edge_vector new_edge_table;
            node_vector new_node_table;
            // segment_queue mimics a min queue of segments w.r.to
            // segment::left.
            std::vector<segment> segment_queue;
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
                segment_queue.clear();
                for (; edge_ptr < edge_end && edge_ptr->parent == u;
                     ++edge_ptr)
                    {
                        // For each edge corresponding to this parent,
                        // we look at all segments from the child.
                        // If the two segments overlap, we add the
                        // minimal
                        // overlap to our queue.
                        // This is Step S3.
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
                // Sort for processing via the overlapper
                std::sort(segment_queue.begin(), segment_queue.end(),
                          [](const segment& a, const segment& b) {
                              return a.left < b.left;
                          });
                // Add sentinel
                segment_queue.emplace_back(
                    segment{ std::numeric_limits<double>::max()
                                 / 2., //TODO fix this ugly hack
                             std::numeric_limits<double>::max(), -1 });
                return edge_ptr;
            }

            void
            buffer_edge(
                //std::map<std::int32_t, std::vector<edge>>& buffered_edges,
                std::vector<edge>& buffered_edges, const double left,
                const double right, const std::int32_t parent,
                const std::int32_t child)
            {
                auto itr = std::find_if(
                    buffered_edges.rbegin(), buffered_edges.rend(),
                    [child](const edge& e) { return e.child == child; });
                if (itr == buffered_edges.rend())
                    {
                        buffered_edges.emplace_back(
                            edge{ left, right, parent, child });
                    }
                else
                    {
                        if (itr->right == left)
                            {
                                itr->right = right;
                            }
                        else
                            {
                                buffered_edges.emplace_back(
                                    edge{ left, right, parent, child });
                            }
                    }
                //auto itr = buffered_edges.find(child);
                //if (itr == buffered_edges.end())
                //    {
                //        buffered_edges[child].emplace_back(
                //            edge{ left, right, parent, child });
                //    }
                //else
                //    {
                //        auto& last = itr->second.back();
                //        if (last.right == left)
                //            {
                //                last.right = right;
                //            }
                //        else
                //            {
                //                itr->second.emplace_back(
                //                    edge{ left, right, parent, child });
                //            }
                //    }
            }

            void
            add_ancestry(std::int32_t input_id, double left, double right,
                         std::int32_t node)
            {
                if (Ancestry[input_id].empty())
                    {
                        //assert(left < right);
                        Ancestry[input_id].emplace_back(left, right, node);
                    }
                else
                    {
                        auto& last = Ancestry[input_id].back();
                        if (last.right == left && last.node == node)
                            {
                                last.right = right;
                            }
                        else
                            {
                                //assert(left < right);
                                Ancestry[input_id].emplace_back(left, right,
                                                                node);
                            }
                    }
            }

            void
            merge_ancestors(const node_vector& input_node_table,
                            const std::int32_t parent_input_id,
                            std::vector<std::int32_t>& idmap)
            // TODO: will have to be made aware of sample labels.
            // in order to handle ancient samples.
            {
                auto output_id = idmap[parent_input_id];
                bool is_sample = (output_id != -1);
                if (is_sample == true)
                    {
                        Ancestry[parent_input_id].clear();
                    }
                double previous_right = 0.0;
                segment_overlapper o(segment_queue);
                std::int32_t ancestry_node = -1;
                //std::map<std::int32_t, std::vector<edge>> buffered_edges;
                //std::vector<edge> buffered_edges;
                E.clear();
                while (o() == true)
                    {
                        if (o.num_overlaps() == 1)
                            {
                                ancestry_node = o.overlapping[0].node;
                                if (is_sample)
                                    {
                                        buffer_edge(E, o.left, o.right,
                                                    output_id, ancestry_node);
                                        ancestry_node = output_id;
                                    }
                            }
                        else
                            {
                                if (output_id == -1)
                                    {
                                        new_node_table.emplace_back(node{
                                            input_node_table[parent_input_id]
                                                .population,
                                            input_node_table[parent_input_id]
                                                .generation });
                                        output_id = new_node_table.size() - 1;
                                        // update sample map
                                        idmap[parent_input_id] = output_id;
                                    }
                                ancestry_node = output_id;
                                for (auto x = o.overlapping.begin();
                                     x < o.overlapping_end; ++x)
                                    {
                                        buffer_edge(E, o.left, o.right,
                                                    output_id, x->node);
                                    }
                            }
                        if (is_sample && o.left != previous_right)
                            {
                                //assert(previous_right < o.left);
                                add_ancestry(parent_input_id, previous_right,
                                             o.left, ancestry_node);
                            }
                        if (!(o.left < o.right))
                            {
                                std::cerr << o.left << ' ' << o.right
                                          << std::endl;
                            }
                        //assert(o.left < o.right);
                        add_ancestry(parent_input_id, o.left, o.right,
                                     ancestry_node);
                        previous_right = o.right;
                    }
                if (is_sample && previous_right != L)
                    {
                        add_ancestry(parent_input_id, previous_right, L,
                                     output_id);
                    }
                if (output_id != -1)
                    {
                        auto n = squash_and_flush_edges(E);
                        if (!n && !is_sample)
                            {
                                new_node_table.erase(new_node_table.begin()
                                                         + output_id,
                                                     new_node_table.end());
                                idmap[parent_input_id] = -1;
                            }
                    }
            }

            int
            squash_and_flush_edges(
                //const std::map<std::int32_t, std::vector<edge>>&
                std::vector<edge>& buffered_edges)
            // Implementation copied from msprime.
            // Squashes identical edges on a per-parent
            // basis and adds them to the output list of edges.
            // Based on implementation found in msprime/lib/msprime.c
            // by Jerome Kelleher, version 0.5.0 of that code.
            // http://github.com/jeromekelleher/msprime.
            {
                std::stable_sort(buffered_edges.begin(), buffered_edges.end(),
                                 [](const edge& a, const edge& b) {
                                     return a.child < b.child;
                                 });
                new_edge_table.insert(new_edge_table.end(),
                                      buffered_edges.begin(),
                                      buffered_edges.end());
                return buffered_edges.size();
                //for (auto& i : buffered_edges)
                //    {
                //        for (auto& e : i.second)
                //            {
                //                new_edge_table.emplace_back(std::move(e));
                //                ++n;
                //            }
                //    }
                //return n;
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
                if (tables.mutation_table.empty())
                    // This skips index building, which
                    // is expensive an un-necessary if there's
                    // nothing to simplify...
                    {
                        return;
                    }

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

                tables.build_indexes(); // Expensive!!!
                std::vector<std::int32_t> remapped_samples(samples.size());
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        remapped_samples[i] = idmap[samples[i]];
                    }
                algorithmL(tables.input_left, tables.output_right,
                           remapped_samples, tables.node_table.size(),
                           tables.L, mutation_counter);
            }

          public:
            table_simplifier(const double maxpos)
                : new_edge_table{}, new_node_table{},
                  segment_queue{}, Ancestry{}, E{}, L{ maxpos }
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
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
                        new_node_table.emplace_back(
                            node{ tables.node_table[s].population,
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
    } // namespace ts
} // namespace fwdpp

#endif
