#ifndef FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP__
#define FWDPP_ANCESTRY_TABLE_SIMPLIFIER_HPP__

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
    namespace ancestry
    {
        class table_simplifier
        {
          private:
            struct segment
            {
                double left, right;
                std::int32_t node;
                segment(double l, double r, std::int32_t n) noexcept
                    : left{ l }, right{ r }, node{ n }
                {
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
                                        //tables_.push_back_node(
                                        //    static_cast<std::int32_t>(
                                        //        tables_.node_table.size()),
                                        //    tables.node_table[parent_input_id]
                                        //        .generation,
                                        //    tables.node_table[parent_input_id]
                                        //        .population);
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
            // by Jerome Kelleher.
            // http://github.com/jeromekelleher/msprime.
            //
            // Text of a rambling email sent to David Lawrie
            // on 19 April 2018:
            // Your comments yesterday prompted me to look at the code again.
            // I don't have time to actually do anything until June, but I
            // realized the following:
            //
            // The number of nodes that (immediately) descend from a node is
            // related to the final value of 'l' in squash_and_flush_edges.
            //
            // In other words, if you have a part of an edge table looking
            // like:
            //
            // 200 9 0 1
            // 200 119 0 1
            //
            // The final value of 'l' for all squashed edges with parent 200 is
            // 1, meaning that 200 leaves two records in the edge table.
            //
            // It is better than that, actually.  Below is a squashed record
            // for node 389 in a test sim.  These are edges (p/c/l/r), and
            // squashing reduces them to non-overlapping intervals.  Thus, the
            // data below can be processed to a lookup table of the form parent
            // -> left -> # immediate descendants, or parent -> (left,right) ->
            // # immediate descendants.  Building that info must help
            // simultaneously simplify and count mutations.
            //
            // 389 53 0.900531 0.990035
            // 389 375 0.517643 0.609341
            // 389 375 0.823912 1
            // 389 380 0.990035 1
            // 389 382 0 0.0121917
            // 389 382 0.823912 0.900531
            // 389 384 0 0.609341
            // 389 388 0.0121917 0.517643
            // 389 384 0 0.0121917
            // 389 384 0.0121917 0.517643
            // 389 384 0.517643 0.609341
            // 389 388 0.0121917 0.517643
            //
            //
            // That information means that you can skip a LOT of counting,
            // etc., which must simplify mutation counting tremendously.
            //
            // There are multiple approaches possible.
            //
            // 1. Convert the data here into (double,(parent,[children])),
            // where double is the left.  These entries can be placed
            // in vectors and sorted.  Essentially, this builds the
            // trees for all marginals as we go along, and it can be searched
            // quickly in the usual ways.
            //
            // 2. Implement algorithms L (and T?) from the msprime paper.
            // Jerome pointed out that they work perfectly well with
            // incompletely-coalesced trees.  I have to reread his
            // paper to see if that should be done here, or after
            // all of simplify is done with.
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

          public:
            table_simplifier(const double region_length)
                : new_edge_table{}, new_node_table{},
                  segment_queue{}, X{}, Ancestry{}, E{}, L{ region_length }
            {
            }

            //template <typename TC>
            //table_simplifier(TC&& initial_table_collection,
            //                 const double region_length)
            //    : tables{ std::forward<TC>(initial_table_collection) },
            //      tables_{}, segment_queue{}, X{}, Ancestry{}, E{},
            //      edge_offset{ static_cast<std::ptrdiff_t>(
            //          tables.edge_table.size()) },
            //      L{ region_length }
            //{
            //    if (!tables.edges_are_sorted())
            //        {
            //            throw std::invalid_argument("edges are not sorted");
            //        }
            //}

            template <typename mutation_container>
            //std::vector<std::int32_t>
            std::pair<std::vector<std::int32_t>, std::vector<std::uint32_t>>
            simplify(table_collection& tables,
                     const std::vector<std::int32_t>& samples,
                     const mutation_container& mutations)
            /// Set theoretic simplify.
            /// TODO: shorten via additional function calls
            /// for readability
            /// TODO: compare against implementation more
            /// closely matching what msprime is doing.
            /// TODO: make a template member so that mutations are
            /// also simplified
            {
                Ancestry.resize(tables.node_table.size());

                // Set some things up for later mutation simplification
                std::unordered_map<std::int32_t, std::vector<std::size_t>>
                    mutation_map; //maps input nodes to locations in input mut table
                std::vector<std::int32_t> mutation_node_map(
                    tables.node_table.size(), -1);
                for (auto& mr : tables.mutation_table)
                    {
                        std::cerr << "mut table record: " << mr.node << ' '
                                  << mutations[mr.key].pos << '\n';
                        mutation_map[mr.node].push_back(mr.key);
                    }

                for (auto& mm : mutation_map)
                    {
                        std::sort(mm.second.begin(), mm.second.end(),
                                  [&mutations, &tables](std::size_t a,
                                                        std::size_t b) {
                                      return mutations[a].pos
                                             < mutations[b].pos;
                                  });
                        std::cerr << "mut map data: " << mm.first << ' ';
                        for (auto& m : mm.second)
                            {
                                std::cerr << mutations[m].pos << ' ';
                            }
                        std::cerr << '\n';
                    }

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

                // Simplify mutations

                // 0. Remap mutation input node ids.  To do this, we use
                // the data stored in Ancestry, which allows us to "push"
                // mutation nodes down the tree.

                //for (std::size_t i = 0; i < tables.mutation_table.size(); ++i)
                for (auto& mm : mutation_map)
                    {
                        for (auto& e : Ancestry[mm.first])
                            {
                                std::cout << "Ancestry of " << mm.first << ' '
                                          << e.node << ' ' << e.left << ' '
                                          << e.right << '\n';
                            }
                        auto seg = Ancestry[mm.first].cbegin();
                        const auto seg_e = Ancestry[mm.first].cend();
                        auto mut = mm.second.cbegin();
                        const auto mute = mm.second.cend();
                        while (seg < seg_e && mut < mute)
                            {
                                //TODO: indirect access alert!
                                auto pos = mutations[*mut].pos;
                                if (seg->left <= pos && pos < seg->right)
                                    {
                                        assert(mut < mute);
                                        assert(*mut < mutations.size());
                                        assert(seg->node
                                               < new_edge_table.size());
                                        std::cout << "found it " << mm.first
                                                  << ' ' << pos << ' '
                                                  << seg->node << '\n';
                                        mutation_node_map.at(mm.first)
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

                // 2. Map input mutation node IDs to output IDs
                // This is fast O(n).
                for (std::size_t i = 0; i < tables.mutation_table.size(); ++i)
                    {
                        std::cout
                            << "mut remapping: "
                            << tables.mutation_table[i].node << " -> "
                            << mutation_node_map[tables.mutation_table[i].node]
                            << ' '
                            << idmap[mutation_node_map[tables.mutation_table[i]
                                                           .node]]
                            << '\n';
                        tables.mutation_table[i].node
                            = mutation_node_map[tables.mutation_table[i].node];
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
                for (auto& mr : tables.mutation_table)
                    {
                        std::cout << "output node and pos: " << mr.node << ' '
                                  << mutations[mr.key].pos << '\n';
                    }
                assert(std::is_sorted(
                    tables.mutation_table.begin(), tables.mutation_table.end(),
                    [&mutations](const mutation_record& a,
                                 const mutation_record& b) {
                        return mutations[a.key].pos < mutations[b.key].pos;
                    }));

                // 3. Now, we use Kelleher et al. (2016)'s Algorithm L
                // to march through each marginal tree and its leaf
                // counts. At the same time, we march through our mutation
                // table, which is sorted by position.
                std::vector<std::uint32_t> mcounts(
                    mutations.size(), 0); //TODO: update the real mcounts!!!

                auto mtable_itr = tables.mutation_table.begin();
                auto mtable_end = tables.mutation_table.end();
                auto mutation_simplifier = [&mutations, &mtable_itr,
                                            mtable_end, &mcounts](
                                               const marginal_tree& marginal) {
                    std::cerr << std::accumulate(marginal.leaf_counts.begin(),
                                                 marginal.leaf_counts.end(), 0)
                              << ' ' << marginal.left << ' ' << marginal.right
                              << ' ' << std::distance(mtable_itr, mtable_end)
                              << ' ';
                    //TODO: do we need to check left here?
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.left)
                        {
                            ++mtable_itr;
                        }
                    while (mtable_itr < mtable_end
                           && mutations[mtable_itr->key].pos < marginal.right)
                        {
                            std::cerr << mutations[mtable_itr->key].pos << " ("
                                      << mtable_itr->node << ','
                                      << marginal.leaf_counts[mtable_itr->node]
                                      << ','
                                      << marginal.parents[mtable_itr->node]
                                      << ") ";
                            assert(mutations[mtable_itr->key].pos
                                   >= marginal.left);
                            assert(mutations[mtable_itr->key].pos
                                   < marginal.right);
                            mcounts[mtable_itr->key]
                                = marginal.leaf_counts[mtable_itr->node];
                            ++mtable_itr;
                        }
                    std::cerr << std::distance(mtable_itr, mtable_end)
                              << std::endl;
                };

                tables.build_indexes();
                std::vector<std::int32_t> remapped_samples(samples.size());
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        remapped_samples[i] = idmap[samples[i]];
                    }
                algorithmL(tables.input_left, tables.output_right,
                           remapped_samples, tables.node_table.size(),
                           tables.L, mutation_simplifier);

                // 4. Remove any elements from table with count == 0
                tables.mutation_table.erase(
                    std::remove_if(tables.mutation_table.begin(),
                                   tables.mutation_table.end(),
                                   [&mcounts](const mutation_record& mr) {
                                       return mcounts[mr.key] == 0;
                                   }),
                    tables.mutation_table.end());

                cleanup();
                return std::make_pair(std::move(idmap), std::move(mcounts));
            }

            //void
            //sort_tables() noexcept
            ////TODO: this can be removed
            //{
            //    tables.sort_edges(edge_offset, E);
            //}

            //table_collection
            //dump_tables()
            ///// Returns the tables via a move-constructed object.
            ///// The table_simplifier instance is now in an inconsistent state.
            //{
            //    table_collection rv(std::move(tables));
            //    return rv;
            //}

            //const node_vector&
            //nodes() const
            //{
            //    return tables.node_table;
            //}

            //std::size_t
            //num_nodes() const
            //{
            //    return tables.node_table.size();
            //}

            //std::size_t
            //num_edges() const
            //{
            //    return tables.edge_table.size();
            //}

            //template <typename visitor>
            //void
            //algorithmL(visitor v, const std::vector<std::int32_t>& samples)
            //{
            //    tables.build_indexes();
            //    ancestry::algorithmL(tables.input_left, tables.output_right,
            //                         samples, tables.node_table.size(), 1.0,
            //                         v);
            //}

            //void
            //algorithmT()
            //{
            //    tables.build_indexes();
            //    ancestry::algorithmT(tables.input_left, tables.output_right,
            //                         tables.node_table.size(), 1.0);
            //}
        };
    } // namespace ancestry
} // namespace fwdpp

#endif
