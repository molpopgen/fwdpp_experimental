#ifndef FWDPP_ANCESTRY_ANCESTRY_TRACKER_HPP__
#define FWDPP_ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <vector>
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include "node.hpp"
#include "edge.hpp"
#include "table_collection.hpp"

namespace fwdpp
{
    namespace ancestry
    {
        class ancestry_tracker
        {
          private:
            struct segment
            {
                double left, right;
                std::int32_t node;
                segment(double l, double r, std::int32_t n) noexcept
                    : left{ l },
                      right{ r },
                      node{ n }
                {
                }
                //inline bool
                //operator>(const segment& rhs) const noexcept
                //{
                //    return left > rhs.left;
                //}
            };

            // tables is the current data set.
            // tables_ is used as temp
            // space during simplification.
            table_collection tables, tables_;
            // Q mimics a min queue of segments w.r.to
            // segment::left. X is a temporary container
            // for storing segments during ancestry merging.
            std::vector<segment> Q, X;
            std::vector<std::vector<segment>> Ancestry;
            /// Temp container used for compacting edges
            edge_vector E;
            // This reflects the length of
            // tables.edge_table after last simplification.
            // It can be used to make sure we only sort
            // newly-added nodes.
            std::ptrdiff_t edge_offset;
            // region length
            const double L;

            // std::pair<std::vector<std::pair<double, double>>,
            //          std::vector<std::pair<double, double>>>
            void
            split_breakpoints(
                const std::vector<double>& breakpoints,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const std::int32_t next_index)
            // TODO add data directly rather than create more intermediate
            // containers
            {
                std::vector<std::pair<double, double>> r1, r2;
                if (breakpoints.empty())
                    {
                        // r1.emplace_back(start, stop);
                        tables.push_back_edge(0., L, std::get<0>(parents),
                                              next_index);
                        // tables.emplace_back_edge(
                        //    start, stop, std::get<0>(parents), next_index);
                        goto out;
                    }
                if (breakpoints.front() != 0.0)
                    {
                        tables.push_back_edge(0., breakpoints.front(),
                                              std::get<0>(parents),
                                              next_index);
                        // tables.emplace_back_edge(start, breakpoints.front(),
                        //                         std::get<0>(parents),
                        //                         next_index);
                        // r1.emplace_back(
                        //    std::make_pair(start, breakpoints.front()));
                    }
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : L;
                        if (j % 2 == 0.)
                            {
                                tables.push_back_edge(
                                    a, b, std::get<0>(parents), next_index);
                                // r1.emplace_back(a, b);
                            }
                        else
                            {
                                tables.push_back_edge(
                                    a, b, std::get<1>(parents), next_index);
                                // r2.emplace_back(a, b);
                            }
                    }
            out:
                return; // std::make_pair(std::move(r1), std::move(r2));
            }

            bool
            add_to_queue(const double left, const double right,
                         const std::int32_t node, const bool added)
            {
                bool rv = added;
                if (!Q.empty() && left > Q.back().left)
                    {
                        return false;
                    }
                Q.emplace_back(left, right, node);
                return rv;
            }

            void
            sort_queue(std::size_t beg) noexcept
            {
                std::sort(Q.begin() + beg, Q.end(),
                          [](const segment& a, const segment& b) {
                              return a.left > b.left;
                          });
                assert(std::is_sorted(Q.begin(), Q.end(),
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
                tables_.clear();
                E.clear();
                for (auto& ai : Ancestry)
                    {
                        ai.clear();
                    }
                edge_offset = tables.edge_table.size();
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
                                        Q.emplace_back(
                                            std::max(seg.left, edge_ptr->left),
                                            std::min(seg.right,
                                                     edge_ptr->right),
                                            seg.node);
                                    }
                            }
                    }
                if (!Q.empty())
                    {
                        sort_queue(0);
                    }
                return edge_ptr;
            }

            void
            compact_tables()
            {
                std::size_t start = 0;
                E.swap(tables_.edge_table);
                assert(tables_.edge_table.empty());

                std::sort(
                    E.begin(), E.end(), [](const edge& a, const edge& b) {
                        return std::tie(a.parent, a.child, a.left, a.right)
                               < std::tie(b.parent, b.child, b.left, b.right);
                    });

                for (std::size_t j = 1; j < E.size(); ++j)
                    {
                        bool condition = E[j - 1].right != E[j].left
                                         || E[j - 1].parent != E[j].parent
                                         || E[j - 1].child != E[j].child;
                        if (condition)
                            {
                                tables_.push_back_edge(
                                    E[start].left, E[j - 1].right,
                                    E[j - 1].parent, E[j - 1].child);
                                start = j;
                            }
                    }
                auto j = E.size();
                tables_.push_back_edge(E[start].left, E[j - 1].right,
                                       E[j - 1].parent, E[j - 1].child);
            }

            void
            defragment(std::vector<segment>& ancestry_segment) noexcept
            {
                assert(std::is_sorted(ancestry_segment.begin(),
                                      ancestry_segment.end(),
                                      [](const segment& a, const segment& b) {
                                          return std::tie(a.left, a.right)
                                                 < std::tie(b.left, b.right);
                                      }));
                auto prev_seg = ancestry_segment.begin();
                auto next_seg = prev_seg + 1;
                //for(;next_seg<ancestry_segment.end();++next_seg)
                while (next_seg < ancestry_segment.end())
                    {
                        if (prev_seg->node == next_seg->node
                            && prev_seg->right == next_seg->left)
                            {
                                prev_seg->right = next_seg->right;
                                next_seg->node = -1;
                                next_seg += 2;
                                prev_seg += 2;
                            }
                        else
                            {
                                ++prev_seg;
                                ++next_seg;
                            }
                    }
                ancestry_segment.erase(
                    std::remove_if(
                        ancestry_segment.begin(), ancestry_segment.end(),
                        [](const segment& seg) { return seg.node == -1; }),
                    ancestry_segment.end());
            }

          public:
            ancestry_tracker(const std::int32_t num_initial_nodes,
                             const double initial_time, std::int32_t pop,
                             const double region_length = 1.0)
                : tables{ num_initial_nodes, initial_time, pop }, tables_{},
                  Q{}, X{}, Ancestry{}, E{}, edge_offset{ 0 },
                  L{ region_length }
            {
            }

            template <typename TC>
            ancestry_tracker(TC&& initial_table_collection,
                             const double region_length)
                : tables{ std::forward<TC>(initial_table_collection) },
                  tables_{}, Q{}, X{}, Ancestry{},
                  edge_offset{ static_cast<std::ptrdiff_t>(
                      tables.edge_table.size()) },
                  L{ region_length }
            {
                if (!tables.edges_are_sorted())
                    {
                        throw std::invalid_argument("edges are not sorted");
                    }
            }
            std::vector<std::int32_t>
            simplify(const std::vector<std::int32_t>& samples)
            /// Set theoretic simplify.
            /// TODO: shorten via additional function calls
            /// for readability
            /// TODO: compare against implementation more
            /// closely matching what msprime is doing.
            {
                Ancestry.resize(tables.node_table.size());

                // Relates input node ids to output node ids
                std::vector<std::int32_t> idmap(tables.node_table.size(), -1);

                // We take our samples and add them to both the output
                // node list and initialize their ancestry with
                // a segment on [0,L).
                for (auto& s : samples)
                    {
                        tables_.push_back_node(
                            tables_.node_table.size(),
                            tables.node_table[s].generation,
                            tables.node_table[s].population);
                        Ancestry[s].emplace_back(
                            0, L, static_cast<std::int32_t>(
                                      tables_.node_table.size() - 1));
                        idmap[s] = static_cast<std::int32_t>(
                            tables_.node_table.size() - 1);
                    }

                std::int32_t anode, znode = -1;
                double aleft, aright, zleft = 0, zright = 0;

                // At this point, our edges are sorted by birth
                // order of parents, from present to past.
                // We can now work our way up the pedigree.
                // This outer loop differs from how we describe it in the
                // paper, but the strict sorting of edges means that this
                // equivalent.
                auto edge_ptr = tables.edge_table.cbegin();
				const auto edge_end = tables.edge_table.cend();
                while (edge_ptr < edge_end )
                    {
                        auto u = edge_ptr->parent;
                        edge_ptr = step_S3(edge_ptr, edge_end, u);
                        std::int32_t v = -1;
                        //This is "stolen" straigh out of Jerome's
                        //code.  GPL ftw.
                        bool defrag_required = false;
                        std::size_t initial_Q_size = Q.size();
                        bool added_to_queue = false;
                        while (!Q.empty())
                            // Steps S4 through S8 of the algorithm.
                            {
                                X.clear();
                                auto l = Q.back().left;
                                double r = L;
                                while (!Q.empty() && Q.back().left == l)
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
                                        r = std::min(r, Q.back().right);
                                        X.emplace_back(std::move(Q.back()));
                                        assert(initial_Q_size > 0);
                                        Q.pop_back();
                                        --initial_Q_size;
                                    }
                                double next_left = 0.0;
                                if (!Q.empty())
                                    {
                                        next_left = Q.back().left;
                                        r = std::min(r, next_left);
                                    }
                                if (X.size() == 1)
                                    {
                                        if (!Q.empty()
                                            && next_left < X[0].right)
                                            {
                                                aleft = X[0].left;
                                                aright = next_left;
                                                anode = X[0].node;
                                                X[0].left = next_left;
                                                added_to_queue = add_to_queue(
                                                    next_left, X[0].right,
                                                    X[0].node, added_to_queue);
                                                ++initial_Q_size;
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
                                                tables_.push_back_node(
                                                    static_cast<std::int32_t>(
                                                        tables_.node_table
                                                            .size()),
                                                    tables.node_table[u]
                                                        .generation,
                                                    tables.node_table[u]
                                                        .population);
                                                v = tables_.node_table.size()
                                                    - 1;
                                                // update sample map
                                                idmap[u] = v;
                                            }
                                        aleft = l;
                                        aright = r;
                                        anode = v;
                                        for (auto& x : X)
                                            {
                                                tables_.push_back_edge(l, r, v,
                                                                       x.node);
                                                if (x.right > r)
                                                    {
                                                        x.left = r;
                                                        added_to_queue
                                                            = add_to_queue(
                                                                x.left,
                                                                x.right,
                                                                x.node,
                                                                added_to_queue);
                                                        ++initial_Q_size;
                                                    }
                                            }
                                    }
                                if (added_to_queue)
                                    {
                                        sort_queue(initial_Q_size);
                                        added_to_queue = false;
                                    }
                                Ancestry[u].emplace_back(aleft, aright, anode);
                                defrag_required
                                    |= zright == aleft && znode == anode;
                                zleft = aleft;
                                zright = aright;
                                znode = anode;
                            }
                        assert(initial_Q_size == 0);
                        if (defrag_required)
                            {
                                defragment(Ancestry[u]);
                            }
                    }

                assert(static_cast<std::size_t>(std::count_if(
                           idmap.begin(), idmap.end(),
                           [](const std::int32_t i) { return i != -1; }))
                       == tables_.node_table.size());

                // Now, we compact the edges,
                // which means removing redundant
                // info due to different edges
                // representing the same ancestry.
                // This is the same as ancestry chain
                // defragmenting, but applied to the final edge
                // collection.
                compact_tables();
                // TODO: allow for exception instead of assert
                assert(tables.edges_are_sorted());
                tables.swap(tables_);
                cleanup();
                return idmap;
            }

            void
            add_offspring_data(
                const std::int32_t next_index,
                const std::vector<double>& breakpoints,
                const std::vector<std::uint32_t>& new_mutations,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const double generation)
            {
                // TODO document why this is generation + 1
                tables.emplace_back_node(next_index, generation + 1, 0);
                // auto split =
                split_breakpoints(breakpoints, parents, next_index);
                // Add the edges
                // for (auto&& brk : split.first)
                //    {
                //        tables.emplace_back_edge(brk.first, brk.second,
                //                                 std::get<0>(parents),
                //                                 next_index);
                //    }
                // for (auto&& brk : split.second)
                //    {
                //        tables.emplace_back_edge(brk.first, brk.second,
                //                                 std::get<1>(parents),
                //                                 next_index);
                //    }
                // TODO: add mutation to tables
            }

            void
            sort_tables() noexcept
            {
                tables.sort_edges();
            }

            table_collection
            dump_tables()
            /// Returns the tables via a move-constructed objec.
            /// The ancestry_tracker instance is now in an inconsistent state.
            {
                table_collection rv(std::move(tables));
                return rv;
            }

            const node_vector&
            nodes() const
            {
                return tables.node_table;
            }

            std::size_t
            num_nodes() const
            {
                return tables.node_table.size();
            }

            std::size_t
            num_edges() const
            {
                return tables.edge_table.size();
            }
        };
    }
}

#endif
