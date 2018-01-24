#ifndef FWDPP_ANCESTRY_ANCESTRY_TRACKER_HPP__
#define FWDPP_ANCESTRY_ANCESTRY_TRACKER_HPP__

#include <vector>
#include <algorithm>
#include <cstddef>
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
            };

            // tables is the current data set.
            // tables_ is used as temp
            // space during simplification.
            table_collection tables, tables_;
            // Q mimics a min-queue, and X
            // is a temp vector for segments
            // while processing Q.  The sorting
			// of Q is handled by 
			// ancestry_tracker::sort_queue.
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

            bool
            sort_queue(bool added2Q) noexcept
            // Sorts the priority queue during
            // simplify as needed.
            {
                if (added2Q)
                    {
                        std::sort(Q.begin(), Q.end(),
                                  [](const segment& a, const segment& b) {
                                      return a.left > b.left;
                                  });
                    }
                return false;
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

          public:
            ancestry_tracker(const std::int32_t num_initial_nodes,
                             const double initial_time, std::int32_t pop,
                             const double region_length = 1.0)
                : tables{ num_initial_nodes, initial_time, pop }, tables_{},
                  Q{}, X{}, Ancestry{}, E{}, edge_offset{ 0 },
                  L{ region_length }
            {
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

                std::int32_t anode;
                double aleft, aright;
                bool added2Q = false;

                // At this point, our edges are sorted by birth
                // order of parents, from present to past.
                // We can now work our way up the pedigree.
                // This outer loop differs from how we describe it in the
                // paper, but the strict sorting of edges means that this
                // equivalent.
                auto edge_ptr = tables.edge_table.cbegin();
                while (edge_ptr < tables.edge_table.cend())
                    {
                        auto u = edge_ptr->parent;
                        for (; edge_ptr < tables.edge_table.end()
                               && edge_ptr->parent == u;
                             ++edge_ptr)
                            {
                                //For each edge corresponding to this parent,
                                //we look at all segments from the child.
                                //If the two segments overlap, we add the minimal
                                //overlap to our queue.
                                //This is Step S3.
								//TODO: can we do binary lookup here?
                                for (auto& seg : Ancestry[edge_ptr->child])
                                    {
                                        if (seg.right > edge_ptr->left
                                            && edge_ptr->right > seg.left)
                                            {
                                                Q.emplace_back(
                                                    std::max(seg.left,
                                                             edge_ptr->left),
                                                    std::min(seg.right,
                                                             edge_ptr->right),
                                                    seg.node);
                                                added2Q = true;
                                            }
                                    }
                            }
                        added2Q = sort_queue(added2Q);
                        std::int32_t v = -1;
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
                                    // TODO: this can be done with reverse iteration,
                                    // but testing on 0.5e9 edges didn't seem to
                                    // make it worthwhile.
                                    {
                                        r = std::min(r, Q.back().right);
                                        X.emplace_back(std::move(Q.back()));
                                        Q.pop_back();
                                    }
                                if (!Q.empty())
                                    {
                                        r = std::min(r, Q.back().left);
                                    }
                                if (X.size() == 1)
                                    {
                                        if (!Q.empty()
                                            && Q.back().left < X[0].right)
                                            {
                                                aleft = X[0].left;
                                                aright = Q.back().left;
                                                anode = X[0].node;
                                                Q.emplace_back(Q.back().left,
                                                               X[0].right,
                                                               X[0].node);
                                                added2Q = true;
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
                                                // Overlap/coalescence, and thus
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
                                                        Q.emplace_back(x.left,
                                                                       x.right,
                                                                       x.node);
                                                        added2Q = true;
                                                    }
                                            }
                                    }
                                added2Q = sort_queue(added2Q);
                                Ancestry[u].emplace_back(aleft, aright, anode);
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
                tables.swap(tables);
                //TODO: allow for exception instead of assert
                assert(tables.edges_are_sorted());
                cleanup();
                return idmap;
            }

            void
            sort_tables() noexcept
            {
                tables.sort_edges(edge_offset);
            }
        };
    }
}

#endif
