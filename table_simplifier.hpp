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

            struct mutation_node_map_entry
            {
                std::int32_t node;
                std::size_t key, location;
                mutation_node_map_entry(std::int32_t n, std::size_t k,
                                        std::size_t l)
                    : node(n), key(k), location(l)
                {
                }
            };

            class segment_overlapper
            /// This class is an iterable object
            /// over [left, right) -> segment
            /// mappings, where the segments
            /// are the genomic intervals in
            /// child nodes overlapping with
            /// the current parent in the
            /// current genomic interval.
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
                }

              public:
                std::vector<segment> overlapping;
                std::vector<segment>::iterator overlapping_end;
                double left, right;
                segment_overlapper()
                    : sbeg(), send(), overlapping{},
                      overlapping_end(overlapping.end()), left(0),
                      right(std::numeric_limits<double>::max())
                {
                }

                void
                init(std::vector<segment>& segs)
                {
                    sbeg = segs.begin();
                    // The - 1 for send assumes a "cap"/sentinel value.
                    send = segs.end() - 1;
                    overlapping.clear();
                    overlapping_end = overlapping.end();
                    left = 0.0;
                    right = std::numeric_limits<double>::max();
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
            segment_overlapper o;
            std::vector<mutation_node_map_entry> mutation_map;

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
                segment_queue.emplace_back(segment{ L, L + 1.0, -1 });
                return edge_ptr;
            }

            void
            buffer_edge(std::vector<edge>& buffered_edges, const double left,
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
            }

            void
            add_ancestry(std::int32_t input_id, double left, double right,
                         std::int32_t node)
            {
                if (Ancestry[input_id].empty())
                    {
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
                o.init(segment_queue);
                std::int32_t ancestry_node = -1;
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
                                add_ancestry(parent_input_id, previous_right,
                                             o.left, output_id);
                            }
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
                        auto n = output_buffered_edges(E);
                        if (!n && !is_sample)
                            {
                                new_node_table.erase(new_node_table.begin()
                                                         + output_id,
                                                     new_node_table.end());
                                idmap[parent_input_id] = -1;
                            }
                    }
            }

            std::size_t
            output_buffered_edges(std::vector<edge>& buffered_edges)
            /// Take our buffered edges and add them to the output edge table
            {
                std::stable_sort(buffered_edges.begin(), buffered_edges.end(),
                                 [](const edge& a, const edge& b) {
                                     return a.child < b.child;
                                 });
                new_edge_table.insert(new_edge_table.end(),
                                      buffered_edges.begin(),
                                      buffered_edges.end());
                return buffered_edges.size();
            }

            template <typename mcont_t>
            void
            prep_mutation_simplification(
                const mcont_t& mutations,
                const mutation_key_vector& mutation_table)
            {
                mutation_map.clear();
                mutation_map.reserve(mutation_table.size());
                for (std::size_t i = 0; i < mutation_table.size(); ++i)
                    {
                        mutation_map.emplace_back(mutation_table[i].node,
                                                  mutation_table[i].key, i);
                    }

                std::sort(mutation_map.begin(), mutation_map.end(),
                          [&mutations](const mutation_node_map_entry& a,
                                       const mutation_node_map_entry& b) {
                              return std::tie(a.node, mutations[a.key].pos)
                                     < std::tie(b.node, mutations[b.key].pos);
                          });
            }

            template <typename mcont_t>
            void
            simplify_mutations(const mcont_t& mutations,
                               mutation_key_vector& mt) const
            // Remove all mutations that do not map to nodes
            // in the simplified tree.  The key here is
            // that Ancestry contains the history of
            // each node, which we use for the remapping.
            {
                // Set all output nodes to null for now.
                for (auto& mr : mt)
                    {
                        mr.node = -1;
                    }

                // Map the input node id of a mutation to
                // its output node id.  If no output ID exists,
                // then the mutation will be removed by the
                // call to erase below.
                auto map_itr = mutation_map.begin();
                const auto map_end = mutation_map.end();

                while (map_itr < map_end)
                    {
                        auto n = map_itr->node;
                        auto seg = Ancestry[n].cbegin();
                        const auto seg_e = Ancestry[n].cend();
                        for (; map_itr < map_end
                               && map_itr->node == n;) //++map_itr)
                            {
                                if (seg == seg_e)
                                    {
                                        ++map_itr;
                                        break;
                                    }
                                while (seg < seg_e && map_itr < map_end
                                       && map_itr->node == n)
                                    {
                                        auto pos = mutations[map_itr->key].pos;
                                        if (seg->left <= pos
                                            && pos < seg->right)
                                            {
                                                mt[map_itr->location].node
                                                    = seg->node;
                                                ++map_itr;
                                            }
                                        else if (pos >= seg->right)
                                            {
                                                ++seg;
                                            }
                                        else
                                            {
                                                ++map_itr;
                                            }
                                    }
                            }
                    }

                // Any mutations with null node values do not have
                // ancestry and may be removed.
                mt.erase(std::remove_if(mt.begin(), mt.end(),
                                        [](const mutation_record& mr) {
                                            return mr.node == -1;
                                        }),
                         mt.end());
                //TODO: replace assert with exception
                assert(std::is_sorted(mt.begin(), mt.end(),
                                      [&mutations](const mutation_record& a,
                                                   const mutation_record& b) {
                                          return mutations[a.key].pos
                                                 < mutations[b.key].pos;
                                      }));
            }

            void
            record_sample_nodes(const std::vector<std::int32_t>& samples,
                                const table_collection& tables,
                                std::vector<std::int32_t>& idmap)
            {
                for (const auto& s : samples)
                    {
                        new_node_table.emplace_back(
                            node{ tables.node_table[s].population,
                                  tables.node_table[s].generation });
                        add_ancestry(s, 0, L,
                                     static_cast<std::int32_t>(
                                         new_node_table.size() - 1));
                        idmap[s] = static_cast<std::int32_t>(
                            new_node_table.size() - 1);
                    }
            }

          public:
            table_simplifier(const double maxpos)
                : new_edge_table{}, new_node_table{}, segment_queue{},
                  Ancestry{}, E{}, L{ maxpos }, o{}, mutation_map{}
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
                     const mutation_container& mutations)
            /// Simplify algorithm is approximately the same
            /// logic as used in msprime 0.6.0
            ///
            /// \param tables A table_collection
            /// \param samples A list of sample (node) ids.
            /// \param mutations A container of mutations
            {
                Ancestry.resize(tables.node_table.size());

                // Set some things up for later mutation simplification
                prep_mutation_simplification(mutations, tables.mutation_table);

                // Relates input node ids to output node ids
                std::vector<std::int32_t> idmap(tables.node_table.size(), -1);

                // We take our samples and add them to both the output
                // node list and initialize their ancestry with
                // a segment on [0,L).
                record_sample_nodes(samples, tables, idmap);
                // Add samples for any preserved nodes in the tables:
                record_sample_nodes(tables.preserved_nodes, tables, idmap);

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

                // When there are preserved nodes, we need to re map
                // their input ids to output ids
                for (auto& p : tables.preserved_nodes)
                    {
                        if (idmap[p] == -1)
                            {
                                throw std::runtime_error(
                                    "preserved node output id maps to null");
                            }
                        p = idmap[p];
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
                simplify_mutations(mutations, tables.mutation_table);

                cleanup();
                return idmap;
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
