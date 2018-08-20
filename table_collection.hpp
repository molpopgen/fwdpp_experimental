#ifndef FWDPP_ANCESTRY_TABLE_COLLECTION_HPP
#define FWDPP_ANCESTRY_TABLE_COLLECTION_HPP

#include <vector>
#include <utility>
#include <cstdint>
#include <cmath>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include "node.hpp"
#include "edge.hpp"
#include "mutation_record.hpp"
#include "indexed_edge.hpp" //TODO: create fewer header dependencies

namespace fwdpp
{
    namespace ts
    {
        using edge_vector = std::vector<edge>;
        using node_vector = std::vector<node>;
        using mutation_key_vector = std::vector<mutation_record>;
        struct table_collection
        {
          private:
            edge_vector temp_edges; //used for sorting
            void
            split_breakpoints(
                const std::vector<double>& breakpoints,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const std::int32_t next_index)
            {
                std::vector<std::pair<double, double>> r1, r2;
                if (breakpoints.empty())
                    {
                        this->push_back_edge(0., L, std::get<0>(parents),
                                             next_index);
                        return;
                    }
                if (breakpoints.front() != 0.0)
                    {
                        this->push_back_edge(0., breakpoints.front(),
                                             std::get<0>(parents), next_index);
                    }
                // TODO: decide how to handle identical positions.
                // If a breakpoint is repeated 2x, if has no effect
                // on the genealogy, as it is a double x-over.  But,
                // if there are an odd number, then it does have an effect.
                // The iteration here needs to be updated to handle this.
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : L;
                        if (b <= a)
                            {
                                throw std::runtime_error(
                                    "right must be > left");
                            }
                        if (j % 2 == 0.)
                            {
                                this->push_back_edge(
                                    a, b, std::get<0>(parents), next_index);
                            }
                        else
                            {
                                this->push_back_edge(
                                    a, b, std::get<1>(parents), next_index);
                            }
                    }
            }

          public:
            node_vector node_table;
            edge_vector edge_table;
            mutation_key_vector mutation_table;
            indexed_edge_container input_left, output_right;
            // This reflects the length of
            // tables.edge_table after last simplification.
            // It can be used to make sure we only sort
            // newly-added nodes.
            std::ptrdiff_t edge_offset;
            const double L;
            table_collection(const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos }
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
                    }
            }

            table_collection(const std::int32_t num_initial_nodes,
                             const double initial_time, std::int32_t pop,
                             const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos }
            {
                if (maxpos < 0 || !std::isfinite(maxpos))
                    {
                        throw std::invalid_argument(
                            "maxpos must be > 0 and finite");
                    }
                for (std::int32_t i = 0; i < num_initial_nodes; ++i)
                    {
                        node_table.push_back(node{ pop, initial_time });
                    }
            }

            void
            sort_edges() noexcept
            /// Sort the edge table.  On PARENT birth times.
            /// The sorting differs from msprime here. The difference
            /// is that we  assume that birth times are recorded forward in
            /// time rather than backwards.
            {
                std::sort(edge_table.begin() + edge_offset, edge_table.end(),
                          [this](const edge& a, const edge& b) {
                              auto ga = this->node_table[a.parent].generation;
                              auto gb = this->node_table[b.parent].generation;
                              if (ga == gb)
                                  {
                                      if (a.parent == b.parent)
                                          {
                                              if (a.child == b.child)
                                                  {
                                                      return a.left < b.left;
                                                  }
                                              return a.child < b.child;
                                          }
                                      return a.parent < b.parent;
                                  }
                              return ga > gb;
                          });
                if (edge_offset > 0)
                    {
                        temp_edges.reserve(edge_table.size());
#ifndef NDEBUG
                        auto size = edge_table.size();
#endif
                        temp_edges.clear();
                        temp_edges.insert(
                            temp_edges.end(),
                            std::make_move_iterator(edge_table.begin()
                                                    + edge_offset),
                            std::make_move_iterator(edge_table.end()));
                        temp_edges.insert(
                            temp_edges.end(),
                            std::make_move_iterator(edge_table.begin()),
                            std::make_move_iterator(edge_table.begin()
                                                    + edge_offset));
                        assert(temp_edges.size() == size);
                        temp_edges.swap(edge_table);
                    }
                temp_edges.clear();
                // TODO: allow for exceptions
                // rather than assertions.
                assert(edges_are_sorted());
            }

            template <typename mutation_container>
            void
            sort_mutations(const mutation_container& mutations)
            {
                //mutations are sorted by increasing position
                std::sort(mutation_table.begin(), mutation_table.end(),
                          [&mutations](const mutation_record& a,
                                       const mutation_record& b) {
                              return mutations[a.key].pos
                                     < mutations[b.key].pos;
                          });
            }

            template <typename mutation_container>
            void
            sort_tables(const mutation_container& mutations)
            /// Sorts the tables
            /// Note that mutations can be mocked via any struct
            /// containing double pos
            {
                sort_edges();
                sort_mutations(mutations);
            }

            bool
            edges_are_sorted() const noexcept
            {
                return std::is_sorted(
                    edge_table.begin(), edge_table.end(),
                    [this](const edge& a, const edge& b) {
                        auto ga = this->node_table[a.parent].generation;
                        auto gb = this->node_table[b.parent].generation;
                        return ga > gb
                               || (ga == gb
                                   && std::tie(a.parent, a.child, a.left)
                                          < std::tie(b.parent, b.child,
                                                     b.left));
                    });
            }

            void
            clear() noexcept
            /// Clears internal vectors.
            /// Mostly used during simplification
            /// where a table_collection is
            /// used as a temp object.
            {
                node_table.clear();
                edge_table.clear();
                mutation_table.clear();
            }

            void
            push_back_node(double generation, std::int32_t pop)
            {
                node_table.push_back(node{ pop, generation });
            }

            template <typename... args>
            void
            emplace_back_node(args&&... Args)
            {
                node_table.emplace_back(node{ std::forward<args>(Args)... });
            }

            void
            push_back_edge(double l, double r, std::int32_t parent,
                           std::int32_t child)
            {
                edge_table.push_back(edge{ l, r, parent, child });
            }

            template <typename... args>
            void
            emplace_back_edge(args&&... Args)
            {
                edge_table.emplace_back(edge{ std::forward<args>(Args)... });
            }

            void
            build_indexes()
            /// Generates the index vectors referred to
            /// as I and O in Kelleher et al. (2016)
            {
                input_left.reserve(edge_table.size());
                output_right.reserve(edge_table.size());
                input_left.clear();
                output_right.clear();
                for (auto& e : edge_table)
                    {
                        assert(e.left < e.right);
                        input_left.emplace_back(
                            e.left, -node_table[e.parent].generation, e.parent,
                            e.child);
                        output_right.emplace_back(
                            e.right, node_table[e.parent].generation, e.parent,
                            e.child);
                    }
                std::sort(input_left.begin(), input_left.end());
                std::sort(output_right.begin(), output_right.end());
            }

            void
            add_offspring_data(
                const std::int32_t next_index,
                const std::vector<double>& breakpoints,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const double generation)
            {
                emplace_back_node(0, generation + 1);
                split_breakpoints(breakpoints, parents, next_index);
            }

            void
            add_offspring_data(
                const std::int32_t next_index,
                const std::vector<double>& breakpoints,
                const std::vector<std::uint32_t>& new_mutations,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const double generation)
            {
                add_offspring_data(next_index, breakpoints, parents,
                                   generation);
                for (auto& m : new_mutations)
                    {
                        mutation_table.emplace_back(
                            mutation_record{ next_index, m });
                        assert(mutation_table.back().node == next_index);
                        assert(mutation_table.back().key == m);
                    }
            }

            std::size_t
            num_nodes() const
            {
                return node_table.size();
            }

            void
            update_offset()
            {
                edge_offset = edge_table.size();
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
