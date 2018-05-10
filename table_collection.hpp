#ifndef FWDPP_ANCESTRY_TABLE_COLLECTION_HPP__
#define FWDPP_ANCESTRY_TABLE_COLLECTION_HPP__

#include <vector>
#include <utility>
#include <cstdint>
#include <tuple>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include "node.hpp"
#include "edge.hpp"
#include "mutation_record.hpp"
#include "msprime_algo.hpp" //TODO: create fewer header dependencies

namespace fwdpp
{
    namespace ancestry
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
                        goto out;
                    }
                if (breakpoints.front() != 0.0)
                    {
                        this->push_back_edge(0., breakpoints.front(),
                                             std::get<0>(parents), next_index);
                    }
                for (unsigned j = 1; j < breakpoints.size(); ++j)
                    {
                        double a = breakpoints[j - 1];
                        double b = (j < breakpoints.size() - 1)
                                       ? breakpoints[j]
                                       : L;
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
            out:
                return;
            }

          public:
            node_vector node_table;
            edge_vector edge_table;
            mutation_key_vector mutation_table;
            index_vector input_left, output_right;
            // This reflects the length of
            // tables.edge_table after last simplification.
            // It can be used to make sure we only sort
            // newly-added nodes.
            // TODO: move to table_collection
            std::ptrdiff_t edge_offset;
            const double L;
            table_collection(const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos }
            {
                //TODO assert maxpos is > 0 and finite
            }

            table_collection(const std::int32_t num_initial_nodes,
                             const double initial_time, std::int32_t pop,
                             const double maxpos)
                : temp_edges{}, node_table{}, edge_table{}, mutation_table{},
                  input_left{}, output_right{}, edge_offset{ 0 }, L{ maxpos }
            {
                //TODO assert maxpos is > 0 and finite
                for (std::int32_t i = 0; i < num_initial_nodes; ++i)
                    {
                        node_table.push_back(node{ i, pop, initial_time });
                    }
            }

            //void
            //sort_edges()
            //{
            //    edge_vector temp_edges;
            //    if (edge_offset > 0)
            //        {
            //            temp_edges.reserve(edge_table.size());
            //        }
            //    sort_edges(edge_offset, temp_edges);
            //}

            void
            sort_edges() noexcept
            /// Sort the edge table.  On PARENT birth times.
            /// The sorting differs from msprime here. The difference
            /// is that we  assume that birth times are recorded forward in
            /// time rather than backwards.
            /// TODO: need offset
            {
                std::sort(edge_table.begin() + edge_offset, edge_table.end(),
                          [this](const edge& a, const edge& b) {
                              auto ga = this->node_table[a.parent].generation;
                              auto gb = this->node_table[b.parent].generation;
                              //
                              //return ga > gb
                              //       || (ga == gb
                              //           && std::tie(a.parent, a.child, a.left)
                              //                  < std::tie(b.parent, b.child,
                              //                             b.left));
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
                        auto size = edge_table.size();
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
                        // std::move(edge_table.begin() + offset,
                        //          edge_table.end(),
                        //          std::back_inserter(temp_edges));
                        // std::move(edge_table.begin(),
                        // edge_table.begin()+offset,
                        //          std::back_inserter(temp_edges));
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
            sort_tables(const mutation_container& mutations)
            //TODO: static_assert that mutations contains the right stuff
            {
                sort_edges();
                //mutations are sorted by increasing position
                std::sort(mutation_table.begin(), mutation_table.end(),
                          [&mutations](const mutation_record& a,
                                       const mutation_record& b) {
                              return mutations[a.key].pos
                                     < mutations[b.key].pos;
                          });
            }

            //void
            //swap(table_collection& t) noexcept
            ///// Swaps out data members.
            ///// Primary use is after simplification,
            ///// where a table_collection is used
            ///// as a temp object.
            //{
            //    edge_table.swap(t.edge_table);
            //    node_table.swap(t.node_table);
            //    mutation_table.swap(t.mutation_table);
            //}

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
            push_back_node(std::int32_t id, double generation,
                           std::int32_t pop)
            {
                node_table.push_back(node{ id, pop, generation });
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
                const std::vector<std::uint32_t>& new_mutations,
                const std::tuple<std::int32_t, std::int32_t>& parents,
                const double generation)
            //TODO: this must move to table_collection
            {
                // TODO document why this is generation + 1
                emplace_back_node(next_index, 0, generation + 1);
                // auto split =
                split_breakpoints(breakpoints, parents, next_index);
                for (auto& m : new_mutations)
                    {
                        mutation_table.emplace_back(
                            mutation_record{ next_index, m });
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
    }
}

#endif
