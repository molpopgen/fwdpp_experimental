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

namespace fwdpp
{
    namespace ancestry
    {
        using edge_vector = std::vector<edge>;
        using node_vector = std::vector<node>;
        using mutation_key_vector
            = std::vector<std::pair<std::int32_t, std::size_t>>;
        struct table_collection
        {
          public:
            node_vector node_table;
            edge_vector edge_table;
            mutation_key_vector mutation_table;

            table_collection() : node_table{}, edge_table{}, mutation_table{}
            {
            }

            table_collection(const std::int32_t num_initial_nodes,
                             const double initial_time, std::int32_t pop)
                : node_table{}, edge_table{}, mutation_table{}
            {
                for (std::int32_t i = 0; i < num_initial_nodes; ++i)
                    {
                        node_table.push_back(node(i, initial_time, pop));
                    }
            }

            void
            sort_edges(const std::size_t offset)
            {
                edge_vector temp_edges;
                if (offset > 0)
                    {
                        temp_edges.reserve(edge_table.size());
                    }
                sort_edges(offset, temp_edges);
            }
            void
            sort_edges(const std::size_t offset,
                       edge_vector& temp_edges) noexcept
            /// Sort the edge table.  On PARENT birth times.
            /// The sorting differs from msprime here. The difference
            /// is that we  assume that birth times are recorded forward in
            /// time rather than backwards.
            /// TODO: need offset
            {
                std::sort(edge_table.begin() + offset, edge_table.end(),
                          [this](const edge& a, const edge& b) {
                              auto ga = this->node_table[a.parent].generation;
                              auto gb = this->node_table[b.parent].generation;
                              return ga > gb
                                     || (ga == gb
                                         && std::tie(a.parent, a.child, a.left)
                                                < std::tie(b.parent, b.child,
                                                           b.left));
                          });
                if (offset > 0)
                    {
                        auto size = edge_table.size();
                        temp_edges.clear();
                        temp_edges.insert(
                            temp_edges.end(), std::make_move_iterator(
                                                  edge_table.begin() + offset),
                            std::make_move_iterator(edge_table.end()));
                        temp_edges.insert(
                            temp_edges.end(),
                            std::make_move_iterator(edge_table.begin()),
                            std::make_move_iterator(edge_table.begin()
                                                    + offset));
                        //std::move(edge_table.begin() + offset,
                        //          edge_table.end(),
                        //          std::back_inserter(temp_edges));
                        //std::move(edge_table.begin(), edge_table.begin()+offset,
                        //          std::back_inserter(temp_edges));
                        assert(temp_edges.size() == size);
                        temp_edges.swap(edge_table);
						temp_edges.clear();
                    }
                // TODO: allow for exceptions
                // rather than assertions.
                assert(edges_are_sorted());
            }

            void
            swap(table_collection& t) noexcept
            /// Swaps out data members.
            /// Primary use is after simplification,
            /// where a table_collection is used
            /// as a temp object.
            {
                edge_table.swap(t.edge_table);
                node_table.swap(t.node_table);
                mutation_table.swap(t.mutation_table);
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
            push_back_node(std::int32_t id, double generation,
                           std::int32_t pop)
            {
                node_table.push_back(node(id, generation, pop));
            }

            template <typename... args>
            void
            emplace_back_node(args&&... Args)
            {
                node_table.emplace_back(std::forward<args>(Args)...);
            }

            void
            push_back_edge(double l, double r, std::int32_t parent,
                           std::int32_t child)
            {
                edge_table.push_back(edge(l, r, parent, child));
            }

            template <typename... args>
            void
            emplace_back_edge(args&&... Args)
            {
                edge_table.emplace_back(std::forward<args>(Args)...);
            }
        };
    }
}

#endif
