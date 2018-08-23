#ifndef FWDPP_TS_DATA_MATRIX_GENERATOR_HPP
#define FWDPP_TS_DATA_MATRIX_GENERATOR_HPP

#include <vector>
#include <cstdint>
#include <stdexcept>
#include "table_collection.hpp"
#include "variant_filler.hpp"

namespace fwdpp
{
    namespace ts
    {
        // TODO: make this conform to
        // fwdpp::data_matrix specs
        struct genotype_matrix
        {
            std::vector<std::int8_t> genotypes;
            std::vector<double> positions;
            std::vector<std::size_t> mut_indexes;
            template <typename G, typename P, typename I>
            genotype_matrix(G&& g, P&& p, I&& i)
                : genotypes(std::forward<G>(g)), positions(std::forward<P>(p)),
                  mut_indexes(std::forward<I>(i))
            {
            }
        };

        //TODO: hide in internal namespace + separate
        struct create_data_matrix_details
        {
            std::vector<std::int8_t> neutral_genotypes, selected_genotypes;
            std::vector<double> neutral_positions, selected_positions;
            std::vector<std::size_t> neutral_indexes, selected_indexes;
            variant_filler vf;
            mutation_key_vector::const_iterator beg, end;
            create_data_matrix_details(
                mutation_key_vector::const_iterator b,
                mutation_key_vector::const_iterator e,
                const std::vector<std::int32_t>& samples,
                const std::int32_t nnodes)
                : neutral_genotypes(), selected_genotypes(),
                  neutral_positions(), selected_positions(), neutral_indexes(),
                  selected_indexes(), vf(nnodes, samples), beg(b), end(e)
            {
                neutral_positions.reserve(e - b);
                selected_positions.reserve(e - b);
                neutral_indexes.reserve(e - b);
                selected_indexes.reserve(e - b);
                neutral_genotypes.reserve(samples.size());
                selected_genotypes.reserve(samples.size());
            }

            template <typename mcont_t>
            inline void
            operator()(const marginal_tree& marginal, const mcont_t& mutations,
                       const bool record_neutral, const bool record_selected)
            {
                while (beg < end && mutations[beg->key].pos < marginal.left)
                    {
                        ++beg;
                    }
                while (beg < end && mutations[beg->key].pos < marginal.right)
                    {
                        auto nv = vf(marginal, beg->node);
                        if (nv)
                            {
                                auto d = vf.view_genotypes();
                                //assert(nv == std::count(d.first,d.first+d.second,1));
                                bool n = mutations[beg->key].neutral;
                                if (n && record_neutral)
                                    {
                                        neutral_genotypes.insert(
                                            neutral_genotypes.end(), d.first,
                                            d.second);
                                        neutral_positions.push_back(
                                            mutations[beg->key].pos);
                                        neutral_indexes.push_back(beg->key);
                                    }
                                else if (!n && record_selected)
                                    {
                                        selected_genotypes.insert(
                                            selected_genotypes.end(), d.first,
                                            d.second);
                                        selected_positions.push_back(
                                            mutations[beg->key].pos);
                                        selected_indexes.push_back(beg->key);
                                    }
                            }
                        ++beg;
                    }
            }
        };

        //TODO: allow for better return value
        template <typename mcont_t>
        std::pair<genotype_matrix, genotype_matrix>
        create_data_matrix(const mcont_t& mutations,
                           const table_collection& tables,
                           const std::vector<std::int32_t>& samples,
                           const bool record_neutral,
                           const bool record_selected)
        {
            if (!record_neutral && !record_selected)
                {
                    throw std::invalid_argument(
                        "must record neutral and/or selected genotypes");
                }
            create_data_matrix_details d(tables.mutation_table.begin(),
                                         tables.mutation_table.end(), samples,
                                         tables.node_table.size());
            auto visitor = [&d, &mutations, record_neutral,
                            record_selected](const marginal_tree& marginal) {
                d(marginal, mutations, record_neutral, record_selected);
            };
            algorithmS(tables.input_left, tables.output_right, samples,
                       tables.node_table.size(), tables.L, visitor);
            return std::make_pair(
                genotype_matrix(std::move(d.neutral_genotypes),
                                std::move(d.neutral_positions),
                                std::move(d.neutral_indexes)),
                genotype_matrix(std::move(d.selected_genotypes),
                                std::move(d.selected_positions),
                                std::move(d.selected_indexes)));
        }
    } // namespace ts
} // namespace fwdpp

#endif
