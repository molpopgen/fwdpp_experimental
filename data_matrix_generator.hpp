#ifndef FWDPP_TS_DATA_MATRIX_GENERATOR_HPP
#define FWDPP_TS_DATA_MATRIX_GENERATOR_HPP

#include "table_collection.hpp"
#include "variant_filler.hpp"

namespace fwdpp
{
    namespace ts
    {
        //TODO: hide in internal namespace + separate
        struct create_data_matrix_details
        {
            std::vector<std::int8_t> genotypes;
            std::vector<double> positions;
            variant_filler vf;
            mutation_key_vector::const_iterator beg, end;
            create_data_matrix_details(
                mutation_key_vector::const_iterator b,
                mutation_key_vector::const_iterator e,
                const std::vector<std::int32_t>& samples,
                const std::int32_t nnodes)
                : genotypes(), positions(), vf(nnodes, samples), beg(b), end(e)
            {
                positions.reserve(e - b);
                genotypes.reserve((e - b) * samples.size());
            }
            template <typename mcont_t>
            inline void
            operator()(const marginal_tree& marginal, const mcont_t& mutations)
            {
                while (beg < end && mutations[beg->key].pos < marginal.left)
                    {
                        ++beg;
                    }
                while (beg < end && mutations[beg->key].pos < marginal.right)
                    {
                        if (mutations[beg->key].neutral == true)
                            {
                                auto nv = vf(marginal, beg->node);
                                if (nv)
                                    {
                                        auto d = vf.view_genotypes();
                                        //assert(nv == std::count(d.first,d.first+d.second,1));
                                        genotypes.insert(genotypes.end(),
                                                         d.first,
                                                         d.first + d.second);
                                        positions.push_back(
                                            mutations[beg->key].pos);
                                    }
                            }
                        ++beg;
                    }
            }
        };

        //TODO: allow for better return value
        template <typename mcont_t>
        std::pair<std::vector<std::int8_t>, std::vector<double>>
        create_data_matrix(const mcont_t& mutations,
                            const table_collection& tables,
                            const std::vector<std::int32_t>& samples)
        {
            create_data_matrix_details d(tables.mutation_table.begin(),
                                          tables.mutation_table.end(), samples,
                                          tables.node_table.size());
            auto visitor = [&d, &mutations](const marginal_tree& marginal) {
                d(marginal, mutations);
            };
            algorithmS(tables.input_left, tables.output_right, samples,
                       tables.node_table.size(), tables.L, visitor);
            return std::make_pair(std::move(d.genotypes),
                                  std::move(d.positions));
        }
    } // namespace ts
} // namespace fwdpp

#endif
