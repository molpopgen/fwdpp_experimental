#ifndef FWDPP_TS_VARIANT_FILLER_HPP
#define FWDPP_TS_VARIANT_FILLER_HPP

#include <cassert>
#include <vector>
#include <cstdint>
#include <algorithm>
#include "msprime_algo.hpp"

namespace fwdpp
{
    namespace ts
    {
        //TODO: rethink name.
        //variant_filler is a misnomer b/c what is really happening
        //is that samples are being tracked via top-2-bottom
        //traversal.
        //Best use of this type is probably composition
        //of a "visitor" for algorithmT.
        class variant_filler
        /// \brief Fill 0/1 genotype array for a mutation.
        {
          private:
            std::vector<std::int32_t> sample_indexes;
            std::vector<std::int32_t> node_stack;
            std::vector<std::int8_t> genotypes;
            std::int32_t stack_top;

          public:
            variant_filler(const std::int32_t nnodes,
                           const std::vector<std::int32_t>& samples)
                : sample_indexes(nnodes, -1), node_stack(nnodes + 1, -1),
                  genotypes(samples.size(), 0), stack_top(-1)
            {
                for (std::size_t i = 0; i < samples.size(); ++i)
                    {
                        auto s = samples[i];
                        if (s < 0 || s >= nnodes)
                            {
                                throw std::out_of_range("sample list contains "
                                                        "invalid node labels");
                            }
                        if (sample_indexes[s] != -1)
                            {
                                throw std::invalid_argument(
                                    "sample labels must be unique");
                            }
                        sample_indexes[static_cast<std::size_t>(s)]
                            = static_cast<std::int32_t>(i);
                    }
            }

            inline std::uint32_t
            operator()(const marginal_tree& marginal, const std::int32_t node)
            /// Returns true if node leads to any samples, false otherwise
            {
                if (node < 0
                    || static_cast<std::size_t>(node) >= sample_indexes.size())
                    {
                        throw std::out_of_range("node value out of range");
                    }
                std::uint32_t variant_in_samples = 0;
                std::fill(std::begin(genotypes), std::end(genotypes), 0);
                stack_top = 0;
                node_stack[stack_top] = node;
                while (stack_top >= 0)
                    {
                        const auto top_node = node_stack[stack_top];
                        assert(top_node!=-1);
                        const auto sample_index = sample_indexes[top_node];
                        if (sample_index != -1)
                            {
                                ++variant_in_samples;
                                assert(genotypes[sample_index]==0);
                                genotypes[sample_index] = 1;
                            }
                        stack_top--;
                        std::int32_t last = -1;
                        for (auto c = marginal.left_child[top_node]; c != -1;
                             c = marginal.right_sib[c])
                            {
                                ++stack_top;
                                node_stack[stack_top] = c;
                                last = c;
                            }
                        assert(last == marginal.right_child[top_node]);
                    }
                return variant_in_samples;
            }

            std::pair<const std::int8_t*, std::size_t>
            view_genotypes() const
            {
                return std::make_pair(genotypes.data(), genotypes.size());
            }
        };
    } // namespace ts
} // namespace fwdpp

#endif
