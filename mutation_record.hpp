#ifndef FWDPP_ANCESTRY_MUTATION_RECORD_HPP__
#define FWDPP_ANCESTRY_MUTATION_RECORD_HPP__

#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct mutation_record
        {
            std::int32_t node;
            std::size_t key;
            double pos;
        };
    } // namespace ts
} // namespace fwdpp

#endif
