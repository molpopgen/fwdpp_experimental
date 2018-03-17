// Copied from fwdpy11_arg_example.
// Author: KRT
// License: GPL3+
#ifndef FWDPP_ANCESTRY_NODE_HPP__
#define FWDPP_ANCESTRY_NODE_HPP__

#include <cstdint>

namespace fwdpp
{
    namespace ancestry
    {
        struct node
        {
            std::int32_t id;
            std::int32_t population;
            double generation;
        };
    }
}
#endif
