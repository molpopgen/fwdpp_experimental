//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef FWDPP_ANCESTRY_EDGE_HPP
#define FWDPP_ANCESTRY_EDGE_HPP

#include <cstdint>

namespace fwdpp
{
    namespace ts
    {
        struct edge
        {
            double left, right;
            std::int32_t parent, child;
        };
    }
}

#endif
