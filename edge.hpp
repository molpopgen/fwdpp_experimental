//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef FWDPP_ANCESTRY_EDGE_HPP__
#define FWDPP_ANCESTRY_EDGE_HPP__

#include <cstdint>

namespace fwdpp
{
    namespace ancestry
    {
        struct edge
        {
            double left, right;
            std::int32_t parent, child;
            edge(double l, double r, std::int32_t p, std::int32_t c) noexcept
                : left{ l },
                  right{ r },
                  parent{ p },
                  child{ c }
            {
            }
        };
    }
}

#endif
