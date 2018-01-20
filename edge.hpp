//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef ANCESTRY_EDGE_HPP__
#define ANCESTRY_EDGE_HPP__

#include <cstdint>

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

//inline edge
//make_edge(double left, double right, std::int32_t parent, std::int32_t child)
//{
//    edge e;
//    e.left = left;
//    e.right = right;
//    e.parent = parent;
//    e.child = child;
//    return e;
//}
#endif
