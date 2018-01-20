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
};

inline edge
make_edge(double left, double right, std::int32_t parent, std::int32_t child)
{
    edge e;
    e.left = left;
    e.right = right;
    e.parent = parent;
    e.child = child;
    return e;
}
#endif
