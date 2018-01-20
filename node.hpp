//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#ifndef ANCESTRY_NODE_HPP__
#define ANCESTRY_NODE_HPP__

#include <cstdint>

struct node
{
    std::int32_t id;
    std::int32_t population;
    double generation;

    node(std::int32_t id_, double generation_, std::int32_t pop_) noexcept
        : id{ id_ },
          population{ pop_ },
          generation{ generation_ }
    {
    }
};

#endif
