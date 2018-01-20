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
};

inline node
make_node(std::uint32_t id, double generation, std::int32_t population)
{
    node n;
    n.id = id;
    n.generation = generation;
    n.population = population;
    return n;
}

#endif
