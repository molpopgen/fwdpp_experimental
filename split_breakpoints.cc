//Copied from fwdpy11_arg_example.
//Author: KRT
//License: GPL3+
#include "split_breakpoints.hpp"

std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints, const double start,
                  const double stop)
{
    std::vector<std::pair<double, double>> r1, r2;
    if (breakpoints.front() != 0.0)
        {
            r1.emplace_back(std::make_pair(start, breakpoints.front()));
        }
    for (unsigned j = 1; j < breakpoints.size(); ++j)
        {
            double a = breakpoints[j - 1];
            double b = (j < breakpoints.size() - 1) ? breakpoints[j] : stop;
            if (j % 2 == 0.)
                {
                    r1.emplace_back(a, b);
                }
            else
                {
                    r2.emplace_back(a, b);
                }
        }
    return std::make_pair(std::move(r1), std::move(r2));
}

