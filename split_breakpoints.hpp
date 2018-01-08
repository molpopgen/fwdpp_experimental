#ifndef SPLIT_BREAKPOINTS_HPP__
#define SPLIT_BREAKPOINTS_HPP__

#include <utility>
#include <vector>

std::pair<std::vector<std::pair<double, double>>,
          std::vector<std::pair<double, double>>>
split_breakpoints(const std::vector<double>& breakpoints, const double start,
                  const double stop);

#endif
