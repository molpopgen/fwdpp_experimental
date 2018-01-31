// Experimental integration of fwdpp engine
// with simplification algorithm from Ralph et al.
// This code starts by copy/pasting several data types
// and functions from the fwdpy11-based example for that paper.

#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <algorithm>
#include <queue>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/infsites.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/util.hpp>
#include "node.hpp"
#include "edge.hpp"
#include "ancestry_tracker.hpp"
//#include "split_breakpoints.hpp"

using namespace fwdpp::ancestry;

// static constexpr std::int32_t ROOTNODE
//     = std::numeric_limits<std::int32_t>::min();
//
// struct table_collection
// {
//     struct segment
//     {
//         double left, right;
//         std::int32_t node;
//         segment() : left{}, right{}, node{} {}
//         segment(double l, double r, std::int32_t n)
//             : left{ l }, right{ r }, node{ n }
//         {
//         }
//     };
//
//     std::vector<node> node_table;
//     std::vector<edge> edge_table;
//     std::vector<std::pair<std::int32_t, std::size_t>> mutation_table;
//
//     table_collection() : node_table{}, edge_table{}, mutation_table{} {}
//
//     table_collection(const std::int32_t num_initial_nodes,
//                      const double initial_time)
//         : node_table{}, edge_table{}, mutation_table{}
//     {
//         for (std::int32_t i = 0; i < num_initial_nodes; ++i)
//             {
//                 node_table.push_back(node(i, initial_time, 0));
//             }
//     }
//
//     std::int32_t
//     add_offspring_data(const std::int32_t next_index,
//                        const std::vector<double>& breakpoints,
//                        const std::vector<fwdpp::uint_t>& new_mutations,
//                        const std::tuple<std::int32_t, std::int32_t>&
//                        parents,
//                        const double generation)
//     {
//         node_table.push_back(
//             node(next_index, generation + 1, 0)); // MUSTDOC
//         auto split = split_breakpoints(breakpoints, 0., 1.);
//         // Add the edges
//         for (auto&& brk : split.first)
//             {
//                 edge_table.push_back(edge(brk.first, brk.second,
//                                           std::get<0>(parents),
//                                           next_index));
//             }
//         for (auto&& brk : split.second)
//             {
//                 edge_table.push_back(edge(brk.first, brk.second,
//                                           std::get<1>(parents),
//                                           next_index));
//             }
//         for (auto&& m : new_mutations)
//             mutation_table.emplace_back(next_index, m);
//         return next_index + 1;
//     }
//
//     void
//     simplify(const std::vector<std::int32_t>& samples)
//     {
//         // reverse time
//         auto maxtime = node_table.back().generation;
//         for (auto& n : node_table)
//             {
//                 n.generation -= maxtime;
//                 // Note: leads to 0 being -0.  SHOULDFIX
//                 n.generation *= -1.0;
//             }
//
//         // Sort the edge table
//         std::sort(edge_table.begin(), edge_table.end(),
//                   [this](const edge& a, const edge& b) {
//                       return std::tie(this->node_table[a.child].generation,
//                                       a.parent, a.child, a.left)
//                              <
//                              std::tie(this->node_table[b.child].generation,
//                                         b.parent, b.child, b.left);
//                   });
//
//         std::vector<edge> Eo;
//         std::vector<node> No;
//         std::vector<std::vector<segment>> Ancestry(node_table.size());
//         // The algorithm using a min queue.  The default C++ queue
//         // is a max queue.  Thus, we must use > rather than <
//         // to generate a min queue;
//         const auto segment_sorter_q = [](const segment& a, const segment& b)
//         {
//             return a.left > b.left;
//         };
//         std::priority_queue<segment, std::vector<segment>,
//                             decltype(segment_sorter_q)>
//             Q(segment_sorter_q);
//
//         for (auto& s : samples)
//             {
//                 No.push_back(node(s, node_table[s].generation, 0));
//                 Ancestry[s].push_back(
//                     segment(0, 1, static_cast<std::int32_t>(No.size() -
//                     1)));
//             }
//
//         auto last_edge = edge_table.begin();
//         segment alpha;
//         std::int32_t u;
//         while (last_edge < edge_table.end())
//             {
//                 u = last_edge->parent;
//                 for (; last_edge < edge_table.end() && last_edge->parent ==
//                 u;
//                      ++last_edge)
//                     {
//                         for (auto& seg : Ancestry[last_edge->child])
//                             {
//                                 if (seg.right > last_edge->left
//                                     && seg.right > last_edge->left)
//                                     {
//                                         Q.emplace(std::max(seg.left,
//                                                            last_edge->left),
//                                                   std::min(seg.right,
//                                                            last_edge->right),
//                                                   seg.node);
//                                     }
//                             }
//                     }
//
//                 std::int32_t v = -1;
//                 while (!Q.empty())
//                     {
//                         auto l = Q.top().left;
//                         double r = 1.0;
//                         std::vector<segment> X;
//                         while (!Q.empty() && Q.top().left == l)
//                             {
//                                 // Can be done w/fewer lines of code.
//                                 auto seg = Q.top();
//                                 Q.pop();
//                                 r = std::min(r, seg.right);
//                                 X.push_back(std::move(seg));
//                             }
//                         if (!Q.empty())
//                             {
//                                 r = std::min(r, Q.top().left);
//                             }
//                         if (X.size() == 1)
//                             {
//                                 alpha = X[0];
//                                 auto x = X[0];
//                                 if (!Q.empty() && Q.top().left < x.right)
//                                     {
//                                         alpha = segment(x.left,
//                                         Q.top().left,
//                                                         x.node);
//                                         x.left = Q.top().left;
//                                         Q.push(x);
//                                     }
//                             }
//                         else
//                             {
//                                 if (v == -1)
//                                     {
//                                         No.push_back(node(
//                                             static_cast<std::int32_t>(
//                                                 No.size()),
//                                             node_table[u].generation, 0));
//                                         v = No.size() - 1;
//                                     }
//                                 alpha = segment(l, r, v);
//                                 for (auto& x : X)
//                                     {
//                                         Eo.push_back(edge(l, r, v, x.node));
//                                         if (x.right > r)
//                                             {
//                                                 x.left = r;
//                                                 Q.emplace(x);
//                                             }
//                                     }
//                             }
//                         Ancestry[u].push_back(alpha);
//                     }
//             }
//
//         std::size_t start = 0;
//         std::vector<edge> compacted_edges;
//         for (std::size_t j = 1; j < Eo.size(); ++j)
//             {
//                 bool condition = Eo[j - 1].right != Eo[j].left
//                                  || Eo[j - 1].parent != Eo[j].parent
//                                  || Eo[j - 1].child != Eo[j].child;
//                 if (condition)
//                     {
//                         compacted_edges.push_back(
//                             edge(Eo[j - 1].left, Eo[j - 1].right,
//                                  Eo[j - 1].parent, Eo[j - 1].child));
//                         start = j;
//                     }
//             }
//         // This is probably really close to the above
//         Eo.erase(std::unique(Eo.begin(), Eo.end(),
//                              [](const edge& a, const edge& b) {
//                                  return a.parent == b.parent
//                                         && a.child == b.child
//                                         && a.right == b.left;
//                              }),
//                  Eo.end());
//         edge_table.swap(compacted_edges);
//         node_table.swap(No);
//     }
// };

using singlepop_t = fwdpp::singlepop<fwdpp::popgenmut>;
using GSLrng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

template <typename fitness_function>
inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
w(singlepop_t& pop, const fitness_function& ff)
{
    auto N_curr = pop.diploids.size();
    std::vector<double> fitnesses(N_curr);
    for (size_t i = 0; i < N_curr; ++i)
        {
            fitnesses[i] = ff(pop.diploids[i], pop.gametes, pop.mutations);
            pop.gametes[pop.diploids[i].first].n = 0;
            pop.gametes[pop.diploids[i].second].n = 0;
        }
    auto lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
    return lookup;
}

std::tuple<std::int32_t, std::int32_t>
get_parent_ids(const std::int32_t first_parental_index,
               const std::uint32_t parent, const int did_swap)
{
    // if (first_parental_index == ROOTNODE)
    //    {
    //        return std::make_tuple(ROOTNODE, ROOTNODE);
    //    }
    return std::make_tuple(
        first_parental_index + 2 * static_cast<std::int32_t>(parent)
            + did_swap,
        first_parental_index + 2 * static_cast<std::int32_t>(parent)
            + !did_swap);
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(const GSLrng_t& rng, singlepop_t& pop,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel,
                  const fwdpp::uint_t generation, ancestry_tracker& ancestry,
                  std::int32_t first_parental_index, std::int32_t next_index)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);

    auto lookup = w(pop, fwdpp::additive_diploid());
    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    auto next_index_local = next_index;
    for (auto& dip : offspring)
        {
            auto p1 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p2 = gsl_ran_discrete(rng.get(), lookup.get());
            auto p1g1 = pop.diploids[p1].first;
            auto p1g2 = pop.diploids[p1].second;
            auto p2g1 = pop.diploids[p2].first;
            auto p2g2 = pop.diploids[p2].second;

            // Mendel
            int swap1 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            int swap2 = (gsl_rng_uniform(rng.get()) < 0.5) ? 1 : 0;
            if (swap1)
                std::swap(p1g1, p1g2);
            if (swap2)
                std::swap(p2g1, p2g2);

            auto p1id = get_parent_ids(first_parental_index, p1, swap1);
            auto p2id = get_parent_ids(first_parental_index, p2, swap2);

            assert(std::get<0>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<0>(p2id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p2id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<0>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p1id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<0>(p2id) < 2 * static_cast<std::int32_t>(N_next));
            assert(std::get<1>(p2id) < 2 * static_cast<std::int32_t>(N_next));
            // auto breakpoints = fwdpp::generate_breakpoints(
            //    pop.diploids[p1], p1g1, p1g2, pop.gametes, pop.mutations,
            //    recmodel);
            auto breakpoints = recmodel();
            auto new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p1],
                pop.gametes, pop.mutations, p1g1, mmodel);
            dip.first = fwdpp::mutate_recombine(
                new_mutations, breakpoints, p1g1, p1g2, pop.gametes,
                pop.mutations, gamete_recycling_bin, pop.neutral,
                pop.selected);

            ancestry.add_offspring_data(next_index_local, breakpoints,
                                        new_mutations, p1id, generation);
            next_index_local++;
            breakpoints
                = recmodel(); // fwdpp::generate_breakpoints(pop.diploids[p2],
                              // p2g1,
            //                 p2g2, pop.gametes,
            //               pop.mutations, recmodel);
            new_mutations = fwdpp::generate_new_mutations(
                mutation_recycling_bin, rng.get(), mu, pop.diploids[p2],
                pop.gametes, pop.mutations, p2g1, mmodel);
            dip.second = fwdpp::mutate_recombine(
                new_mutations, breakpoints, p2g1, p2g2, pop.gametes,
                pop.mutations, gamete_recycling_bin, pop.neutral,
                pop.selected);
            ancestry.add_offspring_data(next_index_local, breakpoints,
                                        new_mutations, p2id, generation);
            next_index_local++;
            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;
        }
    assert(next_index_local == next_index + 2 * N_next);
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    fwdpp::fwdpp_internal::gamete_cleaner(
        pop.gametes, pop.mutations, pop.mcounts, 2 * N_next, std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
}

ancestry_tracker
evolve(const GSLrng_t& rng, singlepop_t& pop,
       const std::vector<std::uint32_t>& popsizes, const double mu_neutral,
       const double mu_selected, const double recrate)
{
    const auto generations = popsizes.size();
    if (!generations)
        throw std::runtime_error("empty list of population sizes");
    if (mu_selected < 0.)
        {
            throw std::runtime_error("negative selected mutation rate: "
                                     + std::to_string(mu_selected));
        }
    if (recrate < 0.)
        {
            throw std::runtime_error("negative recombination rate: "
                                     + std::to_string(recrate));
        }
    pop.mutations.reserve(std::ceil(
        std::log(2 * pop.diploids.size())
        * (4. * double(pop.diploids.size()) * (mu_selected)
           + 0.667 * (4. * double(pop.diploids.size()) * (mu_selected)))));

    fwdpp::poisson_xover recmap(rng.get(), recrate, 0., 1.);
    unsigned generation = 0;
    fwdpp::infsites inf;
    const auto mmodel
        = [&pop, &rng, &inf, &generation, mu_neutral, mu_selected](
            fwdpp::traits::recycling_bin_t<singlepop_t::mcont_t>& recbin,
            singlepop_t::mcont_t& mutations) {
              return inf(recbin, mutations, rng.get(), pop.mut_lookup,
                         generation, mu_neutral, mu_selected,
                         [&rng]() { return gsl_rng_uniform(rng.get()); },
                         []() { return 0.; }, []() { return 0.; });
          };

    ancestry_tracker ancestry(2 * pop.diploids.size(), 0, 0);
    std::int32_t first_parental_index = 0,
                 next_index = 2 * pop.diploids.size();
    double sort_time = 0.0;
    for (; generation < generations; ++generation)
        {
            const auto N_next = popsizes.at(generation);
            evolve_generation(rng, pop, N_next, mu_neutral + mu_selected,
                              mmodel, recmap, generation, ancestry,
                              first_parental_index, next_index);
            // if (first_parental_index == ROOTNODE)
            //    {
            //        first_parental_index = 0;
            //    }
            // else
            //    {
            //        first_parental_index += 2 * pop.diploids.size();
            //    }
            // std::cout<<next_index<<' '<<first_parental_index<<"->";

            std::vector<std::int32_t> samples;
            for (auto i = ancestry.num_nodes() - 2 * pop.diploids.size();
                 i < ancestry.num_nodes(); ++i)
                {
                    assert(ancestry.nodes()[i].generation == generation+1);
                    samples.push_back(i);
                }
            ancestry.sort_tables();
            ancestry.simplify(samples);
            next_index = ancestry.num_nodes();
            first_parental_index = 0;
            // std::cout<<next_index<<' '<<first_parental_index<<"\n";
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation,
                                    2 * pop.diploids.size());
        }
    std::cout << "sort time = " << sort_time << '\n';
    // std::vector<std::int32_t> samples;
    // for (auto i = next_index - 2 * pop.diploids.size(); i < next_index; ++i)
    //    {
    //        samples.push_back(i);
    //    }
    // ancestry.sort_tables();
    // ancestry.simplify(samples);
    return ancestry;
}

int
main(int argc, char** argv)
{
    int argn = 1;
    fwdpp::uint_t N = atoi(argv[argn++]);
    double theta = atof(argv[argn++]);
    double rho = atof(argv[argn++]);
    double pdel = atof(argv[argn++]);
    unsigned seed = atoi(argv[argn++]);

    singlepop_t pop(N);
    std::vector<fwdpp::uint_t> popsizes(10*N, N);
    GSLrng_t rng(seed);
    double mu = theta / (4. * static_cast<double>(N));
    double recrate = rho / (4. * static_cast<double>(N));
    double mudel = mu * pdel;

    auto ancestry = evolve(rng, pop, popsizes, mu, mudel, recrate);
    auto tables = ancestry.dump_tables();
    std::cout << pop.mutations.size() << ' ' << tables.node_table.size() << ' '
              << tables.edge_table.size() << ' '
              << tables.mutation_table.size() << '\n';
    // for (auto& a : tables.edge_table)
    //    {
    //        std::cout << tables.node_table[a.parent].generation << ' '
    //                  << a.parent << ' ' << a.child << ' ' << a.left << ' '
    //                  << a.right << '\n';
    //    }
    // for (auto& n : tables.node_table)
    //    {
    //        std::cout << n.id << ' ' << n.generation << '\n';
    //    }
}
