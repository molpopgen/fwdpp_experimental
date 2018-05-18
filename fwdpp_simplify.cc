// Experimental integration of fwdpp engine
// with simplification algorithm from Ralph et al.
// This code starts by copy/pasting several data types
// and functions from the fwdpy11-based example for that paper.

#include <cmath>
#include <stdexcept>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <gsl/gsl_randist.h>
#include <fwdpp/forward_types.hpp>
#include <fwdpp/fitness_models.hpp>
#include <fwdpp/poisson_xover.hpp>
#include <fwdpp/recbinder.hpp>
#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/slocuspop.hpp>
#include <fwdpp/sugar/GSLrng_t.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>
#include <fwdpp/internal/gamete_cleaner.hpp>
#include <fwdpp/util.hpp>
#include "node.hpp"
#include "edge.hpp"
#include "table_simplifier.hpp"
//#include "split_breakpoints.hpp"

using namespace fwdpp::ancestry;

using slocuspop_t = fwdpp::slocuspop<fwdpp::popgenmut>;
using GSLrng_t = fwdpp::GSLrng_t<fwdpp::GSL_RNG_MT19937>;

template <typename fitness_function>
inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
w(slocuspop_t& pop, const fitness_function& ff)
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
    return std::make_tuple(
        first_parental_index + 2 * static_cast<std::int32_t>(parent)
            + did_swap,
        first_parental_index + 2 * static_cast<std::int32_t>(parent)
            + !did_swap);
}

void
debug_new_edges(std::vector<fwdpp::uint_t>& new_mutations,
                const std::vector<double>& breakpoints,
                const std::size_t nedges,
                const std::tuple<std::int32_t, std::int32_t>& parent_nodes,
                const std::int32_t offspring_node,
                const std::size_t offspring_gamete,
                const std::size_t parent_g1, const std::size_t parent_g2,
                const slocuspop_t& pop, const table_collection& tables)
{
    // All new mutations must be found in the new offspring
    for (auto m : new_mutations)
        {
            if (pop.mutations[m].neutral)
                {
                    if (std::find(
                            pop.gametes[offspring_gamete].mutations.begin(),
                            pop.gametes[offspring_gamete].mutations.end(), m)
                        == pop.gametes[offspring_gamete].mutations.end())
                        {
                            throw std::runtime_error(
                                "mutation not found in offspring");
                        }
                }
            else
                {
                    if (std::find(
                            pop.gametes[offspring_gamete].smutations.begin(),
                            pop.gametes[offspring_gamete].smutations.end(), m)
                        == pop.gametes[offspring_gamete].smutations.end())
                        {
                            throw std::runtime_error(
                                "mutation not found in offspring");
                        }
                }
        }
    // Enforce that all new mutations are found
    // in an edge
    for (auto m : new_mutations)
        {
            const auto pos = pop.mutations[m].pos;
            bool valid = false;
            for (std::size_t i = nedges;
                 !valid && nedges < tables.edge_table.size(); ++i)
                {
                    if (pos >= tables.edge_table[i].left
                        && pos < tables.edge_table[i].right)
                        {
                            valid = true;
                        }
                }
            if (!valid)
                {
                    throw std::runtime_error(
                        "a new mutation is not contained in new edges");
                }
        }
    // All new mutations must be associated with
    // the offspring node
    for (auto m : new_mutations)
        {
            bool valid = false;
            const auto pos = pop.mutations[m].pos;
            for (std::size_t i = 0; !valid && i < tables.mutation_table.size();
                 ++i)
                {
                    if (pop.mutations[tables.mutation_table[i].key].pos == pos)
                        {
                            if (tables.mutation_table[i].node
                                == offspring_node)
                                {
                                    valid = true;
                                }
                        }
                }
            if (!valid)
                {
                    throw std::runtime_error(
                        "new mutation node != offspring node");
                }
        }
    if (!breakpoints.empty())
        {
            std::vector<double> b2(breakpoints);
            auto p1 = std::get<0>(parent_nodes),
                 p2 = std::get<1>(parent_nodes);
            auto pg1 = parent_g1, pg2 = parent_g2;
            if (b2.front() != 0.0)
                {
                    b2.insert(b2.begin(), 0.0);
                }
            else
                {
                    std::swap(p1, p2);
                    std::swap(pg1, pg2);
                }
            if (!std::is_sorted(b2.begin(), b2.end()))
                {
                    throw std::runtime_error("breakpoints not sorted");
                }
            // Delete double x-overs, which is potential source of error
            b2.erase(std::unique(b2.begin(), b2.end()), b2.end());
            if (b2.size() != (tables.edge_table.size() - nedges + 1))
                {
                    throw std::runtime_error("invalid size of new edges");
                }
            for (std::size_t i = 1, j = nedges; i < b2.size(); ++i, ++j)
                {
                    double start = b2[i - 1],
                           stop = (b2[i] == std::numeric_limits<double>::max())
                                      ? tables.L
                                      : b2[i];
                    if (tables.edge_table[j].left
                        == tables.edge_table[j].right)
                        {
                            throw std::runtime_error(
                                "edge left == edge right");
                        }
                    if (start != tables.edge_table[j].left)
                        {
                            throw std::runtime_error("left disagrees");
                        }
                    if (stop != tables.edge_table[j].right)
                        {
                            throw std::runtime_error("right disagrees");
                        }
                    if (tables.edge_table[j].parent != p1)
                        {
                            throw std::runtime_error("parent mismatch");
                        }
                    if (tables.edge_table[j].child != offspring_node)
                        {
                            throw std::runtime_error("child mismatch");
                        }
                    std::swap(p1, p2);

                    //Finally, the offspring must have
                    //all mutations from the parent in the interval [start,stop).
                    auto nneutral
                        = pop.gametes[offspring_gamete].mutations.size();
                    auto nselected
                        = pop.gametes[offspring_gamete].smutations.size();
                    // need to adjust above for new mutations
                    unsigned new_neutral = 0, new_selected = 0;
                    for (auto m : new_mutations)
                        {
                            if (pop.mutations[m].neutral)
                                {
                                    ++new_neutral;
                                }
                            else
                                {
                                    ++new_selected;
                                }
                        }
                    auto itr_i = std::lower_bound(
                        pop.gametes[pg1].mutations.begin(),
                        pop.gametes[pg1].mutations.end(), start,
                        [&pop](const fwdpp::uint_t key, const double d) {
                            return pop.mutations[key].pos < d;
                        });
                    if (itr_i != pop.gametes[pg1].mutations.end())
                        {
                            auto itr_j = std::upper_bound(
                                pop.gametes[pg1].mutations.begin(),
                                pop.gametes[pg1].mutations.end(), start,
                                [&pop](double d, const fwdpp::uint_t key) {
                                    return d < pop.mutations[key].pos;
                                });
                            auto itr_c = itr_i;
                            for (; itr_i < itr_j; ++itr_i)
                                {
                                    if (std::find(pop.gametes[offspring_gamete]
                                                      .mutations.begin(),
                                                  pop.gametes[offspring_gamete]
                                                      .mutations.end(),
                                                  *itr_i)
                                        == pop.gametes[offspring_gamete]
                                               .mutations.end())
                                        {
                                            std::cerr << std::setprecision(15)
                                                      << start << ' ' << stop
                                                      << ' ' << pg1 << '\n';
                                            for (; itr_c < itr_j; ++itr_c)
                                                {
                                                    std::cerr
                                                        << std::setprecision(
                                                               15)
                                                        << "parental mutation "
                                                        << pop.mutations
                                                               [*itr_c]
                                                                   .pos
                                                        << '\n';
                                                }
                                            throw std::runtime_error(
                                                "neutral parental mut not in "
                                                "offspring");
                                        }
                                }
                        }
                    std::swap(pg1, pg2);
                }
        }
}

// Wow, that's a lot of stuff needed:
template <typename breakpoint_function, typename mutation_model,
          typename mrecbin, typename grecbin>
std::int32_t
generate_offspring(const GSLrng_t& rng, const breakpoint_function& recmodel,
                   const mutation_model& mmodel, const double mu,
                   const std::size_t parent, const fwdpp::uint_t parent_g1,
                   const fwdpp::uint_t parent_g2,
                   const std::tuple<std::int32_t, std::int32_t>& parent_nodes,
                   const std::int32_t generation,
                   const std::int32_t next_index, slocuspop_t& pop,
                   std::size_t& offspring_gamete, table_collection& tables,
                   mrecbin& mutation_recycling_bin,
                   grecbin& gamete_recycling_bin)
{
    auto breakpoints = recmodel();
    auto new_mutations = fwdpp::generate_new_mutations(
        mutation_recycling_bin, rng.get(), mu, pop.diploids[parent],
        pop.gametes, pop.mutations, parent_g1, mmodel);
    for (auto& m : new_mutations)
        {
            auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            assert(std::distance(itr.first, itr.second) == 1);
        }
    offspring_gamete = fwdpp::mutate_recombine(
        new_mutations, breakpoints, parent_g1, parent_g2, pop.gametes,
        pop.mutations, gamete_recycling_bin, pop.neutral, pop.selected);
    if (!new_mutations.empty() || !breakpoints.empty())
        {
            assert(offspring_gamete != parent_g1);
        }
    auto nedges = tables.edge_table.size();
    tables.add_offspring_data(next_index, breakpoints, new_mutations,
                              parent_nodes, generation);
#ifndef NDEBUG
    debug_new_edges(new_mutations, breakpoints, nedges, parent_nodes,
                    next_index, offspring_gamete, parent_g1, parent_g2, pop,
                    tables);
#endif
    return next_index + 1;
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(const GSLrng_t& rng, slocuspop_t& pop,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel,
                  const fwdpp::uint_t generation, table_collection& tables,
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

            next_index_local = generate_offspring(
                rng, recmodel, mmodel, mu, p1, p1g1, p1g2, p1id, generation,
                next_index_local, pop, dip.first, tables,
                mutation_recycling_bin, gamete_recycling_bin);
            next_index_local = generate_offspring(
                rng, recmodel, mmodel, mu, p2, p2g1, p2g2, p2id, generation,
                next_index_local, pop, dip.second, tables,
                mutation_recycling_bin, gamete_recycling_bin);
            //auto breakpoints = recmodel();
            //auto new_mutations = fwdpp::generate_new_mutations(
            //    mutation_recycling_bin, rng.get(), mu, pop.diploids[p1],
            //    pop.gametes, pop.mutations, p1g1, mmodel);
            //for (auto& m : new_mutations)
            //    {
            //        auto itr
            //            = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            //        assert(std::distance(itr.first, itr.second) == 1);
            //    }
            //dip.first = fwdpp::mutate_recombine(
            //    new_mutations, breakpoints, p1g1, p1g2, pop.gametes,
            //    pop.mutations, gamete_recycling_bin, pop.neutral,
            //    pop.selected);
            //if (!new_mutations.empty() || !breakpoints.empty())
            //    {
            //        assert(dip.first != p1g1);
            //    }

            //tables.add_offspring_data(next_index_local, breakpoints,
            //                          new_mutations, p1id, generation);
            //next_index_local++;
            //breakpoints = recmodel();
            //new_mutations = fwdpp::generate_new_mutations(
            //    mutation_recycling_bin, rng.get(), mu, pop.diploids[p2],
            //    pop.gametes, pop.mutations, p2g1, mmodel);
            //for (auto& m : new_mutations)
            //    {
            //        auto itr
            //            = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            //        assert(std::distance(itr.first, itr.second) == 1);
            //    }
            //dip.second = fwdpp::mutate_recombine(
            //    new_mutations, breakpoints, p2g1, p2g2, pop.gametes,
            //    pop.mutations, gamete_recycling_bin, pop.neutral,
            //    pop.selected);
            //if (!new_mutations.empty() || !breakpoints.empty())
            //    {
            //        assert(dip.second != p2g1);
            //    }
            //tables.add_offspring_data(next_index_local, breakpoints,
            //                          new_mutations, p2id, generation);
            //next_index_local++;
            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;
        }
    assert(next_index_local
           == next_index + 2 * static_cast<std::int32_t>(N_next));
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    fwdpp::fwdpp_internal::gamete_cleaner(
        pop.gametes, pop.mutations, pop.mcounts, 2 * N_next, std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
}

table_collection
evolve(const GSLrng_t& rng, slocuspop_t& pop,
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

    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());
    unsigned generation = 0;

    const auto mmodel
        = [&pop, &rng, &generation](std::queue<std::size_t>& recbin,
                                    slocuspop_t::mcont_t& mutations) {
              return fwdpp::infsites_popgenmut(
                  recbin, mutations, rng.get(), pop.mut_lookup, generation,
                  0.0, [&rng]() { return gsl_rng_uniform(rng.get()); },
                  []() { return 0.0; }, []() { return 0.0; });
          };

    table_simplifier ancestry(1.0);
    table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    std::int32_t first_parental_index = 0,
                 next_index = 2 * pop.diploids.size();
    for (; generation < generations; ++generation)
        {
            const auto N_next = popsizes.at(generation);
            evolve_generation(rng, pop, N_next, mu_neutral + mu_selected,
                              mmodel, recmap, generation, tables,
                              first_parental_index, next_index);

            std::vector<std::int32_t> samples;
            for (auto i = tables.num_nodes() - 2 * pop.diploids.size();
                 i < tables.num_nodes(); ++i)
                {
                    assert(tables.node_table[i].generation == generation + 1);
                    samples.push_back(i);
                }
            // auto mt = tables.mutation_table;
            // auto nt = tables.node_table;
            // auto et = tables.edge_table;
            // std::ostringstream ofn;
            // ofn << "simoutput/last_edges." << generation << ".bin";
            // std::ofstream out;
            // out.open(ofn.str().c_str());
            // for (auto& e : et)
            //     {
            //         out.write(reinterpret_cast<char*>(&e.parent),
            //                   sizeof(decltype(e.parent)));
            //         out.write(reinterpret_cast<char*>(&e.child),
            //                   sizeof(decltype(e.child)));
            //         out.write(reinterpret_cast<char*>(&e.left),
            //                   sizeof(decltype(e.left)));
            //         out.write(reinterpret_cast<char*>(&e.right),
            //                   sizeof(decltype(e.right)));
            //     }
            // std::int32_t done = -1;
            // out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            // out.close();
            // ofn.str(std::string());
            // ofn << "simoutput/last_nodes." << generation << ".bin";
            // out.open(ofn.str().c_str());
            // for (auto n : nt)
            //     {
            //         out.write(reinterpret_cast<char*>(&n.id), 4);
            //         out.write(reinterpret_cast<char*>(&n.generation),
            //                   sizeof(double));
            //     }
            // out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            // out.close();
            // ofn.str(std::string());
            // ofn << "simoutput/last_mutations." << generation << ".bin";
            // out.open(ofn.str().c_str());
            // for (auto& m : tables.mutation_table)
            //     {
            //         out.write(reinterpret_cast<char*>(&m.node),
            //                   sizeof(std::int32_t));
            //         out.write(
            //             reinterpret_cast<char*>(&pop.mutations[m.key].pos),
            //             sizeof(double));
            //     }
            // out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            // out.close();
            tables.sort_tables(pop.mutations);
            auto xx = ancestry.simplify(tables, samples, pop.mutations);
            if (xx.second != pop.mcounts)
                {
                    for (std::size_t i = 0; i < xx.second.size(); ++i)
                        {
                            std::cout << xx.second[i] << ' ' << pop.mcounts[i];
                            if (xx.second[i] != pop.mcounts[i])
                                std::cout << " *";
                            std::cout << '\n';
                        }
                }
            assert(xx.second == pop.mcounts);

            // TODO: decide how to handle fixations.
            // Ideally, this would be done during simplification,
            // via a policy passed into the simplifier.
            // Until then, we do it manually
            tables.mutation_table.erase(
                std::remove_if(
                    tables.mutation_table.begin(), tables.mutation_table.end(),
                    [&xx, &pop](const fwdpp::ancestry::mutation_record& mr) {
                        return xx.second[mr.key] == 2 * pop.diploids.size();
                    }),
                tables.mutation_table.end());

            assert(pop.mcounts.size() == xx.second.size());
            // ofn.str(std::string());
            // ofn << "simoutput/edges." << generation << ".txt";

            // out.open(ofn.str().c_str());
            // for (auto e : tables.edge_table)
            //     {
            //         out << e.parent << ' ' << e.child << ' ' << e.left << ' '
            //             << e.right << '\n';
            //     }
            // out.close();
            // ofn.str(std::string());
            // ofn << "simoutput/nodes." << generation << ".txt";
            // out.open(ofn.str().c_str());
            // for (auto n : tables.node_table)
            //     {
            //         out << n.id << ' ' << n.generation << '\n';
            //     }
            // out.close();
            // ofn.str(std::string());
            // ofn << "simoutput/idmap." << generation << ".txt";
            // out.open(ofn.str().c_str());
            // for (std::size_t id = 0; id < xx.first.size(); ++id)
            //     {
            //         out << id << ' ' << xx.first[id] << '\n';
            //     }
            // out.close();

            // if (pop.mcounts != xx.second)
            //     {
            //         std::vector<std::size_t> failures;
            //         for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
            //             {
            //                 std::cout << generation << ' ' << pop.mcounts[i]
            //                           << ' ' << xx.second[i] << ' '
            //                           << pop.mutations[i].pos << ' '
            //                           << pop.mutations[i].g << '\n';
            //                 if (pop.mcounts[i] != xx.second[i])
            //                     {
            //                         failures.push_back(i);
            //                     }
            //             }
            //         //Find diploids with failures
            //         std::vector<std::size_t> gfails;
            //         for (std::size_t g = 0; g < pop.gametes.size(); ++g)
            //             {
            //                 auto& gam = pop.gametes[g];
            //                 if (gam.n)
            //                     {
            //                         for (auto f : failures)
            //                             {
            //                                 auto itr = std::find(
            //                                     gam.mutations.begin(),
            //                                     gam.mutations.end(), f);
            //                                 auto itr2 = std::find(
            //                                     gam.smutations.begin(),
            //                                     gam.smutations.end(), f);
            //                                 if (itr != gam.mutations.end()
            //                                     || itr2
            //                                            != gam.smutations.end())
            //                                     {
            //                                         gfails.push_back(g);
            //                                     }
            //                             }
            //                     }
            //             }
            //         std::sort(gfails.begin(), gfails.end());
            //         gfails.erase(std::unique(gfails.begin(), gfails.end()),
            //                      gfails.end());
            //         for (std::size_t dip = 0; dip < pop.diploids.size(); ++dip)
            //             {
            //                 if (std::find(gfails.begin(), gfails.end(),
            //                               pop.diploids[dip].first)
            //                     != gfails.end())
            //                     {
            //                         std::cout << "Node with mut: " << 2 * dip
            //                                   << '\n';
            //                     }
            //                 if (std::find(gfails.begin(), gfails.end(),
            //                               pop.diploids[dip].second)
            //                     != gfails.end())
            //                     {
            //                         std::cout
            //                             << "Node with mut: " << 2 * dip + 1
            //                             << '\n';
            //                     }
            //             }
            //         for (unsigned i = 0; i < mt.size(); ++i)
            //             {
            //                 std::cout << mt[i].node << ' '
            //                           << xx.first[mt[i].node] << ' '
            //                           << pop.mutations[mt[i].key].pos << '\n';
            //             }
            //         std::exit(0);
            //     }
            next_index = tables.num_nodes();
            first_parental_index = 0;
            fwdpp::update_mutations(pop.mutations, pop.fixations,
                                    pop.fixation_times, pop.mut_lookup,
                                    pop.mcounts, generation,
                                    2 * pop.diploids.size());
        }
    return tables;
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

    slocuspop_t pop(N);
    std::vector<fwdpp::uint_t> popsizes(10 * N, N);
    GSLrng_t rng(seed);
    double mu = theta / (4. * static_cast<double>(N));
    double recrate = rho / (4. * static_cast<double>(N));
    double mudel = mu * pdel;

    auto tables = evolve(rng, pop, popsizes, mu, mudel, recrate);
    std::cout << "finished without error " << pop.mutations.size() << ' '
              << tables.node_table.size() << ' ' << tables.edge_table.size()
              << ' ' << tables.mutation_table.size() << '\n';
}
