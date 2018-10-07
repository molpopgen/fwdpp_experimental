// Experimental integration of fwdpp engine
// with simplification algorithm from Kelleher et al.
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
#include <gsl/gsl_matrix_char.h>
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
#include <fwdpp/algorithm/compact_mutations.hpp>
#include "node.hpp"
#include "edge.hpp"
#include "table_simplifier.hpp"
#include "data_matrix_generator.hpp"
//#include "split_breakpoints.hpp"

using namespace fwdpp::ts;

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
#ifndef NDEBUG
    for (auto& m : new_mutations)
        {
            auto itr = pop.mut_lookup.equal_range(pop.mutations[m].pos);
            assert(std::distance(itr.first, itr.second) == 1);
        }
#endif
    // We will only add selected mutations into offspring gametes.
    auto end_of_neutral
        = std::stable_partition(new_mutations.begin(), new_mutations.end(),
                                [&pop](const fwdpp::uint_t key) {
                                    return pop.mutations[key].neutral == true;
                                });
    //if (!std::is_sorted(end_of_neutral, new_mutations.end(),
    //                    [&pop](const std::size_t a, const std::size_t b) {
    //                        return pop.mutations[a].pos < pop.mutations[b].pos;
    //                    }))
    //    {
    //        throw std::runtime_error("bad");
    //    }
    offspring_gamete = fwdpp::mutate_recombine(
        decltype(new_mutations)(end_of_neutral, new_mutations.end()),
        breakpoints, parent_g1, parent_g2, pop.gametes, pop.mutations,
        gamete_recycling_bin, pop.neutral, pop.selected);
    //if (!new_mutations.empty() || !breakpoints.empty())
    //    {
    //        assert(offspring_gamete != parent_g1);
    //    }
    tables.add_offspring_data(next_index, breakpoints, new_mutations,
                              parent_nodes, generation);
    return next_index + 1;
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(const GSLrng_t& rng, slocuspop_t& pop,
                  std::vector<std::int32_t>& ancient_nodes,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel,
                  const fwdpp::uint_t generation, table_collection& tables,
                  table_simplifier& simplifier,
                  std::int32_t first_parental_index, std::int32_t next_index)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);

    auto lookup = w(pop, fwdpp::multiplicative_diploid());
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
            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;
        }
    assert(next_index_local
           == next_index + 2 * static_cast<std::int32_t>(N_next));
    // This is constant-time
    pop.diploids.swap(offspring);
    tables.sort_tables(pop.mutations);
    std::vector<std::int32_t> samples(2 * pop.diploids.size());
    std::iota(samples.begin(), samples.end(),
              tables.num_nodes() - 2 * pop.diploids.size());
    samples.insert(samples.end(), ancient_nodes.begin(), ancient_nodes.end());
    auto idmap
        = simplifier.simplify(tables, samples, pop.mutations);
    // remap ancient nodes
    for (auto& a : ancient_nodes)
        {
            assert(idmap[a] != -1);
            a = idmap[a];
        }
    if (generation > 0 && generation % 100 == 0.0)
        {
            //record every 50th node of the current generation
            //as an ancient sample
            for (std::size_t i = 0; i < 2 * pop.diploids.size(); ++i)
                {
                    if (i && i % 50 == 0.0)
                        {
                            ancient_nodes.push_back(i);
                            assert(tables.node_table[i].generation
                                   == generation + 1);
                        }
                }
        }
    //TODO: deal with the next three calls.
    //The cleansing of fixations from the 
    //mutation table will not be correct if
    //an ancient sample contains a variant.
    //Same with the calls to gamete_cleaner
    //and update_mutations.
    tables.mutation_table.erase(
        std::remove_if(
            tables.mutation_table.begin(), tables.mutation_table.end(),
            [&pop](const fwdpp::ts::mutation_record& mr) {
                return pop.mcounts[mr.key] == 2 * pop.diploids.size();
            }),
        tables.mutation_table.end());
    fwdpp::fwdpp_internal::gamete_cleaner(
        pop.gametes, pop.mutations, pop.mcounts, 2 * N_next, std::true_type());
    fwdpp::update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                            pop.mut_lookup, pop.mcounts, generation,
                            2 * pop.diploids.size());
    if (generation && generation % 100 == 0.0)
        {
            auto new_mut_indexes = fwdpp::compact_mutations(pop);
            // Mutation compacting re-orders the data, so we need to
            // re-index our mutation table
            for (auto& mr : tables.mutation_table)
                {
                    mr.key = new_mut_indexes[mr.key];
                }
        }
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

    // With selection, 2Ns ~ Gamma with mean -5 and shape = 1.
    // Fitness will be 1, 1+sh, 1+s, and we will model additive effects
    // at a site, and multiplicative across
    const double twoN = 2 * pop.diploids.size();
    const double pdeleterious = mu_selected / (mu_selected + mu_neutral);
    const auto mmodel = [&pop, &rng, &generation, pdeleterious,
                         twoN](std::queue<std::size_t>& recbin,
                               slocuspop_t::mcont_t& mutations) {
        return fwdpp::infsites_popgenmut(
            recbin, mutations, rng.get(), pop.mut_lookup, generation,
            pdeleterious, [&rng]() { return gsl_rng_uniform(rng.get()); },
            [&rng, twoN]() {
                return gsl_ran_gamma(rng.get(), 1., -5.) / twoN;
            },
            []() { return 0.5; });
    };

    table_simplifier simplifier(1.0);
    table_collection tables(2 * pop.diploids.size(), 0, 0, 1.0);
    std::int32_t first_parental_index = 0,
                 next_index = 2 * pop.diploids.size();
    std::vector<std::int32_t> ancient_nodes;
    for (; generation < generations; ++generation)
        {
            const auto N_next = popsizes[generation];
            evolve_generation(rng, pop, ancient_nodes, N_next,
                              mu_neutral + mu_selected, mmodel, recmap,
                              generation, tables, simplifier,
                              first_parental_index, next_index);
            if (generation && generation % 10 == 0.0)
                {
                    std::vector<std::int32_t> nodes(2 * pop.diploids.size());
                    std::iota(nodes.begin(), nodes.end(), 0);
                    //std::vector<std::int32_t> samples(2 * pop.diploids.size()
                    //                                  / 10);
                    std::vector<std::int32_t> samples(2 * pop.diploids.size());

                    gsl_ran_choose(rng.get(), samples.data(), samples.size(),
                                   nodes.data(), nodes.size(),
                                   sizeof(int32_t));
                    auto gm = fwdpp::ts::create_data_matrix(
                        pop.mutations, tables, samples, false, true);
                    auto x = std::move(gm.second);
                    if (!x.positions.empty())
                        {
                            assert(std::is_sorted(x.positions.begin(),
                                                  x.positions.end()));
                            auto nrow = x.positions.size();
                            auto ncol = x.genotypes.size() / nrow;
                            assert(ncol == samples.size());
                            decltype(pop.mcounts) mc;
                            fwdpp::fwdpp_internal::process_gametes(
                                pop.gametes, pop.mutations, mc);
                            for (std::size_t i = 0; i < nrow; ++i)
                                {
                                    unsigned nd = 0;
                                    for (std::size_t j = i * ncol;
                                         j < i * ncol + ncol; ++j)
                                        {
                                            if (x.genotypes[j] == 1)
                                                {
                                                    ++nd;
                                                }
                                        }
                                    if (mc[x.mut_indexes[i]] != nd)
                                        {
                                            throw std::runtime_error(
                                                "bad counts from matrix");
                                        }
                                }
                        }
                }
            next_index = tables.num_nodes();
            first_parental_index = 0;
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
    const double mudel = pdel * mu;
    auto tables = evolve(rng, pop, popsizes, mu, mudel, recrate);
    auto nselected = 0;
    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
        {
            if (pop.mcounts[i] && !pop.mutations[i].neutral)
                {
                    ++nselected;
                }
        }
    std::cout << "finished without error " << pop.mutations.size() << ' '
              << tables.node_table.size() << ' ' << tables.edge_table.size()
              << ' ' << tables.mutation_table.size() << ' ' << nselected
              << '\n';
}
