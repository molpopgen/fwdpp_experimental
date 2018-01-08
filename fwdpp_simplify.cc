// Experimental integration of fwdpp engine
// with simplification algorithm from Ralph et al.
// This code starts by copy/pasting several data types
// and functions from the fwdpy11-based example for that paper.

#include <cmath>
#include <stdexcept>

#include <fwdpp/internal/sample_diploid_helpers.hpp>
#include <fwdpp/mutate_recombine.hpp>
#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/sugar/singlepop.hpp>
#include <fwdpp/sugar/GSLrng_t.h>

#include "node.hpp"
#include "edge.hpp"

using singlepop_t = fwdpp::singlepop<fwdpp::popgenmut>;
using GSLrng_t = fwdpp::GSLrng_t<GSL_RNG_MT19937>;

template <fitness_function>
inline fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr
w(singlepop_t& pop, const single_locus_fitness_fxn& ff)
{
    auto N_curr = pop.diploids.size();
    std::vector<double> fitnesses(N_curr);
    for (size_t i = 0; i < N_curr; ++i)
        {
            fitnesses[i] = ff(pop.diploids[i], pop.gametes, pop.mutations);
            pop.gametes[pop.diploids[i].first].n = 0;
            pop.gametes[pop.diploids[i].second].n = 0;
        }
    lookup = fwdpp::fwdpp_internal::gsl_ran_discrete_t_ptr(
        gsl_ran_discrete_preproc(N_curr, &fitnesses[0]));
    return lookup;
}

template <typename breakpoint_function, typename mutation_model>
void
evolve_generation(const GSLrng_t& rng, singlepop_t& pop,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel)
{

    auto gamete_recycling_bin
        = fwdpp::fwdpp_internal::make_gamete_queue(pop.gametes);
    auto mutation_recycling_bin
        = fwdpp::fwdpp_internal::make_mut_queue(pop.mcounts);

    auto lookup = w(pop, ...);
    decltype(pop.diploids) offspring(N_next);

    // Generate the offspring
    std::size_t label = 0;
    for (auto& dip : offspring)
        {
            auto p1 = gsl_ran_discrete_lookup(rng.get(), lookup.get());
            auto p2 = gsl_ran_discrete_lookup(rng.get(), lookup.get());
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

            auto breakpoints = recmodel(pop.gametes[p1g1], pop.gametes[p1g2],
                                        pop.mutations);
            dip.first = ancestry_recombination_details(
                pop, ancestry, gamete_recycling_bin, p1g1, p1g2, breakpoints,
                pid, std::get<0>(offspring_indexes));
            breakpoints = recmodel(pop.gametes[p2g1], pop.gametes[p2g2],
                                   pop.mutations);
            pid = ancestry.get_parent_ids(p2, swap2);

            dip.second = ancestry_recombination_details(
                pop, ancestry, gamete_recycling_bin, p2g1, p2g2, breakpoints,
                pid, std::get<1>(offspring_indexes));

            pop.gametes[dip.first].n++;
            pop.gametes[dip.second].n++;

            // now, add new mutations
            dip.first = fwdpp::mutate_gamete_recycle(
                mutation_recycling_bin, gamete_recycling_bin, rng.get(), mu,
                pop.gametes, pop.mutations, dip.first, mmodel,
                fwdpp::emplace_back());
            dip.second = fwdpp::mutate_gamete_recycle(
                mutation_recycling_bin, gamete_recycling_bin, rng.get(), mu,
                pop.gametes, pop.mutations, dip.second, mmodel,
                fwdpp::emplace_back());
        }
    fwdpp::fwdpp_internal::process_gametes(pop.gametes, pop.mutations,
                                           pop.mcounts);
    fwdpp::fwdpp_internal::gamete_cleaner(
        pop.gametes, pop.mutations, pop.mcounts, 2 * N_next, std::true_type());
    // This is constant-time
    pop.diploids.swap(offspring);
}

double
evolve_singlepop_regions_track_ancestry(
    const GSLrng_t& rng, singlepop_t& pop,
    const std::vector<std::uint32_t>& popsizes, const double mu_selected,
    const double recrate, const fwdpp::extensions::discrete_mut_model& mmodel,
    const fwdpp::extensions::discrete_rec_model& rmodel,
    const double selfing_rate)
{
    if (pop.generation > 0)
        {
            throw std::runtime_error(
                "this population has already been evolved.");
        }
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
    pop.mutations.reserve(
        std::ceil(std::log(2 * pop.N)
                  * (4. * double(pop.N) * (mu_selected)
                     + 0.667 * (4. * double(pop.N) * (mu_selected)))));

    for (unsigned generation = 0; generation < generations;
         ++generation, ++pop.generation)
        {
            const auto N_next = popsizes.at(generation);
            evolve_generation(rng, pop, N_next, mu_selected, mmodels, recmap);
            update_mutations(pop.mutations, pop.fixations, pop.fixation_times,
                             pop.mut_lookup, pop.mcounts, pop.generation,
                             2 * pop.N);
        }
}
