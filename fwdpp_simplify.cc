// Experimental integration of fwdpp engine
// with simplification algorithm from Ralph et al.
// This code starts by copy/pasting several data types
// and functions from the fwdpy11-based example for that paper.

#include <cmath>
#include <stdexcept>
#include <vector>
#include <fstream>
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
evolve_generation(const GSLrng_t& rng, slocuspop_t& pop,
                  const fwdpp::uint_t N_next, const double mu,
                  const mutation_model& mmodel,
                  const breakpoint_function& recmodel,
                  const fwdpp::uint_t generation, table_collection& tables,
                  table_simplifier& ancestry,
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

            tables.add_offspring_data(next_index_local, breakpoints,
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
            tables.add_offspring_data(next_index_local, breakpoints,
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

    //fwdpp::poisson_xover recmap(recrate, 0., 1.);
    //const auto bound_recmodel =[&rng,&recmap](){return recmap(rng.get());};
    const auto recmap
        = fwdpp::recbinder(fwdpp::poisson_xover(recrate, 0., 1.), rng.get());
    unsigned generation = 0;
    //fwdpp::infsites inf;
    // const auto mmodel
    //    = [&pop, &rng, &inf, &generation, mu_neutral, mu_selected](
    //        fwdpp::traits::recycling_bin_t<slocuspop_t::mcont_t>& recbin,
    //        slocuspop_t::mcont_t& mutations) {
    //          return inf(recbin, mutations, rng.get(), pop.mut_lookup,
    //                     generation, mu_neutral, mu_selected,
    //                     [&rng]() { return gsl_rng_uniform(rng.get()); },
    //                     []() { return 0.; }, []() { return 0.; });
    //      };

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
    double sort_time = 0.0;
    for (; generation < generations; ++generation)
        {
            const auto N_next = popsizes.at(generation);
            evolve_generation(rng, pop, N_next, mu_neutral + mu_selected,
                              mmodel, recmap, generation, tables, ancestry,
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
            for (auto i = tables.num_nodes() - 2 * pop.diploids.size();
                 i < tables.num_nodes(); ++i)
                {
                    assert(tables.node_table[i].generation == generation + 1);
                    samples.push_back(i);
                }
            auto mt = tables.mutation_table;
            auto nt = tables.node_table;
            auto et = tables.edge_table;
            std::ostringstream ofn;
            ofn << "simoutput/last_edges." << generation << ".bin";
            std::ofstream out;
            out.open(ofn.str().c_str());
            for (auto& e : et)
                {
                    out.write(reinterpret_cast<char*>(&e.parent),
                              sizeof(decltype(e.parent)));
                    out.write(reinterpret_cast<char*>(&e.child),
                              sizeof(decltype(e.child)));
                    out.write(reinterpret_cast<char*>(&e.left),
                              sizeof(decltype(e.left)));
                    out.write(reinterpret_cast<char*>(&e.right),
                              sizeof(decltype(e.right)));
                }
            std::int32_t done = -1;
            out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            out.close();
            ofn.str(std::string());
            ofn << "simoutput/last_nodes." << generation << ".bin";
            out.open(ofn.str().c_str());
            for (auto n : nt)
                {
                    out.write(reinterpret_cast<char*>(&n.id), 4);
                    out.write(reinterpret_cast<char*>(&n.generation),
                              sizeof(double));
                }
            out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            out.close();
            ofn.str(std::string());
            ofn << "simoutput/last_mutations." << generation << ".bin";
            out.open(ofn.str().c_str());
            for (auto& m : tables.mutation_table)
                {
                    out.write(reinterpret_cast<char*>(&m.node),
                              sizeof(std::int32_t));
                    out.write(
                        reinterpret_cast<char*>(&pop.mutations[m.key].pos),
                        sizeof(double));
                }
            out.write(reinterpret_cast<char*>(&done), sizeof(std::int32_t));
            out.close();
            tables.sort_tables(pop.mutations);
            auto xx = ancestry.simplify(tables, samples, pop.mutations);

            // TODO: decide how to handle fixations.
            // Ideally, this would be done during simplification,
            // via a policy passed into the simplifier.
            // Until then, we do it manually
            tables.mutation_table.erase(
                std::remove_if(tables.mutation_table.begin(),
                               tables.mutation_table.end(),
                               [&xx,&pop](const fwdpp::ancestry::mutation_record& mr) {
                               return xx.second[mr.key] == 2*pop.diploids.size();
                               }),
                tables.mutation_table.end());

            unsigned n = 0;
            for (std::size_t i = 0; i < pop.mutations.size(); ++i)
                {
                    if (pop.mcounts[i])
                        {
                            ++n;
                        }
                }
            std::cerr << "result "
                      << int(std::accumulate(pop.mcounts.begin(),
                                             pop.mcounts.end(), 0))
                      << ' '
                      << int(std::accumulate(xx.second.begin(),
                                             xx.second.end(), 0))
                      << ' ' << n << ' ' << tables.mutation_table.size() << ' '
                      << pop.fixations.size() << '\n';
            assert(pop.mcounts.size() == xx.second.size());
            //for (auto& mr : tables.mutation_table)
            //    {
            //        assert(tables.node_table.at(mr.node).generation
            //               == 1. + pop.mutations[mr.key].g);
            //    }
            //output last set of nodes/edges
            ofn.str(std::string());
            ofn << "simoutput/edges." << generation << ".txt";

            out.open(ofn.str().c_str());
            for (auto e : tables.edge_table)
                {
                    out << e.parent << ' ' << e.child << ' ' << e.left << ' '
                        << e.right << '\n';
                }
            out.close();
            ofn.str(std::string());
            ofn << "simoutput/nodes." << generation << ".txt";
            out.open(ofn.str().c_str());
            for (auto n : tables.node_table)
                {
                    out << n.id << ' ' << n.generation << '\n';
                }
            out.close();
            ofn.str(std::string());
            ofn << "simoutput/idmap." << generation << ".txt";
            out.open(ofn.str().c_str());
            for (std::size_t id = 0; id < xx.first.size(); ++id)
                {
                    out << id << ' ' << xx.first[id] << '\n';
                }
            out.close();

            if (pop.mcounts != xx.second)
                {
                    std::vector<std::size_t> failures;
                    for (std::size_t i = 0; i < pop.mcounts.size(); ++i)
                        {
                            std::cout << generation << ' ' << pop.mcounts[i]
                                      << ' ' << xx.second[i] << ' '
                                      << pop.mutations[i].pos << ' '
                                      << pop.mutations[i].g << '\n';
                            if (pop.mcounts[i] != xx.second[i])
                                {
                                    failures.push_back(i);
                                }
                        }
                    //Find diploids with failures
                    std::vector<std::size_t> gfails;
                    for (std::size_t g = 0; g < pop.gametes.size(); ++g)
                        {
                            auto& gam = pop.gametes[g];
                            if (gam.n)
                                {
                                    for (auto f : failures)
                                        {
                                            auto itr = std::find(
                                                gam.mutations.begin(),
                                                gam.mutations.end(), f);
                                            auto itr2 = std::find(
                                                gam.smutations.begin(),
                                                gam.smutations.end(), f);
                                            if (itr != gam.mutations.end()
                                                || itr2
                                                       != gam.smutations.end())
                                                {
                                                    gfails.push_back(g);
                                                }
                                        }
                                }
                        }
                    std::sort(gfails.begin(), gfails.end());
                    gfails.erase(std::unique(gfails.begin(), gfails.end()),
                                 gfails.end());
                    for (std::size_t dip = 0; dip < pop.diploids.size(); ++dip)
                        {
                            if (std::find(gfails.begin(), gfails.end(),
                                          pop.diploids[dip].first)
                                != gfails.end())
                                {
                                    std::cout << "Node with mut: " << 2 * dip
                                              << '\n';
                                }
                            if (std::find(gfails.begin(), gfails.end(),
                                          pop.diploids[dip].second)
                                != gfails.end())
                                {
                                    std::cout
                                        << "Node with mut: " << 2 * dip + 1
                                        << '\n';
                                }
                        }
                    for (unsigned i = 0; i < mt.size(); ++i)
                        {
                            std::cout << mt[i].node << ' '
                                      << xx.first[mt[i].node] << ' '
                                      << pop.mutations[mt[i].key].pos << '\n';
                        }
                    std::exit(0);
                }
            std::cerr << (pop.mcounts == xx.second) << '\n';
            next_index = tables.num_nodes();
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
