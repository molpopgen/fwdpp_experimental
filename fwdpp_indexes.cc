#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/fwd_functional.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unordered_set>
#include <type_traits>
#include <numeric>
#include <queue>
#include <cassert>
#include <map>
using namespace std;

/*
  Use exising mutation type
*/
using mtype = KTfwd::popgenmut;

/*
  Gamete and population types change
*/
struct gamete_t
{
  unsigned n;
  using mvec_t = vector<size_t>;
  mvec_t neutral,selected;
  gamete_t( unsigned n_ ) noexcept : n(n_),
				     neutral(mvec_t()),
				     selected(mvec_t())
  {
  }
  
};

struct singlepop_t
/*
  Note: this is now trivially copyable w/o any fancy serialization!
*/
{
  using mvec_t = vector<mtype>;
  using gvec_t = vector<gamete_t>;
  using dip_t = pair<size_t,size_t>;
  using dipvec_t = vector<dip_t>;
  using lookup_t = std::unordered_set<double,std::hash<double>,KTfwd::equal_eps>;
  mvec_t mutations;
  gvec_t gametes;
  dipvec_t diploids;
  lookup_t mut_lookup;
  vector<size_t> neutral,selected;
  singlepop_t( unsigned N ) noexcept : mutations(mvec_t()),
				       gametes(gvec_t(1,gamete_t(2*N))),
				       diploids(dipvec_t(2*N, dip_t(0,0))),
				       mut_lookup(lookup_t())
  {
    gametes.reserve(2*N);
    neutral.reserve(100);
    selected.reserve(100);
  }
};

// Fitness scheme.  This is a re-working of KTfwd::site_dependent_fitness
// This code is a LOT simpler...
struct site_dep_w
{
  using result_type = double;
  template<typename het_,
	   typename hom_,
	   typename mvec_t>
  inline result_type operator()(const gamete_t & g1,
				const gamete_t & g2,
				const mvec_t & mutations,
				const het_ & het,
				const hom_ & hom,
				const result_type starting_fitness = 1.) const
  {
    result_type w = starting_fitness;
    if( g1.selected.empty() && g2.selected.empty() ) return w;
    else if (g1.selected.empty()) {
      for( const auto & i : g1.selected ) het( w,mutations[i] );
      return w;
    } else if (g2.selected.empty()) {
      for( const auto & i : g2.selected ) het( w,mutations[i] );
      return w;
    }

    std::vector<size_t>::size_type b1 = 0, b2 = 0;
    bool found = false;
    for( ; b1 < g1.selected.size() ; ++b1 )
      {
	for( ; !found && b2 < g2.selected.size() && !( mutations[g2.selected[b2]].pos > mutations[g1.selected[b1]].pos ) ; ++b2 )
	  {
	    if (g1.selected[b1]==g2.selected[b2])
	      {
		hom(w,mutations[g1.selected[b1]]);
		found=true;
	      }
	    else
	      {
		het(w,g2.selected[b2]);
	      }
	  }
	if(!found) het(w,g1.selected[b1]);
      }
    for( ; b2 < g2.selected.size() ; ++b2 ) het(w,g2.selected[b2]);
    
    return w;
  }
};

//reworking of KTfwd::multiplicative_fitness
struct mult_w
{
  using result_type = double;
  template<typename mvec_t>
  inline result_type operator()( const gamete_t & g1,
				 const gamete_t & g2,
				 const mvec_t & mutations,
				 const double scaling = 1.) const
  {
    using mtype = typename mvec_t::value_type;
    return std::max(0.,site_dep_w()(g1,g2,mutations,
				    [&scaling](result_type & w,const mtype & m)
				    {
				      w *= (1.+scaling*m.s);
				    },
				    [&scaling](result_type & w,const mtype & m)
				    {
				      w *= (1 + m.h*m.s);
				    }));
  }
};

//Mutation policy.  Returns a new mtype
//This is NOT libary quality!!!
struct mpol
{
  using result_type = size_t;
  template<typename pos_t,
	   typename queue_t,
	   typename mvec_t,
	   typename lookup_table_t,
	   typename position_t,
	   typename sdist_t,
	   typename hdist_t>
  inline result_type operator()(queue_t & recycling_bin,
				mvec_t & mutations,
				gsl_rng * r, lookup_table_t * lookup,
				const unsigned & generation,
				const double & neutral_mutation_rate,
				const double & selected_mutation_rate,
				const position_t & posmaker,
				const sdist_t & smaker,
				const hdist_t & hmaker) const
  {
    auto pos = posmaker();
    while(lookup->find(pos) != lookup->end())
      {
	pos = posmaker();
      }
    lookup->insert(pos);
    bool selected = (gsl_rng_uniform(r) < selected_mutation_rate/(neutral_mutation_rate + selected_mutation_rate));
    if ( selected )
      {
	if ( ! recycling_bin.empty() )
	  {
	    auto idx = recycling_bin.front();
	    recycling_bin.pop();
	    mutations[idx] = typename mvec_t::value_type(pos,smaker(),hmaker(),generation,1u);
	    return idx;
	  }
	mutations.emplace_back(pos,smaker(),hmaker(),generation,1u);
	return mutations.size()-1;
      }
    //neutral mutation
    if ( ! recycling_bin.empty() )
      {
	auto idx = recycling_bin.front();
	recycling_bin.pop();
	mutations[idx] = typename mvec_t::value_type(pos,0.,0.,generation,1u);
	return idx;
      }
    mutations.emplace_back(pos,0.,0.,generation,1u);
    return mutations.size()-1;
  }
};

template<typename queue_type,
	 typename mutation_model,
	 typename mvec_t>
void add_N_muts( queue_type & recycling_bin,
		 const mutation_model & mmodel,
		 const unsigned & n,
		 mvec_t & mutations,
		 gamete_t & g)
{
  for(unsigned i=0;i<n;++i)
    {
      size_t idx = mmodel(recycling_bin,mutations);
      if( mutations[idx].neutral )
	{
	  g.neutral.emplace(upper_bound(g.neutral.begin(),g.neutral.end(),mutations[idx].pos,
					[&mutations](const double & __value,const size_t & i) {
					  return __value < mutations[i].pos;
					}),
			    idx);
	}
      else
	{
	  g.selected.emplace(upper_bound(g.selected.begin(),g.selected.end(),mutations[idx].pos,
					 [&mutations](const double & __value,const size_t & i) {
					   return __value < mutations[i].pos;
					 }),
			     idx);
	}
    }
}

//The mutation fxn
template< typename queue_type,
	  typename queue_type2,
	  typename mutation_model,
	  typename gamete_insertion_policy,
	  typename gvec_t,
	  typename mvec_t >
size_t mut_recycle( queue_type & recycling_bin,
		    queue_type2 & gamete_recycling_bin,
		    gsl_rng * r,
		    const double & mu,
		    gvec_t & gametes,
		    mvec_t & mutations,
		    size_t & gamete_index,
		    const mutation_model & mmodel,
		    const gamete_insertion_policy & gpolicy)
{
  unsigned nm = gsl_ran_poisson(r,mu);
  if(!nm) return nm;
  gametes[gamete_index].n--;
  if( ! gamete_recycling_bin.empty() )
    {
      auto idx = gamete_recycling_bin.front();
      gamete_recycling_bin.pop();
      gametes[idx].n=1;
      gametes[idx].neutral=gametes[gamete_index].neutral;
      gametes[idx].selected=gametes[gamete_index].selected;
      add_N_muts(recycling_bin,mmodel,nm,mutations,gametes[idx]);
      return idx;
    }
  typename gvec_t::value_type ng(1,gametes[gamete_index].neutral,gametes[gamete_index].selected);
  add_N_muts(recycling_bin,mmodel,nm,mutations,ng);
  return gpolicy(gametes,ng);
}

//create recycling bins
queue<size_t> make_mut_recycling_bin( const std::vector<mtype> & mutations )
{
  queue<size_t> rv;
  for(size_t i=0;i<mutations.size();++i)
    {
      if(!mutations[i].n && !mutations[i].checked) rv.push(i);
    }
  return rv;
}

queue<size_t> make_gamete_recycling_bin( const std::vector<gamete_t> & gametes )
{
  queue<size_t> rv;
  for(size_t i=0;i<gametes.size();++i)
    {
      if(!gametes[i].n) rv.push(i);
    }
  return rv;
}

//Crossing-over
namespace internal
{
  template<typename itr_type,
	   typename mvec_t>
  inline itr_type rec_update_itr( itr_type __first,
				  itr_type __last,
				  const mvec_t & mutations,
				  const double & val)
  {
    if(__first==__last) return __first;
    return std::upper_bound(__first,__last,
			    std::cref(val),
			    [&mutations](const double __val,const size_t __mut) {
			      return __val < mutations[__mut].pos;
			    });
  }
    
  template< typename itr_type,
	    typename mvec_t>
  itr_type rec_gam_updater( itr_type __first, itr_type __last,
			    const mvec_t & mutations,
			    vector<size_t> & muts,
			    const double & val )
  {
    //O(log_2) comparisons of double plus at most __last - __first copies
    itr_type __ub = rec_update_itr(__first,__last,mutations,val);
    /*
      NOTE: the use of insert here
      instead of std::copy(__first,__ub,std::back_inserter(muts));
      Reduced peak RAM use on GCC 4.9.2/Ubuntu Linux.
    */
    muts.insert(muts.end(),__first,__ub);
    return __ub;
  }
}
template<typename gvec_t,
	 typename mvec_t>
void recombine_gametes_details( const vector<double> & pos,
				const size_t ibeg,
				const size_t jbeg,
				const gvec_t & gametes,
				const mvec_t & mutations,
				vector<size_t> & neutral,
				vector<size_t> & selected)
{
  assert( std::is_sorted(pos.cbegin(),pos.cend()) );
  short SWITCH = 0;

  auto itr = gametes[ibeg].mutations.cbegin(),
    jtr = gametes[jbeg].mutations.cbegin(),
    itr_s = gametes[ibeg].smutations.cbegin(),
    jtr_s = gametes[jbeg].smutations.cbegin(),
    itr_e = gametes[ibeg].mutations.cend(),
    itr_s_e = gametes[ibeg].smutations.cend(),
    jtr_e = gametes[jbeg].mutations.cend(),
    jtr_s_e = gametes[jbeg].smutations.cend();

  for(const double dummy : pos )
    {
      if(!SWITCH)
	{
	  itr = internal::rec_gam_updater(itr,itr_e,
					  mutations,
					  neutral,dummy);
	  itr_s = internal::rec_gam_updater(itr_s,itr_s_e,
					    mutations,
					    selected,dummy);
	  jtr = internal::rec_update_itr(jtr,jtr_e,mutations,dummy);
	  jtr_s = internal::rec_update_itr(jtr_s,jtr_s_e,mutations,dummy);
	}
      else
	{
	  jtr = internal::rec_gam_updater(jtr,jtr_e,
					  mutations,
					  neutral,dummy);
	  jtr_s = internal::rec_gam_updater(jtr_s,jtr_s_e,
					    mutations,
					    selected,dummy);
	  itr = internal::rec_update_itr(itr,itr_e,mutations,dummy);
	  itr_s = internal::rec_update_itr(itr_s,itr_s_e,mutations,dummy);
	}
      SWITCH=!SWITCH;
    }
}
      
template< typename iterator_type,
	  typename recombination_map,
	  typename gvec_t,
	  typename mvec_t,
	  typename glookup_t,
	  typename queue_t>
unsigned recombine_gametes( gsl_rng * r,
			    const double & littler,
			    gvec_t & gametes,
			    mvec_t & mutations,
			    size_t & g1,
			    size_t & g2,
			    glookup_t & gamete_lookup,
			    queue_t & gamete_recycling_bin,
			    vector<size_t> & neutral,
			    vector<size_t> & selected,
			    const recombination_map & mf)
{
  unsigned nbreaks = (littler>0.)?gsl_ran_poisson(r,littler):0;
  if( nbreaks )
    {
      std::vector<double> pos;
      pos.reserve(nbreaks+1);
      for(unsigned i = 0 ; i < nbreaks ; ++i)
	{
	  pos.emplace_back(mf());
	}
      std::sort(pos.begin(),pos.end());
      pos.emplace_back(std::numeric_limits<double>::max());
      neutral.clear();
      selected.clear();
      recombine_gametes_details(pos,g1,g2,gametes,mutations,neutral,selected);

      //Gotta do lookup table thing here
    }
  
  return nbreaks;
}

template<typename mutation_t>
struct glookup_t
{
  using lookup_table_t = std::multimap<double,size_t>;
  using result_type = std::pair< lookup_table_t::iterator, lookup_table_t::iterator>;
  lookup_table_t lookup_table;

  inline double keyit( const std::vector<size_t> & mc,
		       const std::vector<mutation_t> & mutations ) const
  {
    return (mc.empty()) ? -std::numeric_limits<double>::max() : mutations[mc[0]].pos;
  }

  inline void update_details( size_t g,
			      const std::vector<gamete_t> & gametes,
			      const std::vector<mutation_t> & mutations) 
  {
    lookup_table.emplace( std::make_pair( keyit(gametes[g].neutral,mutations)*double(gametes[g].neutral.size()) +
					  keyit(gametes[g].selected,mutations)*double(gametes[g].selected.size()), g) );
  }

  inline result_type lookup( const std::vector<size_t> & n,
			     const std::vector<size_t> & s,
			     const std::vector<mutation_t> & mutations )
  {
    return lookup_table.equal_range(  keyit(n,mutations)*double(n.size()) + keyit(s,mutations)*double(s.size()) );
  }
  
  explicit glookup_t( const std::vector<gamete_t> & gametes,
		      const std::vector<mutation_t> & mutations )
  {
    //for(auto g = gametes.begin();g!=gametes.end();++g)
    for(size_t g=0;g<gametes.size();++g)
      {
	if(gametes[g].n)
	  {
	    update_details(g,gametes,mutations);
	  }
      }
  }
};

template<typename mutation_t>
glookup_t<mutation_t> make_gamete_lookup( const std::vector<gamete_t> & gametes,
					  const std::vector<mutation_t> & mutations )
{
  return glookup_t<mutation_t>(gametes,mutations);
}

template<typename poptype,
	 typename mmodel,
	 typename fmodel,
	 typename recmodel>
void sample(gsl_rng * r,
	    poptype & p,
	    const unsigned N,
	    const double mutrate,
	    const mmodel & m,
	    const fmodel & f,
	    const recmodel & rec)
{
  auto mrec = make_mut_recycling_bin(p.mutations);
  auto grec = make_gamete_recycling_bin(p.gametes);
  auto glookup = make_gamete_lookup(p.gametes,p.mutations);
}


  
int main(int argc, char ** argv)
{
  unsigned N=10000;
  double mu = 0.01;
  singlepop_t pop(N);
  gsl_rng * r;

  for(unsigned gen=0;gen<10*N;++gen)
    {
      sample(r,pop,N,mu,
	     bind(mpol(),placeholders::_1,placeholders::_2,
		  gen,mu,0.,[&r](){return gsl_rng_uniform(r);},
		  [](){return 0.;},[](){return 0.;}),
	     bind(mult_w(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
	     [&r](){return gsl_rng_uniform(r);});
    }	 
}
