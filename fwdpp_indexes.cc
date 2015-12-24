#include <fwdpp/sugar/popgenmut.hpp>
#include <fwdpp/fwd_functional.hpp>
#include <fwdpp/internal/gsl_discrete.hpp>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unordered_set>
#include <type_traits>
#include <numeric>
#include <queue>
#include <functional>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
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
    assert(neutral.empty());
    assert(selected.empty());
  }
  gamete_t(unsigned n_,
	   const mvec_t & neut_,
	   const mvec_t & sel_) noexcept :
				 n(n_),neutral(neut_),selected(sel_)
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
				       diploids(dipvec_t(N, dip_t(0,0))),
				       mut_lookup(lookup_t())
  {
    gametes.reserve(2*N);
    assert(gametes[0].n==2*N);
#ifndef NDEBUG
    for( const auto & d : diploids )
      {
	assert( gametes[d.first].n==2*N );
	assert( gametes[d.second].n==2*N );
      }
#endif
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
		het(w,mutations[g2.selected[b2]]);
	      }
	  }
	if(!found) het(w,mutations[g1.selected[b1]]);
      }
    for( ; b2 < g2.selected.size() ; ++b2 ) het(w,mutations[g2.selected[b2]]);
    
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
  template<typename queue_t,
	   typename mvec_t,
	   typename lookup_table_t,
	   typename position_t,
	   typename sdist_t,
	   typename hdist_t>
  inline result_type operator()(queue_t & recycling_bin,
				mvec_t & mutations,
				gsl_rng * r,
				lookup_table_t * lookup,
				const unsigned & generation,
				const double & neutral_mutation_rate,
				const double & selected_mutation_rate,
				const position_t & posmaker,
				const sdist_t & smaker,
				const hdist_t & hmaker) const
  {
    static_assert( typename std::is_same<typename mvec_t::value_type,mtype>::type(),"foo");
    static_assert( typename std::is_same<typename queue_t::value_type,size_t>::type(),"foo");
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
	    assert(idx < mutations.size());
	    assert(!mutations[idx].n);
	    recycling_bin.pop();
	    mutations[idx].pos=pos;
	    mutations[idx].s=smaker();
	    mutations[idx].h=hmaker();
	    mutations[idx].n=1u;
	    mutations[idx].checked=false;
	    return idx;
	  }
	mutations.emplace_back(pos,smaker(),hmaker(),generation,1u);
	assert(mutations[mutations.size()-1].pos==pos);
	return mutations.size()-1;
      }
    //neutral mutation
    if ( ! recycling_bin.empty() )
      {
	auto idx = recycling_bin.front();
	assert(!mutations[idx].n);
	recycling_bin.pop();
	mutations[idx].pos=pos;
	mutations[idx].s=0.;
	mutations[idx].h=0.;
	mutations[idx].n=1u;
	mutations[idx].checked=false;
	return idx;
      }
    mutations.emplace_back(pos,0.,0.,generation,1u);
    assert(mutations[mutations.size()-1].pos==pos);
    return mutations.size()-1;
  }
};

template<typename queue_type,
	 typename mutation_model,
	 typename mvec_t>
void add_N_muts( queue_type & mut_recycling_bin,
		 const mutation_model & mmodel,
		 const unsigned & n,
		 mvec_t & mutations,
		 gamete_t & g)
{
  for(unsigned i=0;i<n;++i)
    {
      size_t idx = mmodel(mut_recycling_bin,mutations);
      std::cerr << "new mutation index = " << idx << ' ' << mutations.size() << ' ' << mutations.capacity() << '\n';
      if( mutations[idx].neutral )
	{
	  g.neutral.emplace(upper_bound(g.neutral.begin(),g.neutral.end(),mutations[idx].pos,
					[&mutations](const double & __value,const size_t & index) {
					  return __value < mutations[index].pos;
					}),
			    idx);
	}
      else
	{
	  g.selected.emplace(upper_bound(g.selected.begin(),g.selected.end(),mutations[idx].pos,
					 [&mutations](const double & __value,const size_t & index) {
					   return __value < mutations[index].pos;
					 }),
			     idx);
	}
    }
  assert( is_sorted( g.neutral.begin(),g.neutral.end(),[&mutations](size_t i,size_t j) {
 	return mutations[i].pos < mutations[j].pos;
     }) );
  assert( is_sorted( g.selected.begin(),g.selected.end(),[&mutations](size_t i,size_t j) {
 	return mutations[i].pos < mutations[j].pos;
      }) );
}

//The mutation fxn
template< typename queue_type,
	  typename queue_type2,
	  typename mutation_model,
	  typename gamete_insertion_policy,
	  typename gvec_t,
	  typename mvec_t >
size_t mut_recycle( queue_type & mut_recycling_bin,
		    queue_type2 & gamete_recycling_bin,
		    gsl_rng * r,
		    const double & mu,
		    gvec_t & gametes,
		    mvec_t & mutations,
		    const size_t & gamete_index,
		    const mutation_model & mmodel,
		    const gamete_insertion_policy & gpolicy)
{
  static_assert( typename std::is_same<typename queue_type::value_type,size_t>::type(),"foo");
  unsigned nm = gsl_ran_poisson(r,mu);
  if(!nm) return nm;
  assert(gametes[gamete_index].n);
  gametes[gamete_index].n--;
  if( ! gamete_recycling_bin.empty() )
    {
      auto idx = gamete_recycling_bin.front();
      assert(idx != gamete_index);
      assert(!gametes[idx].n);
      gamete_recycling_bin.pop();
      gametes[idx].n=1;
      gametes[idx].neutral=gametes[gamete_index].neutral;
      gametes[idx].selected=gametes[gamete_index].selected;
      add_N_muts(mut_recycling_bin,mmodel,nm,mutations,gametes[idx]);
      return idx;
    }
  typename gvec_t::value_type ng(1,gametes[gamete_index].neutral,gametes[gamete_index].selected);
  add_N_muts(mut_recycling_bin,mmodel,nm,mutations,ng);
  return gpolicy(gametes,std::move(ng));
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

  auto itr = gametes[ibeg].neutral.cbegin(),
    jtr = gametes[jbeg].neutral.cbegin(),
    itr_s = gametes[ibeg].selected.cbegin(),
    jtr_s = gametes[jbeg].selected.cbegin(),
    itr_e = gametes[ibeg].neutral.cend(),
    itr_s_e = gametes[ibeg].selected.cend(),
    jtr_e = gametes[jbeg].neutral.cend(),
    jtr_s_e = gametes[jbeg].selected.cend();

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
  assert( is_sorted( gametes[ibeg].neutral.begin(),gametes[ibeg].neutral.end(),[&mutations](size_t i,size_t j) {
 	return mutations[i].pos < mutations[j].pos;
      }) );
  assert( is_sorted( gametes[jbeg].selected.begin(),gametes[jbeg].selected.end(),[&mutations](size_t i,size_t j) {
 	return mutations[i].pos < mutations[j].pos;
      }) );
}
      
template< typename recombination_map,
	  typename gvec_t,
	  typename mvec_t,
	  typename glookup_t,
	  typename queue_t>
unsigned recombine_gametes( gsl_rng * r,
			    const double & littler,
			    gvec_t & gametes,
			    mvec_t & mutations,
			    size_t & g1,
			    const size_t & g2,
			    glookup_t & gamete_lookup,
			    queue_t & gamete_recycling_bin,
			    vector<size_t> & neutral,
			    vector<size_t> & selected,
			    const recombination_map & mf)
{
  unsigned nbreaks = (littler>0.)?gsl_ran_poisson(r,littler):0;
  cerr << "nb="<< nbreaks << ' ';
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
      auto lookup = gamete_lookup.lookup(neutral,selected,mutations);
      if(lookup.first!=lookup.second)
	{
	  auto itr = find_if(lookup.first,lookup.second,
			     [&gametes,&neutral,&selected](typename glookup_t::inner_t & __p )
			     {
			       return gametes[__p.second].neutral==neutral && gametes[__p.second].selected==selected;
			     });
	  if(itr==lookup.second)
	    {
	      if(!gamete_recycling_bin.empty())
		{
		  auto idx = gamete_recycling_bin.front();
		  assert(g1!=idx);
		  assert(!gametes[idx].n);
		  gamete_recycling_bin.pop();
		  gametes[idx].n=0u;
		  gametes[idx].neutral.swap(neutral);
		  gametes[idx].selected.swap(selected);
		  g1=idx;
		}
	      else
		{
		  gametes.emplace_back(0u,std::move(neutral),std::move(selected));
		  g1=gametes.size()-1;
		}
	      gamete_lookup.update(g1,gametes,mutations);
	    }
	  else
	    {
	      assert(gametes[itr->second].neutral==neutral);
	      assert(gametes[itr->second].selected==selected);
	      g1 = itr->second;
	    }
	}
      else //attempt gamete recycling
	{
	  if(!gamete_recycling_bin.empty())
	    {
	      auto idx = gamete_recycling_bin.front();
	      assert(g1!=idx);
	      assert(!gametes[idx].n);
	      gamete_recycling_bin.pop();
	      gametes[idx].n=0u;
	      gametes[idx].neutral.swap(neutral);
	      gametes[idx].selected.swap(selected);
	      g1=idx;
	    }
	  else
	    {
	      gametes.emplace_back(0u,std::move(neutral),std::move(selected));
	      g1=gametes.size()-1;
	    }
	  gamete_lookup.update(g1,gametes,mutations);
	}
    }
  
  return nbreaks;
}

template<typename mutation_t>
struct glookup_t
{
  using lookup_table_t = std::multimap<double,size_t>;
  using result_type = std::pair< lookup_table_t::iterator, lookup_table_t::iterator>;
  using inner_t = typename lookup_table_t::value_type;
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

  inline void update(size_t idx,
		     const std::vector<gamete_t> & gametes,
		     const std::vector<mutation_t> & mutations)
  {
    update_details(idx,gametes,mutations);
  };
  
  explicit glookup_t( const std::vector<gamete_t> & gametes,
		      const std::vector<mutation_t> & mutations )
  {
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
	    const double littler,
	    const mmodel & m,
	    const fmodel & f,
	    const recmodel & rec)
{
  assert(p.diploids.size()==N);
  auto mrec = make_mut_recycling_bin(p.mutations);
  auto grec = make_gamete_recycling_bin(p.gametes);
  auto glookup = make_gamete_lookup(p.gametes,p.mutations);
  std::cerr << p.mutations.size() << ' ' << mrec.size() << ' ' << p.gametes.size() << ' ' << grec.size() << '\n';
  double wbar = 0.;
  std::vector<double> fitnesses(p.diploids.size());
  for(unsigned i=0;i<N;++i)
    {
      fitnesses[i] = f(p.gametes[p.diploids[i].first],
		       p.gametes[p.diploids[i].second],
		       p.mutations);
      p.gametes[p.diploids[i].first].n = 0;
      p.gametes[p.diploids[i].second].n = 0;
      wbar += fitnesses[i];
    }
#ifndef NDEBUG
  set<size_t> gams_in_dips1;
  for(const auto & d : p.diploids )
    {
      gams_in_dips1.insert(d.first);
      gams_in_dips1.insert(d.second);
    }
  for(unsigned i=0;i<p.gametes.size();++i)
    {
      if(gams_in_dips1.find(i)==gams_in_dips1.end()) assert(!p.gametes[i].n);
      if(!p.gametes[i].n) cerr << i << ' ' << p.gametes[i].n << ' ' << p.gametes[i].neutral.size() << '\n';
      assert(!p.gametes[i].n);
    }
#endif
  wbar /= double(N);
  KTfwd::fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));

  auto  parents=p.diploids;
  static_assert(typename std::is_same<decltype(parents),decltype(p.diploids)>::type(),"foo");

  unsigned i = 0;
  //UPDATE THE DIPLOIDS
  for( ;i<N;++i )
    {
      size_t p1 = gsl_ran_discrete(r,lookup.get());
      size_t p2 = gsl_ran_discrete(r,lookup.get());
      assert(p1 < N);
      assert(p2 < N);
      size_t p1g1 = parents[p1].first;
      size_t p1g2 = parents[p1].second;
      size_t p2g1 = parents[p2].first;
      size_t p2g2 = parents[p2].second;
      assert(p1g1<2*N);
      assert(p1g2<2*N);
      assert(p2g1<2*N);
      assert(p2g2<2*N);
      if(gsl_rng_uniform(r)<0.5) swap(p1g1,p1g2);
      if(gsl_rng_uniform(r)<0.5) swap(p2g1,p2g2);

      std::cerr << "rec1: " << p.gametes[p1g1].neutral.size() << ' ' << p.gametes[p1g2].neutral.size() << ' ' << p1g1 << ' ' << p1g2;
      recombine_gametes(r,littler,p.gametes,p.mutations,
			p1g1,p1g2,glookup,grec,
			p.neutral,p.selected,rec);
      std::cerr << p1g1 << ' ' << p.gametes[p1g1].neutral.size() << '\n';
      assert(p.gametes.size()<=2*N);
      std::cerr << "rec2: " << p.gametes[p2g1].neutral.size() << ' ' << p2g1 << ' ';
      recombine_gametes(r,littler,p.gametes,p.mutations,
			p2g1,p2g2,glookup,grec,
			p.neutral,p.selected,rec);
      cerr << p2g1 << '\n';
      assert(p.gametes.size()<=2*N);
      p.diploids[i].first = p1g1;
      p.diploids[i].second = p2g1;
      p.gametes[p.diploids[i].first].n++;
      p.gametes[p.diploids[i].second].n++;
      assert(p.gametes[p.diploids[i].first].n <= 2*N);
      assert(p.gametes[p.diploids[i].second].n <= 2*N); 

      p.diploids[i].first = mut_recycle(mrec,grec,r,mutrate,p.gametes,p.mutations,p.diploids[i].first,m,
      					[](vector<gamete_t> & gams, gamete_t && g ) {
					  assert(!g.neutral.empty());
#ifndef NDEBUG
					  auto ss=gams.size();
					  auto g2 = g;
#endif
      					  gams.emplace_back(std::forward<gamete_t>(g));
#ifndef NDEBUG
					  assert(gams.size()-1==ss);
					  assert(gams[gams.size()-1].neutral==g2.neutral);
#endif
      					  return size_t(gams.size()-1);
      					});
      assert(p.gametes[p.diploids[i].first].n);
      assert(p.gametes.size()<=2*N);
      p.diploids[i].second = mut_recycle(mrec,grec,r,mutrate,p.gametes,p.mutations,p.diploids[i].second,m,
					 [](vector<gamete_t> & gams, gamete_t && g ) {
#ifndef NDEBUG
    auto ss=gams.size();
#endif
					   gams.emplace_back(std::forward<gamete_t>(g));
#ifndef NDEBUG
					   assert(gams.size()-1==ss);
#endif
					   return size_t(gams.size()-1);
					 });
      assert(p.gametes[p.diploids[i].second].n);
      assert(p.gametes.size()<=2*N);
    }
  assert(i==p.diploids.size());
#ifndef NDEBUG
  std::set<size_t> dgams;
  for( const auto & d : p.diploids )
    {
    dgams.insert(d.first);
    dgams.insert(d.second);
    assert( p.gametes[d.first].n>0 );
    assert( p.gametes[d.second].n>0 );
  }
  for( unsigned i = 0 ; i < p.gametes.size() ; ++i )
    {
    //AHA
    if( dgams.find(i) == dgams.end() )
      {
	if(p.gametes[i].n) cerr << "I = " << i << ' ' << p.gametes[i].neutral.size() << ' '
				<< p.mutations[p.gametes[i].neutral[0]].pos << ' '
				<< p.mutations[p.gametes[i].neutral[0]].g << ' '
				<< p.mutations[p.gametes[i].neutral[0]].n << ' '
				<< '(' << p.mutations[p.gametes[i].neutral[0]].checked << ") "
				<< (p.mut_lookup.find(p.mutations[p.gametes[i].neutral[0]].pos)==p.mut_lookup.end()) << '\n';
	assert (!p.gametes[i].n);
      }
    else assert(p.gametes[i].n);
  }
#endif
  assert(p.gametes.size()<=2*N);
  unsigned NN=0;
  for( auto & g : p.gametes )
    {
      NN+=g.n;
      if(g.n) {
      for( auto & m : g.neutral )
	{
	  if (!p.mutations[m].checked)
	    {
	      p.mutations[m].n=g.n;
	      p.mutations[m].checked=true;
	    }
	  else
	    {
	      p.mutations[m].n+=g.n;
	    }
	}
      for( auto & m : g.selected )
	{
	  if (!p.mutations[m].checked)
	    {
	      p.mutations[m].n=g.n;
	      p.mutations[m].checked=true;
	    }
	  else
	    {
	      p.mutations[m].n+=g.n;
	    }
	}
      }
    }
#ifndef NDEBUG
  if(NN!=2*N)
    {
      cerr << NN << ' ' << 2*N << '\n';
      for(const auto & g : p.gametes) cerr << g.n << ' ' << g.neutral.size() << ' ' << g.selected.size() << '\n';
    }
#endif
  assert(NN==2*N);
#ifndef NDEBUG
  NN=0;
  for(const auto & m : p.mutations) assert(m.n<=2*N);
  for( const auto & g : p.gametes )
    {
      if(g.n)
	{
	  for(const auto & m : g.neutral ) assert(p.mutations[m].n);
	  for(const auto & m : g.selected ) assert(p.mutations[m].n);
	}
    }
#endif
}

template<typename lookup_t>
void update_mutations( vector<mtype> & mutations, unsigned twoN,
		       lookup_t & lookup)
{
  for( auto & m : mutations )
    {
      if(m.n==twoN)
	{
	  //candidate for recycling
	  m.n=0;
	  lookup.erase(m.pos);
	}
      if(!m.checked) {
	lookup.erase(m.pos);
	m.n=0;
      }
      else
	m.checked=0;
    }
}
  
int main(int argc, char ** argv)
{
  unsigned N=10000;
  double mu = 0.01;
  double littler=mu;
  singlepop_t pop(N);
  pop.mutations.reserve( log(2*N)*4.*double(N)*mu + (2./3.)*4.*double(N)*mu );
  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, atoi(argv[1]) );

  for(unsigned gen=0;gen<10*N;++gen)
    {
      cout << gen << '\n';
      sample(r,pop,N,mu,littler,
	     std::bind(mpol(),
		       placeholders::_1,
		       placeholders::_2,
		       r,
		       &pop.mut_lookup,
		       gen,
		       mu,
		       0.0,
		       [&r](){return gsl_rng_uniform(r);},
		       [](){return 0.;},
		       [](){return 0.;}),
	     bind(mult_w(),std::placeholders::_1,std::placeholders::_2,std::placeholders::_3),
	     [&r](){return gsl_rng_uniform(r);});
#ifndef NDEBUG
      set<size_t> foo;
      for( const auto & d : pop.diploids ) {
	foo.insert(d.first);
	foo.insert(d.second);
	assert(pop.gametes[d.first].n);
	assert(pop.gametes[d.second].n);
      }
      for(size_t i=0;i<pop.gametes.size();++i)
	{
	  if(foo.find(i)==foo.end()) assert(!pop.gametes[i].n);
	}
      for( const auto & m : pop.mutations ) assert( pop.mut_lookup.find(m.pos) != pop.mut_lookup.end() );
#endif
      update_mutations(pop.mutations,2*N,pop.mut_lookup);
#ifndef NDEBUG
      for( const auto & m : pop.mutations )
	{
	  if(m.n)
	    assert( pop.mut_lookup.find(m.pos) != pop.mut_lookup.end() );
	  else
	    assert( pop.mut_lookup.find(m.pos) == pop.mut_lookup.end() );
	}
#endif
    }	 
}
