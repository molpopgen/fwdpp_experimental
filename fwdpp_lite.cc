#include <memory>
#include <vector>
#include <list>
#include <functional>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <unordered_set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

/*
  Base types
*/
struct mutation
{
  double pos,s,h;
  unsigned gen;
  bool neutral,fixed;
  mutation(const double & __p,
	   const double & __s,
	   const double & __h,
	   const unsigned & __g) : pos(__p),s(__s),h(__h),gen(__g),neutral(s==0.),fixed(false)
  {
  }

  //Taming compiler errors
  mutation(mutation * m) : pos(m->pos),s(m->s),h(m->h),gen(m->gen)
  {
  }
};

using mptr_t = shared_ptr<mutation>;

struct gamete
{
  using mvector_t = vector<mptr_t>;
  mvector_t mutations,smutations;
  gamete() : mutations(mvector_t()),smutations(mvector_t())
  {
    mutations.reserve(1000);
    smutations.reserve(1000);
  }
  gamete( const mvector_t & n,
	  const mvector_t & s ) : mutations(n), smutations(s)
  {
    mutations.reserve(1000);
    smutations.reserve(1000);
  }

  //constructors req'd to tame compiler errors
  gamete( gamete * g) :  mutations(std::move(g->smutations)),
			 smutations(std::move(g->smutations))
  {
    mutations.reserve(1000);
    smutations.reserve(1000);
  }
			 
};

using gptr_t = shared_ptr<gamete>;

//Below, stuff is basically copied right from fwdpp, but is hacked to avoid tag/dispatch stuff

//Policy copyies
// template<typename T, typename cT>
// inline void push_at_end(  T && t,  cT * ct )
// {
//   ct->emplace_back(std::forward<T>(t));
// }

  struct mutation_remover
  {
    using result_type = bool;
    template<typename iterator_type> 
    inline result_type operator()(const iterator_type & i,
				  const unsigned & x1 ) const
    {
      return i.use_count() == x1;
      //return i->n == x1;
    }
    template<typename iterator_type>
    inline result_type operator()(const iterator_type & i,
				  const unsigned & x1,
				  const unsigned & x2) const
    {
      //PROBLEMO
      return (i.use_count() == x1) ||( i->fixed);//i.use_count() == x2;
      //return i->n == x1 || i->n == x2;
    }
  };


//Mutation-related stuff

template< typename gamete_type,
	  typename mutation_iterator>
void add_new_mutation( mutation_iterator mitr,
		       gamete_type & new_gamete )
{
  if(mitr->neutral)
    {
       new_gamete->mutations.insert(std::lower_bound(new_gamete->mutations.begin(),
						     new_gamete->mutations.end(),mitr->pos,
						     [](const mutation_iterator & __mut,const double & __value){ return __mut->pos < __value;}),
				    std::move(mitr));
				   //std::forward<mutation_iterator>(mitr) );
       //new_gamete->mutations.emplace_back(std::move(mitr));
    }
  else
    {
      new_gamete->smutations.insert(std::lower_bound(new_gamete->smutations.begin(),
						     new_gamete->smutations.end(),mitr->pos,
						     [](const mutation_iterator & __mut,const double & __value){ return __mut->pos < __value;}),
       				    std::move(mitr));
				    //std::forward<mutation_iterator>(mitr) );
      //new_gamete->smutations.emplace_back(std::move(mitr));
    }
}

template<typename mutation_model,
	 typename mutation_insertion_policy,
	 typename mlist_type,
	 typename gamete_type>
void add_N_mutations( const mutation_model & mmodel,
		      const mutation_insertion_policy & mpolicy,
		      const unsigned & n,
		      mlist_type * mutations,
		      gamete_type & g)
{
  for( unsigned i = 0 ; i < n ; ++i )
    {
      //ALLOCATE
      //auto nmut_p = typename mlist_type::value_type(new mutation(mmodel()));
      auto nmut_p = std::allocate_shared<decltype(mmodel())>(std::allocator<decltype(mmodel())>(),mmodel());
      auto mitr = mpolicy(std::move(nmut_p),mutations);
      add_new_mutation(*mitr,g);
    }
}

template< typename gamete_t,
	  typename mutation_model,
	  typename mvector_t,
	  typename gvector_t,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy>
gamete_t mutate_gamete( gsl_rng * r,
			const double & mu, 
			gvector_t * gametes,
			mvector_t * mutations, 
			gamete_t & g,
			const mutation_model &mmodel,
			const mutation_insertion_policy & mpolicy,
			const gamete_insertion_policy & gpolicy)
{
  unsigned nm = gsl_ran_poisson(r,mu);
  //PROBLEM
  //nm = 5;
  if ( nm )
    {
      //ALLOCATE
      //gamete_t ng(new gamete(*g));
      gamete_t ng = std::allocate_shared<gamete>(std::allocator<gamete>(),gamete(*g));

      //"decrement the count"
      g.reset();
      add_N_mutations(mmodel,mpolicy,nm,mutations,ng);
      //Forget insert, swap! -- SLOWER!!!
      // std::sort( ng->mutations.begin(),ng->mutations.end(),
      // 		 []( mptr_t & a, mptr_t & b ) { return a->pos < b->pos; } );
      // std::sort( ng->smutations.begin(),ng->smutations.end(),
      // 		 []( mptr_t & a, mptr_t & b ) { return a->pos < b->pos; } );
      return gpolicy(std::move(ng),gametes);
    }
  return g;
}

//How to recombine


namespace fwdpp_internal {

  template< typename itr_type,
	    typename cont_type >
  itr_type rec_gam_updater( itr_type & __first, itr_type & __last,
			    cont_type & m1, cont_type & m2,
			    const short & SWITCH, const double & val )
  {
    //O(log_2) comparisons of double plus at most __last - __first copies
    itr_type __ub = std::lower_bound(__first,__last,
				     std::cref(val),
				     [](const typename itr_type::value_type & __mut,const double & __val) {
				       return __mut->pos < __val;
				     });
    if (SWITCH)
      std::copy(__first,__ub,std::back_inserter(m1));
    else
      std::copy(__first,__ub,std::back_inserter(m2));
    return __ub;
  }

  template<typename double_vec_type,
	   typename gamete_type,
	   typename gamete_cont_iterator >
  void recombine_gametes( const double_vec_type & pos,
			  gamete_cont_iterator & ibeg,
			  gamete_cont_iterator & jbeg,
			  gamete_type & new_gamete1,
			  gamete_type & new_gamete2 )
  {
    assert( std::is_sorted(pos.cbegin(),pos.cend()) );
    short SWITCH = 0;

    auto itr = ibeg->mutations.cbegin(),
      jtr = jbeg->mutations.cbegin(),
      itr_s = ibeg->smutations.cbegin(),
      jtr_s = jbeg->smutations.cbegin(),
      itr_e = ibeg->mutations.cend(),
      itr_s_e = ibeg->smutations.cend(),
      jtr_e = jbeg->mutations.cend(),
      jtr_s_e = jbeg->smutations.cend();
#ifndef NDEBUG
    decltype(ibeg->mutations.size()) nm1=ibeg->mutations.size()+ibeg->smutations.size(),
      nm2=jbeg->mutations.size()+jbeg->smutations.size();
#endif
    for( const auto & dummy : pos )
      {
	//cout << dummy << ' ';
	itr = fwdpp_internal::rec_gam_updater(itr,itr_e,
					      new_gamete2->mutations,new_gamete1->mutations,SWITCH,dummy);
	itr_s = fwdpp_internal::rec_gam_updater(itr_s,itr_s_e,
						new_gamete2->smutations,new_gamete1->smutations,SWITCH,dummy);
	jtr = fwdpp_internal::rec_gam_updater(jtr,jtr_e,
					      new_gamete1->mutations,new_gamete2->mutations,SWITCH,dummy);
	jtr_s = fwdpp_internal::rec_gam_updater(jtr_s,jtr_s_e,
						new_gamete1->smutations,new_gamete2->smutations,SWITCH,dummy);
	SWITCH=!SWITCH;
      }
#ifndef NDEBUG
    decltype(new_gamete1->mutations.size()) __nm1 = new_gamete1->mutations.size()+new_gamete1->smutations.size(),
      __nm2 = new_gamete2->mutations.size()+new_gamete2->smutations.size();
    /*
    if (!(__nm1+__nm2 == nm1+nm2) )
      {
	std::cerr << nm1 << " -> " << __nm1 
		  << ", "
		  << nm2 << " -> " << __nm2 << '\n';
      }
    */
    assert(__nm1+__nm2 == nm1+nm2);
#endif

    //Through fwdpp 0.2.4, we did a sort here, but it is not necessary
    /*
      std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),
      [](typename gamete_type::mutation_list_type_iterator lhs,
      typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
      std::sort(new_gamete1.smutations.begin(),new_gamete1.smutations.end(),
      [](typename gamete_type::mutation_list_type_iterator lhs,
      typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
      std::sort(new_gamete2.mutations.begin(),new_gamete2.mutations.end(),
      [](typename gamete_type::mutation_list_type_iterator lhs,
      typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
      std::sort(new_gamete2.smutations.begin(),new_gamete2.smutations.end(),
      [](typename gamete_type::mutation_list_type_iterator lhs,
      typename gamete_type::mutation_list_type_iterator rhs) { return lhs->pos < rhs->pos; });
    */
#ifndef NDEBUG
    //using mlist_itr = typename gamete_type::mutation_list_type_iterator;
    //PROBLEM!
    auto am_I_sorted = [](mptr_t & lhs, mptr_t & rhs){return lhs->pos < rhs->pos;};
    assert( std::is_sorted(new_gamete1->mutations.begin(),new_gamete1->mutations.end(),std::cref(am_I_sorted)) );
    assert( std::is_sorted(new_gamete1->smutations.begin(),new_gamete1->smutations.end(),std::cref(am_I_sorted)) );
    assert( std::is_sorted(new_gamete2->mutations.begin(),new_gamete2->mutations.end(),std::cref(am_I_sorted)) );
    assert( std::is_sorted(new_gamete2->smutations.begin(),new_gamete2->smutations.end(),std::cref(am_I_sorted)) );
#endif
  }
    
}

template< typename iterator_type,
	  //typename list_type_allocator,
	  typename vector_type_allocator,
	  template<typename,typename> class vector_type,
	  typename gvector_t>
	  //template<typename,typename> class list_type>
unsigned recombine_gametes_details( const vector_type< double, vector_type_allocator > & pos,
			    //list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			    gvector_t * gametes,
			    iterator_type & g1,
			    iterator_type & g2)
{
  //assert( g1 != gametes->end() );
  //assert( g2 != gametes->end() );
  assert( std::is_sorted(pos.begin(),pos.end()) );
  assert( *(pos.end()-1) == std::numeric_limits<double>::max() );

  //using gtype = typename iterator_type::value_type;
  //using gtype_mcont = typename gtype::mutation_container;
    
  //ALLOCATE
  //iterator_type new_gamete1( new gamete() ), new_gamete2( new gamete() );
  iterator_type new_gamete1 = std::allocate_shared<gamete>(std::allocator<gamete>(),gamete()),
    new_gamete2 = std::allocate_shared<gamete>(std::allocator<gamete>(),gamete());
  //gtype new_gamete1(0u,gtype_mcont(),gtype_mcont()),
  //  new_gamete2(new_gamete1);
    
  new_gamete1->mutations.reserve(g1->mutations.size()+g2->mutations.size());
  new_gamete1->smutations.reserve(g1->smutations.size()+g2->smutations.size());
  new_gamete2->mutations.reserve(g1->mutations.size()+g2->mutations.size());
  new_gamete2->smutations.reserve(g1->smutations.size()+g2->smutations.size());
	
  fwdpp_internal::recombine_gametes(pos,g1,g2,new_gamete1,new_gamete2);
  //cout << '\n';
  auto current_end = gametes->end();
  bool f1 = false, f2 = false;
  for( auto itr = gametes->begin() ;
       (!f1||!f2)&&itr != current_end ; ++itr )
    {
      if(!f1&&*itr == new_gamete1)
	{
	  g1=*itr;
	  f1=true;
	}
      if(!f2&&*itr == new_gamete2)
	{
	  g2=*itr;
	  f2=true;
	}
    }
  if(!f1)
    {
      g1=*gametes->insert(gametes->end(),std::move(new_gamete1));
    }
  if(!f2)
    {
      g2=*gametes->insert(gametes->end(),std::move(new_gamete2));
    }
  return pos.size()-1;
}

//recombination for individual-based simulation
template< typename iterator_type,
	  typename recombination_map,
	  //typename list_type_allocator,
	  typename gvector_t>
	  //template<typename,typename> class list_type>
unsigned recombine_gametes( gsl_rng * r,
			    const double & littler,
			    //list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
			    gvector_t * gametes,
			    iterator_type & g1,
			    iterator_type & g2,
			    const recombination_map & mf)
{
  //assert( g1 != gametes->end() );
  //assert( g2 != gametes->end() );
    
  //Identify cases where recombination cannot result in changed gametes, and get out quick
  if(g1 == g2 ) return 0;
  auto nm1=g1->mutations.size()+g1->smutations.size();
  auto nm2=g2->mutations.size()+g2->smutations.size();
  if((std::min(nm1,nm2)==0 && std::max(nm1,nm2)==1)) return 0;
    
  unsigned nbreaks = (littler > 0) ? gsl_ran_poisson(r,littler) : 0u;
  //PROBLEM
  //nbreaks = 3;
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
      return recombine_gametes_details(pos,gametes,g1,g2);
    }
  return 0;
}

struct genetics101
/*! Genetics 101: simple model of recombination.  r is the probability that the two gametes recombine
 */
{
  using result_type = unsigned;
  template<typename gamete_iterator_type,
	   typename gvector_t,
	   //typename gamete_list_type_allocator,
	   //template<typename,typename> class gamete_list_type,
	   typename rec_pos_generator>
  unsigned operator()( gamete_iterator_type & g1,
		       gamete_iterator_type & g2,
		       //gamete_list_type< typename gamete_iterator_type::value_type, gamete_list_type_allocator > * gametes,
		       gvector_t * gametes,
		       const double & littler,
		       gsl_rng * r,
		       const rec_pos_generator & rp) const
  {
    unsigned NREC = 0;
    if( g1 != g2 )
      //then a non-parental type is inherited from p1 and p1 has two different gametes
      {
	NREC += recombine_gametes(r,littler,gametes,g1,g2,rp);
      }
    return NREC;
  }	   
};
//Fitness

struct no_selection
{
  /*!
    \brief Method for standard diploid simulations of a single locus.
  */
  using result_type = double;
  template<typename iterator_type >
  inline result_type operator()(const iterator_type &, const iterator_type &) const
  {
    /*
      static_assert( std::is_base_of<mutation_base,
      typename iterator_type::value_type::mutation_type>::value,
      "iterator_type::value_type::mutation_type must be derived from KTfwd::mutation_base" );
    */
    return 1.;
  }
  //! \brief Naive implementation for non-standard cases
  template<typename T >
  inline result_type operator()(const T &) const
  {
    return 1.;
  }
};

struct site_dependent_fitness
{
  using result_type = double;
  ///\example diploid_fixed_sh_ind.cc
  template<typename iterator_type,
	   typename fitness_updating_policy_hom,
	   typename fitness_updating_policy_het>
  inline result_type operator()( const iterator_type & g1,
				 const iterator_type & g2,
				 const fitness_updating_policy_hom & fpol_hom,
				 const fitness_updating_policy_het & fpol_het,
				 const double & starting_fitness  = 1. ) const
  {
    /*      static_assert( std::is_base_of<mutation_base,
	    typename iterator_type::value_type::mutation_type>::value,
	    "iterator_type::value_type::mutation_type must be derived from KTfwd::mutation_base" );
    */
    result_type fitness=starting_fitness;
    if( g1->smutations.empty() && g2->smutations.empty() ) return fitness;
    if( !g1->smutations.empty() && g2->smutations.empty() ) 
      {
	std::for_each( g1->smutations.begin(),
		       g1->smutations.end(),
		       std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
	return fitness;
      }
    if( g1->smutations.empty() && !g2->smutations.empty() ) 
      {
	std::for_each( g2->smutations.begin(),
		       g2->smutations.end(),
		       std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
	return fitness;
      }
    typename iterator_type::value_type::mutation_list_type_iterator ib1,ib2;
    typename iterator_type::value_type::mutation_container::const_iterator b1=g1->smutations.cbegin(),
      e1=g1->smutations.cend(),
      b2=g2->smutations.cbegin(),
      e2=g2->smutations.cend();
    //This is a fast way to calculate fitnesses,
    //as it just compares addresses in memory, and 
    //does little in the way of dereferencing and storing
    bool found = false;
    for( ; b1 < e1 ; ++b1 )
      {
	found = false;
	ib1 = *b1;
	for( ; !found && b2 < e2 && !((*b2)->pos > (ib1)->pos) ; ++b2 )
	  {
	    ib2 = *b2;
	    if ( ib2 == ib1 ) //homozygote
	      {
		assert(ib1->s == ib2->s); //just making sure
		assert(ib1->pos == ib2->pos);
		fpol_hom(fitness,ib1);
		found=true;
	      }
	    else 
	      //b2 points to a unique mutation that comes before b1
	      {
		assert(ib2->pos != ib1->pos);
		assert(ib2->pos < ib1->pos);
		fpol_het(fitness,ib2);
	      }
	  }
	if(!found) //het
	  {
	    fpol_het(fitness,ib1);
	  }
      }
    std::for_each( b2,e2,
		   std::bind(fpol_het,std::ref(fitness),std::placeholders::_1) );
    return fitness;
  }
};

namespace fwdpp_internal {
  /*!
    Custom deleter for std::uniq_ptr
  */
  struct gsl_ran_discrete_t_deleter
  {
    void operator()( gsl_ran_discrete_t * l ) noexcept 
    {
      gsl_ran_discrete_free(l);
    }
  };
    
  /*!
    \warning Can only be put into a vector by push_back/emplace back
    because of constraints on std::unique_ptr assignment
  */
  using gsl_ran_discrete_t_ptr = std::unique_ptr< gsl_ran_discrete_t, 
						  gsl_ran_discrete_t_deleter >;
}

//now, sample_diploid, single-deme, N changing.
template< //typename gamete_type,
	  //typename gamete_list_type_allocator,
	  //typename mutation_list_type_allocator,
	  //typename diploid_geno_t,
	  //typename diploid_vector_type_allocator,
	  typename diploid_fitness_function,
	  typename mutation_removal_policy,
	  typename mutation_model,
	  typename recombination_policy,
	  typename mutation_insertion_policy,
	  typename gamete_insertion_policy,
	  typename gvector_t,
	  typename mvector_t,
	  typename dipvector_t>
	  //template<typename,typename> class gamete_list_type,
	  //template<typename,typename> class mutation_list_type,
	  //template<typename,typename> class diploid_vector_type>
double
sample_diploid(gsl_rng * r,
	       gvector_t * gametes,
	       dipvector_t * diploids,
	       mvector_t * mutations,
	       //gamete_list_type<gamete_type,gamete_list_type_allocator > * gametes,
	       //diploid_vector_type<diploid_geno_t,diploid_vector_type_allocator> * diploids,
	       //mutation_list_type<typename gamete_type::mutation_type,mutation_list_type_allocator > * mutations, 
	       const unsigned & N_curr, 
	       const unsigned & N_next, 
	       const double & mu,
	       const mutation_model & mmodel,
	       const recombination_policy & rec_pol,
	       const mutation_insertion_policy & mpolicy,
	       const gamete_insertion_policy & gpolicy_mut,
	       const diploid_fitness_function & ff,
	       const mutation_removal_policy & mp,
	       const double & f)
{
  assert(N_curr == diploids->size());

  //std::for_each( mutations->begin(),mutations->end(),[](typename gamete_type::mutation_type & __m){__m.n=0;});

  std::vector<double> fitnesses(diploids->size());
  double wbar = 0.;
    
  auto dptr = diploids->begin();
  for( unsigned i = 0 ; i < N_curr ; ++i )
    {
      //(dptr+i)->first->n = 0;
      //(dptr+i)->second->n = 0;
      //fitnesses[i] = fwdpp_internal::diploid_fitness_dispatch(ff,(dptr+i),
      //								typename traits::is_custom_diploid_t<diploid_geno_t>::type());
      fitnesses[i] = ff( (dptr+i)->first,(dptr+i)->second );
      wbar += fitnesses[i];
    }
  wbar /= double(diploids->size());

  /*
#ifndef NDEBUG
  std::for_each(gametes->cbegin(),gametes->cend(),[](decltype((*gametes->cbegin())) __g) {
      assert( !__g.n ); } );
#endif
*/
  fwdpp_internal::gsl_ran_discrete_t_ptr lookup(gsl_ran_discrete_preproc(fitnesses.size(),&fitnesses[0]));
  auto parents(std::move(*diploids)); //copy the parents
  auto pptr = parents.begin();
    
  //Change the population size
  if( diploids->size() != N_next)
    {
      diploids->resize(N_next);
      dptr = diploids->begin();
    }
  unsigned NREC=0;
  assert(diploids->size()==N_next);
  //decltype( gametes->begin() ) p1g1,p1g2,p2g1,p2g2;
  decltype( diploids->begin()->first ) p1g1,p1g2,p2g1,p2g2;
  for( unsigned i = 0 ; i < N_next ; ++i )
    {
      assert(dptr==diploids->begin());
      assert( (dptr+i) < diploids->end() );
      size_t p1 = gsl_ran_discrete(r,lookup.get());
#ifdef FWDPP_COMPAT_0_3_0
      size_t p2 = (gsl_rng_uniform(r) <= f) ? p1 : gsl_ran_discrete(r,lookup.get());
#else
      size_t p2 = (f==1. || (f>0. && gsl_rng_uniform(r)<=f)) ? p1 : gsl_ran_discrete(r,lookup.get());
#endif
      assert(p1<parents.size());
      assert(p2<parents.size());
	
      p1g1 = (pptr+typename decltype(pptr)::difference_type(p1))->first;
      p1g2 = (pptr+typename decltype(pptr)::difference_type(p1))->second;
      p2g1 = (pptr+typename decltype(pptr)::difference_type(p2))->first;
      p2g2 = (pptr+typename decltype(pptr)::difference_type(p2))->second;
	
      NREC += rec_pol(p1g1,p1g2);
      NREC += rec_pol(p2g1,p2g2);
	
      (dptr+i)->first = (gsl_rng_uniform(r) <= 0.5) ? p1g1 : p1g2;
      (dptr+i)->second = (gsl_rng_uniform(r) <= 0.5) ? p2g1 : p2g2;
	
      //(dptr+i)->first->n++;
      //assert( (dptr+i)->first->n > 0 );
      //assert( (dptr+i)->first->n <= 2*N_next );
      //(dptr+i)->second->n++;
      //assert( (dptr+i)->second->n > 0 );
      //assert( (dptr+i)->second->n <= 2*N_next );

      //now, add new mutations
      (dptr+i)->first = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->first,mmodel,mpolicy,gpolicy_mut);
      (dptr+i)->second = mutate_gamete(r,mu,gametes,mutations,(dptr+i)->second,mmodel,mpolicy,gpolicy_mut);
    }
// #ifndef NDEBUG
//   for( unsigned i = 0 ; i < diploids->size() ; ++i )
//     {
//       assert( (dptr+i)->first->n > 0 );
//       assert( (dptr+i)->first->n <= 2*N_next );
//       assert( (dptr+i)->second->n > 0 );
//       assert( (dptr+i)->second->n <= 2*N_next );
//     }
// #endif

  parents.clear();

  gametes->erase( std::remove_if(gametes->begin(),
				 gametes->end(),
				 []( typename gvector_t::value_type & __g ) { return __g.use_count() == 1; }),
		  gametes->end() );
  //uh-oh
  std::for_each(mutations->begin(),mutations->end(),
		[&diploids](mptr_t & __m) { if(unsigned(__m.use_count()) == 2u*diploids->size() + 1u) __m->fixed = true; });
  std::for_each( gametes->begin(),
		 gametes->end(),
		 [&mp]( decltype( *(gametes->begin()) ) & __g ) {
		   __g->mutations.erase( std::remove_if(__g->mutations.begin(),__g->mutations.end(),std::cref(mp)), __g->mutations.end() );
		   __g->smutations.erase( std::remove_if(__g->smutations.begin(),__g->smutations.end(),std::cref(mp)), __g->smutations.end() );
		 });
  // decltype(gametes->begin()) temp;
  // for( auto itr = gametes->begin() ; itr != gametes->end() ;  )
  //   {
  //     if( itr->n == 0 ) //this gamete is extinct and need erasing from the list
  // 	{
  // 	  temp = itr;
  // 	  ++itr;
  // 	  gametes->erase(temp);
  // 	}
  //     else //gamete remains extant and we adjust mut counts
  // 	{
  // 	  adjust_mutation_counts(itr,itr->n);
  // 	  ++itr;
  // 	}
  //   }
  // std::for_each( gametes->begin(),
  // 		 gametes->end(),
  // 		 [&mp](decltype( *(gametes->begin()) ) & __g ) {
  // 		   __g.mutations.erase( std::remove_if(__g.mutations.begin(),__g.mutations.end(),std::cref(mp)),__g.mutations.end() );
  // 		   __g.smutations.erase( std::remove_if(__g.smutations.begin(),__g.smutations.end(),std::cref(mp)),__g.smutations.end() );
  // 		 });
  // assert(check_sum(gametes,2*N_next));
  return wbar;
}

struct equal_eps
  {
    /*! \brief Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
      Returns true if std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon()
    */
    using result_type = bool;
    template<typename T>
    inline bool operator()(const T & lhs, const T & rhs) const
    {
      return( std::max(lhs,rhs)-std::min(lhs,rhs) <= std::numeric_limits<T>::epsilon() );
    }
};

//We now have new sets of types

using mvector = vector< mptr_t >;
using gvector = vector< gptr_t >;
using diploid_t = std::pair< gptr_t, gptr_t >;
using dvector = vector<diploid_t>;
using lookup_t = std::unordered_set<double,std::hash<double>,equal_eps>;

//we now need a mutation policy
struct mpolicy
{
  using result_type = mutation;
  inline result_type operator()( gsl_rng * r, lookup_t * lookup, const unsigned & generation ) const
  {
    double pos = gsl_rng_uniform(r);
    while( lookup->find(pos) != lookup->end() )
      {
	pos = gsl_rng_uniform(r);
      }
    return result_type(pos,0.,0.,generation);
  }
};


int main(int argc, char ** argv)
{
  int argn=1;
  unsigned seed = atoi(argv[argn++]);
  unsigned simlen = atoi(argv[argn]++);

  gsl_rng * r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r,seed);

  const unsigned N = 1000;
  //My latest favorit test params
  const double theta = 5000, rho = 5000;
  const double mu = theta/double(4*N),littler=rho/double(4*N);

  mvector mutations; //monomorphic pop
  gvector gametes(1,make_shared<gamete>(gamete()));   //monomorphic pop
  dvector diploids(N, make_pair(*gametes.begin(),*gametes.begin()));
  lookup_t lookup;
  
  //Unlike list-based fwdpp, we can finally reserve!
  mutations.reserve(100000);
  gametes.reserve(2*N);

  const std::function<double(void)> recmap = std::bind(gsl_rng_uniform,r);
  //Sanity check (looks good)
  /*
  for_each( gametes.cbegin(), gametes.cend(),
	    [](const gptr_t & __g ) { cerr << __g.use_count() << '\n'; });

  for_each(diploids.begin(),diploids.end(),
	   [](const diploid_t & __d) {
	     cerr << __d.first->smutations.size() << ' ' << __d.first->smutations.size() << ' '
		  << __d.second->smutations.size() << ' ' << __d.second->smutations.size() << ' '
		  << __d.first.use_count() << ' ' << __d.second.use_count() << '\n';
	   });
  */
  //Here is the minor down-side: the shared_ptr is 2x the size of an interator, presumable b/c it stores other info.
  //cerr << sizeof(mptr_t) << ' ' << sizeof(gptr_t) << ' ' << sizeof( vector<mutation>::iterator ) << ' ' << sizeof(list<mutation>::iterator) << '\n';

  //unsigned simlen = 10*N;

  // //OK, let's make sure that things can even compile...

  // //Mutation
  // /*
  //   iterator_type mutate_gamete( gsl_rng * r,
  //   const double & mu, list_type< typename iterator_type::value_type,list_type_allocator > * gametes,
  //   list_type2<typename iterator_type::value_type::mutation_type,list_type_allocator2 > * mutations, 
  //   iterator_type & g,
  //   const mutation_model &mmodel,
  //   const mutation_insertion_policy & mpolicy,
  //   const gamete_insertion_policy & gpolicy)
  // */
  // unsigned gen = 0;
  // diploids[0].first = mutate_gamete(r,
  // 				    //Mean of 3 mutations
  // 				    3.,
  // 				    &gametes,
  // 				    &mutations,
  // 				    diploids[0].first,
  // 				    std::bind(mpolicy(),r,&lookup,gen),
  // 				    []( mptr_t && m, mvector * mv) { return mv->insert(mv->end(),forward<mptr_t>(m)); },
  // 				    []( gptr_t && g, gvector * gv) { return *(gv->insert(gv->end(),forward<gptr_t>(g))); }
  // 				    );

  // if(! mutations.empty() )
  //   {
  //     for( const auto & m : mutations )
  //     //for(auto m = mutations.begin() ; m != mutations.end() ;++m )
  // 	cout << m.use_count() << ' ' << m->pos << '\n';
  //   }
  // for(auto g = gametes.begin() ; g != gametes.end() ; ++ g)
  //   {
  //     cout << g->use_count() << ' ' << (*g)->mutations.size() << ' ' << (*g)->smutations.size() << '\n';
  //     for( const auto & __m : (*g)->mutations )
  // 	{
  // 	  cout << '\t' << ' ' << __m.use_count() << ' ' << __m->pos << '\n';
  // 	}
  //   }

  // cout << diploids[0].first->mutations.size() << '\n';

  // //OK, bitchin'.  Now, let's make that individual a homozygote:
  // diploids[0].second.reset();
  // diploids[0].second = diploids[0].first;

  // //Rad, that works
  // for(auto g = gametes.begin() ; g != gametes.end() ; ++ g)
  //   {
  //     cout << g->use_count() << ' ' << (*g)->mutations.size() << ' ' << (*g)->smutations.size() << '\n';
  //     for( const auto & __m : (*g)->mutations )
  // 	{
  // 	  cout << '\t' << ' ' << __m.use_count() << ' ' << __m->pos << '\n';
  // 	}
  //   }

  // //Now, lets go for crossing-over

  // //Reset that thing
  // diploids[0].second->mutations.clear();
  // //and mutate it
  // diploids[0].second = mutate_gamete(r,
  // 				    //Mean of 3 mutations
  // 				    3.,
  // 				    &gametes,
  // 				    &mutations,
  // 				    diploids[0].second,
  // 				    std::bind(mpolicy(),r,&lookup,gen),
  // 				    []( mptr_t && m, mvector * mv) { return mv->insert(mv->end(),forward<mptr_t>(m)); },
  // 				    []( gptr_t && g, gvector * gv) { return *(gv->insert(gv->end(),forward<gptr_t>(g))); }
  // 				    );

  // for(auto g = gametes.begin() ; g != gametes.end() ; ++ g)
  //   {
  //     cout << g->use_count() << ' ' << (*g)->mutations.size() << ' ' << (*g)->smutations.size() << '\n';
  //     for( const auto & __m : (*g)->mutations )
  // 	{
  // 	  cout << '\t' << ' ' << __m.use_count() << ' ' << __m->pos << '\n';
  // 	}
  //   }

  // //Cross those bad boyz over

  // //Print
  // cout << "b4:\n";
  // for( const auto & m : diploids[0].first->mutations )
  //   {
  //     cout << m.use_count() << ' ' << m->pos << '\n';
  //   }
  // cout << "//\n";
  // for( const auto & m : diploids[0].second->mutations )
  //   {
  //     cout << m.use_count() << ' ' << m->pos << '\n';
  //   }

  // unsigned NREC = genetics101()(diploids[0].first,diploids[0].second,&gametes,littler,r,recmap);

  // cout << "after:\n";
  // for( const auto & m : diploids[0].first->mutations )
  //   {
  //     cout << m.use_count() << ' ' << m->pos << '\n';
  //   }
  // cout << "//\n";
  // for( const auto & m : diploids[0].second->mutations )
  //   {
  //     cout << m.use_count() << ' ' << m->pos << '\n';
  //   }

  //cout << diploids[0].second->mutations.size() << '\n';
  //[]( mutation && __m, mvector * mv ) { mv->emplace_back(mptr_t(new mutation(std::forward<mutation>(__m)))); return mv->end()-1; },
  //[]( gamete && __g, gvector * gv ) { gv->emplace_back(mptr_t(new gamete(std::forward<gamete>(__g)))); return gv->end()-1; } );

  for( unsigned gen = 0 ; gen < simlen ; ++gen )
    {
      //if(gen % 100 == 0.) cout << gen << endl;
      double wbar = sample_diploid(r,
  				   &gametes,
  				   &diploids,
  				   &mutations,
  				   //constant size...
  				   N,
  				   N,   
  				   //Mutation
  				   mu,
  				   std::bind(mpolicy(),r,&lookup,gen),
  				   //crossing-over
  				   std::bind(genetics101(),std::placeholders::_1,std::placeholders::_2,
  					     &gametes,
  					     littler,
  					     r,
  					     recmap),
  				   //insertion/update policies--these gotta be custom, as the existing ones are NOT applicab
				   []( mptr_t && m, mvector * mv) { return mv->insert(mv->end(),forward<mptr_t>(m)); },
				   []( gptr_t && g, gvector * gv) { return *(gv->insert(gv->end(),forward<gptr_t>(g))); },
  				   //Fitness
  				   std::bind(no_selection(),std::placeholders::_1,std::placeholders::_2),
  				   //mut removal policy
  				   std::bind(mutation_remover(),std::placeholders::_1,1,2*N+1),
  				   //[&N](mptr_t & __g){ return __g.use_count() == 1 || __g.use_count() == (2*N+1); },
  				   0.
  				   );
      mutations.erase( remove_if(mutations.begin(),mutations.end(),
				 [](mptr_t & m) { return m.use_count() == 1 || m.use_count() == 2*N+1; }),
		       mutations.end() );
    }
  std::cerr << mutations.size() << ' ' << gametes.size() << ' ' << diploids.size() << '\n';
}
