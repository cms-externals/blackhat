/*
 * cached_EP.h
 *
 *  Created on: FEB 9, 2011
 *      Author: harald
 */

#ifndef CACHED_EP_H_
#define CACHED_EP_H_

#include <map>
#include <vector>
#include "settings.h"
#include "BH_typedefs.h"
#include "mom_conf.h"
#include "eval_param.h"
#ifndef BH_PUBLIC
#include "BH_grass.h"
#endif
#include "mode_dependent_typedefs.h"

#if BH_USE_OMP
#include <omp.h>
#include "BH_omp.h"
#endif

namespace BH {

//class TreeHelAmpl;
//template <class T> class SeriesC;
template <class T> class momentum_configuration;
template <class T> class eval_param;
#ifndef BH_PUBLIC
template <class T> class BH_grass;
#endif
namespace CachedTHA {
	
class Cached_EP 
#if BH_USE_OMP
: private Tools::HasOMPLock
#endif
{
	std::vector<std::vector<int> > d_index_vectors;
	std::vector<eval_param<R>* > d_values;
	std::vector<eval_param<RHP>* > d_values_HP;
	std::vector<eval_param<RVHP>* > d_values_VHP;
#if BH_USE_OMP
	Tools::FSArray<volatile long,20> d_mcIDs;
	Tools::FSArray<volatile long,20> d_mcIDs_HP;
	Tools::FSArray<volatile long,20> d_mcIDs_VHP;
#else
	std::vector<long > d_mcIDs;
	std::vector<long > d_mcIDs_HP;
	std::vector<long > d_mcIDs_VHP;
#endif
	//! keeps track of the number of Cached_EP_user using each entry
	std::vector<int> d_load;
    //template <class T> eval_param<T>* eval_fn(size_t n, momentum_configuration<T>& mc);
    /**/
    // for grassman tree cache 
	std::vector<std::vector<int> > d_index_vectors_gr;
#ifndef BH_PUBLIC
	std::vector<BH_grass<R>* > d_values_gr;
	std::vector<BH_grass<RHP>* > d_values_gr_HP;
	std::vector<BH_grass<RVHP>* > d_values_gr_VHP;
#if BH_USE_OMP
	Tools::FSArray<volatile long,20> d_mcIDs_gr;
	Tools::FSArray<volatile long,20> d_mcIDs_gr_HP;
	Tools::FSArray<volatile long,20> d_mcIDs_gr_VHP;



#else
	std::vector<long > d_mcIDs_gr;
	std::vector<long > d_mcIDs_gr_HP;
	std::vector<long > d_mcIDs_gr_VHP;
#endif
//! keeps track of the number of Cached_EP_user using each entry
	std::vector<int> d_load_gr;
    //template <class T> eval_param<T>* eval_fn(size_t n, momentum_configuration<T>& mc);
#endif
public:
	//! Constructor
	Cached_EP();
	//! Return the eval_param with the n th index vector in the list
	/** n is zero based */
	//! Evaluate the tree with the n th index vector in the list
	/** n is zero based */
	size_t add(const std::vector<int>& indices);
	size_t add_gr(const std::vector<int>& indices);
	eval_param<R>* eval(size_t n, momentum_configuration<R>& mc);
	eval_param<RHP>* eval(size_t n, momentum_configuration<RHP>& mc);
	eval_param<RVHP>* eval(size_t n, momentum_configuration<RVHP>& mc);
    //
#ifndef BH_PUBLIC

	BH_grass<R>* eval_gr(size_t n, momentum_configuration<R>& mc);
	BH_grass<RHP>* eval_gr(size_t n, momentum_configuration<RHP>& mc);
	BH_grass<RVHP>* eval_gr(size_t n, momentum_configuration<RVHP>& mc);
#endif
	//
    //! prints some statistics
	void print_stat();
	const std::vector<int>& get_index_vector(int n) const {return d_index_vectors[n];} ;
	const std::vector<int>& get_index_vector_gr(int n) const {return d_index_vectors[n];} ;
	//void dry_run();
	void refresh(momentum_configuration<R>& mc);
	template <class T> inline void set_value(int n, eval_param<T>* );
	template <class T> inline eval_param<T>* get_value(int n);
#ifndef BH_PUBLIC
	template <class T> inline void set_value_gr(int n, BH_grass<T>* );
	template <class T> inline BH_grass<T>* get_value_gr(int n);
	template <class T> inline void set_mcID_gr(int n,long id);
	template <class T> inline long get_mcID_gr(int n);
#endif
	template <class T> inline void set_mcID(int n,long id);
	template <class T> inline long get_mcID(int n);
	//! destructor
	virtual ~Cached_EP();
};


template <> inline void Cached_EP::set_value(int n, eval_param<R>* val){ if(d_values[n]!=0) {delete d_values[n];};  d_values[n]=val;}
template <> inline void Cached_EP::set_value(int n, eval_param<RHP>* val){ if(d_values_HP[n]!=0) {delete d_values_HP[n];};  d_values_HP[n]=val;}
template <> inline void Cached_EP::set_value(int n, eval_param<RVHP>* val){ if(d_values_VHP[n]!=0) {delete d_values_VHP[n];};  d_values_VHP[n]=val;}

template <> inline eval_param<R>* Cached_EP::get_value(int n){ return d_values[n];}
template <> inline eval_param<RHP>* Cached_EP::get_value(int n){ return d_values_HP[n];}
template <> inline eval_param<RVHP>* Cached_EP::get_value(int n){ return d_values_VHP[n];}

#if BH_USE_OMP
template <> inline void Cached_EP::set_mcID<R>(int n,long id){ const_cast<volatile long&>(d_mcIDs[n])=id;}
template <> inline void Cached_EP::set_mcID<RHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_HP[n])=id;}
template <> inline void Cached_EP::set_mcID<RVHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_VHP[n])=id;}
#else 
template <> inline void Cached_EP::set_mcID<R>(int n,long id){ d_mcIDs[n]=id;}
template <> inline void Cached_EP::set_mcID<RHP>(int n,long id){ d_mcIDs_HP[n]=id;}
template <> inline void Cached_EP::set_mcID<RVHP>(int n,long id){ d_mcIDs_VHP[n]=id;}
#endif

template <> inline long Cached_EP::get_mcID<R>(int n){ return d_mcIDs[n];}
template <> inline long Cached_EP::get_mcID<RHP>(int n){ return d_mcIDs_HP[n];}
template <> inline long Cached_EP::get_mcID<RVHP>(int n){ return d_mcIDs_VHP[n];}

/**/
#ifndef BH_PUBLIC
template <> inline void Cached_EP::set_value_gr(int n, BH_grass<R>* val){ if(d_values_gr[n]!=0) {delete d_values_gr[n];};  d_values_gr[n]=val;}
template <> inline void Cached_EP::set_value_gr(int n, BH_grass<RHP>* val){ if(d_values_gr_HP[n]!=0) {delete d_values_gr_HP[n];};  d_values_gr_HP[n]=val;}
template <> inline void Cached_EP::set_value_gr(int n, BH_grass<RVHP>* val){ if(d_values_gr_VHP[n]!=0) {delete d_values_gr_VHP[n];};  d_values_gr_VHP[n]=val;}

template <> inline BH_grass<R>* Cached_EP::get_value_gr(int n){ return d_values_gr[n];}
template <> inline BH_grass<RHP>* Cached_EP::get_value_gr(int n){ return d_values_gr_HP[n];}
template <> inline BH_grass<RVHP>* Cached_EP::get_value_gr(int n){ return d_values_gr_VHP[n];}

#if BH_USE_OMP
template <> inline void Cached_EP::set_mcID_gr<R>(int n,long id){ const_cast<volatile long&>(d_mcIDs_gr[n])=id;}
template <> inline void Cached_EP::set_mcID_gr<RHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_gr_HP[n])=id;}
template <> inline void Cached_EP::set_mcID_gr<RVHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_gr_VHP[n])=id;}
#else 
template <> inline void Cached_EP::set_mcID_gr<R>(int n,long id){ d_mcIDs_gr[n]=id;}
template <> inline void Cached_EP::set_mcID_gr<RHP>(int n,long id){ d_mcIDs_gr_HP[n]=id;}
template <> inline void Cached_EP::set_mcID_gr<RVHP>(int n,long id){ d_mcIDs_gr_VHP[n]=id;}
#endif
template <> inline long Cached_EP::get_mcID_gr<R>(int n){ return d_mcIDs_gr[n];}
template <> inline long Cached_EP::get_mcID_gr<RHP>(int n){ return d_mcIDs_gr_HP[n];}
template <> inline long Cached_EP::get_mcID_gr<RVHP>(int n){ return d_mcIDs_gr_VHP[n];}

#endif /*BH_PUBLIC*/
class Cached_EP_factory {
    std::map<size_t,Cached_EP*> d_Cached_EP;
public :
	Cached_EP_factory() {};
	// integer n_mom_conf is meant to be useful for dipoles when many distinct 
    // momentum configurations are floating around at the same time.
    Cached_EP*  new_CEP(size_t n_mom_conf);
	void print_state();
	void refresh(momentum_configuration<R>& mc);
	static Cached_EP_factory* default_CEP;
	~Cached_EP_factory();

};


}
}
#endif /* CACHED_EP_H_ */
