/*
 * cached_THA.h
 *
 *  Created on: Sep 30, 2009
 *      Author: harald
 */

#ifndef CACHED_THA_H_
#define CACHED_THA_H_

#include <map>
#include <vector>
#include "settings.h"
#include "BH_typedefs.h"
#include "process.h"
#include "mom_conf.h"
#include "mode_dependent_typedefs.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "amplitudes.h"
#include "OneLoopHelAmpl.h"
#include "cached_EP.h"
#ifndef BH_PUBLIC
#include "grass_tree.h"
#endif

#if BH_USE_OMP
#include <omp.h>
#include "BH_omp.h"
#endif

namespace BH {

class TreeHelAmpl;
//template <class T> class SeriesC;
template <class T> class momentum_configuration;
class process;


namespace CachedTHA {
class Cached_EP;
	
//! THA objects that cache the eval results for more than one index vector
class Cached_THA 
#if BH_USE_OMP
: private Tools::HasOMPLock
#endif
{
//	TreeHelAmpl* d_THA_p;
	std::vector<std::complex<R> > d_values;
	std::vector<std::complex<RHP> > d_values_HP;
	std::vector<std::complex<RVHP> > d_values_VHP;
#if BH_USE_OMP
	Tools::FSArray<volatile long,20> d_mcIDs;
	Tools::FSArray<volatile long,20> d_mcIDs_HP;
	Tools::FSArray<volatile long,20> d_mcIDs_VHP;
#else
	std::vector<long > d_mcIDs;
	std::vector<long > d_mcIDs_HP;
	std::vector<long > d_mcIDs_VHP;
#endif
public:
	std::vector<std::vector<int> > d_index_vectors;
	std::vector<int> d_load;
	//Cached_THA(TreeHelAmpl* THAb);
	Cached_THA(){};
    //template <class T> std::complex<T> eval(int n, momentum_configuration<T>& mc);
	//template <class T> std::complex<T> eval(int n, momentum_configuration<T>& mc, Cached_EP* CEP, size_t ep_ind);
    virtual std::complex<R> eval(int n, momentum_configuration<R>& mc, Cached_EP* CEP, size_t ep_ind){return std::complex<R>(0,0);};
    virtual std::complex<RHP> eval(int n, momentum_configuration<RHP>& mc,Cached_EP* CEP, size_t ep_ind){return std::complex<RHP>(0,0);};
    virtual std::complex<RVHP> eval(int n, momentum_configuration<RVHP>& mc,Cached_EP* CEP, size_t ep_ind){return std::complex<RVHP>(0,0);};
	//! Evaluate the tree with the n th index vector in the list
	/** n is zero based */
	size_t add(const std::vector<int>& indices);
	//! destructor
	//! prints some statistics
	virtual void print_stat(){};
	//const process& get_process() const ;
	const std::vector<int>& get_index_vector(int n) const {return d_index_vectors[n];} ;
	//void dry_run();
	//void refresh(momentum_configuration<R>& mc);
	template <class T> inline void set_value(int n,const std::complex<T>& val);
	template <class T> inline std::complex<T> get_value(int n);
	template <class T> inline void set_mcID(int n,long id);
	template <class T> inline long get_mcID(int n);
	virtual ~Cached_THA(){};
};

template <> inline void Cached_THA::set_value(int n,const std::complex<R>& val){d_values[n]=val;}
template <> inline void Cached_THA::set_value(int n,const std::complex<RHP>& val){d_values[n]=to_double(val); d_values_HP[n]=val;}
template <> inline void Cached_THA::set_value(int n,const std::complex<RVHP>& val){d_values[n]=to_double(val); d_values_HP[n]=to_HP(val); d_values_VHP[n]=val;}

template <> inline std::complex<R> Cached_THA::get_value(int n){ return d_values[n];}
template <> inline std::complex<RHP> Cached_THA::get_value(int n){ return d_values_HP[n];}
template <> inline std::complex<RVHP> Cached_THA::get_value(int n){ return d_values_VHP[n];}

#if BH_USE_OMP
template <> inline void Cached_THA::set_mcID<R>(int n,long id){ const_cast<volatile long&>(d_mcIDs[n])=id;}
template <> inline void Cached_THA::set_mcID<RHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_HP[n])=id;}
template <> inline void Cached_THA::set_mcID<RVHP>(int n,long id){ const_cast<volatile long&>(d_mcIDs_VHP[n])=id;}
#else 
template <> inline void Cached_THA::set_mcID<R>(int n,long id){ d_mcIDs[n]=id;}
template <> inline void Cached_THA::set_mcID<RHP>(int n,long id){ d_mcIDs_HP[n]=id;}
template <> inline void Cached_THA::set_mcID<RVHP>(int n,long id){ d_mcIDs_VHP[n]=id;}
#endif
template <> inline long Cached_THA::get_mcID<R>(int n){ return d_mcIDs[n];}
template <> inline long Cached_THA::get_mcID<RHP>(int n){ return d_mcIDs_HP[n];}
template <> inline long Cached_THA::get_mcID<RVHP>(int n){ return d_mcIDs_VHP[n];}


	
class Cached_THA_ep: public Cached_THA
{
	TreeHelAmpl* d_THA_p;
public:
	Cached_THA_ep(TreeHelAmpl* THAb): d_THA_p(THAb) {};
    std::complex<R> eval(int n, momentum_configuration<R>& mc,Cached_EP* CEP, size_t ep_ind){return eval_fn(n,mc,CEP,ep_ind);}
    std::complex<RHP> eval(int n, momentum_configuration<RHP>& mc,Cached_EP* CEP, size_t ep_ind){return eval_fn(n,mc,CEP,ep_ind);}
    std::complex<RVHP> eval(int n, momentum_configuration<RVHP>& mc,Cached_EP* CEP, size_t ep_ind){return eval_fn(n,mc,CEP,ep_ind);}
	void print_stat();
	const process& get_process() const ;
	virtual ~Cached_THA_ep();
private:
    template <class T> std::complex<T> eval_fn(int n, momentum_configuration<T>& mc, Cached_EP* CEP, size_t ep_ind);

};

#ifndef BH_PUBLIC
class Cached_THA_grass: public Cached_THA
{
	//TreeHelAmpl* d_THA_p;
    Tree_Fn_Ptr_eval_grass d_tree_eval;
    Tree_Fn_Ptr_eval_grass_HP d_tree_eval_HP; 
    Tree_Fn_Ptr_eval_grass_VHP d_tree_eval_VHP;

public:
	Cached_THA_grass(Tree_Fn_Ptr_eval_grass tree_eval, Tree_Fn_Ptr_eval_grass_HP tree_eval_HP, Tree_Fn_Ptr_eval_grass_VHP tree_eval_VHP ): d_tree_eval(tree_eval), d_tree_eval_HP(tree_eval_HP), d_tree_eval_VHP(tree_eval_VHP)   {};
    std::complex<R> eval(int n, momentum_configuration<R>& mc,Cached_EP* CEP, size_t ep_ind);//{return eval_fn(n,mc,CEP,ep_ind);}
    std::complex<RHP> eval(int n, momentum_configuration<RHP>& mc,Cached_EP* CEP, size_t ep_ind);//{return eval_fn(n,mc,CEP,ep_ind);}
    std::complex<RVHP> eval(int n, momentum_configuration<RVHP>& mc,Cached_EP* CEP, size_t ep_ind);//{return eval_fn(n,mc,CEP,ep_ind);}
    virtual ~Cached_THA_grass();
//private:
//    template <class T> std::complex<T> eval_fn(int n, momentum_configuration<T>& mc, Cached_EP* CEP, size_t ep_ind);
};


#endif
class Cached_THA_user  {
protected:
    Cached_THA* d_CTHA_p;
	size_t d_n;
    Cached_EP* d_Cached_EP;
    size_t d_ep_ind;
public:
    //Cached_THA_user(Cached_THA& CTHA,size_t n): d_CTHA_p(&CTHA), d_n(n) {};
	Cached_THA_user(Cached_THA* CTHA,size_t n): d_CTHA_p(CTHA), d_n(n) {};
	Cached_THA_user(Cached_THA* CTHA,size_t n,Cached_EP* CEP,size_t ep_ind): d_CTHA_p(CTHA), d_n(n), d_Cached_EP(CEP), d_ep_ind(ep_ind) {};
    virtual std::complex<R> eval(momentum_configuration<R>& mc){return std::complex<R>(0,0);};
	virtual std::complex<RHP> eval(momentum_configuration<RHP>& mc){return std::complex<RHP>(0,0);};
	virtual std::complex<RVHP> eval(momentum_configuration<RVHP>& mc){return std::complex<RVHP>(0,0);};
    //void dry_run(){ d_CTHA_p->dry_run();}
    //const process& get_process(){ return d_CTHA_p->get_process();};
    //const std::vector<int>& get_index_vector(){ return d_CTHA_p->get_index_vector(d_n); }
    virtual ~Cached_THA_user(){};
};


class Cached_THA_user_normal: public Cached_THA_user  {
public:
	Cached_THA_user_normal(Cached_THA* CTHA,size_t n): Cached_THA_user(CTHA,n) {};
	Cached_THA_user_normal(Cached_THA* CTHA,size_t n,Cached_EP* CEP,size_t ep_ind): Cached_THA_user(CTHA,n,CEP,ep_ind) {};
    std::complex<R> eval(momentum_configuration<R>& mc){return d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind);};
	std::complex<RHP> eval(momentum_configuration<RHP>& mc){return d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind);};
	std::complex<RVHP> eval(momentum_configuration<RVHP>& mc){return d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind);};
    virtual ~Cached_THA_user_normal();
};


class Cached_THA_user_conjugate: public Cached_THA_user  {
protected:
    std::complex<R> d_sign_R;
    std::complex<RHP> d_sign_RHP;
    std::complex<RVHP> d_sign_RVHP;
public:
	Cached_THA_user_conjugate(Cached_THA* CTHA,size_t n,int sign): Cached_THA_user(CTHA,n), d_sign_R(sign,0), d_sign_RHP(sign,0), d_sign_RVHP(sign,0) {};
	Cached_THA_user_conjugate(Cached_THA* CTHA,size_t n, Cached_EP* CEP,size_t ep_ind,  int sign): Cached_THA_user(CTHA,n,CEP,ep_ind), d_sign_R(sign,0), d_sign_RHP(sign,0), d_sign_RVHP(sign,0) {};
    virtual std::complex<R> eval(momentum_configuration<R>& mc){return d_sign_R*conj(d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind));};
	virtual std::complex<RHP> eval(momentum_configuration<RHP>& mc){return d_sign_RHP*conj(d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind));};
	virtual std::complex<RVHP> eval(momentum_configuration<RVHP>& mc){return d_sign_RVHP*conj(d_CTHA_p->eval(d_n,mc,d_Cached_EP,d_ep_ind));};
    virtual ~Cached_THA_user_conjugate();
};




class Cached_THA_factory {
	std::map<const process,Cached_THA*> d_amplitudes;
public :
	Cached_THA_factory() {};
//	Cached_THA_user new_THA(process pro,const std::vector<int>& ind, short conj=0);
	Cached_THA_user* new_THA(process pro,const std::vector<int>& ind, short conj=0);
	void print_state();
	//void refresh(momentum_configuration<R>& mc);
	//void fill_process_list(std::vector<process >&);
	static Cached_THA_factory* default_CTHA;
	~Cached_THA_factory();

};

class Cached_TA;
class Cached_TA_factory;


class Cached_TA_factory {
    std::map<const process, TreeHelAmpl* > d_amplitudes;
    std::map<const std::vector<int>, eval_param<R>* > d_eval_params;
    std::map<eval_param<R>*, size_t> d_ep_n;
    public:
    std::vector<TreeHelAmpl* > d_TreeHelAmpl;
    
    std::vector< Cached_TA* > d_Cached_TA;
    std::vector< eval_param<R>* > d_eval_param;
    std::vector< std::vector<int>* > d_ind;
     

    //size_t add_Cached_TA(){ d_Cached_TA.push_back(new Cached_TA(this));
    size_t add_Cached_TA();
    Cached_TA* get_Cached_TA(size_t n){return d_Cached_TA[n];};

    void delete_Cached_TA(size_t n){ d_Cached_TA.erase(d_Cached_TA.begin()+n); }
    Cached_TA_factory();
    void new_TreeAmpl(const process  pro, const std::vector<int>&  ind, short conj=0); 
	static Cached_TA_factory* default_CTA;
    ~Cached_TA_factory();
};

class Cached_TA {
    private:
    bool d_complete;
    std::map< std::pair<TreeHelAmpl*, size_t >, std::pair<size_t,size_t> > d_map_vals;
    //vector over thas; map< n_ep, pos in d_p_ep >
    //std::vector< std::map< size_t,size_t > > d_map_eps;
    std::map<size_t,size_t > d_known_eps; 
    std::map<TreeHelAmpl*, size_t > d_known_tha; 

    std::vector< TreeHelAmpl* > d_p_tree_hel_ampls;
    std::vector< std::vector< size_t > > d_p_ep;
    
    std::vector< eval_param<R>* > d_eval_param;
    std::vector< eval_param<RHP>* > d_eval_param_HP;
    std::vector< eval_param<RVHP>* > d_eval_param_VHP;
    std::vector< std::vector<int> * > d_ind;
    
    std::vector< TreeHelAmpl* >::const_iterator it_tha;
    std::vector< eval_param<R>* >::const_iterator it_ep;
    std::vector< size_t >::const_iterator it_n_ep;
    std::vector< std::vector<int> * >::const_iterator it_ind;
    
    std::complex<R> ** d_p_vals_c;
    std::complex<R> ** d_p_vals_mc;
   
    std::complex<R> ** a_vals;
    std::complex<R> * a_vals_c;
    std::complex<R> * a_vals_mc;
    size_t d_nbr_c;
    size_t d_nbr_mc;
    //complex<R> ** p_vals;
    std::vector<pair<size_t,int> > d_poll;
    std::vector<pair<size_t,size_t> > d_poll_c;
    std::vector<pair<size_t,size_t> > d_poll_mc;


    public:
    std::complex<R> ** p_vals;
    //Cached_TA(Cached_TA_factory* CTAf): d_Cached_TA_factory(CTAf) {};
    Cached_TA();
    void add(TreeHelAmpl*,size_t,short conj=0);
    size_t add(size_t n_ep);
    void complete_construction();

    //std::vector< std::vector< complex<R> > > d_vals;
    //std::vector< complex<R> > d_vals_c;
    //std::vector< complex<R> > d_p_vals_mc;
    void eval(momentum_configuration<R>& mc);
    void eval(momentum_configuration<RHP>& mc);
    void eval(momentum_configuration<RVHP>& mc);
    static Cached_TA_factory* default_CTA;
    ~Cached_TA();

};




}
}
#endif /* CACHED_THA_H_ */
