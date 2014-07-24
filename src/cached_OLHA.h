/*
 * cached_OLHA.h
 *
 *  Created on: 12-Nov-2008
 *      Author: daniel
 */

#ifndef CACHED_OLHA_H_
#define CACHED_OLHA_H_

#include <map>
#include <vector>
#include "BH_typedefs.h"
#include "process.h"   // needed for member of pro_cs
//#include "assembly.h"   // needed to derive from partial_amplitude
#include "IR_checked.h"

#if BH_USE_GMP

#include "gmp_r.h"
#endif



namespace BH {

class OneLoopAmplitude_base;
template <class T> class SeriesC;
template <class T> class momentum_configuration;
class process;
template <class CUT> class cut_part_factory;
class Cut_Part_base;
class multi_precision_fraction;


namespace CachedOLHA {



//! OLHA objects that cache the eval results for more than one index vector
class Cached_OLHA {
	OneLoopAmplitude_base* d_OLHA_p;
	std::vector<std::vector<int> > d_index_vectors;
	std::vector<SeriesC<R> > d_values;
	std::vector<SeriesC<RHP> > d_values_HP;
	std::vector<SeriesC<RVHP> > d_values_VHP;
	std::vector<std::complex<R> > d_tree_values;
	std::vector<std::complex<RHP> > d_tree_values_HP;
	std::vector<std::complex<RVHP> > d_tree_values_VHP;
    //charge conj	
	std::vector<SeriesC<R> > d_conj_values;
	std::vector<SeriesC<RHP> > d_conj_values_HP;
	std::vector<SeriesC<RVHP> > d_conj_values_VHP;
    
	std::vector<double > d_accuracy;
	std::vector<long > d_mcIDs;
	std::vector<long > d_mcIDs_HP;
	std::vector<long > d_mcIDs_VHP;
	std::vector<long > d_mu_index;
	std::vector<long > d_mu_index_HP;
	std::vector<long > d_mu_index_VHP;
	//! keeps track of the number of Cached_OLHA_user using each entry
	std::vector<int> d_load;

#if BH_USE_GMP
	std::vector<long > d_mcIDs_GMP;
	std::vector<long > d_mu_index_GMP;
	std::vector<std::complex<RGMP> > d_tree_values_GMP;
	std::vector<SeriesC<RGMP> > d_values_GMP;
	std::vector<SeriesC<RGMP> > d_conj_values_GMP;
#endif

public:
	//! Constructor
//	Cached_OLHA(const process&,color_structure cs);
	//! Constructor
	Cached_OLHA(OneLoopAmplitude_base* OLHAb);
	//! Evaluate the OLHA with the n th index vector in the list
	/** n is zero based */
	//template <class T> SeriesC<T> eval(size_t n, momentum_configuration<T>& mc,int mu_index);
	SeriesC<R> eval(size_t n, momentum_configuration<R>& mc,int mu_index);
	SeriesC<RHP> eval(size_t n, momentum_configuration<RHP>& mc,int mu_index);
	SeriesC<RVHP> eval(size_t n, momentum_configuration<RVHP>& mc,int mu_index);
#ifdef BH_USE_GMP
	SeriesC<RGMP> eval(size_t n, momentum_configuration<RGMP>& mc,int mu_index);
#endif
	//template <class T> SeriesC<T> eval_conj(size_t n, momentum_configuration<T>& mc,int mu_index);
	SeriesC<R> eval_conj(size_t n, momentum_configuration<R>& mc,int mu_index);
	SeriesC<RHP> eval_conj(size_t n, momentum_configuration<RHP>& mc,int mu_index);
	SeriesC<RVHP> eval_conj(size_t n, momentum_configuration<RVHP>& mc,int mu_index);
#ifdef BH_USE_GMP
	SeriesC<RGMP> eval_conj(size_t n, momentum_configuration<RGMP>& mc,int mu_index);
#endif
	//! Evaluate the tree with the n th index vector in the list
	/** n is zero based */
	template <class T> std::complex<T> eval_tree(size_t n, momentum_configuration<T>& mc,int mu_index);
	//! Checks whether the index vector indices is already in the list, if yes returns its position, if not inserts it and returns its position.
	/** the index n is zero based */
	size_t add(const std::vector<int>& indices);
	//! destructor
	virtual ~Cached_OLHA();
	//! prints some statistics
	void print_stat();
	const process& get_process() const ;
	const std::vector<int>& get_index_vector(int n) const {return d_index_vectors[n];} ;
	color_structure get_color_structure() const;
	void dry_run();
	void refresh(momentum_configuration<R>& mc, int mu_index);
	template <class T> inline void set_value(int n,const SeriesC<T>& val);
	template <class T> inline void set_conj_value(int n,const SeriesC<T>& val);
	template <class T> inline SeriesC<T> get_value(int n);
	template <class T> inline SeriesC<T> get_conj_value(int n);
	template <class T> inline void set_tree_value(int n,const std::complex<T>& val);
	template <class T> inline std::complex<T> get_tree_value(int n);
	void set_accuracy(int n,const double accuracy){d_accuracy[n]=accuracy;};
	double get_accuracy(int n){return d_accuracy[n];};
    template <class T> inline void set_mcID(int n,long id);
	template <class T> inline long get_mcID(int n);
	template <class T> inline void set_mu_index(int n,int id);
	template <class T> inline int get_mu_index(int n);
};

template <> inline void Cached_OLHA::set_value(int n,const SeriesC<R>& val){d_values[n]=val;}
template <> inline void Cached_OLHA::set_value(int n,const SeriesC<RHP>& val){d_values[n]=to_double(val);d_values_HP[n]=val;}
template <> inline void Cached_OLHA::set_value(int n,const SeriesC<RVHP>& val){d_values[n]=to_double(val);d_values_HP[n]=to_HP(val);d_values_VHP[n]=val;}


template <> inline void Cached_OLHA::set_conj_value(int n,const SeriesC<R>& val){d_conj_values[n]=val;}
template <> inline void Cached_OLHA::set_conj_value(int n,const SeriesC<RHP>& val){d_conj_values[n]=to_double(val);d_conj_values_HP[n]=val;}
template <> inline void Cached_OLHA::set_conj_value(int n,const SeriesC<RVHP>& val){d_conj_values[n]=to_double(val);d_conj_values_HP[n]=to_HP(val);d_conj_values_VHP[n]=val;}


template <> inline void Cached_OLHA::set_tree_value(int n,const std::complex<R>& val){d_tree_values[n]=val;}
template <> inline void Cached_OLHA::set_tree_value(int n,const std::complex<RHP>& val){d_tree_values[n]=to_double(val);d_tree_values_HP[n]=val;}
template <> inline void Cached_OLHA::set_tree_value(int n,const std::complex<RVHP>& val){d_tree_values[n]=to_double(val);d_tree_values_HP[n]=to_HP(val);d_tree_values_VHP[n]=val;}

template <> inline SeriesC<R> Cached_OLHA::get_value(int n){ return d_values[n];}
template <> inline SeriesC<RHP> Cached_OLHA::get_value(int n){ return d_values_HP[n];}
template <> inline SeriesC<RVHP> Cached_OLHA::get_value(int n){ return d_values_VHP[n];}

template <> inline SeriesC<R> Cached_OLHA::get_conj_value(int n){ return d_conj_values[n];}
template <> inline SeriesC<RHP> Cached_OLHA::get_conj_value(int n){ return d_conj_values_HP[n];}
template <> inline SeriesC<RVHP> Cached_OLHA::get_conj_value(int n){ return d_conj_values_VHP[n];}

template <> inline std::complex<R> Cached_OLHA::get_tree_value(int n){ return d_tree_values[n];}
template <> inline std::complex<RHP> Cached_OLHA::get_tree_value(int n){ return d_tree_values_HP[n];}
template <> inline std::complex<RVHP> Cached_OLHA::get_tree_value(int n){ return d_tree_values_VHP[n];}

template <> inline void Cached_OLHA::set_mcID<R>(int n,long id){ d_mcIDs[n]=id;}
template <> inline void Cached_OLHA::set_mcID<RHP>(int n,long id){ d_mcIDs_HP[n]=id;}
template <> inline void Cached_OLHA::set_mcID<RVHP>(int n,long id){ d_mcIDs_VHP[n]=id;}

template <> inline long Cached_OLHA::get_mcID<R>(int n){ return d_mcIDs[n];}
template <> inline long Cached_OLHA::get_mcID<RHP>(int n){ return d_mcIDs_HP[n];}
template <> inline long Cached_OLHA::get_mcID<RVHP>(int n){ return d_mcIDs_VHP[n];}

template <> inline void Cached_OLHA::set_mu_index<R>(int n,int id){ d_mu_index[n]=id;}
template <> inline void Cached_OLHA::set_mu_index<RHP>(int n,int id){ d_mu_index_HP[n]=id;}
template <> inline void Cached_OLHA::set_mu_index<RVHP>(int n,int id){ d_mu_index_VHP[n]=id;}

template <> inline int Cached_OLHA::get_mu_index<R>(int n){ return d_mu_index[n];}
template <> inline int Cached_OLHA::get_mu_index<RHP>(int n){ return d_mu_index_HP[n];}
template <> inline int Cached_OLHA::get_mu_index<RVHP>(int n){ return d_mu_index_VHP[n];}

#if BH_USE_GMP

template <> inline void Cached_OLHA::set_value(int n,const SeriesC<RGMP>& val){d_values[n]=to_double(val);d_values_HP[n]=to_HP(val);d_values_VHP[n]=to_VHP(val);d_values_GMP[n]=val;}
template <> inline void Cached_OLHA::set_conj_value(int n,const SeriesC<RGMP>& val){d_conj_values[n]=to_double(val);d_conj_values_HP[n]=to_HP(val);d_conj_values_VHP[n]=to_VHP(val);d_conj_values_GMP[n]=val;}
template <> inline void Cached_OLHA::set_tree_value(int n,const std::complex<RGMP>& val){d_tree_values[n]=to_double(val);d_tree_values_HP[n]=to_HP(val);d_tree_values_VHP[n]=to_VHP(val);d_tree_values_GMP[n]=val;}
template <> inline SeriesC<RGMP> Cached_OLHA::get_value(int n){ return d_values_GMP[n];}
template <> inline SeriesC<RGMP> Cached_OLHA::get_conj_value(int n){ return d_conj_values_GMP[n];}
template <> inline std::complex<RGMP> Cached_OLHA::get_tree_value(int n){ return d_tree_values_GMP[n];}
template <> inline void Cached_OLHA::set_mcID<RGMP>(int n,long id){ d_mcIDs_GMP[n]=id;}
template <> inline long Cached_OLHA::get_mcID<RGMP>(int n){ return d_mcIDs_GMP[n];}
template <> inline void Cached_OLHA::set_mu_index<RGMP>(int n,int id){ d_mu_index_GMP[n]=id;}
template <> inline int Cached_OLHA::get_mu_index<RGMP>(int n){ return d_mu_index_GMP[n];}

#endif



class Cached_OLHA_user  {
protected:
	Cached_OLHA* d_COLHA_p;
	size_t d_n;
public:
	//Cached_OLHA_user(Cached_OLHA& COLHA,size_t n): d_COLHA_p(&COLHA), d_n(n) {};
	Cached_OLHA_user(Cached_OLHA* COLHA_p,size_t n): d_COLHA_p(COLHA_p), d_n(n) {};
	virtual SeriesC<R> eval(momentum_configuration<R>& mc,int mu_index){return SeriesC<R>(-2,0);};
	virtual SeriesC<RHP> eval(momentum_configuration<RHP>& mc,int mu_index){return SeriesC<RHP>(-2,0);};
	virtual SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,int mu_index){return SeriesC<RVHP>(-2,0);};

	std::complex<R> get_tree(momentum_configuration<R>& mc,int mu_index){return std::complex<R>(0,0);}; // needs mu so if we need to compute the OLHA, we know it
	std::complex<RHP> get_tree(momentum_configuration<RHP>& mc,int mu_index){return std::complex<RHP>(0,0);}; // needs mu so if we need to compute the OLHA, we know it
	std::complex<RVHP> get_tree(momentum_configuration<RVHP>& mc,int mu_index){return std::complex<RVHP>(0,0);}; // needs mu so if we need to compute the OLHA, we know it
	
	void dry_run(){ d_COLHA_p->dry_run();}
	//const process& get_process(){return d_COLHA_p->get_process();};
	//const std::vector<int>& get_index_vector(){return d_COLHA_p->get_index_vector(d_n); }
	//color_structure color_struct(){return d_COLHA_p->get_color_structure();};
	virtual ~Cached_OLHA_user(){};
};


class Cached_OLHA_user_normal: public Cached_OLHA_user  {
public:
	Cached_OLHA_user_normal(Cached_OLHA* COLHA_p,size_t n): Cached_OLHA_user(COLHA_p, n) {};
	virtual SeriesC<R> eval(momentum_configuration<R>& mc,int mu_index){ return d_COLHA_p->eval(d_n,mc,mu_index);};
	virtual SeriesC<RHP> eval(momentum_configuration<RHP>& mc,int mu_index){ return d_COLHA_p->eval(d_n,mc,mu_index);};
	virtual SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,int mu_index){ return d_COLHA_p->eval(d_n,mc,mu_index);};

	virtual std::complex<R> get_tree(momentum_configuration<R>& mc,int mu_index){ return d_COLHA_p->eval_tree(d_n,mc,mu_index);};
	virtual std::complex<RHP> get_tree(momentum_configuration<RHP>& mc,int mu_index){ return d_COLHA_p->eval_tree(d_n,mc,mu_index);};
	virtual std::complex<RVHP> get_tree(momentum_configuration<RVHP>& mc,int mu_index){ return d_COLHA_p->eval_tree(d_n,mc,mu_index);};
	virtual ~Cached_OLHA_user_normal();
};

class Cached_OLHA_user_conjugate: public Cached_OLHA_user  {
protected:
    std::complex<R> d_sign_R;
    std::complex<RHP> d_sign_RHP;
    std::complex<RVHP> d_sign_RVHP;
#if BH_USE_GMP
    std::complex<RGMP> d_sign_RGMP;
#endif
public:
#if BH_USE_GMP
	Cached_OLHA_user_conjugate(Cached_OLHA* COLHA_p,size_t n, int sign): Cached_OLHA_user(COLHA_p, n), d_sign_R(sign,0), d_sign_RHP(sign,0), d_sign_RVHP(sign,0), d_sign_RGMP(sign,0)  {};
#else
	Cached_OLHA_user_conjugate(Cached_OLHA* COLHA_p,size_t n, int sign): Cached_OLHA_user(COLHA_p, n), d_sign_R(sign,0), d_sign_RHP(sign,0), d_sign_RVHP(sign,0)  {};
#endif

	virtual SeriesC<R> eval(momentum_configuration<R>& mc,int mu_index){ return d_sign_R*(d_COLHA_p->eval_conj(d_n,mc,mu_index));};
	virtual SeriesC<RHP> eval(momentum_configuration<RHP>& mc,int mu_index){ return d_sign_RHP*(d_COLHA_p->eval_conj(d_n,mc,mu_index));};
	virtual SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,int mu_index){ return d_sign_RVHP*(d_COLHA_p->eval_conj(d_n,mc,mu_index));};

	virtual std::complex<R> get_tree(momentum_configuration<R>& mc,int mu_index){ return d_sign_R*conj(d_COLHA_p->eval_tree(d_n,mc,mu_index));};
	virtual std::complex<RHP> get_tree(momentum_configuration<RHP>& mc,int mu_index){ return d_sign_RHP*conj(d_COLHA_p->eval_tree(d_n,mc,mu_index));};
	virtual std::complex<RVHP> get_tree(momentum_configuration<RVHP>& mc,int mu_index){ return d_sign_RVHP*conj(d_COLHA_p->eval_tree(d_n,mc,mu_index));};
#ifdef BH_USE_GMP
	virtual std::complex<RGMP> get_tree(momentum_configuration<RGMP>& mc,int mu_index){ return d_sign_RGMP*conj(d_COLHA_p->eval_tree(d_n,mc,mu_index));};
#endif    
	virtual ~Cached_OLHA_user_conjugate();
};



class Cached_OLHA_factory {
public :
	enum rational_type { recursive, ratext };
	enum cut_type { FHZ, Darren, FHZ_wS, Darren_wS };
	Cached_OLHA_factory() {}
	virtual Cached_OLHA_user* new_OLHA(const process& pro,color_structure cs,const std::vector<int>& ind,int conj=0)=0;
	virtual void print_state()=0;
	virtual void refresh(momentum_configuration<R>& mc, int mu_index)=0;
	virtual ~Cached_OLHA_factory(){};
	static Cached_OLHA_factory* default_COLHA;

};

class pro_cs {
	private:
		process d_pro;
		color_structure d_cs;
		friend bool operator<(const pro_cs&,const pro_cs&);
	public:
		pro_cs(){};
		pro_cs(const process& pro,color_structure cs): d_pro(pro), d_cs(cs){};
		const process& get_process() const {return d_pro;};
		color_structure get_color_structure() const {return d_cs;};
		friend class Cached_OLHA_factory;
};

template <class OLHA> class Cached_OLHA_factory_impl : public Cached_OLHA_factory  {

	typename std::map<pro_cs,Cached_OLHA*> d_amplitudes;

public :
	Cached_OLHA_factory_impl() {}
	virtual Cached_OLHA_user* new_OLHA(const process& pro,color_structure cs,const std::vector<int>& ind,int conj=0);
	virtual ~Cached_OLHA_factory_impl();
	virtual void print_state();
	virtual void fill_process_list(std::vector<std::pair<process,color_structure> >&);
	virtual void refresh(momentum_configuration<R>& mc, int mu_index);
	static Cached_OLHA_factory* d_global_p;
};

template <> Cached_OLHA_user* Cached_OLHA_factory_impl<OneLoopHelAmpl>::new_OLHA(const process& pro,color_structure cs,const std::vector<int>& ind, int conj);
template <> Cached_OLHA_user* Cached_OLHA_factory_impl<IR_checked_OLHA>::new_OLHA(const process& pro,color_structure cs,const std::vector<int>& ind, int conj);

bool operator<(const pro_cs& pc1,const pro_cs& pc2);
/*
class partial_amplitude_cached : public partial_amplitude_base {
	std::vector<Cached_OLHA_user> d_ampls;
	std::vector<multi_precision_fraction> d_factors;
	std::vector<R> d_factors_R;
	std::vector<std::vector<int> > d_indices;
	std::vector<std::vector<int> > d_subs_indices;
	std::vector<subtraction*> d_subs;
public:
	partial_amplitude_cached(QCDorder lo_or_nlo=nlo) : partial_amplitude_base(lo_or_nlo) {};
	virtual ~partial_amplitude_cached();

	virtual void add(const process& pro, color_structure cs, const std::vector<int>& ind,int num, int den);
	virtual void add(const process& pro, color_structure cs, const std::vector<int>& ind, R weight);
	virtual void add_subtraction(const process& pro, const std::vector<int>& ind, multi_precision_fraction mpf,int order);
	virtual SeriesC<R> eval(momentum_configuration<R>& mc, const std::vector<int>& ind,int mu_index){return eval_fn(mc,ind,mu_index);};
	virtual SeriesC<RHP> eval(momentum_configuration<RHP>& mc, const std::vector<int>& ind,int mu_index){return eval_fn(mc,ind,mu_index);};
	virtual SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc, const std::vector<int>& ind,int mu_index){return eval_fn(mc,ind,mu_index);};
	virtual void dry_run(const std::vector<int>&);
private:
	//the copy and assignments are expensive, we restrict their use.
	partial_amplitude_cached(const partial_amplitude_cached&);
	partial_amplitude_cached& operator=(const partial_amplitude_cached&);
	template <class T> SeriesC<T> eval_fn(momentum_configuration<T>& mc, const std::vector<int>& ind,int mu_index);

};
*/

void use_IR_checked();
void use_OneLoopHelAmpl();

}
}
#endif /* CACHED_OLHA_H_ */
