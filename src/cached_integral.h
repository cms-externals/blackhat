/*
 * cached_integral.h
 *
 *  Created on: 28-Nov-2008
 *      Author: daniel
 */

#ifndef CACHED_INTEGRAL_H_
#define CACHED_INTEGRAL_H_

#include <vector>
#include <map>
#include "BH_typedefs.h"
#include "Series.h"
#include "cut_part.h"
#include "index_vector.h"
#if BH_USE_GMP
#include "gmp_r.h"
#endif


namespace BH {

template <class T> class momentum_configuration;
typedef momentum_configuration<double> mom_conf;

//! Class for index vectors
/** an Index vector is basically a vector<int> but allows to carry more inforamtion, like the permutation code. */


namespace CachedIntegral {

/*
class scalar_integral{
public:
	enum scalar_integral_type { box_integral, triangle_integral, bubble_integral };
	enum mass_type { zero_mass, one_mass, two_mass, two_mass_easy, two_mass_hard, three_mass, four_mass };
	virtual scalar_integral_type type()=0;
};

class box_scalar_integral : public scalar_integral{

public:
	virtual scalar_integral_type type(){return box_integral;};
};
*/



//! abstract class for cached integrals
class Cached_Integral {
public:
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf& mc,int mu)=0;
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_HP& mc,int mu)=0;
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_VHP& mc,int mu)=0;
#if BH_USE_GMP
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(momentum_configuration<RGMP>& mc,int mu)=0;
#endif
	virtual ~Cached_Integral(){};
};


class Cached_Integral_impl : public Cached_Integral {
	SeriesC<double> d_value;
	SeriesC<dd_real> d_value_HP;
	SeriesC<qd_real> d_value_VHP;
	//! keeps track of the number of Cached_Integral_user using this integral
	mutable size_t d_load;
	long d_mcID;
	long d_mcID_HP;
	long d_mcID_VHP;
#if BH_USE_GMP
	SeriesC<RGMP> d_value_GMP;
	long d_mcID_GMP;
#endif
public:
	Cached_Integral_impl() : d_load(0){};
	size_t get_load() const {return d_load;};
	//! gets the value for the current mom_conf
	virtual const SeriesC<double>* get() const ;
	virtual const SeriesC<dd_real>* get_HP() const;
	virtual const SeriesC<qd_real>* get_VHP() const ;
#if BH_USE_GMP
	virtual const SeriesC<RGMP>* get_GMP() const ;
#endif
//	virtual const SeriesC<double>& get_value() const {return d_value;};
//	virtual const SeriesC<dd_real>& get_value_HP() const {return d_value_HP;};
//	virtual const SeriesC<qd_real>& get_value_VHP() const {return d_value_VHP;};
	long get_mcID() const {return d_mcID;};
	long get_mcID_HP() const {return d_mcID_HP;};
	long get_mcID_VHP() const {return d_mcID_VHP;};
	void set_mcID(long mcID) { d_mcID=mcID;};
	void set_mcID_HP(long mcID) { d_mcID_HP=mcID;};
	void set_mcID_VHP(long mcID) { d_mcID_VHP=mcID;};
	void increase_load(){ ++d_load; };
	void set_value(SeriesC<R> val){d_value=val;}
	void set_value(SeriesC<RHP> val){d_value_HP=val;}
	void set_value(SeriesC<RVHP> val){d_value_VHP=val;}
#if BH_USE_GMP
	long get_mcID_GMP() const {return d_mcID_GMP;};
	void set_mcID_GMP(long mcID) { d_mcID_GMP=mcID;};
	void set_value(SeriesC<RGMP> val){d_value_GMP=val;}
#endif

	template <class T> inline long get_mc_ID() const ;
	template <class T> inline const SeriesC<T>* get_value() const ;
	virtual ~Cached_Integral_impl(){};
};

template <> inline long Cached_Integral_impl::get_mc_ID<R>() const { return d_mcID;}
template <> inline long Cached_Integral_impl::get_mc_ID<RHP>() const { return d_mcID_HP;}
template <> inline long Cached_Integral_impl::get_mc_ID<RVHP>() const { return d_mcID_VHP;}


template <> inline const SeriesC<R>* Cached_Integral_impl::get_value<R>() const { return this->get();}
template <> inline const SeriesC<RHP>* Cached_Integral_impl::get_value<RHP>() const { return this->get_HP();}
template <> inline const SeriesC<RVHP>* Cached_Integral_impl::get_value<RVHP>() const { return this->get_VHP();}

#if BH_USE_GMP
template <> inline long Cached_Integral_impl::get_mc_ID<RGMP>() const { return d_mcID_GMP;}
template <> inline const SeriesC<RGMP>* Cached_Integral_impl::get_value<RGMP>() const { return this->get_GMP();}
#endif


class Cached_Bubble_Integral : public Cached_Integral_impl {
	std::vector<int> d_index_vector_1;
	std::vector<int> d_index_vector_2;
	//! keeps track of the number of Cached_Integral_user using this integral
	size_t d_nbr_indices;
	long d_scode;
public:
	friend bool operator<(const Cached_Bubble_Integral& c1,const Cached_Bubble_Integral& c2);
	friend std::ostream& operator<<(std::ostream& s,const Cached_Bubble_Integral& c2);
	Cached_Bubble_Integral(const std::vector<int>&,const std::vector<int>&);
	long get_scode() const {return d_scode;};
	size_t get_nbr() const {return d_nbr_indices;};
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_HP& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_VHP& mc,int mu);
#if BH_USE_GMP
	virtual void prepare(momentum_configuration<RGMP>& mc,int mu);
#endif
	virtual ~Cached_Bubble_Integral(){};
};


class Cached_Triangle_Integral : public Cached_Integral_impl {
	std::vector<int> d_index_vector_1;
	std::vector<int> d_index_vector_2;
	std::vector<int> d_index_vector_3;
	//! keeps track of the number of Cached_Integral_user using this integral
	size_t d_nbr_indices;
	long d_scode1;
	long d_scode2;
	long d_scode3;
public:
	friend bool operator<(const Cached_Triangle_Integral& c1,const Cached_Triangle_Integral& c2);
	friend std::ostream& operator<<(std::ostream& s,const Cached_Triangle_Integral& c2);
	Cached_Triangle_Integral(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	long get_scode(size_t i) const {switch(i){ case 1: return d_scode1;case 2: return d_scode2;case 3: return d_scode3;}};
	size_t get_nbr() const {return d_nbr_indices;};
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_HP& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_VHP& mc,int mu);
#if BH_USE_GMP
	virtual void prepare(momentum_configuration<RGMP>& mc,int mu);
#endif
	virtual ~Cached_Triangle_Integral(){};
};

class Cached_Box_Integral : public Cached_Integral_impl {
	std::vector<int> d_index_vector_1;
	std::vector<int> d_index_vector_2;
	std::vector<int> d_index_vector_3;
	std::vector<int> d_index_vector_4;
	//! keeps track of the number of Cached_Integral_user using this integral
	size_t d_nbr_indices;
	long d_scode1;
	long d_scode2;
	long d_scode3;
	long d_scode4;
public:
	friend bool operator<(const Cached_Box_Integral& c1,const Cached_Box_Integral& c2);
	friend std::ostream& operator<<(std::ostream& s,const Cached_Box_Integral& c2);
	Cached_Box_Integral(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	long get_scode(size_t i) const {switch(i){ case 1: return d_scode1;case 2: return d_scode2;case 3: return d_scode3; case 4: return d_scode4;}};
	size_t get_nbr() const {return d_nbr_indices;};
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_HP& mc,int mu);
	//!computes all the integrals for the specified mom_conf
	virtual void prepare(mom_conf_VHP& mc,int mu);
#if BH_USE_GMP
	virtual void prepare(momentum_configuration<RGMP>& mc,int mu);
#endif
	virtual ~Cached_Box_Integral(){};
};

class Cached_Integral_Factory {
	std::vector<Cached_Box_Integral*> d_box_integrals;
	std::vector<Cached_Triangle_Integral*> d_triangle_integrals;
	std::vector<Cached_Bubble_Integral*> d_bubble_integrals;
public:
	Cached_Bubble_Integral* new_integral(const std::vector<int>&,const std::vector<int>&);
	Cached_Triangle_Integral* new_integral(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	Cached_Box_Integral* new_integral(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	void print_state();
	void prepare(mom_conf& mc,int mu);
	~Cached_Integral_Factory();
	static Cached_Integral_Factory s_default_CIF;
};

class Cached_Integral_User {
public:
	virtual void print_info()=0;
	virtual const SeriesC<R>* get_value(momentum_configuration<R>& mc,const Index_Vector& IV,int mu)=0;
	virtual const SeriesC<RHP>* get_value(momentum_configuration<RHP>& mc,const Index_Vector& IV,int mu)=0;
	virtual const SeriesC<RVHP>* get_value(momentum_configuration<RVHP>& mc,const Index_Vector& IV,int mu)=0;
#if BH_USE_GMP
	virtual const SeriesC<RGMP>* get_value(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu)=0;
#endif
	virtual ~Cached_Integral_User(){};
};

class Cached_Bubble_Integral_User: public Cached_Integral_User {
	std::vector<int> d_corner1;
	std::vector<int> d_corner2;
	std::map<long,Cached_Bubble_Integral*> d_values;
	typedef std::map<long,Cached_Bubble_Integral*>::iterator iter;
public:
	template <class T> const SeriesC<T>* get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu);
	virtual const SeriesC<double>* get_value(mom_conf& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RHP>* get_value(mom_conf_HP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RVHP>* get_value(mom_conf_VHP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#if BH_USE_GMP
	virtual const SeriesC<RGMP>* get_value(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#endif
	Cached_Bubble_Integral_User(const std::vector<int>&,const std::vector<int>&);
	void print_info();
	virtual ~Cached_Bubble_Integral_User(){};

};

class Cached_Triangle_Integral_User : public Cached_Integral_User {
	std::vector<int> d_corner1;
	std::vector<int> d_corner2;
	std::vector<int> d_corner3;
	std::map<long,Cached_Triangle_Integral*> d_values;
	typedef std::map<long,Cached_Triangle_Integral*>::iterator iter;
public:
	template <class T> const SeriesC<T>* get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu);
	virtual const SeriesC<R>* get_value(mom_conf& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RHP>* get_value(mom_conf_HP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RVHP>* get_value(mom_conf_VHP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#if BH_USE_GMP
	virtual const SeriesC<RGMP>* get_value(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#endif
	Cached_Triangle_Integral_User(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	void print_info();
	virtual ~Cached_Triangle_Integral_User(){};
};

class Cached_Box_Integral_User : public Cached_Integral_User {
	std::vector<int> d_corner1;
	std::vector<int> d_corner2;
	std::vector<int> d_corner3;
	std::vector<int> d_corner4;
	std::map<long,Cached_Box_Integral*> d_values;
	typedef std::map<long,Cached_Box_Integral*>::iterator iter;
public:
	template <class T> const SeriesC<T>* get_value_fn(momentum_configuration<T>& mc,const Index_Vector& IV,int mu);
	virtual const SeriesC<R>* get_value(mom_conf& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RHP>* get_value(mom_conf_HP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
	virtual const SeriesC<RVHP>* get_value(mom_conf_VHP& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#if BH_USE_GMP
	virtual const SeriesC<RGMP>* get_value(momentum_configuration<RGMP>& mc,const Index_Vector& IV,int mu){return get_value_fn(mc,IV,mu);};
#endif
	Cached_Box_Integral_User(const std::vector<int>&,const std::vector<int>&,const std::vector<int>&,const std::vector<int>&);
	void print_info();
	virtual ~Cached_Box_Integral_User(){};
};

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};

class Cut_Part_wCI {
protected:
	std::vector<BH::CachedIntegral::Cached_Integral_User*> CI_users;
public:
	virtual SeriesC<R> eval(momentum_configuration<R>& mc,const Index_Vector& ind,int mu)=0;
	virtual SeriesC<RHP> eval(momentum_configuration<RHP>& mc,const Index_Vector& ind,int mu)=0;
	virtual SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int mu)=0;
#if BH_USE_GMP
	virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const Index_Vector& ind,int mu)=0;
#endif
	virtual ~Cut_Part_wCI();
};


class Known_Cut_wCI : public Cut_Part_base {
	CachedIntegral::Cut_Part_wCI* d_CPwCI_p;
public:
			virtual SeriesC<R> eval(mom_conf& mc,const std::vector<int>& ind);
			virtual SeriesC<RHP> eval(mom_conf_HP& mc,const std::vector<int>& ind);
			virtual SeriesC<RVHP> eval(mom_conf_VHP& mc,const std::vector<int>& ind);
			virtual SeriesC<R> eval(const eval_param<R>&);
			virtual SeriesC<RHP> eval(const eval_param<RHP>&);
			virtual SeriesC<RVHP> eval(const eval_param<RVHP>&);
#if BH_USE_GMP
			virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif
			Known_Cut_wCI(const process&,color_structure cs);
			Known_Cut_wCI(const process&,CachedIntegral::Cut_Part_wCI*);
			~Known_Cut_wCI();
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
#if BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){return SeriesC<RGMP>(-2,0);};
#endif
};

class Known_Cut_wCI_offset : public Cut_Part_base {
	size_t d_offset;
	CachedIntegral::Cut_Part_wCI* d_CPwCI_p;
public:
		virtual SeriesC<R> eval(mom_conf& mc,const std::vector<int>& ind);
		virtual SeriesC<RHP> eval(mom_conf_HP& mc,const std::vector<int>& ind);
		virtual SeriesC<RVHP> eval(mom_conf_VHP& mc,const std::vector<int>& ind);
		virtual SeriesC<R> eval(const eval_param<R>&);
		virtual SeriesC<RHP> eval(const eval_param<RHP>&);
		virtual SeriesC<RVHP> eval(const eval_param<RVHP>&);
#if BH_USE_GMP
	virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
#endif
		Known_Cut_wCI_offset(const process&,CachedIntegral::Cut_Part_wCI*,size_t offset);
		Known_Cut_wCI_offset(const process&,color_structure cs,size_t offset);
		size_t get_offset() const {return d_offset;};
		~Known_Cut_wCI_offset();
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
#if BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){return SeriesC<RGMP>(-2,0);};
#endif

};


Cut_Part_wCI* CwCI_Ptr(const process& pro,color_structure cs,const std::vector<int>& ind);


#ifdef USE_OLD_CUT_PART

class Unknown_Cut_Part_wCI : public Cut_Part {
	std::vector<CachedIntegral::Cached_Box_Integral_User*> d_box_integrals;
	std::vector<CachedIntegral::Cached_Triangle_Integral_User*> d_triangle_integrals;
	std::vector<CachedIntegral::Cached_Bubble_Integral_User*> d_bubble_integrals;
public :
	virtual SeriesC<R> eval(mom_conf&,const std::vector<int>&);
	virtual SeriesC<RHP> eval(mom_conf_HP&,const std::vector<int>&);
	virtual SeriesC<RVHP> eval(mom_conf_VHP&,const std::vector<int>&);
#if BH_USE_GMP
	virtual SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const std::vector<int>&);
#endif
	//! Constructor
	Unknown_Cut_Part_wCI(Cut_Part* cp, cutD_factory* CF = cutD_factory::default_CF,option* opt=option::always_true) : Cut_Part(cp,CF,opt){construct_integral_table();};
	//! Constructor
	Unknown_Cut_Part_wCI(const process& p,const std::vector<particle_ID>& possible_props, cutD_factory* CF = cutD_factory::default_CF, option* opt = option::always_true): Cut_Part(p,possible_props, CF, opt){construct_integral_table();} ;
	//! Constructor
	Unknown_Cut_Part_wCI(const process& p,const std::vector<particle_ID>& possible_props, const ordering_constraint& ordered, cutD_factory* CF= cutD_factory::default_CF, option* opt= option::always_true): Cut_Part(p,possible_props,ordered,CF,opt){construct_integral_table();};
	virtual ~Unknown_Cut_Part_wCI();
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
    	SeriesC<R> get_conjugate_cut_part(){return SeriesC<R>(-2,0);};
    	SeriesC<RHP> get_conjugate_cut_part_HP(){return SeriesC<RHP>(-2,0);};
    	SeriesC<RVHP> get_conjugate_cut_part_VHP(){return SeriesC<RVHP>(-2,0);};
#if BH_USE_GMP
    	SeriesC<RGMP> get_conjugate_cut_part_GMP(){return SeriesC<RGMP>(-2,0);};
#endif
private:
	void construct_integral_table();
	template <class T> SeriesC<T> Unknown_Cut_Part_wCI_eval(momentum_configuration<T>& mc,const std::vector<int>& ind);
};

#endif

} /* cachedIntegral*/
} /* BH */

#if BH_USE_GMP
#define GENERATE_CUT_PART_WCI_CLASS_DEFINITION(NAME) \
	class NAME : public Cut_Part_wCI { \
public:\
NAME (const std::vector<int>&);\
	   SeriesC<R> eval(momentum_configuration<R>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
       SeriesC<RGMP> eval(momentum_configuration<RGMP>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
private:\
       template <class T> SeriesC<T> eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,int mu);\
};
#else
#define GENERATE_CUT_PART_WCI_CLASS_DEFINITION(NAME) \
	class NAME : public Cut_Part_wCI { \
public:\
NAME (const std::vector<int>&);\
SeriesC<R> eval(momentum_configuration<R>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
       SeriesC<RHP> eval(momentum_configuration<RHP>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
       SeriesC<RVHP> eval(momentum_configuration<RVHP>& mc,const Index_Vector& ind,int mu){return eval_fn(mc,ind,mu);};\
private:\
       template <class T> SeriesC<T> eval_fn(momentum_configuration<T>& mc,const Index_Vector& ind,int mu);\
};
#endif

#endif /* CACHED_INTEGRAL_H_ */
