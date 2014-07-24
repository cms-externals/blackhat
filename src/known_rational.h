/*
 * known_rational.h
 *
 *  Created on: 6 Jul 2009
 *      Author: daniel
 */

#ifndef KNOWN_RATIONAL_H_
#define KNOWN_RATIONAL_H_

#include "rational.h"
#include "settings.h"



namespace BH {

template <class T> class eval_param;
class process;
template <class T> class momentum_configuration;

bool Rec_Rational_is_zero(const process& pro,color_structure cs);

//! factory class for known rational part
class Known_Rec_Rational_base : public Rational_base {
protected:
#ifdef USE_MC_RAT
	Tree_Fn_Ptr _eval_C_ptr;
	Tree_Fn_Ptr_HP _eval_CHP_ptr;
	Tree_Fn_Ptr_VHP _eval_CVHP_ptr;
#endif
	Tree_Fn_Ptr_eval _eval_C_Ptr_eval;
	Tree_Fn_Ptr_eval_HP _eval_CHP_Ptr_eval;
	Tree_Fn_Ptr_eval_VHP _eval_CVHP_Ptr_eval;
#if BH_USE_GMP
#ifdef USE_MC_RAT
	Tree_Fn_Ptr_GMP _eval_CGMP_ptr;
#endif
	Tree_Fn_Ptr_eval_GMP _eval_CGMP_Ptr_eval;
#endif
	mass_param_coll* _masses;
public:
	Known_Rec_Rational_base(const process& amp_pro): Rational_base(amp_pro) {_masses=new mass_param_coll(amp_pro);} //Default constructor
	virtual ~Known_Rec_Rational_base(){ delete _masses;};
	
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
};


//! factory class for known rational part
class Known_Rec_Rational : public Known_Rec_Rational_base {
//	friend class Rational_factory;
public:
	Known_Rec_Rational(const process& amp_pro,color_structure cs=glue); //Default constructor to set up HelAmpl
//	Known_Rec_Rational(process& amp_pro, vector<ph_type>& possible_props); //Usual constructor but will automatically choose the correct legs to recurse with if needed
	//! constructor
	/** \param amp_pro process \param possible_props list of the possible propagators \param amp_i [ part if the shift \param < part of the shift  */
	Known_Rec_Rational(const process& amp_pro, size_t amp_i, size_t amp_j);
	//! evaluation of the rational term in simple precision
#ifdef USE_MC_RAT
	virtual C eval(mom_conf& mc,const std::vector<int>& ind){if(settings::general::s_use_ep_only){eval_param<R> ep(mc,ind); return (*_eval_C_Ptr_eval)(ep,*_masses);}else{return (*_eval_C_ptr)(mc,ind);}};
	//! evaluation of the rational term in double precision
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){if(settings::general::s_use_ep_only){eval_param<RHP> ep(mc,ind); return (*_eval_CHP_Ptr_eval)(ep,*_masses);}else{return (*_eval_CHP_ptr)(mc,ind);}};
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){if(settings::general::s_use_ep_only){eval_param<RVHP> ep(mc,ind); return (*_eval_CVHP_Ptr_eval)(ep,*_masses);}else{return (*_eval_CVHP_ptr)(mc,ind);}};
#else
	virtual C eval(mom_conf& mc,const std::vector<int>& ind){ eval_param<R> ep(mc,ind); return (*_eval_C_Ptr_eval)(ep,*_masses);};
	//! evaluation of the rational term in double precision
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){eval_param<RHP> ep(mc,ind); return (*_eval_CHP_Ptr_eval)(ep,*_masses);};
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){eval_param<RVHP> ep(mc,ind);return (*_eval_CVHP_Ptr_eval)(ep,*_masses);};

#endif
	virtual C eval(const eval_param<R>& ep){return (*_eval_C_Ptr_eval)(ep,*_masses);};
	//! evaluation of the rational term in double precision
	virtual CHP eval(const eval_param<RHP>& ep){return (*_eval_CHP_Ptr_eval)(ep,*_masses);};
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(const eval_param<RVHP>& ep){return (*_eval_CVHP_Ptr_eval)(ep,*_masses);};

#if BH_USE_GMP
#ifdef USE_MC_RAT
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){if(settings::general::s_use_ep_only){eval_param<RGMP> ep(mc,ind); return (*_eval_CGMP_Ptr_eval)(ep,*_masses);}else{return (*_eval_CGMP_ptr)(mc,ind);}};
#else
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){eval_param<RGMP> ep(mc,ind); return (*_eval_CGMP_Ptr_eval)(ep,*_masses);};
#endif
	virtual CGMP eval(const eval_param<RGMP>& ep){return (*_eval_CGMP_Ptr_eval)(ep,*_masses);};
#endif
	virtual bool is_zero();
	virtual ~Known_Rec_Rational(){};//Destructor to delete all the new rational terms we have computed
	// This gives the estimated accuracy of the computation
	double get_accuracy(){return double(16);};
};


class Known_Rec_Rational_offset : public Known_Rec_Rational_base {
	size_t _offset;
public:
	Known_Rec_Rational_offset(const process& amp_pro,color_structure cs=glue,size_t offset=0);
	//! evaluation of the rational term in simple precision
	virtual C eval(mom_conf& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	//! evaluation of the rational term in double precision
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};

	virtual C eval(const eval_param<R>& ep);
	//! evaluation of the rational term in double precision
	virtual CHP eval(const eval_param<RHP>& ep);
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(const eval_param<RVHP>& ep);

#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CGMP eval(const eval_param<RGMP>& ep);
#endif

	virtual ~Known_Rec_Rational_offset(){};//Destructor to delete all the new rational terms we have computed
private :
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){eval_param<T> ep(mc,ind); return eval(ep);};
};

class Known_Rec_Rational_permutation : public Known_Rec_Rational_base {
	std::vector<int> _perm_ind;
public:
	Known_Rec_Rational_permutation(const process& amp_pro,color_structure cs,const std::vector<int>& per);
	//! evaluation of the rational term in simple precision
	virtual C eval(mom_conf& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	//! evaluation of the rational term in double precision
	virtual CHP eval(mom_conf_HP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(mom_conf_VHP& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};

	virtual C eval(const eval_param<R>& ep);
	//! evaluation of the rational term in double precision
	virtual CHP eval(const eval_param<RHP>& ep);
	//! evaluation of the rational term in quadruple precision
	virtual CVHP eval(const eval_param<RVHP>& ep);

#if BH_USE_GMP
	virtual CGMP eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){return eval_fn(mc,ind);};
	virtual CGMP eval(const eval_param<RGMP>& ep);
#endif

	virtual ~Known_Rec_Rational_permutation(){};//Destructor to delete all the new rational terms we have computed
private :
	template <class T> std::complex<T> eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){eval_param<T> ep(mc,ind); return eval(ep);};

};


}

#endif /* KNOWN_RATIONAL_H_ */
