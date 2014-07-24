/*!\file computable.h
\brief header for the computable objects
 */

#ifndef COMPUTABLE_H_
#define COMPUTABLE_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"

#ifdef BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {
//! class for objects that are evaluated with a mom_conf and an index vector.
/** Objects of this class have three eval members, one for each precision level. Their values are cached for later use.  */

template <class T> class momentum_configuration;
template <class T> class eval_param;

template <template <typename> class containerType> class computable {
protected:

	long int _last_ID;
	long int _last_ID_HP;
	long int _last_ID_VHP;
	containerType<R> _last_value;
	containerType<RHP> _last_value_HP;
	containerType<RVHP> _last_value_VHP;
	std::vector<int> _last_ind;
	std::vector<int> _last_ind_HP;
	std::vector<int> _last_ind_VHP;
#ifdef BH_USE_GMP
	long int _last_ID_GMP;
	containerType<RGMP> _last_value_GMP;
	std::vector<int> _last_ind_GMP;
#endif


public:
	computable() :
		_last_ID(-1),
		_last_ID_HP(-1),
		_last_ID_VHP(-1)
#ifdef BH_USE_GMP
		,_last_ID_GMP(-1)
#endif
		{}  ;
	virtual containerType<double> get_value(momentum_configuration<R>& mc,const std::vector<int>& ind);
	virtual containerType<dd_real> get_value(momentum_configuration<RHP>& mc,const std::vector<int>& ind);
	virtual containerType<qd_real> get_value(momentum_configuration<RVHP>& mc,const std::vector<int>& ind);

	virtual containerType<double> get_value(const eval_param<R>& ep);
	virtual containerType<dd_real> get_value(const eval_param<RHP>& ep);
	virtual containerType<qd_real> get_value(const eval_param<RVHP>& ep);

	//Keep these as non virtual until we add definitions of eval(ep) into cutD, bubbleD etc in partitions.h
#ifndef NOTSWIG
	virtual containerType<double> eval(momentum_configuration<R>& mc,const std::vector<int>& ind);
	virtual containerType<dd_real> eval(momentum_configuration<RHP>&,const std::vector<int>&);
	virtual containerType<qd_real> eval(momentum_configuration<RVHP>&,const std::vector<int>&);
	virtual containerType<double> eval(const eval_param<R>& ep);
	virtual containerType<dd_real> eval(const eval_param<RHP>& ep);
	virtual containerType<qd_real> eval(const eval_param<RVHP>& ep);
#else
	virtual containerType<R> eval(momentum_configuration<R>& mc,const std::vector<int>& ind)=0;
	virtual containerType<RHP> eval(momentum_configuration<RHP>&,const std::vector<int>&)=0;
	virtual containerType<RVHP> eval(momentum_configuration<RVHP>&,const std::vector<int>&)=0;
	virtual containerType<R> eval(const eval_param<R>& ep)=0;
	virtual containerType<RHP> eval(const eval_param<RHP>& ep)=0;
	virtual containerType<RVHP> eval(const eval_param<RVHP>& ep)=0;
#endif

#ifdef BH_USE_GMP
	virtual containerType<RGMP> get_value(momentum_configuration<RGMP>& mc,const std::vector<int>& ind);
	virtual containerType<RGMP> get_value(const eval_param<RGMP>& ep);
	virtual containerType<RGMP> eval(momentum_configuration<RGMP>&,const std::vector<int>&);
	virtual containerType<RGMP> eval(const eval_param<RGMP>& ep);
#endif


	virtual ~computable(){};
};


#ifdef SWIG
%template(Ccomputable) computable<std::complex>;
#endif


//class computable {
//protected:
//
//	long int _last_ID;
//	long int _last_ID_HP;
//	long int _last_ID_VHP;
//	C _last_value;
//	CHP _last_value_HP;
//	CVHP _last_value_VHP;
//	std::vector<int> _last_ind;
//	std::vector<int> _last_ind_HP;
//	std::vector<int> _last_ind_VHP;
//public:
//	virtual C get_value(mom_conf& mc,std::vector<int> ind);
//	virtual CHP get_value(mom_conf_HP& mc,std::vector<int> ind);
//	virtual CVHP get_value(mom_conf_VHP& mc,std::vector<int> ind);
//	virtual C eval(mom_conf& mc,std::vector<int> ind)=0;
//	virtual CHP eval(mom_conf_HP&,std::vector<int>)=0;
//	virtual CVHP eval(mom_conf_VHP&,std::vector<int>)=0;
//	virtual ~computable(){};
//};







//! template class for computable objects that might be zero
/** Objects of the class zero_checked_computable (and those derived from it) are computed and if they are found to be zero often enough, the evaluation will be skipped and zero will be returned immediatly. */
template <template <typename> class container> class zero_checked_computable : public computable<container> {
protected:
	size_t _passed_C, _unknown_C,_failed_C;
	size_t _passed_CHP, _unknown_CHP,_failed_CHP;
	size_t _passed_CVHP, _unknown_CVHP,_failed_CVHP;
	static  R _C_tolerance;
	static  R _C_hard_limit;
	static  R _CHP_tolerance;
	static  R _CHP_hard_limit;
	static  R _CVHP_tolerance;
	static  R _CVHP_hard_limit;
	static  size_t _min_passed_C;
	static  size_t _min_passed_CHP;
	static  size_t _min_passed_CVHP;
	bool _is_zero;
	computable<container>* _true_computable;
public :
	zero_checked_computable(computable<container>* true_computable) : _is_zero(false), _true_computable(true_computable), _passed_C(0), _unknown_C(0), _failed_C(0), _passed_CHP(0), _unknown_CHP(0), _failed_CHP(0), _passed_CVHP(0), _unknown_CVHP(0), _failed_CVHP(0) {};
	container<R> get_value(momentum_configuration<R>& mc, const std::vector<int>& ind);
	container<RHP> get_value(momentum_configuration<RHP>& mc, const std::vector<int>& ind);
	container<RVHP> get_value(momentum_configuration<RVHP>& mc, const std::vector<int>& ind);
	container<R> eval(momentum_configuration<R>& mc, const std::vector<int>& ind);
	container<RHP> eval(momentum_configuration<RHP>& mc, const std::vector<int>& ind);
	container<RVHP> eval(momentum_configuration<RVHP>& mc, const std::vector<int>& ind);

	container<R> get_value(const eval_param<R>& ep);
	container<RHP> get_value(const eval_param<RHP>& ep);
	container<RVHP> get_value(const eval_param<RVHP>& ep);
	container<R> eval(const eval_param<R>& ep);
	container<RHP> eval(const eval_param<RHP>& ep);
	container<RVHP> eval(const eval_param<RVHP>& ep);

	computable<container>* true_computable(){return _true_computable;};
	static void set_defaults(R _tolerance,R _hard_limit, size_t min_passed);
	static void set_HP_defaults(R _tolerance,R _hard_limit, size_t min_passed);
	static void set_VHP_defaults(R _tolerance,R _hard_limit, size_t min_passed);
};

//template <template <typename> class container> class limited_size_ZCC : public zero_checked_computable<container> {
//	size_t nbr_indices;
//public :
//	limited_size_ZCC(computable<container>* true_computable,size_t nbr) : zero_checked_computable<container>(true_computable) , nbr_indices(nbr){};
//	container<R> get_value(mom_conf& mc, const vector<int>& ind);
//	container<RHP> get_value(mom_conf_HP& mc, const vector<int>& ind);
//	container<RVHP> get_value(mom_conf_VHP& mc, const vector<int>& ind);
//
//};
//

}
#endif /*COMPUTABLE_H_*/
