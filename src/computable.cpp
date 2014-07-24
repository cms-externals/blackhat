/*!\file computable.cpp
\brief implementation for the computable objects
 */
#include "computable.h"
#include <vector>
#include <complex>
#include "mom_conf.h"
#include "eval_param.h"
#include "Series.h"
#include "qd_suppl.h"
#include "BH_utilities.h"

using namespace std;
#define _VERBOSE 0

namespace BH {


template <template <typename> class container> container<R> computable<container>::eval(momentum_configuration<R>& mc,const std::vector<int>& ind){		throw BHerror("Not implemented");
};
template <template <typename> class container> container<RHP> computable<container>::eval(momentum_configuration<RHP>&,const std::vector<int>&){		throw BHerror("Not implemented");
};
template <template <typename> class container> container<RVHP> computable<container>::eval(momentum_configuration<RVHP>&,const std::vector<int>&){		throw BHerror("Not implemented");
};
template <template <typename> class container> container<R> computable<container>::eval(const eval_param<R>& ep){		throw BHerror("Not implemented");
};
template <template <typename> class container> container<RHP> computable<container>::eval(const eval_param<RHP>& ep){		throw BHerror("Not implemented");
};
template <template <typename> class container> container<RVHP> computable<container>::eval(const eval_param<RVHP>& ep){		throw BHerror("Not implemented");
};



template <template <typename> class container> container<R> computable<container>::get_value(mom_conf& mc,const vector<int>& ind){
	if (!((mc.get_ID() == _last_ID) && (ind == _last_ind)) ){
		_last_ID = mc.get_ID();
		_last_ind = ind;
		_last_value = this->eval(mc,ind);
//		_MESSAGE4("setting value ",_last_value," for mom_conf ",_last_ID );

	}
//	else {_MESSAGE4("using old value ",_last_value," for mom_conf ",_last_ID);}
	return _last_value;
}

template <template <typename> class container> container<RHP> computable<container>::get_value(mom_conf_HP& mc,const vector<int>& ind){
	if (!((mc.get_ID() == _last_ID_HP) && (ind == _last_ind_HP)) ){
		_last_ID_HP = mc.get_ID();
		_last_ind_HP = ind;
		_last_value_HP = this->eval(mc,ind);
//		_last_ID = _last_ID_HP ;
//		_last_ind = _last_ind_HP;
//		_last_value = to_double(_last_value_HP);
	}
	return _last_value_HP;
}
template <template <typename> class container> container<RVHP>  computable<container>::get_value(mom_conf_VHP& mc,const vector<int>& ind){
	if (!((mc.get_ID() == _last_ID_VHP) && (ind == _last_ind_VHP)) ){
		_last_ID_VHP = mc.get_ID();
		_last_ind_VHP = ind;
		_last_value_VHP = this->eval(mc,ind);
//		_last_ID = _last_ID_VHP ;
//		_last_ind = _last_ind_VHP;
//		_last_value = to_double(_last_value_VHP);
//		_last_ID_HP = _last_ID_VHP ;
//		_last_ind_HP = _last_ind_VHP;
//		_last_value_HP = to_HP(_last_value_VHP);
	}
	return _last_value_VHP;
}

#if BH_USE_GMP
template <template <typename> class container> container<RGMP>  computable<container>::get_value(momentum_configuration<RGMP>& mc,const vector<int>& ind){
	if (!((mc.get_ID() == _last_ID_GMP) && (ind == _last_ind_GMP)) ){
		_last_ID_GMP = mc.get_ID();
		_last_ind_GMP = ind;
		_last_value_GMP = this->eval(mc,ind);
	}
	return _last_value_GMP;
}
#endif


template <template <typename> class container> container<R> computable<container>::get_value(const eval_param<R>& ep){
	if ((ep.get_ID() != _last_ID)){
		_last_ID = ep.get_ID();
		_last_value = this->eval(ep);
//		_MESSAGE4("setting value ",_last_value," for mom_conf ",_last_ID );

	}
//	else {_MESSAGE4("using old value ",_last_value," for mom_conf ",_last_ID);}
	return _last_value;
}

template <template <typename> class container> container<RHP> computable<container>::get_value(const eval_param<RHP>& ep){
	if ((ep.get_ID() != _last_ID_HP)){
		_last_ID_HP = ep.get_ID();
		_last_value_HP = this->eval(ep);
	}
	return _last_value_HP;
}

template <template <typename> class container> container<RVHP> computable<container>::get_value(const eval_param<RVHP>& ep){
	if ((ep.get_ID() != _last_ID_VHP)){
		_last_ID_VHP = ep.get_ID();
		_last_value_VHP = this->eval(ep);
	}
	return _last_value_VHP;
}

#if BH_USE_GMP
//template <template <typename> class container> container<RVHP> computable<container>::get_value(const eval_param<RVHP>& ep){
//	if ((ep.get_ID() != _last_ID_GMP)){
//		_last_ID_GMP = ep.get_ID();
//		_last_value_GMP = this->eval(ep);
//	}
//	return _last_value_GMP;
//}
#endif

template <template <typename> class container> container<R> zero_checked_computable<container>::get_value(mom_conf& mc,const std::vector<int>& ind){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return C(0.0);
	}
//		_MESSAGE("-----------------------------------------------");
	if (!((mc.get_ID() ==computable<container>::_last_ID) && (ind == computable<container>::_last_ind)) ){
		computable<container>::_last_ID = mc.get_ID();
		computable<container>::_last_ind = ind;
		computable<container>::_last_value = _true_computable->eval(mc,ind);
	}
	if ( abs( computable<container>::_last_value) < _C_tolerance ) {
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value," tagged as zero... (passed:",_passed_C ," unknown: ",_unknown_C," failed: ",_failed_C, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		_passed_C++;
		if ( _passed_C >= _min_passed_C && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
		if ( abs( computable<container>::_last_value) > _C_hard_limit ) {
			_failed_C++;
		}
		else _unknown_C++;
	}

	return computable<container>::_last_value;
}

template <template <typename> class container> container<RHP> zero_checked_computable<container>::get_value(mom_conf_HP& mc,const std::vector<int>& ind){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return CHP(0,0);
	}
#if _VERBOSE
		_MESSAGE("-----------------------------------------------");
#endif
		if (!((mc.get_ID() ==computable<container>::_last_ID_HP) && (ind == computable<container>::_last_ind_HP)) ){
		computable<container>::_last_ID_HP = mc.get_ID();
		computable<container>::_last_ind_HP = ind;
		computable<container>::_last_value_HP = _true_computable->eval(mc,ind);
	}
	if ( abs( computable<container>::_last_value_HP) < _CHP_tolerance ) {
		_passed_CHP++;
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value_HP," tagged as zero... (passed:",_passed_CHP ," unknown: ",_unknown_CHP," failed: ",_failed_CHP, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		if ( _passed_CHP >= _min_passed_CHP && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value_HP," not tagged as zero... (passed:",_passed_CHP ," unknown: ",_unknown_CHP," failed: ",_failed_CHP, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		if ( abs( computable<container>::_last_value_HP) > _CHP_hard_limit ) {
			_failed_CHP++;
		}
		else _unknown_CHP++;
	}

	return computable<container>::_last_value_HP;
}

template <template <typename> class container> container<RVHP> zero_checked_computable<container>::get_value(mom_conf_VHP& mc,const std::vector<int>& ind){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return CVHP(0,0);
	}
//		_MESSAGE("-----------------------------------------------");
	if (!((mc.get_ID() ==computable<container>::_last_ID_VHP) && (ind == computable<container>::_last_ind_VHP)) ){
		computable<container>::_last_ID_VHP = mc.get_ID();
		computable<container>::_last_ind_VHP = ind;
		computable<container>::_last_value_VHP = _true_computable->eval(mc,ind);
	}
	if ( abs( computable<container>::_last_value_VHP) < _CVHP_tolerance ) {
//		_MESSAGE7("tagged as zero... (passed:",_passed ," unknown: ",_unknown," failed: ",_failed, ")");
//		_MESSAGE("-----------------------------------------------");
		_passed_CVHP++;
		if ( _passed_CVHP >= _min_passed_CVHP && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
		if ( abs( computable<container>::_last_value_VHP) > _CVHP_hard_limit ) {
			_failed_CVHP++;
		}
		else _unknown_CVHP++;
	}

	return computable<container>::_last_value_VHP;
}



template <template <typename> class container> container<R> zero_checked_computable<container>::get_value(const eval_param<R>& ep){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return C(0.0);
	}
//		_MESSAGE("-----------------------------------------------");
	if (!((ep.get_ID() ==computable<container>::_last_ID)) ){
		computable<container>::_last_ID = ep.get_ID();
		computable<container>::_last_value = _true_computable->eval(ep);
	}
	if ( abs( computable<container>::_last_value) < _C_tolerance ) {
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value," tagged as zero... (passed:",_passed_C ," unknown: ",_unknown_C," failed: ",_failed_C, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		_passed_C++;
		if ( _passed_C >= _min_passed_C && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
		if ( abs( computable<container>::_last_value) > _C_hard_limit ) {
			_failed_C++;
		}
		else _unknown_C++;
	}

	return computable<container>::_last_value;
}

template <template <typename> class container> container<RHP> zero_checked_computable<container>::get_value(const eval_param<RHP>& ep){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return CHP(0,0);
	}
#if _VERBOSE
		_MESSAGE("-----------------------------------------------");
#endif
		if (!((ep.get_ID() ==computable<container>::_last_ID_HP)) ){
		computable<container>::_last_ID_HP = ep.get_ID();
		computable<container>::_last_value_HP = _true_computable->eval(ep);
	}
	if ( abs( computable<container>::_last_value_HP) < _CHP_tolerance ) {
		_passed_CHP++;
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value_HP," tagged as zero... (passed:",_passed_CHP ," unknown: ",_unknown_CHP," failed: ",_failed_CHP, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		if ( _passed_CHP >= _min_passed_CHP && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
#if _VERBOSE
		_MESSAGE8(computable<container>::_last_value_HP," not tagged as zero... (passed:",_passed_CHP ," unknown: ",_unknown_CHP," failed: ",_failed_CHP, ")");
		_MESSAGE("-----------------------------------------------");
#endif
		if ( abs( computable<container>::_last_value_HP) > _CHP_hard_limit ) {
			_failed_CHP++;
		}
		else _unknown_CHP++;
	}

	return computable<container>::_last_value_HP;
}

template <template <typename> class container> container<RVHP> zero_checked_computable<container>::get_value(const eval_param<RVHP>& ep){
	if (_is_zero ) {
//		_MESSAGE2("using zero_checked_computable for object of type: ",typeid(*_true_computable).name());
		return CVHP(0,0);
	}
//		_MESSAGE("-----------------------------------------------");
	if (!((ep.get_ID() ==computable<container>::_last_ID_VHP)) ){
		computable<container>::_last_ID_VHP = ep.get_ID();
		computable<container>::_last_value_VHP = _true_computable->eval(ep);
	}
	if ( abs( computable<container>::_last_value_VHP) < _CVHP_tolerance ) {
//		_MESSAGE7("tagged as zero... (passed:",_passed ," unknown: ",_unknown," failed: ",_failed, ")");
//		_MESSAGE("-----------------------------------------------");
		_passed_CVHP++;
		if ( _passed_CVHP >= _min_passed_CVHP && _failed_C == 0 && _failed_CHP == 0&& _failed_CVHP == 0) {
			_is_zero=true;
		}
	}
	else {
		if ( abs( computable<container>::_last_value_VHP) > _CVHP_hard_limit ) {
			_failed_CVHP++;
		}
		else _unknown_CVHP++;
	}

	return computable<container>::_last_value_VHP;
}



//
//
//
//template <template <typename> class container> container<R> limited_size_ZCC<container>::get_value(mom_conf& mc,const std::vector<int>& ind){
//	if (zero_checked_computable<container>::_is_zero ) {
////		_MESSAGE2("using limited_size_ZCC for object of type: ",typeid(*_true_computable).name());
//		return C(0.0);
//	}
////		_MESSAGE("-----------------------------------------------");
//	bool passed=true;
//	if (zero_checked_computable<container>::_last_ind.size() < nbr_indices || ind.size()< nbr_indices  ) {passed=false;}
//	else {
//		for (size_t k=0;k<nbr_indices;k++){
//			if (zero_checked_computable<container>::_last_ind[k]!=ind[k]){ passed=false; }
//		}
//	}
//	if (!((passed && mc.get_ID() ==computable<container>::_last_ID) ) ){
//		computable<container>::_last_ID = mc.get_ID();
//		computable<container>::_last_ind = ind;
//		computable<container>::_last_value = zero_checked_computable<container>::_true_computable->eval(mc,ind);
//	}
//	if ( abs( computable<container>::_last_value) < zero_checked_computable<container>::_C_tolerance ) {
//#if _VERBOSE
//		_MESSAGE8(computable<container>::_last_value," tagged as zero... (passed:",_passed_C ," unknown: ",_unknown_C," failed: ",_failed_C, ")");
//		_MESSAGE("-----------------------------------------------");
//#endif
//		zero_checked_computable<container>::_passed_C++;
//		if (zero_checked_computable<container>::_passed_C >= zero_checked_computable<container>::_min_passed_C && zero_checked_computable<container>::_failed_C == 0 && zero_checked_computable<container>::_failed_CHP == 0&& zero_checked_computable<container>::_failed_CVHP == 0) {
//			zero_checked_computable<container>::_is_zero=true;
//		}
//	}
//	else {
//		if ( abs( computable<container>::_last_value) > zero_checked_computable<container>::_C_hard_limit ) {
//			zero_checked_computable<container>::_failed_C++;
//		}
//		else zero_checked_computable<container>::_unknown_C++;
//	}
//
//	return computable<container>::_last_value;
//}
//
//template <template <typename> class container> container<RHP> limited_size_ZCC<container>::get_value(mom_conf_HP& mc,const std::vector<int>& ind){
//	if (zero_checked_computable<container>::_is_zero ) {
////		_MESSAGE2("using limited_size_ZCC for object of type: ",typeid(*_true_computable).name());
//		return CHP(0,0);
//	}
//#if _VERBOSE
//		_MESSAGE("-----------------------------------------------");
//#endif
//		bool passed=true;
//		if (zero_checked_computable<container>::_last_ind_HP.size() < nbr_indices || ind.size()< nbr_indices  ) {passed=false;}
//		else {
//		for (size_t k=0;k<nbr_indices;k++){
//			if (zero_checked_computable<container>::_last_ind_HP[k]!=ind[k]){ passed=false; }
//		}
//		}
//		if (!(passed && (mc.get_ID() ==computable<container>::_last_ID_HP) ) ){
//		computable<container>::_last_ID_HP = mc.get_ID();
//		computable<container>::_last_ind_HP = ind;
//		computable<container>::_last_value_HP = zero_checked_computable<container>::_true_computable->eval(mc,ind);
//	}
//	if ( abs( computable<container>::_last_value_HP) < zero_checked_computable<container>::_CHP_tolerance ) {
//		zero_checked_computable<container>::_passed_CHP++;
//#if _VERBOSE
//		_MESSAGE8(computable<container>::_last_value_HP," tagged as zero... (passed:",_passed_CHP ," unknown: ",_unknown_CHP," failed: ",_failed_CHP, ")");
//		_MESSAGE("-----------------------------------------------");
//#endif
//		if ( zero_checked_computable<container>::_passed_CHP >= zero_checked_computable<container>::_min_passed_CHP && zero_checked_computable<container>::_failed_C == 0 && zero_checked_computable<container>::_failed_CHP == 0&& zero_checked_computable<container>::_failed_CVHP == 0) {
//			zero_checked_computable<container>::_is_zero=true;
//		}
//	}
//	else {
//#if _VERBOSE
//		_MESSAGE8(computable<container>::_last_value_HP," not tagged as zero... (passed:",zero_checked_computable<container>::_passed_CHP ," unknown: ",zero_checked_computable<container>::_unknown_CHP," failed: ",zero_checked_computable<container>::_failed_CHP, ")");
//		_MESSAGE("-----------------------------------------------");
//#endif
//		if ( abs( computable<container>::_last_value_HP) > zero_checked_computable<container>::_CHP_hard_limit ) {
//			zero_checked_computable<container>::_failed_CHP++;
//		}
//		else zero_checked_computable<container>::_unknown_CHP++;
//	}
//
//	return computable<container>::_last_value_HP;
//}
//
//template <template <typename> class container> container<RVHP> limited_size_ZCC<container>::get_value(mom_conf_VHP& mc,const std::vector<int>& ind){
//	if (zero_checked_computable<container>::_is_zero ) {
////		_MESSAGE2("using limited_size_ZCC for object of type: ",typeid(*_true_computable).name());
//		return CVHP(0,0);
//	}
////		_MESSAGE("-----------------------------------------------");
//	bool passed=true;
//	if (zero_checked_computable<container>::_last_ind_VHP.size() < nbr_indices || ind.size()< nbr_indices  ) {passed=false;}
//	else {
//for (size_t k=0;k<nbr_indices;k++){
//		if (zero_checked_computable<container>::_last_ind_VHP[k]!=ind[k]){ passed=false; }
//	}}
//	long getID=mc.get_ID();
//	long lID=computable<container>::_last_ID_VHP;
//	bool test1= (getID == lID);
//	bool test2=(passed && test1 );
//	bool test = (!(test2));
//	if (test ){
//		computable<container>::_last_ID_VHP = mc.get_ID();
//		computable<container>::_last_ind_VHP = ind;
//		computable<container>::_last_value_VHP = zero_checked_computable<container>::_true_computable->eval(mc,ind);
//	}
//	if ( abs( computable<container>::_last_value_VHP) < zero_checked_computable<container>::_CVHP_tolerance ) {
////		_MESSAGE7("tagged as zero... (passed:",_passed ," unknown: ",_unknown," failed: ",_failed, ")");
////		_MESSAGE("-----------------------------------------------");
//		zero_checked_computable<container>::_passed_CVHP++;
//		if ( zero_checked_computable<container>::_passed_CVHP >= zero_checked_computable<container>::_min_passed_CVHP && zero_checked_computable<container>::_failed_C == 0 && zero_checked_computable<container>::_failed_CHP == 0&& zero_checked_computable<container>::_failed_CVHP == 0) {
//			zero_checked_computable<container>::_is_zero=true;
//		}
//	}
//	else {
//		if ( abs( computable<container>::_last_value_VHP) > zero_checked_computable<container>::_CVHP_hard_limit ) {
//			zero_checked_computable<container>::_failed_CVHP++;
//		}
//		else zero_checked_computable<container>::_unknown_CVHP++;
//	}
//
//	return computable<container>::_last_value_VHP;
//}
//
//
//
//




template <template <typename> class container> R zero_checked_computable<container>::_C_tolerance(1e-7);
template <template <typename> class container> R zero_checked_computable<container>::_C_hard_limit(1e-5);
template <template <typename> class container> R zero_checked_computable<container>::_CHP_tolerance(1e-16);
template <template <typename> class container> R zero_checked_computable<container>::_CHP_hard_limit(1e-12);
template <template <typename> class container> R zero_checked_computable<container>::_CVHP_tolerance(1e-32);
template <template <typename> class container> R zero_checked_computable<container>::_CVHP_hard_limit(1e-24);
template <template <typename> class container> size_t zero_checked_computable<container>::_min_passed_C(5);
template <template <typename> class container> size_t zero_checked_computable<container>::_min_passed_CHP(2);
template <template <typename> class container> size_t zero_checked_computable<container>::_min_passed_CVHP(1);

template <template <typename> class container> void zero_checked_computable<container>::set_defaults(R tolerance,R hard_limit, size_t min_passed){
	zero_checked_computable<container>::_C_tolerance=tolerance;
	zero_checked_computable<container>::_C_hard_limit=hard_limit;
	zero_checked_computable<container>::_min_passed_C=min_passed;
}
template <template <typename> class container> void zero_checked_computable<container>::set_HP_defaults(R tolerance,R hard_limit, size_t min_passed){
	zero_checked_computable<container>::_CHP_tolerance=tolerance;
	zero_checked_computable<container>::_CHP_hard_limit=hard_limit;
	zero_checked_computable<container>::_min_passed_CHP=min_passed;
}
template <template <typename> class container> void zero_checked_computable<container>::set_VHP_defaults(R tolerance,R hard_limit, size_t min_passed){
	zero_checked_computable<container>::_CVHP_tolerance=tolerance;
	zero_checked_computable<container>::_CVHP_hard_limit=hard_limit;
	zero_checked_computable<container>::_min_passed_CVHP=min_passed;
}

template <template <typename> class container> container<R> zero_checked_computable<container>::eval(mom_conf& mc,const std::vector<int>& ind){
	return _true_computable->eval(mc,ind);
}
template <template <typename> class container> container<RHP> zero_checked_computable<container>::eval(mom_conf_HP& mc,const std::vector<int>& ind){
	return _true_computable->eval(mc,ind);
}
template <template <typename> class container> container<RVHP> zero_checked_computable<container>::eval(mom_conf_VHP& mc,const std::vector<int>& ind){
	return _true_computable->eval(mc,ind);
}

template <template <typename> class container> container<R> zero_checked_computable<container>::eval(const eval_param<R>& ep){
	return _true_computable->eval(ep);
}
template <template <typename> class container> container<RHP> zero_checked_computable<container>::eval(const eval_param<RHP>& ep){
	return _true_computable->eval(ep);
}
template <template <typename> class container> container<RVHP> zero_checked_computable<container>::eval(const eval_param<RVHP>& ep){
	return _true_computable->eval(ep);
}


#if BH_USE_GMP
template <template <typename> class container> container<RGMP> computable<container>::get_value(const eval_param<RGMP>& ep){
	_WARNING3("RGMP no yet implemented for (eval_param) for objects of the type ",typeid(*this).name(),".");
	throw BHerror("Not implemented");
};
template <template <typename> class container> container<RGMP> computable<container>::eval(momentum_configuration<RGMP>&,const std::vector<int>&){
	_WARNING3("RGMP no yet implemented (for mom_conf) for objects of the type ",typeid(*this).name(),".");
	throw BHerror("Not implemented");
};
template <template <typename> class container> container<RGMP> computable<container>::eval(const eval_param<RGMP>& ep){ _WARNING3("eval not yet implemented for RGMP for type: ",typeid(*this).name(),".");};
#endif


// EXPLICIT INSTANCIATION
template class computable<std::complex>;
template class computable<SeriesC>;
template class zero_checked_computable<std::complex > ;
//template class limited_size_ZCC<std::complex > ;

}
