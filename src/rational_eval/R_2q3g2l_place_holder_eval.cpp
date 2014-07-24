/*
 * R_2q3g2l.cpp
 *
 *  Created on: Sep 18, 2008
 *      Author: daniel
 */



#include "R_2q3g2l_eval.h"
#include "eval_param.h"


using namespace std;

namespace BH  {


#define _VERBOSE 0

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)


template <class T> complex<T> R2q3g2l_place_holder(const eval_param<T>& ep,const mass_param_coll& mpc){
	_WARNING("rational term set to zero.");
	 return complex<T>(0,0);

}






template <class T> complex<T> (*R2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {


return &R2q3g2l_place_holder;

}


template complex<R> (*R2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
template complex<RHP> (*R2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
template complex<RVHP> (*R2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

template complex<RGMP> (*R2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif



}

