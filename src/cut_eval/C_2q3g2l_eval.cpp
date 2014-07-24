/*
 * C_2q3g2l.cpp
 *
 *  Created on: Sep 18, 2008
 *      Author: daniel
 */



#include "C_2q3g2l_eval.h"
#include "tree_amp.h"
#include "integrals_ep.h"

using namespace std;

namespace BH  {


#define _VERBOSE 0

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)


template <class T> SeriesC<T> C2q3g2l_place_holder(const eval_param<T>& ep, const T& mu){
	 return SeriesC<T>(-2,0);

}






template <class T> SeriesC<T> (*C2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<T>&, const T&) {


return &C2q3g2l_place_holder;

}


template SeriesC<R> (*C2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<R>&, const R&) ;
template SeriesC<RHP> (*C2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<RHP>&, const RHP&) ;
template SeriesC<RVHP> (*C2q3g2l_place_holder_Ptr_eval(int hc))(const eval_param<RVHP>&, const RVHP&) ;



}
