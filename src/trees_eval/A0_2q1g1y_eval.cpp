/*
 * A0_2q1g1y_eval.cpp
 *
 *  Created on: Sep 22, 2008
 *      Author_eval: daniel
 */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include <iostream>

using namespace std;

namespace BH {

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)


template <class T> complex<T> A_qbgqy_zero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }

template <class T> complex<T> A_qbgqy_mmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	-complex<T>(0,-1)*SPA(1,2)*SPB(4,3)/(SPA(1,4)*SPB(2,1));
 }

template <class T> complex<T> A_qbgqy_mppm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	-complex<T>(0,1)*SPA(1,4)*SPB(3,2)/(SPA(1,2)*SPB(4,1));
 }

template <class T> complex<T> A_qbgqy_pmmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	-complex<T>(0,-1)*SPA(2,3)*SPB(4,1)/(SPA(1,4)*SPB(2,1));
 }

template <class T> complex<T> A_qbgqy_ppmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	-complex<T>(0,1)*SPA(3,4)*SPB(2,1)/(SPA(1,2)*SPB(4,1));
 }


template <class T> complex<T> A_qbgyq_mmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,-1)*SPA(1,2)*SPB(3,4)/(SPA(1,3)*SPB(2,1));
 }

template <class T> complex<T> A_qbgyq_mpmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,1)*SPA(1,3)*SPB(4,2)/(SPA(1,2)*SPB(3,1));
 }

template <class T> complex<T> A_qbgyq_pmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,-1)*SPA(2,4)*SPB(3,1)/(SPA(1,3)*SPB(2,1));
 }

template <class T> complex<T> A_qbgyq_ppmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,1)*SPA(4,3)*SPB(2,1)/(SPA(1,2)*SPB(3,1));
 }



template <class T> complex<T> A_qbygq_mmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,1)*SPA(1,2)*SPB(4,3)/(SPA(1,3)*SPB(2,1));
 }

template <class T> complex<T> A_qbygq_mpmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,-1)*SPA(1,3)*SPB(2,4)/(SPA(1,2)*SPB(3,1));
 }

template <class T> complex<T> A_qbygq_pmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,1)*SPA(4,2)*SPB(3,1)/(SPA(1,3)*SPB(2,1));
 }

template <class T> complex<T> A_qbygq_ppmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return	complex<T>(0,-1)*SPA(3,4)*SPB(2,1)/(SPA(1,2)*SPB(3,1));
 }




template <class T> complex<T>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {

   case 1153: return &A_qbgqy_mmpp_eval;
   case 955: return &A_qbgqy_mppm_eval;

   case 613: return &A_qbgyq_mmpp_eval;
   case 595: return &A_qbgyq_mpmp_eval;

   case 1118: return &A_qbgqy_pmmp_eval;
   case 920: return &A_qbgqy_ppmm_eval;

   case 398: return &A_qbgyq_pmpm_eval;
   case 380: return &A_qbgyq_ppmm_eval;

   case 350: return &A_qbygq_pmpm_eval;
   case 248: return &A_qbygq_ppmm_eval;
   case 565: return &A_qbygq_mmpp_eval;
   case 463: return &A_qbygq_mpmp_eval;



   case 937: // qm, m, qp, ym
   case 1171: // qm, p, qp, yp
   case 631: // qm, p, yp, qp
   case 577: // qm, m, ym, qp
   case 457: // qm, ym, m, qp
   case 571: // qm, yp, p, qp
   case 242: // qp, m, ym, qm
   case 356: // qp, p, yp, qm

	   return &A_qbgqy_zero_eval;


   default:return 0;

}
}



template complex<R>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q1g1y_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


