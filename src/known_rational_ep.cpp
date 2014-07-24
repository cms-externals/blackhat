/*
 * known_rational_ep.cpp
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#include "known_rational.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "BH_utilities.h"
#include "BH_A0.h"
#include "trees_eval/amplitudes_tree_eval.h"
#include "rational_eval/amplitudes_rat_eval.h"
#include "process_utils.h"
#include "settings.h"
#include "tree_amp.h"
#include "helcode.h"


#define _VERBOSE 0

using namespace std;

namespace BH {

bool Rec_Rational_is_zero(const process& pro,color_structure cs);//Defined in rec_rational.cpp

template <class T> complex<T> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&,const mass_param_coll&){
	
	if ( Rec_Rational_is_zero(pro,cs) ) {
			return &ZeroF_eval;
	};

	if ( BH::settings::rational_settings::s_set_all_zero ) {
			return &ZeroF_eval;
	};

	switch(pro.pcode()){
		case odd_nbr_q:
		case odd_nbr_q2:
		case odd_nbr_l:
		case odd_nbr_l1:
		case odd_nbr_l2:

		case 200020: /* two photons two quarks should not appear */
			return &ZeroF_eval;
		break;

		case 3:
		// At one-loop the 3-pt vertex is either zero or doen't factorise in which case we
		// want to drop it anyway and calculate an Inf term
			return &ZeroF_eval;
		break;

		case 4: {
			switch (cs) {
				case glue: case LT: case leading_color: return R4g_Ptr_eval<T>(helcode_g(pro));
				case nf:   return R4g_nf_Ptr_eval<T>(helcode_g(pro));
				case RT: return R4g_SLC_Ptr_eval<T>(helcode_g(pro));
			}; break;
		}
		case 5: {
			switch (cs) {
			case glue: case LT: case leading_color: return R5g_Ptr_eval<T>(helcode_g(pro));
			case nf:return R5g_nf_Ptr_eval<T>(helcode_g(pro));
			}; break;
		}
		case 6: {
			switch (cs) {
			case glue: case LT: case leading_color: return R6g_Ptr_eval<T>(helcode_g(pro));
			case nf:return R6g_Ptr_nf_eval<T>(helcode_g(pro));
			}; break;
		}

		case 7: {
			switch (cs) {
			case glue: return R7g_Ptr_eval<T>(helcode_g(pro));
			}; break;

		}
		case 22:
			switch (cs) {
			case nfLT: return R2q2g_nfLT_Ptr_eval<T>(helcode_2q(pro));
			case LT: return R2q2g_LT_Ptr_eval<T>(helcode_2q(pro));
			case nf: return R2g2q_nf_Ptr_eval<T>(helcode_2q(pro));
// 			case LT: return R2g2q_L_Ptr_eval<T>(pro);
// 			case RT: return R2g2q_SLC_Ptr_eval<T>(pro);
			case leading_color: return R_Ptr_eval<T>(pro,right_direction(pro,cs));
				}
			break;
		case 23:
			switch (cs) {
			case LT: return R2q3g_L_Ptr_eval<T>(helcode_2q(pro));
			case RT: return R2q3g_SLC_Ptr_eval<T>(helcode_2q(pro));
			case nf: return R2q3g_nf_Ptr_eval<T>(helcode_2q(pro));
			case nfLT: return R2q3g_nf_Ptr_eval<T>(helcode_2q(pro));
			case leading_color: return R_Ptr_eval<T>(pro,right_direction(pro,cs));
			}
			break;
		case 24:
			switch (cs) {
			case LT: return R2q4g_L_Ptr_eval<T>(helcode_2q(pro));
			case RT: return R2q4g_SLC_Ptr_eval<T>(helcode_2q(pro));
			case nf: return R2q4g_nf_Ptr_eval<T>(helcode_2q(pro));
			case leading_color: return R_Ptr_eval<T>(pro,right_direction(pro,cs));
			}
			break;
		//case 40:
		//	switch (cs) {
		//		case leading_color: return R4q_L_Ptr_eval<T>(helcode_4q(pro));
		//	}
		case 40:
			switch (cs) {
			case LLT: return R2q2Q_LLT_Ptr_eval<T>(helcode_2q2Q(pro));
			case LRT: return R2q2Q_LRT_Ptr_eval<T>(helcode_2q2Q(pro));
			case nfLLT: return R2q2Q_nfLLT_Ptr_eval<T>(helcode_2q2Q(pro));
			case nfLRT: return R2q2Q_nfLRT_Ptr_eval<T>(helcode_2q2Q(pro));
// 			case leading_color: return R2q2Q_L_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
// 			case nf: return R2q2Q_nf_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
			}
			break;
		case 41:
			switch (cs) {
			case LRT: return R2q2Q1g_LRT_Ptr_eval<T>(helcode_2q2Q(pro));
			case LLT: return R2q2Q1g_LLT_Ptr_eval<T>(helcode_2q2Q(pro));
			case nfLLT: return R2q2Q1g_nfLLT_Ptr_eval<T>(helcode_2q2Q(pro));
// 			case leading_color: return R2q2Q1g_PentP_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
			//case nf: return R2q2Q1g_nf_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
			}
			break;
		case 220:
			switch (cs) {
			case leading_color: return R2q2l_L_Ptr_eval<T>(helcode_2q2l(pro));
			case nf: return R2q2l_nf_Ptr_eval<T>(helcode_2q2l(pro));  // all zero
			}
			break;
		case 221:
			switch (cs) {
			case leading_color: return R2q1g2l_L_Ptr_eval<T>(helcode_2q2l(pro));
			case sub_leading_color: return R2q1g2l_SLC_Ptr_eval<T>(helcode_2q2l(pro));
//			case nf: return R2q1g2l_nf_Ptr_eval<T>(helcode_2q2l(pro)); // all zero
			case AX: return R2q1g2l_AX_Ptr_eval<T>(helcode_2q2l(pro)); // all zero
			}
			break;
		case 222:
			switch (cs) {
			case leading_color: return R2q2g2l_L_Ptr_eval<T>(helcode_2q2l(pro));
			case sub_leading_color: return R2q2g2l_SLC_Ptr_eval<T>(helcode_2q2l(pro));
			case nf: return R2q2g2l_nf_Ptr_eval<T>(helcode_2q2l(pro));
			case nf_top: return R2q2g2l_nf_top_Ptr_eval<T>(helcode_2q2l(pro));
			case VECT: return R2q2g2l_VECT_Ptr_eval<T>(helcode_2q2l(pro));
			case AX: return R2q2g2l_AX_Ptr_eval<T>(helcode_2q2l(pro));
			case AXSL: return R2q2g2l_AXSL_Ptr_eval<T>(helcode_2q2l(pro));
			}
			break;
		case 223:
			switch (cs) {
            case leading_color: case LT: return	R2q3g2l_L_Ptr_eval<T>(helcode_2q2l(pro));
            case nf: case sub_leading_color:
				if (BH::settings::rational_settings::s_skip_sub_leading_color){
					return	R2q3g2l_place_holder_Ptr_eval<T>(helcode_2q2l(pro));
				} else {
					return 0;
				}
			}
			break;
		case 224:
			switch (cs) {
			//case leading_color: return R2q4g2l_place_holder_Ptr_eval<T>(helcode_2q2l(pro));
//			case sub_leading_color: return R2q4g2l_place_holder_Ptr_eval<T>(helcode_2q2l(pro));
//			case nf: return R2q4g2l_place_holder_Ptr_eval<T>(helcode_2q2l(pro));
			}
			break;
		case 240:
			switch (cs) {
			case leading_color: return R2q2Q2l_L_Ptr_eval<T>(helcode_2q2l2Q(pro));
			case sub_leading_color: return R2q2Q2l_sl_Ptr_eval<T>(helcode_2q2l2Q(pro));
			case nf: return R2q2Q2l_nf_Ptr_eval<T>(helcode_2q2l2Q(pro));
			case nf_top: return R2q2Q2l_nf_top_Ptr_eval<T>(helcode_2q2l2Q(pro));
			case AX: return R2q2Q2l_AX_Ptr_eval<T>(helcode_2q2l2Q(pro));
			}
			break;
		case 241:
			switch (cs) {
			//case leading_color: return R2q2Q1g2l_place_holder_Ptr_eval<T>(helcode_2q2l2Q(pro));
//			case sub_leading_color: return R2q2Q1g2l_place_holder_Ptr_eval<T>(helcode_2q2l2Q(pro));
//			case nf: return R2q2Q1g2l_place_holder_Ptr_eval<T>(helcode_2q2l2Q(pro));
//			case AX: return R2q2Q1g2l_place_holder_Ptr_eval<T>(helcode_2q2l2Q(pro));
			}
			break;
		case 100021:
			switch (cs) {
			//case leading_color: return R_Ptr_eval<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
			case leading_color: return R2q1g1y_L_Ptr_eval<T>(helcode_2q1y(pro));
			case sub_leading_color: return R2q1g1y_SLC_Ptr_eval<T>(helcode_2q1y(pro));
				}
			break;

		case 100022:
		  	switch (cs) {case LT: case RT: case nfLT: case nfRT: return 0;}
			complex<T> (*RPtr)(const eval_param<T>&, const mass_param_coll&);
			switch (cs) {
				case leading_color: RPtr= R2q2g1y_L_Ptr_eval<T>( helcode_2q1y(pro)); break;
				case sub_leading_color: RPtr = R2q2g1y_SLC_Ptr_eval<T>( helcode_2q1y(pro)); break;
				case nf: RPtr = R2q2g1y_nf_Ptr_eval<T>( helcode_2q1y(pro)); break;
			}
			
			if (RPtr != 0) { return RPtr; } else return R_Ptr_eval<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
			break;
		case 100023:
			switch (cs) {
			case leading_color: return R_Ptr_eval<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
				}
			break;
		case 100040:
			switch (cs) {
				case leading_color: return R2q2Q1y_L_Ptr_eval<T>(helcode_2q2Q(pro));
				case sub_leading_color: return R2q2Q1y_sl_Ptr_eval<T>(helcode_2q2Q(pro));
				case nf: return R2q2Q1y_nf_Ptr_eval<T>(helcode_2q2Q(pro));
			}
			break;
		case 2000002: case 2000003: case 2000004:
			return R_Ptr_eval<T>(replace_gluino_with_quark(pro),cs);
		case 2000020: case 2000021:
			return R_Ptr_eval<T>(replace_gluino_with_quark(pro),cs);
		case 2100020:
			return R_Ptr_eval<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
		case 2000200: case 2000201: case 2000202:
			return R_Ptr_eval<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
		case 2000220:
			return R_Ptr_eval<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
			break;
		case 2000221:
			switch (cs) {
				 case leading_color: return R2q2Q1g2l_L_Ptr_eval<T>(helcode_2q2l2G(pro));
			}
			break;


		case 100000003: {
			switch (cs) {
			case glue: return R3g1ph_G_Ptr_eval<T>(helcode_Ng1ph(pro));
			}; break;
		}
		case 100000004: {
			switch (cs) {
			case glue: return R4g1ph_G_Ptr_eval<T>(helcode_Ng1ph(pro));
			}; break;
		}

		case 100000022: {  // higgs 2 quarks 2 gluons
			switch (cs) {
			case LT: return R2q2g1ph_LT_Ptr_eval<T>(helcode_phi_1q(pro));
			case RT: return R2q2g1ph_RT_Ptr_eval<T>(helcode_phi_1q(pro));
			case nfLT: return R2q2g1ph_nfLT_Ptr_eval<T>(helcode_phi_1q(pro));
			}; break;
		}

		case 100000040: {  // higgs 4 quarks
			switch (cs) {
			case leading_color:  return R2q2Q1ph_lc_Ptr_eval<T>(helcode_phi_2q2Q(pro));
			case sub_leading_color: return R2q2Q1ph_slc_Ptr_eval<T>(helcode_phi_2q2Q(pro));
			case nf:  return R2q2Q1ph_nf_Ptr_eval<T>(helcode_phi_2q2Q(pro));
			}; break;
		}




	}
	return 0;
}


complex<R> Known_Rec_Rational_offset::eval(const eval_param<R>& ep){
	eval_param<R> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%ep.size()));
	}
	return (*_eval_C_Ptr_eval)(rotated_ep,*_masses);
}
complex<RHP> Known_Rec_Rational_offset::eval(const eval_param<RHP>& ep){
	eval_param<RHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%ep.size()));
	}
	return (*_eval_CHP_Ptr_eval)(rotated_ep,*_masses);
}
complex<RVHP> Known_Rec_Rational_offset::eval(const eval_param<RVHP>& ep){
	eval_param<RVHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%ep.size()));
	}
	return (*_eval_CVHP_Ptr_eval)(rotated_ep,*_masses);
}

#if BH_USE_GMP
complex<RGMP> Known_Rec_Rational_offset::eval(const eval_param<RGMP>& ep){
	eval_param<RGMP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%ep.size()));
	}
	return (*_eval_CGMP_Ptr_eval)(rotated_ep,*_masses);
}
#endif

complex<R> Known_Rec_Rational_permutation::eval(const eval_param<R>& ep){
	eval_param<R> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_C_Ptr_eval)(rotated_ep,*_masses);

}
complex<RHP> Known_Rec_Rational_permutation::eval(const eval_param<RHP>& ep){
	eval_param<RHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CHP_Ptr_eval)(rotated_ep,*_masses);

}
complex<RVHP> Known_Rec_Rational_permutation::eval(const eval_param<RVHP>& ep){
	eval_param<RVHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CVHP_Ptr_eval)(rotated_ep,*_masses);

}

#if BH_USE_GMP
complex<RGMP> Known_Rec_Rational_permutation::eval(const eval_param<RGMP>& ep){
	eval_param<RGMP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CGMP_Ptr_eval)(rotated_ep,*_masses);

}
#endif

template complex<R> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<R>&,const mass_param_coll&);
template complex<RHP> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<RGMP>&,const mass_param_coll&);
#endif
}
