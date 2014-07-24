/*
 * rec_rational_factory.cpp
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#include "rec_rational_factory.h"
#include "process_utils.h"
#ifndef BH_PUBLIC
#include "rec_rational.h"
#endif
#include "known_rational.h"
#include "BH_typedefs.h"
#include "settings.h"
#include "ratext/ratext_part_worker.h"
#ifndef BH_PUBLIC 
#include "ratext/ratext_part_normal.h"
#endif
using namespace std;

namespace BH {
#ifdef USE_MC_RAT
template <class T> std::complex<T> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<T>& ,const  std::vector<int>&) ;
#else
template <class T> std::complex<T> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&,const mass_param_coll&) ;
#endif


#ifndef BH_PUBLIC

Rational_base* Rec_Rational_factory::new_rational(const process& pro,color_structure cs){

//	 if Rec_Rational_factory::s_set_all_zero is true , all rational terms will be recognized as zero
	 if (BH::settings::rational_settings::s_set_all_zero){
		 return new Known_Rec_Rational(pro,cs);
	 };

	 // here is the place to put special hacks, like for the photon amplitudes

	switch(pro.pcode()){
		case 100021:
			switch (cs) {
			case leading_color: case nf: return new_rational(replace_photon_with_gluon(pro),right_direction(pro,cs));
				}
			break;
		case 100022:
			if( R_Ptr_eval<R>(pro,cs) != 0){
				return new Known_Rec_Rational(pro,cs);
			}
			switch (cs) {
			case leading_color: case nf: return new_rational(replace_photon_with_gluon(pro),right_direction(pro,cs));
				}
			break;
		case 100023:
			switch (cs) {
			case leading_color: case nf: return new_rational(replace_photon_with_gluon(pro),right_direction(pro,cs));
				}
			break;
		case 22:case 23: case 24:case 25:case 26:case 27: {
			vector<particle_ID> pp;size_t qb_pos;
			for (int i=1;i<=pro.n();i++){if (pro.p(i).is_a(quark) && pro.p(i).is_anti_particle()) qb_pos=i; }
			for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((qb_pos+(i-2))%pro.n()+1)) ; }
			process new_pro(pp);
			if( R_Ptr_eval<R>(new_pro,cs) == 0){
				return new Unknown_Rec_Rational(pro,cs);
			} else {
				return new Known_Rec_Rational_offset(pro,cs,qb_pos);
			}
			 break;}

		case 220: case 221: case 222: case 223: {
			switch (cs) {
			case leading_color: {
				vector<particle_ID> pp;size_t q_pos;
				// for the 2q2l amplitudes, the q comes before the antiquark in the canonical ordering
				for (int i=1;i<=pro.n();i++){if (pro.p(i).is_a(quark) && !(pro.p(i).is_anti_particle())) q_pos=i; }
				for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((q_pos+(i-2))%pro.n()+1)) ; }
				process new_pro(pp);
					if( R_Ptr_eval<R>(new_pro,leading_color) == 0){
						return new Unknown_Rec_Rational(pro,leading_color);
					} else {
						return new Known_Rec_Rational_offset(pro,leading_color,q_pos);
					}
					break;}
			case sub_leading_color: {
				vector<int> new_ind;
				process new_pro=order_2q2e(pro,new_ind);
					if( R_Ptr_eval<R>(new_pro,sub_leading_color) == 0){
						return new Unknown_Rec_Rational(pro,sub_leading_color);
					} else {
						return new Known_Rec_Rational_permutation(pro,sub_leading_color,new_ind);
					}
					break;}
			}
			break;
			}
	}
	if( R_Ptr_eval<R>(pro,cs) != 0){
//		_MESSAGE("USING KNOWN RATIONAL");
		return new Known_Rec_Rational(pro,cs);
	}
	else {
		for (int offset=1;offset<= pro.n();offset++) {
			vector<particle_ID> pp;
			for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((offset+(i-2))%pro.n()+1)) ; }
			process new_pro(pp);
			if( R_Ptr_eval<R>(new_pro,cs) != 0){
//				_MESSAGE("USING KNOWN RATIONAL");

				return new Known_Rec_Rational_offset(pro,cs,offset);
			}
		}
	}
	return new Unknown_Rec_Rational(pro,cs);

}


#endif

}
