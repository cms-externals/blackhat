/*
 * rational_factory.cpp
 *
 *  Created on: 8 Jul 2009
 *      Author: daniel
 */

#include "rational_factory.h"
#include "process_utils.h"
#include "BH_typedefs.h"
#include "settings.h"
#ifndef BH_PUBLIC
#include "rec_rational.h"
#endif
#include "known_rational.h"
#include "ratext/ratext_part_worker.h"
#ifndef BH_PUBLIC
#include "ratext/ratext_part_normal.h"
#endif
#include "rec_rational_factory.h"
#include "BH_error.h"
//#include "ratext/rat_mix.h"

namespace BH {

#ifdef USE_MC_RAT
template <class T> std::complex<T> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<T>& ,const  std::vector<int>&) ;
#else
template <class T> std::complex<T> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&,const mass_param_coll&) ;
#endif
Rational_base* Known_Rational_factory::new_rational(const process& pro,color_structure cs){
	if (settings::rational_settings::s_use_known_formulae_in_ratext){
		if( R_Ptr_eval<R>(pro,cs) != 0){
			//_MESSAGE2("Using known rational amplitude ",pro);
			return new Known_Rec_Rational(pro,cs);
		}
		else {
			for (int offset=1;offset<= pro.n();offset++) {
				std::vector<particle_ID> pp;
				for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((offset+(i-2))%pro.n()+1)) ; }
				process new_pro(pp);
				if( R_Ptr_eval<R>(pro,cs) != 0){
					//_MESSAGE2("Known rat process ",pro);
					return new Known_Rec_Rational_offset(pro,cs,offset);
				}
			}
		}
	}
	return 0;
}



ratext::worker_rational_factory explicit_WRF;
#ifndef BH_PUBLIC
Rec_Rational_factory explicit_RRF;
ratext::normal_ratext_factory explicit_NRF;
#endif
//ratmix::normal_rat_mix_factory explicit_NRM;

Known_Rational_factory explicit_KRF;
Rational_factory<Rational_base>* Known_Rational_factory::s_default_KRF=&explicit_KRF;
#ifndef BH_PUBLIC
Rational_factory<Rational_base>* Rec_Rational_factory::s_default_RRF=&explicit_RRF;
#endif

template <> Rational_factory<Rational_base>* Rational_factory<Rational_base>::default_rational_factory(){
	switch ( settings::general::s_rat_type ){
	case  settings::general::ratext:
#ifndef BH_PUBLIC
		return &explicit_NRF;
#else
		throw BHerror("Not possible for public version");
#endif
	case  settings::general::recursive:
#ifndef BH_PUBLIC
		return &explicit_RRF;
#else
		throw BHerror("Not possible for public version");
#endif
	case  settings::general::ratext_worker:
		return &explicit_WRF;
//	case  settings::general::ratmix:
//		return &explicit_NRM;
	default : {
		_WARNING("Missing case in Rational_factory<Rational_base>::default_rational_factory()");
		throw BHerror("Missing case in Rational_factory");
	}
	}
}

Rational_factory<Rational_base>* ratext::worker_rational_factory::s_default_WRF=&explicit_WRF;


}
