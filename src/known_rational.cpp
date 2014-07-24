/*
 * known_rational.cpp
 *
 *  Created on: 6 Jul 2009
 *      Author: daniel
 */

#include "known_rational.h"
#include "tree_amp.h"
#include "settings.h"
#include "process_utils.h"
#include "helcode.h"

using namespace std;

namespace BH {

bool Rec_Rational_is_zero(const process& pro,color_structure cs){
	size_t nbr_particles=0;
	size_t nbr_anti_particles=0;
	size_t nbr_pos_hel=0;
	size_t nbr_neg_hel=0;
	size_t nbr_leptons=0;
	size_t nbr_photons=0;
	size_t nbr_quarks=0;
	size_t nbr_gluons=0;
	size_t nbr_higgs=0;

	for(int i=1;i<=pro.n();i++)
	{
		if ( (pro.p(i).is_a(quark)) ){
			++nbr_quarks;
			if (pro.p(i).is_anti_particle()) {
				++nbr_anti_particles;
			} else {
				++nbr_particles;
			}
			if (pro.p(i).helicity() == 1) {
				++nbr_pos_hel;
			} else {
				++nbr_neg_hel;
			}
		}

		if ( (pro.p(i).is_a(lepton)) ){
			++nbr_leptons;
		}
		if ( (pro.p(i).is_a(photon)) ){
			++nbr_photons;
		}
		if ( (pro.p(i).is_a(gluon)) ){
			++nbr_gluons;
		}
		if ( (pro.p(i).is_a(higgs)) ){
			++nbr_higgs;
	}
	}


	if (nbr_anti_particles != nbr_particles ) {
		//			_MESSAGE2(pro," is zero because of flavor");
		return true;

	}
	if (nbr_neg_hel != nbr_pos_hel ) {
		//			_MESSAGE2(pro," is zero because of helicity");
		return true;
	}

	//
//	if (nbr_leptons > 0 && nbr_photons > 0 && (nbr_anti_particles+nbr_particles)> 0) return true;
// no photons if there is no leptons to couple to
//	if (nbr_leptons == 0 && nbr_photons > 0 ) return true;

	// no photons if there is no leptons or no quraks to couple to
	// wrong for ng my
	//if (nbr_leptons == 0 && nbr_quarks == 0 && nbr_photons > 0 ) return true;
switch (cs){
	case RT: if (nbr_gluons == pro.n()) return true;
	case sub_leading_color :
		switch ( pro.pcode()) {
			case 220: return true;
		}
}

	return false;
}



  template <class T> complex<T> dummy2(momentum_configuration<T>& mc, const vector<int>& ind){_WARNING("should not be in dummy function...");throw BHerror("Wrong usage");return complex<T>(0,0);}
  template <class T> complex<T> (*dummy(int hc))(momentum_configuration<T>& ,const  vector<int>&) {return &dummy2;}


template <class T> complex<T> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&,const mass_param_coll&); // Located in known_rational_ep.cpp

#ifdef USE_MC_RAT

template <class T> complex<T> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<T>& ,const  vector<int>&){
	if ( Rec_Rational_is_zero(pro,cs) ) {
			return &ZeroF;
	};

	if ( BH::settings::rational_settings::s_set_all_zero ) {
			return &ZeroF;
	};
	switch(pro.pcode()){
		case odd_nbr_q:
		case odd_nbr_q2:
		case odd_nbr_l:
		case odd_nbr_l1:
		case odd_nbr_l2:

		case 200020: /* two photons two quarks should not appear */
			return &ZeroF;
		break;


		case 3:
		// At one-loop the 3-pt vertex is either zero or doen't factorise in which case we
		// want to drop it anyway and calculate an Inf term
			return &ZeroF;
		break;

		case 4: {
			switch (cs) {
				case glue: case LT: case leading_color: return R4g_Ptr<T>(helcode_g(pro));
				case nf:   return R4g_nf_Ptr<T>(helcode_g(pro));
				case RT: return R4g_SLC_Ptr<T>(helcode_g(pro));
			}; break;
		}
		case 5: {
			switch (cs) {
			case glue: case LT: case leading_color: return R5g_Ptr<T>(helcode_g(pro));
			case nf:return R5g_Ptr_nf<T>(helcode_g(pro));
			}; break;
		}
		case 6: {
			switch (cs) {
			case glue: case LT: case leading_color: return R6g_Ptr<T>(helcode_g(pro));
			case nf:return R6g_Ptr_nf<T>(helcode_g(pro));
			}; break;
		}

		case 7: {
			switch (cs) {
			case glue: return R7g_Ptr<T>(helcode_g(pro));
			}; break;

		}
		case 22:
			switch (cs) {
			case nf: return R2g2q_nf_Ptr<T>(helcode_2q(pro));
			case nfLT: return R2g2q_nf_Ptr<T>(helcode_2q(pro));
			case LT: return R2g2q_L_Ptr<T>(pro);
			case RT: return R2g2q_SLC_Ptr<T>(pro);
			case leading_color: return R_Ptr<T>(pro,right_direction(pro,cs));
				}
			break;
		case 23:
			switch (cs) {
			case LT: return R2q3g_L_Ptr<T>(helcode_2q(pro));
			case RT: return R2q3g_SLC_Ptr<T>(helcode_2q(pro));
			case nf: return R2q3g_nf_Ptr<T>(helcode_2q(pro));
			case nfLT: return R2q3g_nf_Ptr<T>(helcode_2q(pro));
			case leading_color: return R_Ptr<T>(pro,right_direction(pro,cs));
			}
			break;
		case 24:
			switch (cs) {
			case LT: return R2q4g_L_Ptr<T>(helcode_2q(pro));
			case RT: return R2q4g_SLC_Ptr<T>(helcode_2q(pro));
			case nf: return R2q4g_nf_Ptr<T>(helcode_2q(pro));
			case leading_color: return R_Ptr<T>(pro,right_direction(pro,cs));
			}
			break;
		//case 40:
		//	switch (cs) {
		//		case leading_color: return R4q_L_Ptr<T>(helcode_4q(pro));
		//	}
		case 40:
			switch (cs) {
#ifndef BH_PUBLIC
			case LLT: return R2q2Q_L_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			case LRT: return R2q2Q_L_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			case nfLLT: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			case nfLRT: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));

			case leading_color: return R2q2Q_L_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			case nf: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
#else
			/*the fixed process should be written in the files, no need to fix again*/
			case LLT: return R2q2Q_L_Ptr<T>(helcode_2q2Q(pro));
			case LRT: return R2q2Q_L_Ptr<T>(helcode_2q2Q(pro));
			case nfLLT: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(pro));
			case nfLRT: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(pro));

			case leading_color: return R2q2Q_L_Ptr<T>(helcode_2q2Q(pro));
			case nf: return R2q2Q_nf_Ptr<T>(helcode_2q2Q(pro));
#endif
			}
			break;
		case 41:
			switch (cs) {
	                case LLT: return dummy<T>(0);  //hack. this is just so as not to return a zero pointer.
	                case nfLLT: return dummy<T>(0);  //hack. this is just so as not to return a zero pointer.
	                case LRT: return dummy<T>(0);  //hack. this is just so as not to return a zero pointer.
// 			case leading_color: return R2q2Q1g_PentP_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			//case nf: return R2q2Q1g_nf_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
			}
			break;
		case 220:
			switch (cs) {
			case leading_color: return R2q2l_L_Ptr<T>(helcode_2q2l(pro));
			case nf: return R2q2l_nf_Ptr<T>(helcode_2q2l(pro));  // all zero
			}
			break;
		case 221:
			switch (cs) {
			case leading_color: return R2q1g2l_L_Ptr<T>(helcode_2q2l(pro));
			case sub_leading_color: return R2q1g2l_SLC_Ptr<T>(helcode_2q2l(pro));
//			case nf: return R2q1g2l_nf_Ptr<T>(helcode_2q2l(pro)); // all zero
			case AX: return R2q1g2l_AX_Ptr<T>(helcode_2q2l(pro)); // all zero
			}
			break;
		case 222:
			switch (cs) {
			case leading_color: return R2q2g2l_L_Ptr<T>(helcode_2q2l(pro));
			case sub_leading_color: return R2q2g2l_SLC_Ptr<T>(helcode_2q2l(pro));
			case nf: return R2q2g2l_nf_Ptr<T>(helcode_2q2l(pro));
			case nf_top: return R2q2g2l_nf_top_Ptr<T>(helcode_2q2l(pro));
			case VECT: return R2q2g2l_VECT_Ptr<T>(helcode_2q2l(pro));
			case AX: return R2q2g2l_AX_Ptr<T>(helcode_2q2l(pro));
			case AXSL: return R2q2g2l_AXSL_Ptr<T>(helcode_2q2l(pro));
			}
			break;
		case 223:
			switch (cs) {
                case leading_color: case LT: return	R2q3g2l_L_Ptr<T>(helcode_2q2l(pro));
			case nf: case sub_leading_color:
				if (BH::settings::rational_settings::s_skip_sub_leading_color){
					return	R2q3g2l_place_holder_Ptr<T>(helcode_2q2l(pro));
				} else {
					return 0;
				}
			}
			break;
		case 224:
			switch (cs) {
			//case leading_color: return R2q4g2l_place_holder_Ptr<T>(helcode_2q2l(pro));
//			case sub_leading_color: return R2q4g2l_place_holder_Ptr<T>(helcode_2q2l(pro));
//			case nf: return R2q4g2l_place_holder_Ptr<T>(helcode_2q2l(pro));
			}
			break;
		case 240:
			switch (cs) {
			case leading_color: return R2q2Q2l_L_Ptr<T>(helcode_2q2l2Q(pro));
			case sub_leading_color: return R2q2Q2l_sl_Ptr<T>(helcode_2q2l2Q(pro));
			case nf: return R2q2Q2l_nf_Ptr<T>(helcode_2q2l2Q(pro));
			case nf_top: return R2q2Q2l_nf_top_Ptr<T>(helcode_2q2l2Q(pro));
			case AX: return R2q2Q2l_AX_Ptr<T>(helcode_2q2l2Q(pro));
			}
			break;
		case 241:
			switch (cs) {
			//case leading_color: return R2q2Q1g2l_place_holder_Ptr<T>(helcode_2q2l2Q(pro));
//			case sub_leading_color: return R2q2Q1g2l_place_holder_Ptr<T>(helcode_2q2l2Q(pro));
//			case nf: return R2q2Q1g2l_place_holder_Ptr<T>(helcode_2q2l2Q(pro));
//			case AX: return R2q2Q1g2l_place_holder_Ptr<T>(helcode_2q2l2Q(pro));
			}
			break;
		case 100021:
			switch (cs) {
			//case leading_color: return R_Ptr<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
			case leading_color: return R2q1g1y_L_Ptr<T>(helcode_2q1y(pro));
			case sub_leading_color: return R2q1g1y_SLC_Ptr<T>(helcode_2q1y(pro));
				}
			break;

		case 100022:
		  switch (cs) {case LT: case RT: case nfLT: case nfRT: return 0;}
			complex<T> (*RPtr)(momentum_configuration<T>& ,const  vector<int>&);
			switch (cs) {
				case leading_color: RPtr= R2q2g1y_L_Ptr<T>( helcode_2q1y(pro)); break;
				case sub_leading_color: RPtr = R2q2g1y_SLC_Ptr<T>( helcode_2q1y(pro)); break;
				case nf: RPtr = R2q2g1y_nf_Ptr<T>( helcode_2q1y(pro)); break;
			}
			if (RPtr != 0) { return RPtr; } else return R_Ptr<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
			break;
		case 100023:
			switch (cs) {
			case leading_color: return R_Ptr<T>(replace_photon_with_gluon(pro),right_direction(pro,cs));
				}
			break;
		case 100040:	
			switch (cs) {
				case leading_color: return R2q2Q1y_L_Ptr<T>(helcode_2q2Q(pro));
				case sub_leading_color: return R2q2Q1y_sl_Ptr<T>(helcode_2q2Q(pro));
				case nf: return R2q2Q1y_nf_Ptr<T>(helcode_2q2Q(pro));
			}
			break;
		case 2000002: case 2000003: case 2000004:
			return R_Ptr<T>(replace_gluino_with_quark(pro),cs);
		case 2000020: case 2000021: case 2000022:
			return R_Ptr<T>(replace_gluino_with_quark(pro),cs);
		case 2100020:
			return R_Ptr<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
		case 2000200: case 2000201: case 2000202:
			return R_Ptr<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
		case 2000220:
			return R_Ptr<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1),cs);
			break;
		case 2000221:
			switch (cs) {
				 case leading_color: return R2q2Q1g2l_L_Ptr<T>(helcode_2q2l2G(pro));
			}
			break;

	case 100000003: {  //phi plus 3 gluons
			switch (cs) {
			case glue: return dummy<T>(helcode_g(pro));  //hack. this is just so as not to return a zero pointer.
			}; break;
		}
	case 100000004: {  //phi plus 4 gluons
			switch (cs) {
			case glue: return dummy<T>(helcode_g(pro));  //hack. this is just so as not to return a zero pointer.
			}; break;
		}

	case 100000022: {  //phi plus 2 quarks 2 gluons
			switch (cs) {
// 			  _PRINT_("inside 100000022 switch in known_rational",cs):
			case LT: case RT: case nfLT: return dummy<T>(helcode_phi_1q(pro));  //hack. this is just so as not to return a zero pointer.
			}; break;
		}

	case 100000040: {  // higgs 4 quarks
	  switch (cs) {
	  case leading_color: case sub_leading_color: case nf: return dummy<T>(helcode_phi_2q2Q(pro));  //hack. this is just so as not to return a zero pointer.
	  }; break;
	}


	}
	return 0;
}
#endif


Known_Rec_Rational_offset::Known_Rec_Rational_offset(const process& pro,color_structure cs,size_t offset): Known_Rec_Rational_base(pro), _offset(offset)
{
	// we use the implicit conversion pointer 0 -> false
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((_offset+(i-2))%pro.n()+1)) ; }
	process new_pro(pp);
//_MESSAGE("using Known_Rec_Rational_offset");
#ifdef USE_MC_RAT
	if ( !(_eval_C_ptr = R_Ptr<R>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_ptr = R_Ptr<RHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_ptr = R_Ptr<RVHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}
#endif
	//Not implemented yet
	if ( !(_eval_C_Ptr_eval = R_Ptr_eval<R>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_Ptr_eval = R_Ptr_eval<RHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_Ptr_eval = R_Ptr_eval<RVHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational_offset)",pro);throw BHerror("Unknown known amplitude");}

}

#ifdef USE_MC_RAT


complex<R> Known_Rec_Rational_offset::eval(momentum_configuration<R>& mc, const vector<int>& ind ){
	vector<int> new_ind; rotate_copy(ind.begin(),ind.begin()+(_offset-1),ind.end(),back_inserter(new_ind)); 
	if(settings::general::s_use_ep_only){
		eval_param<R> ep(mc,new_ind); 
		return (*_eval_C_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_C_ptr)(mc,new_ind);
	}
}
complex<RHP> Known_Rec_Rational_offset::eval(momentum_configuration<RHP>& mc, const vector<int>& ind ){
	vector<int> new_ind; rotate_copy(ind.begin(),ind.begin()+(_offset-1),ind.end(),back_inserter(new_ind));
	if(settings::general::s_use_ep_only){
		eval_param<RHP> ep(mc,new_ind); 
		return (*_eval_CHP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CHP_ptr)(mc,new_ind);
	}
}
complex<RVHP> Known_Rec_Rational_offset::eval(momentum_configuration<RVHP>& mc, const vector<int>& ind ){
	vector<int> new_ind; rotate_copy(ind.begin(),ind.begin()+(_offset-1),ind.end(),back_inserter(new_ind));
	if(settings::general::s_use_ep_only){
		eval_param<RVHP> ep(mc,new_ind); 
		return (*_eval_CVHP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CVHP_ptr)(mc,new_ind);
	}
}
#if BH_USE_GMP
complex<RGMP> Known_Rec_Rational_offset::eval(momentum_configuration<RGMP>& mc, const vector<int>& ind ){
	vector<int> new_ind; rotate_copy(ind.begin(),ind.begin()+(_offset-1),ind.end(),back_inserter(new_ind));
	if(settings::general::s_use_ep_only){
		eval_param<RGMP> ep(mc,new_ind);
		return (*_eval_CGMP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CGMP_ptr)(mc,new_ind);
	}
}
#endif

#endif
Known_Rec_Rational_permutation::Known_Rec_Rational_permutation(const process& pro,color_structure cs,const vector<int>& per): Known_Rec_Rational_base(pro), _perm_ind(per)
{
	// we use the implicit conversion pointer 0 -> false
	vector<particle_ID> pp;
	for (int i=1;i<=pro.n();i++){
		pp.push_back(pro.p(_perm_ind[i-1])) ;
	}
	process new_pro(pp);
//_MESSAGE("using Known_Rec_Rational_offset");
#ifdef USE_MC_RAT

	if ( !(_eval_C_ptr = R_Ptr<R>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_ptr = R_Ptr<RHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_ptr = R_Ptr<RVHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}
#endif
	//Not implemented yet
	if ( !(_eval_C_Ptr_eval = R_Ptr_eval<R>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_Ptr_eval = R_Ptr_eval<RHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_Ptr_eval = R_Ptr_eval<RVHP>(new_pro,cs)) ) {_WARNING2("Unknown known amplitude for (rat_permutation)",pro);throw BHerror("Unknown known amplitude");}

}


#ifdef USE_MC_RAT

complex<R> Known_Rec_Rational_permutation::eval(momentum_configuration<R>& mc, const vector<int>& ind ){
	
	vector<int> new_ind;
	for (int i=0;i<ind.size();i++){
		new_ind.push_back(ind[_perm_ind[i]-1]) ;
	}	
	
	if(settings::general::s_use_ep_only){
		eval_param<R> ep(mc,new_ind); 
		return (*_eval_C_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_C_ptr)(mc,new_ind);
	}
}
complex<RHP> Known_Rec_Rational_permutation::eval(momentum_configuration<RHP>& mc, const vector<int>& ind ){
	vector<int> new_ind;
	for (int i=0;i<ind.size();i++){
		new_ind.push_back(ind[_perm_ind[i]-1]) ;
	}
	if(settings::general::s_use_ep_only){
		eval_param<RHP> ep(mc,new_ind); 
		return (*_eval_CHP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CHP_ptr)(mc,new_ind);
	}
}
complex<RVHP> Known_Rec_Rational_permutation::eval(momentum_configuration<RVHP>& mc, const vector<int>& ind ){
	vector<int> new_ind;
	for (int i=0;i<ind.size();i++){
		new_ind.push_back(ind[_perm_ind[i]-1]) ;
	}
	if(settings::general::s_use_ep_only){
		eval_param<RVHP> ep(mc,new_ind); 
		return (*_eval_CVHP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CVHP_ptr)(mc,new_ind);
	}
}
#if BH_USE_GMP
complex<RGMP> Known_Rec_Rational_permutation::eval(momentum_configuration<RGMP>& mc, const vector<int>& ind ){
	vector<int> new_ind;
	for (int i=0;i<ind.size();i++){
		new_ind.push_back(ind[_perm_ind[i]-1]) ;
	}
	if(settings::general::s_use_ep_only){
		eval_param<RGMP> ep(mc,new_ind);
		return (*_eval_CGMP_Ptr_eval)(ep,*_masses);
	}else{
		return (*_eval_CGMP_ptr)(mc,new_ind);
	}
}
#endif

#endif
Known_Rec_Rational::Known_Rec_Rational(const process& pro,color_structure cs): Known_Rec_Rational_base(pro)
{
	// we use the implicit conversion pointer 0 -> false#ifdef USE_MC_RAT
#ifdef USE_MC_RAT

	if ( !(_eval_C_ptr = R_Ptr<R>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_ptr = R_Ptr<RHP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_ptr = R_Ptr<RVHP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
#endif
	if ( !(_eval_C_Ptr_eval = R_Ptr_eval<R>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CHP_Ptr_eval = R_Ptr_eval<RHP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
	if ( !(_eval_CVHP_Ptr_eval = R_Ptr_eval<RVHP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}

#if BH_USE_GMP
#ifdef USE_MC_RAT
	if ( !(_eval_CGMP_ptr = R_Ptr<RGMP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
#endif
	if ( !(_eval_CGMP_Ptr_eval = R_Ptr_eval<RGMP>(pro,cs)) ) {_WARNING2("Unknown known amplitude for (rational)",pro);throw BHerror("Unknown known amplitude");}
#endif

}



bool Known_Rec_Rational::is_zero(){
switch( d_process.pcode()){
// At one-loop the 3-pt vertex is either zero or doen't factorize in which case we
// want to drop it anyway and calculate an Inf term
case 3:	return true;
default: return false;
}
}

#ifdef USE_MC_RAT

template complex<R> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<R>& ,const  vector<int>&);
template complex<RHP> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<RHP>& ,const  vector<int>&);
template complex<RVHP> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<RVHP>& ,const  vector<int>&);

#if BH_USE_GMP
template complex<RGMP> (*R_Ptr(const process& pro,color_structure cs))(momentum_configuration<RGMP>& ,const  vector<int>&);
#endif

#endif

}
