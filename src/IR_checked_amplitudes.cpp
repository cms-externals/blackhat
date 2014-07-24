//#include "OneLoopHelAmpl.h"
//#include "BH_A0.h"
#include "qd_suppl.h"
#include <typeinfo>
#include "IR_checked.h"
#ifndef BH_PUBLIC
#include "IR_checked_pole.h"
#endif
#include "IR_checked_cut_part.h"
#include "cached_integral.h"
#include <cassert>
#include "ratext/ratext_part.h"
#ifndef BH_PUBLIC
#include "rec_rational.h"
#endif
#include "rec_rational_factory.h"
#include "BH_debug.h"
#define _VERBOSE 0   // there are two levels of verbosity : 1 and 2 (more output)
#define _USE_SMALL_TEST 0 /* if 0 then no test of small distance between poles is performed */
using namespace std;
namespace BH {



void IR_fix_Rational(Rational_base* Rb){
#ifndef BH_PUBLIC
	Rec_Rational* RR=static_cast<Rec_Rational*>(Rb);
	Spurious_Rational* rat=RR->get_spurious_part();

	if ( rat!=0 ){
		for (int k=1;k<=rat->nbr_daughters();k++){
			Rec_BB* old_sp = rat->get_daughter(k);
			rat->set_daughter(k,new IR_checked_Spurious_Pole(dynamic_cast<Spurious_Pole*>(old_sp)));
			// don't delete the old daughters
			// delete old_sp;
		}
		for (int k=1; k<=rat->nbr_daughters(); k++ ) {
			computable<std::complex>* old_eval=rat->get_eval_daughter(k);
			rat->set_eval_daughter(k,new zero_checked_computable<std::complex>(rat->get_daughter(k)));
			delete old_eval;
		}
	}
#endif
}


void IR_checked_OLHA::construct(){
/* Nothing to do in the public version*/
#ifndef BH_PUBLIC
	Unknown_Rec_Rational* URR=dynamic_cast<Unknown_Rec_Rational*>(_rational_part);
	if ( URR != 0){
		IR_fix_Rational(_rational_part);
		for (int ii=1;ii<=URR->nbr_daughters();ii++){
			Rec_Pair* RP = dynamic_cast<Rec_Pair*>(URR->get_daughter(ii));
			if ( RP != 0 ){
				Rec_Rational* RR;
				RR= dynamic_cast<Rec_Rational*>(RP->left());
				if ( RR !=0 ) {
					IR_fix_Rational(RR);
//					_MESSAGE("found a rational term, fixing IR...");
				}
				RR = dynamic_cast<Rec_Rational*>(RP->right());
				if ( RR !=0 ) {
					IR_fix_Rational(RR);
//					_MESSAGE("found a rational term, fixing IR...");
				}
			}
		}
	}
#endif
}


//IR_checked_OLHA::IR_checked_OLHA(const process& p,const vector<particle_ID>& possible_props, cutD_factory* cf, option* op): One_Loop_Helicity_Amplitude<IRC_cut_part_factory,Rec_Rational_factory>(p,possible_props,cf,op), IR_checked(1e-5,1e-15,1e-30) {
//	construct();
//}

IR_checked_OLHA::IR_checked_OLHA(const process& p, color_structure cs): One_Loop_Helicity_Amplitude(p,cs,Rec_Rational_factory::default_rational_factory(),IR_checked_cut_part_factory::s_default_IR_checked_cut_part_factory), IR_checked(1e-5,1e-15,1e-30) {
	construct();
}

//IR_checked_OLHA::IR_checked_OLHA(const process& p, color_structure cs, Rational_factory<Rec_Rational>* RRF , cut_part_factory<IR_checked_Cut_Part>* CPF): One_Loop_Helicity_Amplitude<IRC_cut_part_factory,Rec_Rational_factory>(p,cs,RRF,CPF), IR_checked(1e-5,1e-15,1e-30) {
//	construct();
//}



void IR_checked_OLHA::set_littleIR_tolerances(R t_R,R t_RHP,R t_RVHP){
#ifndef BH_PUBLIC
	Rec_Rational* RR = dynamic_cast<Rec_Rational*>(_rational_part);
	assert(RR);
	Spurious_Rational* rat=RR->get_spurious_part();
	if (rat!=0) {
		for (int k=1;k<=rat->nbr_daughters();k++){
		dynamic_cast<IR_checked_Spurious_Pole*>(rat->get_daughter(k))->set_tolerances(t_R,t_RHP,t_RVHP);
		}
	}
#endif
}

void IR_checked_OLHA::set_bigIR_tolerances(R t_R,R t_RHP,R t_RVHP){
	IR_checked_Cut_Part* cut= dynamic_cast<IR_checked_Cut_Part*>(_cut_part);
		cut->set_tolerances(t_R,t_RHP,t_RVHP);
}


SeriesC<R> IR_checked_OLHA::eval(mom_conf& mc,const vector<int>& ind){
  BH_DEBUG_MESSAGE2("Computing amplitude ",d_process);
  C tree=get_tree(mc,ind);
  BH_DEBUG_MESSAGE2("tree ",tree);

	SeriesC<R> cut_result=_cut_part->get_value(mc,ind);
	C rat_result=_rational_part->get_value(mc,ind);
  BH_DEBUG_MESSAGE2("rat ",rat_result);
 BH_DEBUG_MESSAGE2("cut ",cut_result);
	R tree_abs=abs(tree);
	R cut_result_abs=abs(cut_result[0]);
	R rat_result_abs=abs(rat_result);
	
	bool test;
	if (cut_result_abs == 0.0 ){
	  if (rat_result_abs == 0.0){
          BH_DEBUG_MESSAGE("Confirm change with Daniel.");
          //old: test =!(tree_abs < _IR_tolerance_R);
          test=false; 
	  } else {
		test =( tree_abs/rat_result_abs < _IR_tolerance_R
				&& !(tree_abs < _IR_tolerance_R)
				);
	  }
	}
	else {
		if (rat_result_abs == 0.0){
			test=
					tree_abs/cut_result_abs < _IR_tolerance_R
					&&
					!(tree_abs < _IR_tolerance_R
			);
		} else {
			test= ( abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_R   // big cancellation
					||
					tree_abs/cut_result_abs < _IR_tolerance_R  || tree_abs/rat_result_abs < _IR_tolerance_R
					&&
					!(tree_abs < _IR_tolerance_R
					));
		}
}
	BH_DEBUG(

		 if (test){
		   BH_DEBUG_MESSAGE("Test passed");
		 } else {
		   BH_DEBUG_MESSAGE("Test not passed");
		 }
	);
	if ( test
	){
		multi_precision_reader* mpr=dynamic_cast<multi_precision_reader*>(&mc);
		if ( mpr == 0 ){

			BH_DEBUG_MESSAGE("using higher precision for non non multi_precision_reader");
//			throw BHerror("IR_checked_OLHA::eval called with an argument that is not a multi_precision_reader");


			SeriesC<RHP> HP_result;
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_HP mcHP=mc.extend<RHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_HP;
			new_mu_index_HP = DefineMu<RHP>(mcHP,RHP(mc.mom(_cut_part->get_mu_index<R>()).E().real()));

			_cut_part->set_mu_HP(new_mu_index_HP);
			HP_result=eval(mcHP,new_ind);
			BH_DEBUG_MESSAGE2("HP_result: ",HP_result);
			return to_double(HP_result);
		} else {
			SeriesC<RHP> HP_result=eval(mpr->mc_HP(),ind);
			return to_double(HP_result);
		}
	}
	
	// Collect the accuracy estimate
	double ratacc=_rational_part->get_accuracy();
	double cutacc=_cut_part->get_accuracy();
	ratacc>cutacc?_accuracy=cutacc:_accuracy=ratacc;
   
    // fill in conjugate value of amplitude
    _conjugate_rat_part=conj(rat_result);
    _conjugate_rat_part_HP=to_HP(conj(rat_result));
    _conjugate_rat_part_VHP=to_VHP(conj(rat_result));
    _conjugate_cut_part=_cut_part->get_conjugate_cut_part();
    _conjugate_cut_part_HP=_cut_part->get_conjugate_cut_part_HP();
    _conjugate_cut_part_VHP=_cut_part->get_conjugate_cut_part_VHP();

	return rat_result+cut_result;
}
SeriesC<RHP> IR_checked_OLHA::eval(mom_conf_HP& mc,const vector<int>& ind){
	CHP tree=get_tree(mc,ind);
	SeriesC<RHP> cut_result=_cut_part->get_value(mc,ind);
	CHP rat_result=_rational_part->get_value(mc,ind);

	RHP tree_abs=abs(tree);
	RHP cut_result_abs=abs(cut_result[0]);
	RHP rat_result_abs=abs(rat_result);
	BH_DEBUG_MESSAGE2("rat ",rat_result);
	BH_DEBUG_MESSAGE2("cut ",cut_result);

	bool test;

	if (cut_result_abs == 0.){
		if (rat_result_abs == 0.0) {
			test = !(tree_abs < _IR_tolerance_RHP);
		} else {
		test=
				tree_abs/rat_result_abs < _IR_tolerance_RHP
				&&
				!(tree_abs < _IR_tolerance_RHP);
		}
	} else {
		if (rat_result_abs == 0.0) {
			test=	tree_abs/cut_result_abs < _IR_tolerance_RHP
					&&
					!(tree_abs < _IR_tolerance_RHP);
		} else {
		test=abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RHP   // big cancellation
				||
				tree_abs/cut_result_abs < _IR_tolerance_RHP  || tree_abs/rat_result_abs < _IR_tolerance_RHP
				&&
				!(tree_abs < _IR_tolerance_RHP);
		}
	}
	if ( test ){
		multi_precision_reader_HP* mpr=dynamic_cast<multi_precision_reader_HP*>(&mc);
		if ( mpr ==0 ){
#if _VERBOSE
			_WARNING("pointer is not a multi_precision_reader, we improvise");
//			throw BHerror("IR_checked_OLHA::eval called with an argument that is not a multi_precision_reader");
#endif
			SeriesC<RVHP> VHP_result;
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_VHP mcVHP=mc.extend<RVHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_VHP ;
// no need for all this mcVHP has just been created, and muVHP nevertheless is RVHP(muHP)
//			if (_cut_part->get_mu_index<RVHP>() == 0) {
//				new_mu_index_VHP = DefineMu<RVHP>(mcVHP,_cut_part->get_mu());
//			} else {
//				new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(_cut_part->get_mu_index<RHP>()).E().real()));
//			}
			new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(_cut_part->get_mu_index<RHP>()).E().real()));
			_cut_part->set_mu_VHP(new_mu_index_VHP);
			VHP_result=eval(mcVHP,new_ind);
			return to_HP(VHP_result);
		}
		 else {
			 SeriesC<RVHP> VHP_result=eval(mpr->mc_VHP(),ind);
			return to_HP(VHP_result);
		}

	}
	
	// Collect the accuracy estimate
	double ratacc=_rational_part->get_accuracy();
	double cutacc=_cut_part->get_accuracy();
	ratacc>cutacc?_accuracy=cutacc:_accuracy=ratacc;	

    // fill in conjugate value of amplitude
    _conjugate_rat_part=to_double(conj(rat_result));
    _conjugate_rat_part_HP=to_HP(conj(rat_result));
    _conjugate_rat_part_VHP=to_VHP(conj(rat_result));
    _conjugate_cut_part=_cut_part->get_conjugate_cut_part();
    _conjugate_cut_part_HP=_cut_part->get_conjugate_cut_part_HP();
    _conjugate_cut_part_VHP=_cut_part->get_conjugate_cut_part_VHP();

	return rat_result+cut_result;
}

SeriesC<RVHP> IR_checked_OLHA::eval(mom_conf_VHP& mc,const vector<int>& ind){
	CVHP tree=get_tree(mc,ind);
	SeriesC<RVHP> cut_result=_cut_part->get_value(mc,ind);
	CVHP rat_result=_rational_part->get_value(mc,ind);
	
	RVHP tree_abs=abs(tree);
	RVHP cut_result_abs=abs(cut_result[0]);
	RVHP rat_result_abs=abs(rat_result);
	
	bool test;
	if (cut_result_abs == 0.){
		if (rat_result_abs == 0.0) {
			test = !(tree_abs < _IR_tolerance_RVHP);
		} else {
		test=
				tree_abs/rat_result_abs < _IR_tolerance_RVHP
				&&
				!(tree_abs < _IR_tolerance_RVHP);
		}
	} else {
		if (rat_result_abs == 0.0) {
			test=	tree_abs/cut_result_abs < _IR_tolerance_RVHP
					&&
					!(tree_abs < _IR_tolerance_RVHP);
		} else {
		test=abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RVHP   // big cancellation
				||
				tree_abs/cut_result_abs < _IR_tolerance_RVHP  || tree_abs/rat_result_abs < _IR_tolerance_RVHP
				&&
				!(tree_abs < _IR_tolerance_RVHP);
		}
	}


	if ( test){
#if _VERBOSE
		_MESSAGE("very very large cancellation!");
#endif
	}
	
	// Collect the accuracy estimate
	double ratacc=_rational_part->get_accuracy();
	double cutacc=_cut_part->get_accuracy();
	ratacc>cutacc?_accuracy=cutacc:_accuracy=ratacc;
	
    // fill in conjugate value of amplitude
    _conjugate_rat_part=to_double(conj(rat_result));
    _conjugate_rat_part_HP=to_HP(conj(rat_result));
    _conjugate_rat_part_VHP=to_VHP(conj(rat_result));
    _conjugate_cut_part=_cut_part->get_conjugate_cut_part();
    _conjugate_cut_part_HP=_cut_part->get_conjugate_cut_part_HP();
    _conjugate_cut_part_VHP=_cut_part->get_conjugate_cut_part_VHP();

	return rat_result+cut_result;
}

}


#ifdef BH_USE_GMP
SeriesC<RGMP> IR_checked_OLHA::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind){
	CGMP tree=get_tree(mc,ind);
	SeriesC<RGMP> cut_result=_cut_part->get_value(mc,ind);
	CGMP rat_result=_rational_part->get_value(mc,ind);

	RGMP tree_abs=abs(tree);
	RGMP cut_result_abs=abs(cut_result[0]);
	RGMP rat_result_abs=abs(rat_result);

	bool test;
	if (cut_result_abs == 0.){
		if (rat_result_abs == 0.0) {
			test = !(tree_abs < _IR_tolerance_RGMP);
		} else {
		test=
				tree_abs/rat_result_abs < _IR_tolerance_RGMP
				&&
				!(tree_abs < _IR_tolerance_RGMP);
		}
	} else {
		if (rat_result_abs == 0.0) {
			test=	tree_abs/cut_result_abs < _IR_tolerance_RGMP
					&&
					!(tree_abs < _IR_tolerance_RGMP);
		} else {
		test=abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RGMP   // big cancellation
				||
				tree_abs/cut_result_abs < _IR_tolerance_RGMP  || tree_abs/rat_result_abs < _IR_tolerance_RGMP
				&&
				!(tree_abs < _IR_tolerance_RGMP);
		}
	}


	if ( test){
#if _VERBOSE
		_MESSAGE("very very large cancellation!");
#endif
	}

	// Collect the accuracy estimate
	double ratacc=_rational_part->get_accuracy();
	double cutacc=_cut_part->get_accuracy();
	ratacc>cutacc?_accuracy=cutacc:_accuracy=ratacc;

    // fill in conjugate value of amplitude
    _conjugate_rat_part=to_double(conj(rat_result));
    _conjugate_rat_part_HP=to_HP(conj(rat_result));
    _conjugate_rat_part_VHP=to_VHP(conj(rat_result));
    _conjugate_cut_part=_cut_part->get_conjugate_cut_part();
    _conjugate_cut_part_HP=_cut_part->get_conjugate_cut_part_HP();
    _conjugate_cut_part_VHP=_cut_part->get_conjugate_cut_part_VHP();

	return rat_result+cut_result;
}
#endif

