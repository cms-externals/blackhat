/*
 * IR_checked_amplitudes_ep.cpp
 *
 *  Created on: 21-May-2009
 *      Author: darrenforde
 */
//#include "OneLoopHelAmpl.h"
//#include "BH_A0.h"
#include "qd_suppl.h"
#include <typeinfo>
#include "IR_checked.h"
#include "IR_checked_cut_part.h"
#include "cached_integral.h"
#include <cassert>
#include "BH_debug.h"
#define _VERBOSE 0   // there are two levels of verbosity : 1 and 2 (more output)
using namespace std;
namespace BH {




SeriesC<R> IR_checked_OLHA::eval(const eval_param<R>& ep){
	BH_DEBUG_MESSAGE2("Computing amplitude (ep) ",d_process);
	C tree=get_tree(ep);
	SeriesC<R> cut_result=_cut_part->get_value(ep);
	C rat_result=_rational_part->get_value(ep);
	
	R tree_abs=abs(tree);
	R cut_result_abs=abs(cut_result[0]);
	R rat_result_abs=abs(rat_result);

	if ( abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_R   // big cancellation
		||
		tree_abs/cut_result_abs < _IR_tolerance_R  || tree_abs/rat_result_abs < _IR_tolerance_R
		&&
		!(tree_abs < _IR_tolerance_R)
	){
		SeriesC<RHP> HP_result;
		vector<Cmom<RHP> > *new_mom=extend_momenta<R,RHP>(ep);
		eval_param<RHP> epHP(*new_mom);
		HP_result=eval(epHP);
		cut_result=to_double(HP_result);
		delete new_mom;
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

SeriesC<RHP> IR_checked_OLHA::eval(const eval_param<RHP>& ep){
	CHP tree=get_tree(ep);
	SeriesC<RHP> cut_result=_cut_part->get_value(ep);
	CHP rat_result=_rational_part->get_value(ep);

	RHP tree_abs=abs(tree);
	RHP cut_result_abs=abs(cut_result[0]);
	RHP rat_result_abs=abs(rat_result);
	
	if ( abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RHP   // big cancellation
		||
		tree_abs/cut_result_abs < _IR_tolerance_RHP  || tree_abs/rat_result_abs < _IR_tolerance_RHP
		&&
		!(tree_abs < _IR_tolerance_RHP)
	){
		SeriesC<RVHP> VHP_result;
		vector<Cmom<RVHP> > *new_mom=extend_momenta<RHP,RVHP>(ep);
		eval_param<RVHP> epVHP(*new_mom);
		VHP_result=eval(epVHP);
		cut_result = to_HP(VHP_result);
		delete new_mom;
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

SeriesC<RVHP> IR_checked_OLHA::eval(const eval_param<RVHP>& ep){
	CVHP tree=get_tree(ep);
	SeriesC<RVHP> cut_result=_cut_part->get_value(ep);
	CVHP rat_result=_rational_part->get_value(ep);

	RVHP tree_abs=abs(tree);
	RVHP cut_result_abs=abs(cut_result[0]);
	RVHP rat_result_abs=abs(rat_result);
	
	if ( abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RVHP   // big cancellation
		||
		tree_abs/cut_result_abs < _IR_tolerance_RVHP  || tree_abs/rat_result_abs < _IR_tolerance_RVHP
		&&
		!(tree_abs < _IR_tolerance_RVHP)
	){
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


#ifdef BH_USE_GMP
SeriesC<RGMP> IR_checked_OLHA::eval(const eval_param<RGMP>& ep){
	CGMP tree=get_tree(ep);
	SeriesC<RGMP> cut_result=_cut_part->get_value(ep);
	CGMP rat_result=_rational_part->get_value(ep);

	RGMP tree_abs=abs(tree);
	RGMP cut_result_abs=abs(cut_result[0]);
	RGMP rat_result_abs=abs(rat_result);

	if ( abs(cut_result[0]+rat_result) /(cut_result_abs+rat_result_abs) <  _IR_tolerance_RGMP   // big cancellation
		||
		tree_abs/cut_result_abs < _IR_tolerance_RGMP  || tree_abs/rat_result_abs < _IR_tolerance_RGMP
		&&
		!(tree_abs < _IR_tolerance_RGMP)
	){
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


}

