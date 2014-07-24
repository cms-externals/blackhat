/*
 * OneLoopHelAmpl.cpp
 *
 *  Created on: 22-Apr-2009
 *      Author: daniel
 */

#include "amplitudes.h"
#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "integrals.h"
#include "BH_A0.h"
#include "OneLoopHelAmpl.h"
#include "IR_checked.h"
#include "cut_part_factory.h"
#include "scheme.h"
#include "process_utils.h"
#include "process_utils_massive.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "settings.h"
#include <cassert>
#include "ratext/ratext_part.h"
#include "BH_debug.h"
//#include "BG_tree.h"

//#define __ONE_LOOP_TREE_TYPE BG_TREE_FACTORY_V
#define __ONE_LOOP_TREE_TYPE TREE_FACTORY_TYPE

using namespace std;
namespace BH {

/*
BH_DEBUG(
using CachedIntegral::Known_Cut_wCI;
using CachedIntegral::Known_Cut_wCI_offset;
)
*/
One_Loop_Helicity_Amplitude::One_Loop_Helicity_Amplitude(const process& PRO,color_structure cs, Rational_factory<Rational_base>* RRF,cut_part_factory<Cut_Part_base>* CPF) : OneLoopAmplitude_base(PRO,cs) , _scheme(FDH), _accuracy(0), 
	_conjugate_cut_part(SeriesC<R>(-2,0)), 
	_conjugate_cut_part_HP(SeriesC<RHP>(-2,0)), 
	_conjugate_cut_part_VHP(SeriesC<RVHP>(-2,0)), 
	_conjugate_rat_part(complex<R>(0,0)),
	_conjugate_rat_part_HP(complex<RHP>(0,0)),
	_conjugate_rat_part_VHP(complex<RVHP>(0,0)) {

	init(PRO,cs,RRF,CPF);
}

One_Loop_Helicity_Amplitude::One_Loop_Helicity_Amplitude(const process& PRO,color_structure cs, Rational_factory<Rational_base>* RRF) : OneLoopAmplitude_base(PRO,cs) , _scheme(FDH), _accuracy(0), 
	_conjugate_cut_part(SeriesC<R>(-2,0)), 
	_conjugate_cut_part_HP(SeriesC<RHP>(-2,0)), 
	_conjugate_cut_part_VHP(SeriesC<RVHP>(-2,0)), 
	_conjugate_rat_part(complex<R>(0,0)),
	_conjugate_rat_part_HP(complex<RHP>(0,0)),
	_conjugate_rat_part_VHP(complex<RVHP>(0,0)) {
	cut_part_factory<Cut_Part_base>* CPF = cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(PRO);
    
	init(PRO,cs,RRF,CPF);
}

void One_Loop_Helicity_Amplitude::init(const process& PRO,color_structure cs, Rational_factory<Rational_base>* RRF,cut_part_factory<Cut_Part_base>* CPF) {
	__ONE_LOOP_TREE_TYPE TF;
	d_tree_ptr=TF.new_tree(PRO);

	_cut_part=CPF->new_cut_part(PRO, cs);
/*
BH_DEBUG(
    //CachedIntegral::Known_Cut_Part_offset* KCBO;
    CachedIntegral::Known_Cut_wCI* KCwCI;
    //CachedIntegral::Known_Cut_wCI_offset* KCwCIo;
	if ( KCwCI = dynamic_cast<Known_Cut_wCI*>(_cut_part) ){
		_MESSAGE4("Known Cut_Part with Cached Integral ",PRO," ",cs);
    }
	else {
		_MESSAGE4("Numerical Cut_Part ",PRO," ",cs);
    }; 
)
*/
	_rational_part=RRF->new_rational(PRO,cs);

/*
BH_DEBUG(
    //CachedIntegral::Known_Cut_Part_offset* KCBO;
    Known_Rec_Rational* KRR;
    //CachedIntegral::Known_Cut_wCI_offset* KCwCIo;
	if ( KRR = dynamic_cast<Known_Rec_Rational*>(_rational_part) ){
		_MESSAGE4("Known Rat_Part with Cached Integral ",PRO," ",cs);
    }
	else {
		_MESSAGE4("Numerical Rat_Part ",PRO," ",cs);
    }; 
)
*/


	//makes sure we have a tree
	assert(d_tree_ptr);
	//makes sure we have a cut part
	assert(_cut_part);
	//makes sure we have a rational part
	assert(_rational_part);
}


One_Loop_Helicity_Amplitude::~One_Loop_Helicity_Amplitude(){
	delete d_tree_ptr;
	delete _rational_part;
	delete _cut_part;
}

//
//SeriesC<R> One_Loop_Helicity_Amplitude::eval(mom_conf& mc,const vector<int>& ind){
//	C rat_result=_rational_part->get_value(mc,ind);
//	SeriesC<R> cut_result=_cut_part->get_value(mc,ind);
//	C Atree=_cut_part->get_tree(mc,ind);
////	CSeries cut_result=_cut_part->eval(mc,ind);
//#if _VERBOSE
//	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
//	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
//#endif
//	return rat_result+cut_result+scheme_shift<R>(get_process(),_scheme)*Atree;
//}
//SeriesC<RHP> One_Loop_Helicity_Amplitude::eval(mom_conf_HP& mc,const vector<int>& ind){
//	CHP rat_result=_rational_part->get_value(mc,ind);
//	SeriesC<RHP> cut_result=_cut_part->get_value(mc,ind);
//	CHP Atree=_cut_part->get_tree(mc,ind);
//#if _VERBOSE
//	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
//	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
//#endif
//	return rat_result+cut_result;
//}
//SeriesC<RVHP> One_Loop_Helicity_Amplitude::eval(mom_conf_VHP& mc,const vector<int>& ind){
//	CVHP rat_result=_rational_part->get_value(mc,ind);
//	SeriesC<RVHP> cut_result=_cut_part->get_value(mc,ind);
//	CVHP Atree=_cut_part->get_tree(mc,ind);
//#if _VERBOSE
//	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
//	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
//#endif
//	return rat_result+cut_result;
//}

template <class T> SeriesC<T> One_Loop_Helicity_Amplitude::eval_fn(momentum_configuration<T>& mc,const vector<int>& ind){
	// Evaluate the amplitude
	complex<T> rat_result=_rational_part->get_value(mc,ind);
	SeriesC<T> cut_result=_cut_part->get_value(mc,ind);
	complex<T> Atree=_cut_part->get_tree(mc,ind);
	
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


#if _VERBOSE
	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
#endif
	return rat_result+cut_result;
}



SeriesC<R> One_Loop_Helicity_Amplitude::eval(mom_conf& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
SeriesC<RHP> One_Loop_Helicity_Amplitude::eval(mom_conf_HP& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
SeriesC<RVHP> One_Loop_Helicity_Amplitude::eval(mom_conf_VHP& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};

#if BH_USE_GMP
SeriesC<RGMP> One_Loop_Helicity_Amplitude::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind){ return eval_fn(mc,ind);};
#endif

SeriesC<R> One_Loop_Helicity_Amplitude::eval(const eval_param<R>& ep){
	C rat_result=_rational_part->get_value(ep);
	SeriesC<R> cut_result=_cut_part->get_value(ep);
	C Atree=_cut_part->get_tree(ep);

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


#if _VERBOSE
	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
#endif
	return rat_result+cut_result+scheme_shift<R>(get_process(),_scheme)*Atree;
}

SeriesC<RHP> One_Loop_Helicity_Amplitude::eval(const eval_param<RHP>& ep){
	CHP rat_result=_rational_part->get_value(ep);
	SeriesC<RHP> cut_result=_cut_part->get_value(ep);
	CHP Atree=_cut_part->get_tree(ep);
	
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


#if _VERBOSE
	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
#endif
	return rat_result+cut_result;
}

SeriesC<RVHP> One_Loop_Helicity_Amplitude::eval(const eval_param<RVHP>& ep){
	CVHP rat_result=_rational_part->get_value(ep);
	SeriesC<RVHP> cut_result=_cut_part->get_value(ep);
	CVHP Atree=_cut_part->get_tree(ep);
	
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


#if _VERBOSE
	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
#endif
	return rat_result+cut_result;
}

#if BH_USE_GMP
SeriesC<RGMP> One_Loop_Helicity_Amplitude::eval(const eval_param<RGMP>& ep){
	CGMP rat_result=_rational_part->get_value(ep);
	SeriesC<RGMP> cut_result=_cut_part->get_value(ep);
	CGMP Atree=_cut_part->get_tree(ep);

	// Collect the accuracy estimate
	double ratacc=_rational_part->get_accuracy();
	double cutacc=_cut_part->get_accuracy();
	ratacc>cutacc?_accuracy=cutacc:_accuracy=ratacc;

    // fill in conjugate value of amplitude
    _conjugate_rat_part=to_double(conj(rat_result));
    _conjugate_cut_part=to_double(_cut_part->get_conjugate_cut_part());



#if _VERBOSE
	_MESSAGE4("rational: ", rat_result ," normalized: ", rat_result/Atree);
	_MESSAGE4("cut part: ", cut_result ," normalized: ", cut_result/Atree);
#endif
	return rat_result+cut_result;
}
#endif





}
