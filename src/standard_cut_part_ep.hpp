/*
 * standard_cut_part_ep.hpp
 *
 *  Created on: 6 Aug 2009
 *      Author: daniel
 */

#ifndef STANDARD_CUT_PART_EP_HPP_
#define STANDARD_CUT_PART_EP_HPP_


#include "cut_part.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include<vector>
#include<complex>

#include "BH_typedefs.h"
#include "amplitudes.h"
#include "partitions.h"
#ifndef BH_PUBLIC
#include "options.h"
#endif
#include "integrals_ep.h"
#include "BH_A0.h"
#include "rational.h"
#include "cut_eval/amplitudes_cut_eval.h"
#include "process_utils.h"
#include "settings.h"
#include "cut_part_worker.h"
#ifndef BH_PUBLIC
#include "cut_part_normal.h"
#endif
#include "BH_debug.h"


using namespace std;

#define _VERBOSE 0
#define _CONSERVATIVE_ERROR_ESTIMATE 1 //Set this to 1 if you want to take into account the actual error of the coefficeint divided by the tree

namespace BH {

namespace cut {


template <class BOX,class TRI,class BUB> template <class T> SeriesC<T> standard_cut_part<BOX,TRI,BUB>::eval_without_check(const eval_param<T>& ep){
	SeriesC<T> intval(-2,0);
	// Add up the bubbles
	SeriesC<T> totalbubble(-2,0);
	SeriesC<T> totalbubble_conj(-2,0);
	complex<T> coeffbubble;

	//Get the correct mu scale for the integrals
	multi_precision_constant mmu(get_mu());
	T musqr(mmu*mmu);
	
	double curacc(to_double(MaxDigits<T>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	T treeres=abs(get_tree(ep)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<T>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=T(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(T(10)));
	}
#endif

#if _VERBOSE==1
	complex<T> bubc(0,0), tric(0,0), boxc(0,0);
	complex<T> treeres_full=get_tree(ep);
	if(abs(treeres_full)<BH::DeltaZero<T>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=T(1);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres_full);}
	_MESSAGE("");
	_MESSAGE("**************************** Cut Part NO IR test (ep) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

	 for (size_t i=1;i<=nbr_bubbles();i++){
		 	BUB* thebubble=bubble(i);

		vector<int> c1( Indicesm1(thebubble,1));
		vector<int> c2( Indicesm1(thebubble,2));
		intval = IntM(ep,musqr,c1,c2,_masses);

	       coeffbubble = d_bubbles[i-1]->eval(ep);

		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
#if _VERBOSE==1
//		_MESSAGE7(*thebubble,"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy()+acc_shift);
//		 _MESSAGE5(*thebubble,"=",coeffbubble*intval/treeres_full," estimated accuracy=",thebubble->get_accuracy()+acc_shift);
	    _MESSAGE7(*thebubble,"=",coeffbubble/treeres_full," * ",intval," estimated accuracy=",thisacc+acc_shift);
		 bubc+=coeffbubble;
#endif
	       totalbubble += coeffbubble*intval;
	       totalbubble_conj += conj(coeffbubble)*intval;
	  }
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
	_MESSAGE4("      Sum of all bubbles(n)=",totalbubble/treeres_full," from ",bubc);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
	 SeriesC<T> totaltri(-2,0);
	 SeriesC<T> totaltri_conj(-2,0);
	complex<T> coefftri;

	 for (size_t i=1;i<=nbr_triangles();i++){

	       //  Call Int in integrals.cpp
		vector<int> c1( Indicesm1(triangle(i),1) );
		vector<int> c2( Indicesm1(triangle(i),2) );
		vector<int> c3( Indicesm1(triangle(i),3) );
		intval = IntM(ep,musqr,c1,c2,c3,_masses);

	       // The contribution to the coefficient divided by tree
	       coefftri = d_triangles[i-1]->eval(ep);
#if _VERBOSE==1
//		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
//		_MESSAGE3(*d_triangles[i-1],"=",coefftri*intval/treeres_full);
		_MESSAGE3(*d_triangles[i-1],"=",coefftri/treeres_full);
		 tric+=coefftri;
#endif
	       totaltri += coefftri*intval;
	       totaltri_conj += conj(coefftri)*intval;
	  }

	// Now add up the boxes

	 SeriesC<T> totalbox(-2,0);
	 SeriesC<T> totalbox_conj(-2,0);
	 complex<T> coeffbox;
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all triangles=",totaltri);
	_MESSAGE4("      Sum of all triangles(n)=",totaltri/treeres_full," from ",tric);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_boxes(), " Boxes");
#endif
	 for (unsigned int i=1;i<=nbr_boxes();i++){

	     //  Here we call the integrals .
		vector<int> d1( Indicesm1(box(i),1) );
		vector<int> d2( Indicesm1(box(i),2) );
		vector<int> d3( Indicesm1(box(i),3) );
		vector<int> d4( Indicesm1(box(i),4) );
		intval = IntM(ep,musqr,d1,d2,d3,d4,_masses);

	       coeffbox = d_boxes[i-1]->eval(ep);
#if _VERBOSE==1
//		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
//	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox*intval/treeres_full);
	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox/treeres_full);
		 boxc+=coeffbox;
#endif
	       totalbox += coeffbox*intval;
	       totalbox_conj += conj(coeffbox)*intval;
	  }
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all boxes=",totalbox);
	_MESSAGE4("      Sum of all boxes(n)=",totalbox/treeres_full," from ",boxc);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
//	_MESSAGE2("      Final result=",totalbubble + totaltri + totalbox);
	_MESSAGE2("      Final result(n)=",(totalbubble + totaltri + totalbox)/treeres_full);
#endif
	
	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<T> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    	_conjugate_cut_part=to_double(conjugate_cut_part);
    	_conjugate_cut_part_HP=to_HP(conjugate_cut_part);
    	_conjugate_cut_part_VHP=to_VHP(conjugate_cut_part);

	 // Sum together the bubble triangle and box results
	 return totalbubble + totaltri + totalbox ;
}


template <class BOX,class TRI,class BUB> SeriesC<R> standard_cut_part<BOX,TRI,BUB>::eval_with_check(const eval_param<R>& ep){
	SeriesC<R> intval(-2,0);
	// Add up the bubbles
	SeriesC<R> totalbubble(-2,0);
	SeriesC<R> totalbubble_conj(-2,0);
	complex<R> coeffbubble;

	//Get the correct mu scale for the integrals
	multi_precision_constant mmu(get_mu());
	R musqr(mmu*mmu);
	
	double curacc(MaxDigits<R>()); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	R treeres=abs(get_tree(ep)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<R>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=R(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(10.));
	}
#endif

#if _VERBOSE==1
	complex<R> treeres_full=get_tree(ep);
	if(abs(treeres_full)<BH::DeltaZero<R>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=R(1);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres_full);}
	_MESSAGE("");
	_MESSAGE("**************************** Cut Part WITH IR test (ep) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

		vector<int> c1( Indicesm1(thebubble,1) );
		vector<int> c2( Indicesm1(thebubble,2) );
		intval = IntM(ep,musqr,c1,c2,_masses);

		coeffbubble = d_bubbles[i-1]->eval(ep);

		if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE2("WARNING : Recomputing bubble at HP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
			vector<Cmom<RHP> > *new_mom=extend_momenta<R,RHP>(ep);
			eval_param<RHP> epHP(*new_mom);
			complex<RHP> HP_result=(d_bubbles[i-1]->eval(epHP));
			coeffbubble=to_double(HP_result);
			delete new_mom; //delete the new momentum vector created in the extend function

			//Check if this was accurate enough
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
				// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
				vector<Cmom<RVHP> > *new_mom=extend_momenta<R,RVHP>(ep);
				eval_param<RVHP> epVHP(*new_mom);
				complex<RVHP> VHP_result=(d_bubbles[i-1]->eval(epVHP));
				coeffbubble=to_double(VHP_result);
				delete new_mom; //delete the new momentum vector created in the extend function
				 //If this was to fail then there is not much more we can do
#if _VERBOSE==1
				if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
					_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
				}
#endif
			}
		}

		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
#if _VERBOSE==1
//		_MESSAGE7(d_bubbles[i-1],"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy());
//	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble*intval/treeres_full," estimated accuracy=",d_bubbles[i-1]->get_accuracy());
	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble/treeres_full," estimated accuracy=",thisacc+acc_shift);
#endif

		totalbubble += coeffbubble*intval;
	    totalbubble_conj += conj(coeffbubble)*intval;
	}
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
	_MESSAGE2("      Sum of all bubbles(n)=",totalbubble/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
	SeriesC<R> totaltri(-2,0);
	SeriesC<R> totaltri_conj(-2,0);
	complex<R> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

		//  Call Int in integrals.cpp
		vector<int> c1( Indicesm1(triangle(i),1) );
		vector<int> c2( Indicesm1(triangle(i),2) );
		vector<int> c3( Indicesm1(triangle(i),3) );
		intval = IntM(ep,musqr,c1,c2,c3,_masses);

		// The contribution to the coefficient divided by tree
		coefftri = d_triangles[i-1]->eval(ep);
#if _VERBOSE==1
//		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
//		_MESSAGE3(*d_triangles[i-1],"=",coefftri*intval/treeres_full);
		_MESSAGE3(*d_triangles[i-1],"=",coefftri/treeres_full);
#endif
		totaltri += coefftri*intval;
		totaltri_conj += conj(coefftri)*intval;
	}

#if _VERBOSE==1
//	_MESSAGE2("      Sum of all triangles=",totaltri);
	_MESSAGE2("      Sum of all triangles(n)=",totaltri/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_boxes(), " Boxes");
#endif
	// Now add up the boxes

	SeriesC<R> totalbox(-2,0);
	SeriesC<R> totalbox_conj(-2,0);
	complex<R> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
		vector<int> d1( Indicesm1(box(i),1) );
		vector<int> d2( Indicesm1(box(i),2) );
		vector<int> d3( Indicesm1(box(i),3) );
		vector<int> d4( Indicesm1(box(i),4) );
		intval = IntM(ep,musqr,d1,d2,d3,d4,_masses);

		coeffbox = d_boxes[i-1]->eval(ep);
#if _VERBOSE==1
//		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
//	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox*intval/treeres_full);
	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox/treeres_full);
#endif
		totalbox += coeffbox*intval;
		totalbox_conj += conj(coeffbox)*intval;
	}
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all boxes=",totalbox);
	_MESSAGE2("      Sum of all boxes(n)=",totalbox/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
//	_MESSAGE2("      Final result=",totalbubble + totaltri + totalbox);
	_MESSAGE2("      Final result(n)=",(totalbubble + totaltri + totalbox)/treeres_full);
#endif
	
	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<R> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    	_conjugate_cut_part=to_double(conjugate_cut_part);
    	_conjugate_cut_part_HP=to_HP(conjugate_cut_part);
    	_conjugate_cut_part_VHP=to_VHP(conjugate_cut_part);

	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}

template <class BOX,class TRI,class BUB> SeriesC<RHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(const eval_param<RHP>& ep){
	SeriesC<RHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RHP> totalbubble(-2,0);
	SeriesC<RHP> totalbubble_conj(-2,0);
	complex<RHP> coeffbubble;

	//Get the correct mu scale for the integrals
	multi_precision_constant mmu(get_mu());
	RHP musqr(mmu*mmu);
	
	double curacc(to_double(MaxDigits<RHP>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	RHP treeres=abs(get_tree(ep)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<RHP>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=RHP(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(RHP(10)));
	}
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

		vector<int> c1( Indicesm1(thebubble,1) );
		vector<int> c2( Indicesm1(thebubble,2) );
		intval = IntM(ep,musqr,c1,c2,_masses);

		coeffbubble = d_bubbles[i-1]->eval(ep);
		if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
			// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
			vector<Cmom<RVHP> > *new_mom=extend_momenta<RHP,RVHP>(ep);
			eval_param<RVHP> epVHP(*new_mom);
			complex<RVHP> VHP_result=(d_bubbles[i-1]->eval(epVHP));
			coeffbubble=to_HP(VHP_result);
			delete new_mom; //delete the new momentum vector created in the extend function
			//If this was to fail then there is not much more we can do
#if _VERBOSE==1
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
				_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
			}
#endif
		}

		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
		totalbubble += coeffbubble*intval;
	    totalbubble_conj += conj(coeffbubble)*intval;
	}
#if _VERBOSE
	_PRINT(totalbubble);
#endif
	SeriesC<RHP> totaltri(-2,0);
	SeriesC<RHP> totaltri_conj(-2,0);
	complex<RHP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

		//  Call Int in integrals.cpp
		vector<int> c1( Indicesm1(triangle(i),1) );
		vector<int> c2( Indicesm1(triangle(i),2) );
		vector<int> c3( Indicesm1(triangle(i),3) );
		intval = IntM(ep,musqr,c1,c2,c3,_masses);

		// The contribution to the coefficient divided by tree
		coefftri = d_triangles[i-1]->eval(ep);
		totaltri += coefftri*intval;
		totaltri_conj += conj(coefftri)*intval;
	}

#if _VERBOSE
	_PRINT(totaltri);
#endif
	// Now add up the boxes

	SeriesC<RHP> totalbox(-2,0);
	SeriesC<RHP> totalbox_conj(-2,0);
	complex<RHP> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
		vector<int> d1( Indicesm1(box(i),1) );
		vector<int> d2( Indicesm1(box(i),2) );
		vector<int> d3( Indicesm1(box(i),3) );
		vector<int> d4( Indicesm1(box(i),4) );
		intval = IntM(ep,musqr,d1,d2,d3,d4,_masses);

		coeffbox = d_boxes[i-1]->eval(ep);
		totalbox += coeffbox*intval;
		totalbox_conj += conj(coeffbox)*intval;
	}
#if _VERBOSE
	_PRINT(totalbox);
#endif
	
	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<RHP> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    _conjugate_cut_part=to_double(conjugate_cut_part);
    _conjugate_cut_part_HP=to_HP(conjugate_cut_part);
    _conjugate_cut_part_VHP=to_VHP(conjugate_cut_part);

	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}

template <class BOX,class TRI,class BUB> SeriesC<RVHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(const eval_param<RVHP>& ep){
	SeriesC<RVHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RVHP> totalbubble(-2,0);
	SeriesC<RVHP> totalbubble_conj(-2,0);
	complex<RVHP> coeffbubble;

	//Get the correct mu scale for the integrals
	multi_precision_constant mmu(get_mu());
	RVHP musqr(mmu*mmu);
	
	double curacc(to_double(MaxDigits<RVHP>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	RVHP treeres=abs(get_tree(ep)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<RVHP>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=RVHP(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(RVHP(10)));
	}
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

		vector<int> c1( Indicesm1(thebubble,1) );
		vector<int> c2( Indicesm1(thebubble,2) );
		intval = IntM(ep,musqr,c1,c2,_masses);

		coeffbubble = d_bubbles[i-1]->eval(ep);
#if _VERBOSE==1
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
				_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
			}
#endif
		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
		totalbubble += coeffbubble*intval;
	    totalbubble_conj += conj(coeffbubble)*intval;
	}
#if _VERBOSE
	_PRINT(totalbubble);
#endif
	SeriesC<RVHP> totaltri(-2,0);
	SeriesC<RVHP> totaltri_conj(-2,0);
	complex<RVHP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

		//  Call Int in integrals.cpp
		vector<int> c1( Indicesm1(triangle(i),1) );
		vector<int> c2( Indicesm1(triangle(i),2) );
		vector<int> c3( Indicesm1(triangle(i),3) );
		intval = IntM(ep,musqr,c1,c2,c3,_masses);

		// The contribution to the coefficient divided by tree
		coefftri = d_triangles[i-1]->eval(ep);
		totaltri += coefftri*intval;
		totaltri_conj += conj(coefftri)*intval;
	}

#if _VERBOSE
	_PRINT(totaltri);
#endif
	// Now add up the boxes

	SeriesC<RVHP> totalbox(-2,0);
	SeriesC<RVHP> totalbox_conj(-2,0);
	complex<RVHP> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
		vector<int> d1( Indicesm1(box(i),1) );
		vector<int> d2( Indicesm1(box(i),2) );
		vector<int> d3( Indicesm1(box(i),3) );
		vector<int> d4( Indicesm1(box(i),4) );
		intval = IntM(ep,musqr,d1,d2,d3,d4,_masses);

		coeffbox = d_boxes[i-1]->eval(ep);
		totalbox += coeffbox*intval;
		totalbox_conj += conj(coeffbox)*intval;
	}
#if _VERBOSE
	_PRINT(totalbox);
#endif
	
	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<RVHP> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    _conjugate_cut_part=to_double(conjugate_cut_part);
    _conjugate_cut_part_HP=to_HP(conjugate_cut_part);
    _conjugate_cut_part_VHP=to_VHP(conjugate_cut_part);

	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}

#if BH_USE_GMP
template <class BOX,class TRI,class BUB> SeriesC<RGMP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(const eval_param<RGMP>& ep){
	SeriesC<RGMP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RGMP> totalbubble(-2,0);
	SeriesC<RGMP> totalbubble_conj(-2,0);
	complex<RGMP> coeffbubble;

	//Get the correct mu scale for the integrals
	multi_precision_constant mmu(get_mu());
	RGMP musqr(mmu*mmu);

	double curacc(to_double(MaxDigits<RGMP>())); // keep the lowest recorded accuracy



	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	RGMP treeres=abs(get_tree(ep)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<RGMP>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=RGMP(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(RGMP(10)));
	}
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

		vector<int> c1( Indicesm1(thebubble,1) );
		vector<int> c2( Indicesm1(thebubble,2) );
		intval = IntM(ep,musqr,c1,c2,_masses);

		coeffbubble = d_bubbles[i-1]->eval(ep);
#if _VERBOSE==1
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
				_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at GMP accuracy.");
			}
#endif
		double thisacc=thebubble->get_accuracy();

		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		BH_DEBUG_PRINT(*this);
		BH_DEBUG_PRINT(coeffbubble);
		BH_DEBUG_PRINT(intval);

		totalbubble += coeffbubble*intval;
		totalbubble_conj += conj(coeffbubble)*intval;
	}
#if _VERBOSE
	_PRINT(totalbubble);
#endif
	SeriesC<RGMP> totaltri(-2,0);
	SeriesC<RGMP> totaltri_conj(-2,0);
	complex<RGMP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

		//  Call Int in integrals.cpp
		vector<int> c1( Indicesm1(triangle(i),1) );
		vector<int> c2( Indicesm1(triangle(i),2) );
		vector<int> c3( Indicesm1(triangle(i),3) );
		intval = IntM(ep,musqr,c1,c2,c3,_masses);

		// The contribution to the coefficient divided by tree
		coefftri = d_triangles[i-1]->eval(ep);
		totaltri += coefftri*intval;
		totaltri_conj += conj(coefftri)*intval;
	}

#if _VERBOSE
	_PRINT(totaltri);
#endif
	// Now add up the boxes

	SeriesC<RGMP> totalbox(-2,0);
	SeriesC<RGMP> totalbox_conj(-2,0);
	complex<RGMP> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
		vector<int> d1( Indicesm1(box(i),1) );
		vector<int> d2( Indicesm1(box(i),2) );
		vector<int> d3( Indicesm1(box(i),3) );
		vector<int> d4( Indicesm1(box(i),4) );
		intval = IntM(ep,musqr,d1,d2,d3,d4,_masses);

		BH_DEBUG_MESSAGE9(d1," ",d2," ",d3," ",d4,": ",intval)

		coeffbox = d_boxes[i-1]->eval(ep);
		totalbox += coeffbox*intval;
		totalbox_conj += conj(coeffbox)*intval;
	}
#if _VERBOSE
	_PRINT(totalbox);
#endif

	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<RGMP> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    _conjugate_cut_part=to_double(conjugate_cut_part);
    //_MESSAGE("missing _conjugate_cut_part_GMP");


	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;}
#endif

} /* cut */
} /* BH */

#endif /* STANDARD_CUT_PART_EP_HPP_ */
