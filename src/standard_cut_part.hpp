/*
 * standard_cut_part.hpp
 *
 *  Created on: 6 Aug 2009
 *      Author: daniel
 */

#ifndef STANDARD_CUT_PART_HPP_
#define STANDARD_CUT_PART_HPP_

#include "standard_cut_part.h"
#include <iostream>
#include <cassert>
#ifndef BH_PUBLIC
#include "rec_tree.h"
#include "partial_order.h"
#endif
#include "cut_Darren.h"
#include "settings.h"
#include "settings_reader.h"
#ifndef BH_PUBLIC
#include "cut_part_normal.h"
#endif
#include "cut_part_worker.h"
#include "integrals_ep.h"

#include "BH_debug.h"
#include "timing.h"

#define _VERBOSE 0 // 0 is silent, 1 gives some information
#define _CONSERVATIVE_ERROR_ESTIMATE 1 //Set this to 1 if you want to take into account the actual error of the coefficeint divided by the tree

namespace BH {

namespace cut {

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};
template <class BOX,class TRI,class BUB>  standard_cut_part<BOX,TRI,BUB>::~standard_cut_part(){
        for_each(d_boxes.begin(),d_boxes.end(),do_delete<BOX>());
        for_each(d_triangles.begin(),d_triangles.end(),do_delete<TRI>());
        for_each(d_bubbles.begin(),d_bubbles.end(),do_delete<BUB>());
}

template <class BOX,class TRI,class BUB> template <class T> SeriesC<T> standard_cut_part<BOX,TRI,BUB>::eval_without_check(momentum_configuration<T>& mc,const vector<int>& ind){
        SeriesC<T> intval(-2,0);

        // Add up the bubbles
        SeriesC<T> totalbubble(-2,0);
        complex<T> coeffbubble(0.,0.);
	
	double curacc(to_double(MaxDigits<T>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	std::complex<T> Ctreeres=get_tree(mc,ind); // Get the absolute power of the tree
	T treeres=abs(Ctreeres); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<T>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=T(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(T(10)));
	}
#endif

#if _VERBOSE==1
	complex<T> treeres_full=get_tree(mc,ind);
	if(abs(treeres_full)<BH::DeltaZero<T>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=T(1);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres_full);}

	_MESSAGE("");
	_MESSAGE("**************************** Cut Part NO IR test (mc) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

        int mu_index=get_mu_index<T>();
        if (mu_index == 0){
                set_mu_index<T>( DefineMu<T>(mc,m_MU));
                mu_index=get_mu_index<T>();
        }

         for (size_t i=1;i<=nbr_bubbles();i++){
                 BUB* thebubble=bubble(i);
               intval=*d_bubble_integrals[i-1]->get_value(mc,ind,mu_index);

               coeffbubble = thebubble->eval(mc,ind);
#if _VERBOSE==1
		_MESSAGE7(*thebubble,"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy()+acc_shift);
#endif
			 double thisacc=thebubble->get_accuracy();
			 if(thisacc+acc_shift<curacc){
				 curacc=thisacc+acc_shift;
			 }
			 
			 totalbubble += coeffbubble*intval;
          }
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
    _MESSAGE2("      Sum of all bubbles(n)=",totalbubble/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
         SeriesC<T> totaltri(-2,0);
         complex<T> coefftri(0.,0.);

         for (size_t i=1;i<=nbr_triangles();i++){

               //  Call Int in integrals.cpp

               intval=*d_triangle_integrals[i-1]->get_value(mc,ind,mu_index);

               // The contribution to the coefficient divided by tree
               coefftri = triangle(i)->eval(mc,ind);
#if _VERBOSE==1
		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
#endif
               totaltri += coefftri*intval;
          }

#if _VERBOSE==1
//	_MESSAGE2("      Sum of all triangles=",totaltri);
    _MESSAGE2("      Sum of all triangles(n)=",totaltri/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_boxes(), " Boxes");
#endif
        // Now add up the boxes

         SeriesC<T> totalbox(-2,0);
         complex<T> coeffbox;

         for (unsigned int i=1;i<=nbr_boxes();i++){

               intval=*d_box_integrals[i-1]->get_value(mc,ind,mu_index);

               coeffbox = box(i)->eval(mc,ind);
#if _VERBOSE==1
		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
#endif
               totalbox += coeffbox*intval;
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

	// Sum together the bubble, triangle and box results
	return totalbubble + totaltri + totalbox ;
}



template <class BOX,class TRI,class BUB> SeriesC<R> standard_cut_part<BOX,TRI,BUB>::eval_with_check_wCI(
	momentum_configuration<R>& mc,const std::vector<int>& ind){
	SeriesC<R> intval(-2,0);
	SeriesC<R> intval2(-2,0);
	// Add up the bubbles
	SeriesC<R> totalbubble(-2,0);
	SeriesC<R> totalbubble_conj(-2,0);
	complex<R> coeffbubble;

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<R>();
    if (mu_index == 0 || mc.n()<mu_index /*typically the case with a new mc*/){
            set_mu_index<R>( DefineMu<R>(mc,m_MU));
            mu_index=get_mu_index<R>();
    }
    eval_param<R> ep(mc,ind);

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

		intval=*d_bubble_integrals[i-1]->get_value(mc,ind,mu_index);

		coeffbubble = d_bubbles[i-1]->eval(ep);
		double accuracy=d_bubbles[i-1]->get_accuracy()+acc_shift;
		if(accuracy < settings::general::s_bub_cut_precision && ! isnan(accuracy)){
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
		intval=*d_triangle_integrals[i-1]->get_value(mc,ind,mu_index);

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

		intval=*d_box_integrals[i-1]->get_value(mc,ind,mu_index);

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

	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}




template <class BOX,class TRI,class BUB> SeriesC<RHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check_wCI(
		momentum_configuration<RHP>& mc,const std::vector<int>& ind){
	SeriesC<RHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RHP> totalbubble(-2,0);
	SeriesC<RHP> totalbubble_conj(-2,0);
	complex<RHP> coeffbubble;

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<RHP>();
    if (mu_index == 0 || mc.n()<mu_index /*typically the case with a new mc*/){
            set_mu_index<RHP>( DefineMu<RHP>(mc,m_MU));
            mu_index=get_mu_index<RHP>();
    }
    eval_param<RHP> ep(mc,ind);

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

		acc_shift=to_double(log(treeres)/log(RHP(10.)));
	}
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

		intval=*d_bubble_integrals[i-1]->get_value(mc,ind,mu_index);

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

		intval=*d_triangle_integrals[i-1]->get_value(mc,ind,mu_index);

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

		intval=*d_box_integrals[i-1]->get_value(mc,ind,mu_index);

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


	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}




template <class BOX,class TRI,class BUB> SeriesC<RVHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check_wCI(momentum_configuration<RVHP>& mc,const std::vector<int>& ind){
	SeriesC<RVHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RVHP> totalbubble(-2,0);
	SeriesC<RVHP> totalbubble_conj(-2,0);
	complex<RVHP> coeffbubble;

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<RVHP>();
    if (mu_index == 0 || mc.n()<mu_index /*typically the case with a new mc*/){
            set_mu_index<RVHP>( DefineMu<RVHP>(mc,m_MU));
            mu_index=get_mu_index<RVHP>();
    }
    eval_param<RVHP> ep(mc,ind);


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

		intval=*d_bubble_integrals[i-1]->get_value(mc,ind,mu_index);

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
		intval=*d_triangle_integrals[i-1]->get_value(mc,ind,mu_index);

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
		intval=*d_box_integrals[i-1]->get_value(mc,ind,mu_index);

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

	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}
#if BH_USE_GMP
std::vector<Cmom<RGMP> >* extend_momentaGMP(const eval_param<RGMP>& ep);

template <class BOX,class TRI,class BUB> SeriesC<RGMP> standard_cut_part<BOX,TRI,BUB>::eval_with_check_wCI(
		momentum_configuration<RGMP>& mc,const std::vector<int>& ind){
	SeriesC<RGMP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RGMP> totalbubble(-2,0);
	SeriesC<RGMP> totalbubble_conj(-2,0);
	complex<RGMP> coeffbubble;

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<RGMP>();
    if (mu_index == 0 || mc.n()<mu_index /*typically the case with a new mc*/){
            set_mu_index<RGMP>( DefineMu<RGMP>(mc,m_MU));
            mu_index=get_mu_index<RGMP>();
    }
    eval_param<RGMP> ep(mc,ind);

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
		BH_START_TIMER(bubble_integral)
		intval=*d_bubble_integrals[i-1]->get_value(mc,ind,mu_index);
		BH_STOP_TIMER(bubble_integral)


		BH_START_TIMER(bubble_coefficient)
		coeffbubble = d_bubbles[i-1]->eval(ep);
		BH_STOP_TIMER(bubble_coefficient)

		double accuracy=d_bubbles[i-1]->get_accuracy()+acc_shift;
		bool enough= ( accuracy > settings::general::s_bub_cut_precision) || ( isnan(accuracy));
		int oldPrecision=RGMP::get_current_precision();
		while (!enough){
			BH_DEBUG_MESSAGE2("WARNING : Recomputing bubble at GMP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);\
			RGMP::set_precision(RGMP::get_current_precision()+8);

			// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
			vector<Cmom<RGMP> > *new_mom=extend_momentaGMP(ep);
			eval_param<RGMP> epHP(*new_mom);
			coeffbubble=(d_bubbles[i-1]->eval(epHP));
			delete new_mom; //delete the new momentum vector created in the extend function
			double accuracy=d_bubbles[i-1]->get_accuracy()+acc_shift;
			enough= accuracy > settings::general::s_bub_cut_precision || (isnan(accuracy));
			BH_DEBUG(
			if(enough){
				_MESSAGE5("Managed to get ",accuracy," digits, target was: ",settings::general::s_bub_cut_precision," digits of precision.");
			}
			)
		}
		RGMP::set_precision(oldPrecision);
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
	SeriesC<RGMP> totaltri(-2,0);
	SeriesC<RGMP> totaltri_conj(-2,0);
	complex<RGMP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

		BH_START_TIMER(triangle_integral)
		intval=*d_triangle_integrals[i-1]->get_value(mc,ind,mu_index);
		BH_STOP_TIMER(triangle_integral)

		// The contribution to the coefficient divided by tree
		BH_START_TIMER(triangle_coefficient)
		coefftri = d_triangles[i-1]->eval(ep);
		BH_STOP_TIMER(triangle_coefficient)
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

		BH_START_TIMER(box_integral)
		intval=*d_box_integrals[i-1]->get_value(mc,ind,mu_index);
		BH_STOP_TIMER(box_integral)

		BH_START_TIMER(box_coefficient)
		coeffbox = d_boxes[i-1]->eval(ep);
		BH_STOP_TIMER(box_coefficient)
		totalbox += coeffbox*intval;
		totalbox_conj += conj(coeffbox)*intval;
	}
#if _VERBOSE
	_PRINT(totalbox);
#endif

	_accuracy=curacc; // Store the accuracy for this part of the computation
    SeriesC<RGMP> conjugate_cut_part=totalbubble_conj + totaltri_conj + totalbox_conj;
    _conjugate_cut_part=to_double(conjugate_cut_part);


	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}


#endif



template <class BOX,class TRI,class BUB> SeriesC<R> standard_cut_part<BOX,TRI,BUB>::eval_with_check(momentum_configuration<R>& mc,const vector<int>& ind){
	SeriesC<R> intval(-2,0);
	// Add up the bubbles
	SeriesC<R> totalbubble(-2,0);
	complex<R> coeffbubble;

	double curacc(MaxDigits<R>()); // keep the lowest recorded accuracy

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<R>();
    if (mu_index == 0){
            set_mu_index<R>( DefineMu<R>(mc,m_MU));
            mu_index=get_mu_index<R>();
    }

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	R treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<R>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=R(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(10.));
	}
#endif

#if _VERBOSE==1
	complex<R> treeres_full=get_tree(mc,ind);
	if(abs(treeres_full)<BH::DeltaZero<R>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=R(1);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres_full);}
	_MESSAGE("");
	_MESSAGE("**************************** Cut Part WITH IR test (mc) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

        vector<int> cc1,cc2;
         for (int j=1;j<=thebubble->corner_size(1);j++){
                     cc1.push_back(ind[thebubble->corner_ind(1,j)-1]);
                 }
         for (int j=1;j<=thebubble->corner_size(2);j++){
                     cc2.push_back(ind[thebubble->corner_ind(2,j)-1]);
                 }

         intval = Int(mc,mu_index,cc1,cc2);

		coeffbubble = d_bubbles[i-1]->eval(mc,ind);

		if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE6("WARNING : Recomputing ",*(d_bubbles[i-1])," = ",coeffbubble/treeres_full," bubble at HP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_HP mcHP=mc.extend<RHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_HP;
			int old_mu_index_HP=this->get_mu_index<RHP>();
			new_mu_index_HP = DefineMu<RHP>(mcHP,RHP(mc.mom(this->get_mu_index<R>()).E().real()));
			this->set_mu_HP(new_mu_index_HP);
			complex<RHP> HP_result=d_bubbles[i-1]->eval(mcHP,new_ind);
			coeffbubble= to_double(HP_result);
			this->set_mu_HP(old_mu_index_HP);

			//Check if this was accurate enough
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE6("WARNING : Recomputing ",*(d_bubbles[i-1])," = ",coeffbubble/treeres_full," bubble at VHP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
				// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
				vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
				mom_conf_VHP mcVHP=mc.extend<RVHP>(ind);
				// if d_index_mu is not set, we use the values from the multi_precision_constant,
				// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
				int new_mu_index_VHP;
				int old_mu_index_VHP=this->get_mu_index<RVHP>();
				new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(this->get_mu_index<R>()).E().real()));
				this->set_mu_VHP(new_mu_index_VHP);
				complex<RVHP> VHP_result=d_bubbles[i-1]->eval(mcVHP,new_ind);
				coeffbubble= to_double(VHP_result);
				this->set_mu_VHP(old_mu_index_VHP);

				//If this was to fail then there is not much more we can do
#if _VERBOSE==1
				if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
					_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
				}
#endif
			}
		}

#if _VERBOSE==1
//		_MESSAGE7(d_bubbles[i-1],"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy());
//	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble*intval/treeres," estimated accuracy=",d_bubbles[i-1]->get_accuracy());
	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble/treeres_full," estimated accuracy=",thebubble->get_accuracy());
#endif

		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
		totalbubble += coeffbubble*intval;
	}
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
	_MESSAGE2("      Sum of all bubbles(n)=",totalbubble/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
	SeriesC<R> totaltri(-2,0);
	complex<R> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

        //  Call Int in integrals.cpp
        vector<int> cc1,cc2,cc3;
        for (int j=1;j<=triangle(i)->corner_size(1);j++){
                    cc1.push_back(ind[triangle(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(2);j++){
                    cc2.push_back(ind[triangle(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(3);j++){
                    cc3.push_back(ind[triangle(i)->corner_ind(3,j)-1]);
                }

        intval = Int(mc,mu_index,cc1,cc2,cc3);

        // The contribution to the coefficient divided by tree
        coefftri = triangle(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
//		_MESSAGE3(*d_triangles[i-1],"=",coefftri*intval/treeres);
		_MESSAGE3(*d_triangles[i-1],"=",coefftri/treeres_full);
#endif
		totaltri += coefftri*intval;
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
	complex<R> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
        vector<int> dd1,dd2,dd3,dd4;
        for (int j=1;j<=box(i)->corner_size(1);j++){
                    dd1.push_back(ind[box(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(2);j++){
                    dd2.push_back(ind[box(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(3);j++){
                    dd3.push_back(ind[box(i)->corner_ind(3,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(4);j++){
                    dd4.push_back(ind[box(i)->corner_ind(4,j)-1]);
                }

        intval = Int(mc,mu_index,dd1,dd2,dd3,dd4);

        coeffbox = box(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
//	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox*intval/treeres);
	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox/treeres_full);
#endif
		totalbox += coeffbox*intval;
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
	
	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}

template <class BOX,class TRI,class BUB> SeriesC<RHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(momentum_configuration<RHP>& mc,const vector<int>& ind){
	SeriesC<RHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RHP> totalbubble(-2,0);
	complex<RHP> coeffbubble;
	
	double curacc(to_double(MaxDigits<RHP>())); // keep the lowest recorded accuracy

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<RHP>();
    if (mu_index == 0){
            set_mu_index<RHP>( DefineMu<RHP>(mc,m_MU));
            mu_index=get_mu_index<RHP>();
    }

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	RHP treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<RHP>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=RHP(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(RHP(10)));
	}
#endif

#if _VERBOSE==1
	complex<RHP> treeres_full=get_tree(mc,ind);
	if(abs(treeres_full)<BH::DeltaZero<RHP>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=complex<RHP>(1,0);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres);}
	_MESSAGE("");
	_MESSAGE("**************************** Cut Part WITH IR test HP (mc) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

        vector<int> cc1,cc2;
         for (int j=1;j<=thebubble->corner_size(1);j++){
                     cc1.push_back(ind[thebubble->corner_ind(1,j)-1]);
                 }
         for (int j=1;j<=thebubble->corner_size(2);j++){
                     cc2.push_back(ind[thebubble->corner_ind(2,j)-1]);
                 }

         intval = Int(mc,mu_index,cc1,cc2);

		coeffbubble = d_bubbles[i-1]->eval(mc,ind);

		if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
#if _VERBOSE==1
			_MESSAGE6("WARNING : Recomputing ",*(d_bubbles[i-1])," = ",coeffbubble/treeres_full," bubble at VHP as accuracy is ",d_bubbles[i-1]->get_accuracy()+acc_shift);
#endif
			vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
			mom_conf_VHP mcVHP=mc.extend<RVHP>(ind);
			// if d_index_mu is not set, we use the values from the multi_precision_constant,
			// otherwise, we construct a mu vector in the extended mc corresponding to the energy component of the vector d_mu+index in the NP mc.
			int new_mu_index_VHP;
			int old_mu_index_VHP=this->get_mu_index<RVHP>();
			new_mu_index_VHP = DefineMu<RVHP>(mcVHP,RVHP(mc.mom(this->get_mu_index<RHP>()).E().real()));
			this->set_mu_VHP(new_mu_index_VHP);
			complex<RVHP> VHP_result=d_bubbles[i-1]->eval(mcVHP,new_ind);
			coeffbubble= to_HP(VHP_result);
			this->set_mu_VHP(old_mu_index_VHP);

			//Check if this was accurate enough
			if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
			//If this was to fail then there is not much more we can do
#if _VERBOSE==1
				_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
#endif
			}
		}
#if _VERBOSE==1
//		_MESSAGE7(d_bubbles[i-1],"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy());
//	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble*intval/treeres_full," estimated accuracy=",d_bubbles[i-1]->get_accuracy());
	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble/treeres_full," estimated accuracy=",thebubble->get_accuracy());
#endif
		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
		totalbubble += coeffbubble*intval;
	}
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
	_MESSAGE2("      Sum of all bubbles(n)=",totalbubble/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
	SeriesC<RHP> totaltri(-2,0);
	complex<RHP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

        //  Call Int in integrals.cpp
        vector<int> cc1,cc2,cc3;
        for (int j=1;j<=triangle(i)->corner_size(1);j++){
                    cc1.push_back(ind[triangle(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(2);j++){
                    cc2.push_back(ind[triangle(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(3);j++){
                    cc3.push_back(ind[triangle(i)->corner_ind(3,j)-1]);
                }

        intval = Int(mc,mu_index,cc1,cc2,cc3);

        // The contribution to the coefficient divided by tree
        coefftri = triangle(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
//		_MESSAGE3(*d_triangles[i-1],"=",coefftri*intval/treeres);
		_MESSAGE3(*d_triangles[i-1],"=",coefftri/treeres_full);
#endif
		totaltri += coefftri*intval;
	}

#if _VERBOSE==1
//	_MESSAGE2("      Sum of all triangles=",totaltri);
	_MESSAGE2("      Sum of all triangles(n)=",totaltri/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_boxes(), " Boxes");
#endif
	// Now add up the boxes

	SeriesC<RHP> totalbox(-2,0);
	complex<RHP> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
        vector<int> dd1,dd2,dd3,dd4;
        for (int j=1;j<=box(i)->corner_size(1);j++){
                    dd1.push_back(ind[box(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(2);j++){
                    dd2.push_back(ind[box(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(3);j++){
                    dd3.push_back(ind[box(i)->corner_ind(3,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(4);j++){
                    dd4.push_back(ind[box(i)->corner_ind(4,j)-1]);
                }

        intval = Int(mc,mu_index,dd1,dd2,dd3,dd4);

        coeffbox = box(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
//	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox*intval/treeres);
	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox/treeres_full);
#endif
		totalbox += coeffbox*intval;
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
	
	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}

template <class BOX,class TRI,class BUB> SeriesC<RVHP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(momentum_configuration<RVHP>& mc,const vector<int>& ind){
	SeriesC<RVHP> intval(-2,0);
	// Add up the bubbles
	SeriesC<RVHP> totalbubble(-2,0);
	complex<RVHP> coeffbubble;
	
	double curacc(to_double(MaxDigits<RVHP>())); // keep the lowest recorded accuracy

	//Get the correct mu scale for the integrals
    int mu_index=get_mu_index<RVHP>();
    if (mu_index == 0){
            set_mu_index<RVHP>( DefineMu<RVHP>(mc,m_MU));
            mu_index=get_mu_index<RVHP>();
    }

	// As the result will be proportional to the tree we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	RVHP treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<RVHP>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=RVHP(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(RVHP(10)));
	}
#endif

#if _VERBOSE==1
	complex<RVHP> treeres_full=get_tree(mc,ind);
	if(abs(treeres_full)<BH::DeltaZero<RVHP>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=complex<R>(1,0);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres);}
	_MESSAGE("");
	_MESSAGE("**************************** Cut Part WITH IR test VHP (mc) *****************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_bubbles(), " Bubbles");
#endif

	for (size_t i=1;i<=nbr_bubbles();i++){
		BUB* thebubble=bubble(i);

        vector<int> cc1,cc2;
         for (int j=1;j<=thebubble->corner_size(1);j++){
                     cc1.push_back(ind[thebubble->corner_ind(1,j)-1]);
                 }
         for (int j=1;j<=thebubble->corner_size(2);j++){
                     cc2.push_back(ind[thebubble->corner_ind(2,j)-1]);
                 }

         intval = Int(mc,mu_index,cc1,cc2);

		coeffbubble = d_bubbles[i-1]->eval(mc,ind);

#if _VERBOSE==1
		if(d_bubbles[i-1]->get_accuracy()+acc_shift<settings::general::s_bub_cut_precision){
			//If this was to fail then there is not much more we can do
			_MESSAGE3("WARNING : Failed to get ",settings::general::s_bub_cut_precision," digits of precision even at VHP accuracy.");
		}
#endif

#if _VERBOSE==1
//		_MESSAGE7(d_bubbles[i-1],"=",real(coeffbubble),imag(coeffbubble)>T(0)?"+I ":"-I ",abs(imag(coeffbubble))," estimated accuracy=",thebubble->get_accuracy());
//	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble*intval/treeres_full," estimated accuracy=",d_bubbles[i-1]->get_accuracy());
	    _MESSAGE5(*(d_bubbles[i-1]),"=",coeffbubble/treeres_full," estimated accuracy=",thebubble->get_accuracy()+acc_shift);
#endif
		double thisacc=thebubble->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		BH_DEBUG_PRINT(*this);
		BH_DEBUG_PRINT(coeffbubble);
		BH_DEBUG_PRINT(intval);
		totalbubble += coeffbubble*intval;
	}
#if _VERBOSE==1
//	_MESSAGE2("      Sum of all bubbles=",totalbubble);
	_MESSAGE2("      Sum of all bubbles(n)=",totalbubble/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_triangles(), " Triangles");
#endif
	SeriesC<RVHP> totaltri(-2,0);
	complex<RVHP> coefftri;

	for (size_t i=1;i<=nbr_triangles();i++){

        //  Call Int in integrals.cpp
        vector<int> cc1,cc2,cc3;
        for (int j=1;j<=triangle(i)->corner_size(1);j++){
                    cc1.push_back(ind[triangle(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(2);j++){
                    cc2.push_back(ind[triangle(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=triangle(i)->corner_size(3);j++){
                    cc3.push_back(ind[triangle(i)->corner_ind(3,j)-1]);
                }

        intval = Int(mc,mu_index,cc1,cc2,cc3);

        // The contribution to the coefficient divided by tree
        coefftri = triangle(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_triangles[i-1],"=",real(coefftri),imag(coefftri)>T(0)?"+I ":"-I ",abs(imag(coefftri)));
//		_MESSAGE3(*d_triangles[i-1],"=",coefftri*intval/treeres);
		_MESSAGE3(*d_triangles[i-1],"=",coefftri/treeres_full);
#endif
		totaltri += coefftri*intval;
	}

#if _VERBOSE==1
//	_MESSAGE2("      Sum of all triangles=",totaltri);
	_MESSAGE2("      Sum of all triangles(n)=",totaltri/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",nbr_boxes(), " Boxes");
#endif
	// Now add up the boxes

	SeriesC<RVHP> totalbox(-2,0);
	complex<RVHP> coeffbox;

	for (unsigned int i=1;i<=nbr_boxes();i++){

		//  Here we call the integrals .
        vector<int> dd1,dd2,dd3,dd4;
        for (int j=1;j<=box(i)->corner_size(1);j++){
                    dd1.push_back(ind[box(i)->corner_ind(1,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(2);j++){
                    dd2.push_back(ind[box(i)->corner_ind(2,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(3);j++){
                    dd3.push_back(ind[box(i)->corner_ind(3,j)-1]);
                }
        for (int j=1;j<=box(i)->corner_size(4);j++){
                    dd4.push_back(ind[box(i)->corner_ind(4,j)-1]);
                }

        intval = Int(mc,mu_index,dd1,dd2,dd3,dd4);
        BH_DEBUG_MESSAGE4("intval: ",intval," cached: ",d_box_integrals[i]->get_value(mc,ind,mu_index));

        coeffbox = box(i)->eval(mc,ind);
#if _VERBOSE==1
//		_MESSAGE5(*d_boxes[i-1],"=",real(coeffbox),imag(coeffbox)>T(0)?"+I ":"-I ",abs(imag(coeffbox)));
//	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox*intval/treeres_full);
	    _MESSAGE3(*d_boxes[i-1],"=",coeffbox/treeres_full);
#endif
		totalbox += coeffbox*intval;
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
	
	// Sum together the bubble triangle and box results
	return totalbubble+totaltri+totalbox;
}
#if BH_USE_GMP
template <class BOX,class TRI,class BUB> SeriesC<RGMP> standard_cut_part<BOX,TRI,BUB>::eval_with_check(momentum_configuration<RGMP>& mc,const vector<int>& ind){
	return eval_with_check_wCI(mc,ind);
}
#endif
//template <class BOX,class TRI,class BUB> template <class T> SeriesC<T> standard_cut_part<BOX,TRI,BUB>::eval_fn(momentum_configuration<T>& mc,const vector<int>& ind){
//
//	if(settings::general::s_use_ep_only){// If we only want to use the eval_param code
//		// If the scale is stored as a component in a variable then we need to create
//		//  a multi_precison_constant instead.
//	    int mu_index=get_mu_index<R>();
//	    if(mu_index!=0){
//	    	multi_precision_constant mmu(mc.mom(mu_index).E().real(),mc.mom(get_mu_index<RHP>()).E().real(),mc.mom(get_mu_index<RVHP>()).E().real());
//	    	set_mu(mmu);
//	    }
//
//	    // Create an eval_param from the mom_conf, do this in a dangerous way for now
//		eval_param<T> ep(mc,ind);
//
//		if(settings::general::s_use_check_in_cut_part){
//			return eval_with_check(ep);
//		}else{
//			return eval_without_check(ep);
//		}
//	}
//	else{
//		if(settings::general::s_use_check_in_cut_part){
//			return eval_with_check(mc,ind);
//		}else{
//			return eval_without_check(mc,ind);
//		}
//	}
//}

template <class BOX,class TRI,class BUB> SeriesC<R> standard_cut_part<BOX,TRI,BUB>::eval_fn(momentum_configuration<R>& mc,const vector<int>& ind){

	if(settings::general::s_use_ep_only){// If we only want to use the eval_param code
		// If the scale is stored as a component in a variable then we need to create
		//  a multi_precison_constant instead.
	    int mu_index=this->get_mu_index<R>();
	    if(mu_index!=0){
	    	R muvalue=mc.mom(mu_index).E().real();
	    	multi_precision_constant mmu(muvalue);
	    	set_mu(mmu);
	    }

	    // Create an eval_param from the mom_conf, do this in a dangerous way for now
		eval_param<R> ep(mc,ind);

		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check(ep);
		}else{
			return eval_without_check(ep);
		}
	}
	else{
		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check_wCI(mc,ind);
		}else{
			return eval_without_check(mc,ind);
		}
	}
}

template <class BOX,class TRI,class BUB> SeriesC<RHP> standard_cut_part<BOX,TRI,BUB>::eval_fn(momentum_configuration<RHP>& mc,const vector<int>& ind){

	if(settings::general::s_use_ep_only){// If we only want to use the eval_param code
		// If the scale is stored as a component in a variable then we need to create
		//  a multi_precison_constant instead.
	    int mu_index=this->get_mu_index<RHP>();
	    if(mu_index!=0){
	    	RHP muvalue=mc.mom(mu_index).E().real();
	    	multi_precision_constant mmu(to_double(muvalue));
	    	set_mu(mmu);
	    }

	    // Create an eval_param from the mom_conf, do this in a dangerous way for now
		eval_param<RHP> ep(mc,ind);

		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check(ep);
		}else{
			return eval_without_check(ep);
		}
	}
	else{
		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check_wCI(mc,ind);
		}else{
			return eval_without_check(mc,ind);
		}
	}
}

template <class BOX,class TRI,class BUB> SeriesC<RVHP> standard_cut_part<BOX,TRI,BUB>::eval_fn(momentum_configuration<RVHP>& mc,const vector<int>& ind){

	if(settings::general::s_use_ep_only){// If we only want to use the eval_param code
		// If the scale is stored as a component in a variable then we need to create
		//  a multi_precison_constant instead.
	    int mu_index=this->get_mu_index<RVHP>();
	    if(mu_index!=0){
	    	RVHP muvalue=mc.mom(mu_index).E().real();
	    	multi_precision_constant mmu(muvalue);
	    	set_mu(mmu);
	    }

	    // Create an eval_param from the mom_conf, do this in a dangerous way for now
		eval_param<RVHP> ep(mc,ind);

		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check(ep);
		}else{
			return eval_without_check(ep);
		}
	}
	else{
		if(settings::general::s_use_check_in_cut_part){
			return eval_with_check_wCI(mc,ind);
		}else{
			return eval_without_check(mc,ind);
		}
	}
}
#if BH_USE_GMP
template <class BOX,class TRI,class BUB> SeriesC<RGMP> standard_cut_part<BOX,TRI,BUB>::eval_fn(momentum_configuration<RGMP>& mc,const vector<int>& ind){

	if(settings::general::s_use_ep_only){// If we only want to use the eval_param code
		_WARNING("Sorry, not ready for that..."); throw;
//		// If the scale is stored as a component in a variable then we need to create
//		//  a multi_precison_constant instead.
//	    int mu_index=this->get_mu_index<RVHP>();
//	    if(mu_index!=0){
//	    	RVHP muvalue=mc.mom(mu_index).E().real();
//	    	multi_precision_constant mmu(muvalue);
//	    	set_mu(mmu);
//	    }
//
//	    // Create an eval_param from the mom_conf, do this in a dangerous way for now
//		eval_param<RVHP> ep(mc,ind);
//
//		if(settings::general::s_use_check_in_cut_part){
//			return eval_with_check(ep);
//		}else{
//			return eval_without_check(ep);
//		}
	}
	else{
		if(settings::general::s_use_check_in_cut_part){
//			_WARNING("Sorry, not ready for that..."); throw;
			return eval_with_check(mc,ind);
		}else{
			return eval_without_check(mc,ind);
		}
	}
}
#endif

template <class BOX,class TRI,class BUB> template <class T> SeriesC<T> standard_cut_part<BOX,TRI,BUB>::eval_fn(const eval_param<T>& ep){
	if(settings::general::s_use_check_in_cut_part){
		return eval_with_check(ep);
	}else{
		return eval_without_check(ep);
}
}

#ifndef BH_PUBLIC
template <class BOX,class TRI,class BUB> template <class OrigType> standard_cut_part<BOX,TRI,BUB>::standard_cut_part(OrigType& orig, option* opt,cutD_factory* cf) : Cut_Part_base( orig.get_process() ), _accuracy(0) {
	TREE_FACTORY_TYPE TF;
//	_tree_ptr=TF.new_tree(fix_flavors(orig.get_process()));

	vector<size_t> box_map;
	vector<size_t> tri_map;
	vector<size_t> bub_map;
	map<long,size_t> bub_to_ind;
	map<long,size_t> tri_to_ind;
	map<long,size_t> box_to_ind;

	for (size_t bub=1; bub<=orig.nbr_bubbles();bub++){
		if ( (*opt)(orig.bubble(bub)) ) {
			bubbleD* newbub=cf->new_bubble(*orig.bubble(bub));
			newbub->clear_daughters();
			d_bubbles.push_back(newbub);
			bub_map.push_back(d_bubbles.size());
			bub_to_ind.insert(pair<long,size_t>(orig.bubble(bub)->get_ID(),d_bubbles.size()));
		}
		else{
			bub_map.push_back(0);
			bub_to_ind.insert(pair<long,size_t>(orig.bubble(bub)->get_ID(),0));
		}
	};

	for (size_t tri=1; tri<=orig.nbr_triangles();tri++){
		if ( (*opt)(orig.triangle(tri)) ) {
			triangleD* newtri=cf->new_triangle(*orig.triangle(tri));
			newtri->clear_daughters();
			newtri->clear_parents();
			d_triangles.push_back(newtri);
			tri_map.push_back(d_triangles.size());
			tri_to_ind.insert(pair<long,size_t>(orig.triangle(tri)->get_ID(),d_triangles.size()));
		}
		else{
			tri_to_ind.insert(pair<long,size_t>(orig.triangle(tri)->get_ID(),0));
			tri_map.push_back(0);
		}
	};

	for (size_t box=1; box<=orig.nbr_boxes();box++){
		if ( (*opt)(orig.box(box)) ) {
			boxD* newbox=cf->new_box(*orig.box(box));
			newbox->clear_parents();
			d_boxes.push_back(newbox);
			box_map.push_back(d_boxes.size());
			box_to_ind.insert(pair<long,size_t>(orig.box(box)->get_ID(),d_boxes.size()));
		}
		else{
			box_map.push_back(0);
		}
	}

	for (size_t bub=1;bub<=orig.nbr_bubbles();bub++){
		if (bub_map[bub-1]!=0){
			bubbleD* old_bubble=orig.bubble(bub);
			bubbleD* new_bubble=d_bubbles[bub_map[bub-1]-1];
			for (size_t dau=1;dau<=orig.bubble(bub)->daughters_nbr();dau++){
				triangleD* old_triangle=old_bubble->get_daughter(dau);
				map<long,size_t>::iterator pos=tri_to_ind.find(old_triangle->get_ID());
				if (pos==tri_to_ind.end()){throw;}
				size_t new_triangle_index= (pos->second);
				if (new_triangle_index!=0){
					triangleD* new_triangle=d_triangles[new_triangle_index-1]	;
					new_bubble->add_daughter(new_triangle,old_bubble->get_opened_corner(dau));
					new_triangle->add_parent(new_bubble,old_bubble->get_opened_corner(dau));
				}
			}
		}
	}

	//	boxes[box]->get_parent(par)->add_daughter(boxes[box],boxes[box]->get_closed_corner(par));

	for (size_t tri=1;tri<=orig.nbr_triangles();tri++){
		if (tri_map[tri-1]!=0){
			triangleD* old_triangle=orig.triangle(tri);
			triangleD* new_triangle=d_triangles[tri_map[tri-1]-1];
			for (size_t dau=1;dau<=orig.triangle(tri)->daughters_nbr();dau++){
				boxD* old_box=old_triangle->get_daughter(dau);
				size_t new_box_index= (box_to_ind.find(old_box->get_ID())->second);
				if (new_box_index!=0){
					boxD* new_box=d_boxes[new_box_index-1]	;
					new_triangle->add_daughter(new_box,old_triangle->get_opened_corner(dau));
					new_box->add_parent(new_triangle,old_triangle->get_opened_corner(dau));
				}
			}
		}
	}
	for (size_t box=1;box<=orig.nbr_boxes();box++){
		if (box_map[box-1]!=0){
			boxD* old_box=orig.box(box);
			boxD* new_box=d_boxes[box_map[box-1]-1];
			for (size_t par=1;par<=orig.box(box)->parents_nbr();par++){
				triangleD* old_triangle=old_box->get_parent(par);
				size_t new_triangle_index= (tri_to_ind.find(old_triangle->get_ID())->second);
				if (new_triangle_index!=0){
					triangleD* new_triangle=d_triangles[new_triangle_index-1]	;
					for (int new_par=1;new_par<=new_box->parents_nbr();new_par++){
						if ( new_box->get_parent(new_par)->get_ID() == old_triangle->get_ID() ){
							new_box->set_closed_corner(new_par,old_box->get_closed_corner(par));
//							new_triangle->set_opened_corner(tri_map[tri-1],old_triangle->get_opened_corner(dau));

						}
					}

				}
			}
		}
	}

	for (size_t tri=1;tri<=orig.nbr_triangles();tri++){
		if (tri_map[tri-1]!=0){
			triangleD* old_tri=orig.triangle(tri);
			triangleD* new_tri=d_triangles[tri_map[tri-1]-1];
			for (size_t par=1;par<=orig.triangle(tri)->parents_nbr();par++){
				bubbleD* old_bubble=old_tri->get_parent(par);
				size_t new_bubble_index= (bub_to_ind.find(old_bubble->get_ID())->second);
				if (new_bubble_index!=0){
					bubbleD* new_bubble=d_bubbles[new_bubble_index-1]	;
					for (int new_par=1;new_par<=new_tri->parents_nbr();new_par++){
						if ( new_tri->get_parent(new_par)->get_ID() == old_bubble->get_ID() ){
							new_tri->set_closed_corner(new_par,old_tri->get_closed_corner(par));

						}
					}

				}
			}
			for (size_t dau=1;dau<=orig.triangle(tri)->daughters_nbr();dau++){
				boxD* old_box=old_tri->get_daughter(dau);
				size_t new_box_index= (box_to_ind.find(old_box->get_ID())->second);
				if (new_box_index!=0){
					boxD* new_box=d_boxes[new_box_index-1]	;
					for (int new_dau=1;new_dau<=new_tri->daughters_nbr();new_dau++){
						if ( new_tri->get_daughter(new_dau)->get_ID() == old_box->get_ID() ){
							new_tri->set_opened_corner(new_dau,old_tri->get_opened_corner(dau));

						}
					}

				}
			}

		}
	}

	for (size_t bub=1;bub<=orig.nbr_bubbles();bub++){
		if (bub_map[bub-1]!=0){
			bubbleD* old_bub=orig.bubble(bub);
			bubbleD* new_bub=d_bubbles[bub_map[bub-1]-1];
			for (size_t dau=1;dau<=orig.bubble(bub)->daughters_nbr();dau++){
				triangleD* old_tri=old_bub->get_daughter(dau);
				size_t new_tri_index= (tri_to_ind.find(old_tri->get_ID())->second);
				if (new_tri_index!=0){
					triangleD* new_tri=d_triangles[new_tri_index-1]	;
					for (int new_dau=1;new_dau<=new_bub->daughters_nbr();new_dau++){
						if ( new_bub->get_daughter(new_dau)->get_ID() == old_tri->get_ID() ){
							new_bub->set_opened_corner(new_dau,old_bub->get_opened_corner(dau));

						}
					}

				}
			}

		}
	}

//	for (int i=0;i<d_bubbles.size();i++){
//		eval_bubbles.push_back(new zero_checked_computable<std::complex>(d_bubbles[i]));
//	}
//	for (int i=0;i<d_triangles.size();i++){
//		eval_triangles.push_back(new zero_checked_computable<std::complex>(d_triangles[i]));
//	}
//	for (int i=0;i<d_boxes.size();i++){
//		eval_boxes.push_back(new zero_checked_computable<std::complex>(d_boxes[i]));
//	}

	construct_integral_table();




}

#endif

template <class BOX,class TRI,class BUB> void standard_cut_part<BOX,TRI,BUB>::construct_integral_table(){
	for (typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();++it){
		vector<int> cor[2];
		for (int c=1 ; c <=2 ; c++){
			for (int j=1;j<=(*it)->corner_size(c);j++){
				cor[c-1].push_back((*it)->corner_ind(c,j));
			}
		}
		d_bubble_integrals.push_back(new CachedIntegral::Cached_Bubble_Integral_User(cor[0],cor[1]));
	}
	for (typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();++it){
		vector<int> cor[3];
		for (int c=1 ; c <=3 ; c++){
			for (int j=1;j<=(*it)->corner_size(c);j++){
				cor[c-1].push_back((*it)->corner_ind(c,j));
			}
		}
		d_triangle_integrals.push_back(new CachedIntegral::Cached_Triangle_Integral_User(cor[0],cor[1],cor[2]));
	}
	for (typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();++it){
		vector<int> cor[4];
		for (int c=1 ; c <=4 ; c++){
			for (int j=1;j<=(*it)->corner_size(c);j++){
				cor[c-1].push_back((*it)->corner_ind(c,j));
			}
		}
		d_box_integrals.push_back(new CachedIntegral::Cached_Box_Integral_User(cor[0],cor[1],cor[2],cor[3]));
	}
}

}
}

#endif /* STANDARD_CUT_PART_HPP_ */
