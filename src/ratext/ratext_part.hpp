/*
 * ratext_part.cpp
 *
 *  Created on: 17-May-2009
 *      Author: daniel
 */

#include "ratext_part.h"
#include "worker_utils.h"  // for read_process_from_stream
#include <iostream>
#include <string>
#include <cassert>
#include "OneLoopHelAmpl.h"
#include <fstream>
#include <typeinfo>
#ifndef BH_PUBLIC
#include "ratext/data_files.h"
#endif
#include "ratext/filename.h"
#include "BH_typedefs.h"
#include "known_rational.h"
#include "BH_debug.h"

#define _VERBOSE 0 
#define _CONSERVATIVE_ERROR_ESTIMATE 1 //Set this to 1 if you want to take into account the actual error of the coefficient divided by the tree


using BH::worker::read_process_from_stream;

using namespace std;

namespace BH {

// In rat_ext_fac.cpp
bool process_has_helicity_conserving_quark_pair(const process& pro);


namespace ratext {


template <class PENT,class BOX,class TRI,class BUB> template <class T> complex<T> ratext_part<PENT,BOX,TRI,BUB>::eval_rat(momentum_configuration<T>& mc, const vector<int>& ind)
{
	complex<T> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0), kb_bub(0,0), kb_tri(0,0), kb_box(0,0), ksum, kstep;
	T clr_fac=this->template get_colour_fac<T>();
	
	double curacc(to_double(MaxDigits<T>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
	//  difference of the exponent of the tree result from 1
	double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
	T treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
	if(treeres > BH::DeltaZero<T>()){ // If the tree is zero then there is no error_shift
		if(treeres>1){treeres=T(1)/treeres;} // we turn all exponents into positive powers

		acc_shift=to_double(log(treeres)/log(T(10)));
	}
#endif

#if _VERBOSE==1
	_MESSAGE("");
	_MESSAGE("****************************** Rat Ext (mc) **********************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
#endif

	for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
		kstep=clr_fac*((*it)->eval(mc,ind))-kb_bub;
		ksum=resultbubsum+kstep;
		kb_bub=(ksum-resultbubsum)-kstep;
		resultbubsum=ksum;
		
		double thisacc=(*it)->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}
		
#if _VERBOSE==1
		_MESSAGE8("BUB : ",*(*it),"=",real(kstep),imag(kstep)>T(0)?"+I ":"-I ",abs(imag(kstep))," estimated acc=",thisacc+acc_shift);
#endif
	}

#if _VERBOSE==1
	_MESSAGE2("      Sum of all bubbles=",resultbubsum);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

	for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
		kstep=clr_fac*((*it)->eval(mc,ind))-kb_tri;
		ksum=resulttrisum+kstep;
		kb_tri=(ksum-resulttrisum)-kstep;
		resulttrisum=ksum;
#if _VERBOSE==1
		_MESSAGE6("TRI : ",*(*it),"=",real(kstep),imag(kstep)>T(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
	}

#if _VERBOSE==1
	_MESSAGE2("      Sum of all triangles=",resulttrisum);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

	for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
		kstep=clr_fac*((*it)->eval(mc,ind))-kb_box;
		ksum=resultboxsum+kstep;
		kb_box=(ksum-resultboxsum)-kstep;
		resultboxsum=ksum;
#if _VERBOSE==1
		_MESSAGE6("BOX : ",*(*it),"=",real(kstep),imag(kstep)>T(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
	}
#if _VERBOSE==1
	_MESSAGE2("      Sum of all boxes=",resultboxsum);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
	for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
		complex<T> resultpent=clr_fac*((*it)->eval(mc,ind));

		_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>T(0)?"+I ":"-I ",abs(imag(resultpent)));
	}
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
#endif

	_accuracy=curacc; // Store the accuracy for this part of the computation

	return resultbubsum+resulttrisum+resultboxsum;
}


template <class PENT,class BOX,class TRI,class BUB> complex<R> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_IR_checked(momentum_configuration<R>& mc, const vector<int>& ind)
{

	if ( ! settings::rational_settings::s_use_IR_in_ratext){
		return eval_rat(mc,ind);
	} else {
		complex<R> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0);
		R clr_fac=this->template get_colour_fac<R>();
		
		double curacc(MaxDigits<R>()); // keep the lowest recorded accuracy

		// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
		//  difference of the exponent of the tree result from 1
		double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
		R treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
		if(treeres > BH::DeltaZero<R>()){ // If the tree is zero then there is no error_shift
			if(treeres>1){treeres=R(1)/treeres;} // we turn all exponents into positive powers

			acc_shift=to_double(log(treeres)/log(10.));
		}
#endif

		multi_precision_reader* mpr=dynamic_cast<multi_precision_reader*>(&mc);

#if _VERBOSE==1
		_MESSAGE("");
		_MESSAGE("************************ Rat Ext with Check (mc) ******************************");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
#endif

		for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
			complex<R> resultbub=clr_fac*((*it)->eval(mc,ind));
			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at HP as accuracy is ",(*it)->get_accuracy()+acc_shift);
#endif
				if ( mpr ==0 ){
					complex<RHP> HP_result;
					vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
					mom_conf_HP mcHP=mc.extend<RHP>(ind);
					HP_result=((*it)->eval(mcHP,new_ind));
					resultbub=clr_fac*to_double(HP_result);
				} else {
					complex<RHP> HP_result=((*it)->eval(mpr->mc_HP(),ind));
					resultbub=clr_fac*to_double(HP_result);
				}
			}
			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",(*it)->get_accuracy()+acc_shift);
#endif
				// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
				complex<RVHP> VHP_result;
				vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
				mom_conf_VHP mcVHP=mc.extend<RVHP>(ind);
				VHP_result=((*it)->eval(mcVHP,new_ind));
				resultbub=clr_fac*to_double(VHP_result);
			}
#if _VERBOSE==1
			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
				_MESSAGE3("WARNING : Failed to get ",settings::rational_settings::s_rat_ext_precision," digits of precision even at VHP accuracy.");
			}
#endif
			double thisacc=(*it)->get_accuracy();
			if(thisacc+acc_shift<curacc){
				curacc=thisacc+acc_shift;
			}
			resultbubsum+=resultbub;
#if _VERBOSE==1
			_MESSAGE8("BUB : ",*(*it),"=",real(resultbub),imag(resultbub)>R(0)?"+I ":"-I ",abs(imag(resultbub))," estimated acc=",thisacc+acc_shift);
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all bubbles=",resultbubsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

		for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
#if _VERBOSE==1
			complex<R> resulttri=clr_fac*((*it)->eval(mc,ind));
			resulttrisum+=resulttri;

			_MESSAGE6("TRI : ",*(*it),"=",real(resulttri),imag(resulttri)>R(0)?"+I ":"-I ",abs(imag(resulttri)));
#else
			resulttrisum+=clr_fac*((*it)->eval(mc,ind));
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all triangles=",resulttrisum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

		for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
#if _VERBOSE==1
			complex<R> resultbox=clr_fac*((*it)->eval(mc,ind));
			resultboxsum+=resultbox;

			_MESSAGE6("BOX : ",*(*it),"=",real(resultbox),imag(resultbox)>R(0)?"+I ":"-I ",abs(imag(resultbox)));
#else
			resultboxsum+=clr_fac*((*it)->eval(mc,ind));
#endif
		}
#if _VERBOSE==1
		_MESSAGE2("      Sum of all boxes=",resultboxsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
		for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
			complex<R> resultpent=clr_fac*((*it)->eval(mc,ind));

			_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>R(0)?"+I ":"-I ",abs(imag(resultpent)));
		}
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
#endif

		_accuracy=curacc; // Store the accuracy for this part of the computation

		return resultbubsum+resulttrisum+resultboxsum;
	}
}






template <class PENT,class BOX,class TRI,class BUB> complex<RHP> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_IR_checked(momentum_configuration<RHP>& mc, const vector<int>& ind)
{
	if ( ! settings::rational_settings::s_use_IR_in_ratext){
		return eval_rat(mc,ind);
	} else {
		complex<RHP> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0);
		RHP clr_fac=this->template get_colour_fac<RHP>();
		
		double curacc(to_double(MaxDigits<RHP>())); // keep the lowest recorded accuracy

		// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
		//  difference of the exponent of the tree result from 1
		double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
		RHP treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
		if(treeres > BH::DeltaZero<RHP>()){ // If the tree is zero then there is no error_shift
			if(treeres>1){treeres=RHP(1)/treeres;} // we turn all exponents into positive powers

			acc_shift=to_double(log(treeres)/log(RHP(10.)));
		}
#endif

		multi_precision_reader_HP* mpr=dynamic_cast<multi_precision_reader_HP*>(&mc);

	#if _VERBOSE==1
		_MESSAGE("");
		_MESSAGE("*********************** Rat Ext with Check HP (mc) *****************************");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
	#endif

		for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
			complex<RHP> resultbub=clr_fac*((*it)->eval(mc,ind));
			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
	#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",(*it)->get_accuracy()+acc_shift);
	#endif
				if ( mpr ==0 ){
					complex<RVHP> HP_result;
					vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
					mom_conf_VHP mcHP=mc.extend<RVHP>(ind);
					HP_result=((*it)->eval(mcHP,new_ind));
					resultbub=clr_fac*to_HP(HP_result);
				} else {
					complex<RVHP> HP_result=((*it)->eval(mpr->mc_VHP(),ind));
					resultbub=clr_fac*to_HP(HP_result);
				}
			}
	#if _VERBOSE==1
			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
				_MESSAGE4("WARNING : Failed to get ",settings::rational_settings::s_rat_ext_precision," even at VHP as accuracy, we have ",(*it)->get_accuracy());
			}
	#endif
			double thisacc=(*it)->get_accuracy();
			if(thisacc+acc_shift<curacc){
				curacc=thisacc+acc_shift;
			}
			
			resultbubsum+=resultbub;
	#if _VERBOSE==1
			_MESSAGE8("BUB : ",*(*it),"=",real(resultbub),imag(resultbub)>RHP(0)?"+I ":"-I ",abs(imag(resultbub))," estimated acc=",(*it)->get_accuracy()+acc_shift);
	#endif
		}

	#if _VERBOSE==1
		_MESSAGE2("      Sum of all bubbles=",resultbubsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
	#endif

		for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
	#if _VERBOSE==1
			complex<RHP> resulttri=clr_fac*((*it)->eval(mc,ind));
			resulttrisum+=resulttri;

			_MESSAGE6("TRI : ",*(*it),"=",real(resulttri),imag(resulttri)>RHP(0)?"+I ":"-I ",abs(imag(resulttri)));
	#else
			resulttrisum+=clr_fac*((*it)->eval(mc,ind));
	#endif
		}

	#if _VERBOSE==1
		_MESSAGE2("      Sum of all triangles=",resulttrisum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
	#endif

		for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
	#if _VERBOSE==1
			complex<RHP> resultbox=clr_fac*((*it)->eval(mc,ind));
			resultboxsum+=resultbox;

			_MESSAGE6("BOX : ",*(*it),"=",real(resultbox),imag(resultbox)>RHP(0)?"+I ":"-I ",abs(imag(resultbox)));
	#else
			resultboxsum+=clr_fac*((*it)->eval(mc,ind));
	#endif
		}
	#if _VERBOSE==1
		_MESSAGE2("      Sum of all boxes=",resultboxsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
		for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
			complex<RHP> resultpent=clr_fac*((*it)->eval(mc,ind));

			_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>RHP(0)?"+I ":"-I ",abs(imag(resultpent)));
		}
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
	#endif
		
		_accuracy=curacc; // Store the accuracy for this part of the computation

		return resultbubsum+resulttrisum+resultboxsum;
	}
}


template <class PENT,class BOX,class TRI,class BUB> template <class T> complex<T> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_ep(const eval_param<T>& ep)
{
	complex<T> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0), kb_bub(0,0), kb_tri(0,0), kb_box(0,0), ksum, kstep;
	T clr_fac=this->template get_colour_fac<T>();
	
	double curacc(to_double(MaxDigits<T>())); // keep the lowest recorded accuracy

	// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
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
//	Tree_factory TF;
//	process pro_n(ph,p,m,m,m);
//	Rec_Tree* tree=TF.new_tree(pro_n);
	complex<T> treeres_full=get_tree(ep);
//	complex<T> treeres_full=tree->eval(ep);
	if(abs(treeres_full)<BH::DeltaZero<T>()||!BH::settings::general::s_normalise_debug_coeff_output){treeres_full=T(1);_MESSAGE("Not dividing by the tree");}
	else{_MESSAGE2("Tree=",treeres_full);}
	_MESSAGE("");
	_MESSAGE("****************************** Rat Ext (ep) **********************************");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
#endif

	for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
		kstep=clr_fac*((*it)->eval(ep))-kb_bub;
		ksum=resultbubsum+kstep;
		kb_bub=(ksum-resultbubsum)-kstep;
		resultbubsum=ksum;
		
#if _VERBOSE==1
		_MESSAGE8("BUB : ",*(*it),"(n) =",real(kstep/treeres_full),imag(kstep/treeres_full)>T(0)?"+I ":"-I ",abs(imag(kstep/treeres_full))," estimated acc=",(*it)->get_accuracy()+acc_shift);
#endif
		double thisacc=(*it)->get_accuracy();
		if(thisacc+acc_shift<curacc){
			curacc=thisacc+acc_shift;
		}		
	}

#if _VERBOSE==1
	_MESSAGE2("      Sum of all bubbles (n)=",resultbubsum/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

	for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
		kstep=clr_fac*((*it)->eval(ep))-kb_tri;
		ksum=resulttrisum+kstep;
		kb_tri=(ksum-resulttrisum)-kstep;
		resulttrisum=ksum;		
#if _VERBOSE==1
		_MESSAGE6("TRI : ",*(*it),"(n) =",real(kstep/treeres_full),imag(kstep/treeres_full)>T(0)?"+I ":"-I ",abs(imag(kstep/treeres_full)));
#endif
	}

#if _VERBOSE==1
	_MESSAGE2("      Sum of all triangles (n)=",resulttrisum/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

	for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
		kstep=clr_fac*((*it)->eval(ep))-kb_box;
		ksum=resultboxsum+kstep;
		kb_box=(ksum-resultboxsum)-kstep;
		resultboxsum=ksum;		
#if _VERBOSE==1
		_MESSAGE6("BOX : ",*(*it),"(n) =",real(kstep/treeres_full),imag(kstep/treeres_full)>T(0)?"+I ":"-I ",abs(imag(kstep/treeres_full)));
#endif
	}
#if _VERBOSE==1
	_MESSAGE2("      Sum of all boxes (n)=",resultboxsum/treeres_full);
	_MESSAGE("");
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
	for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
		complex<T> resultpent=clr_fac*((*it)->eval(ep));

		_MESSAGE6("PENT : ",*(*it),"(n) =",real(resultpent/treeres_full),imag(resultpent/treeres_full)>T(0)?"+I ":"-I ",abs(imag(resultpent/treeres_full)));
	}
	_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
	_MESSAGE2("      Final result (n)=",(resultbubsum+resulttrisum+resultboxsum)/treeres_full);
#endif
	
	_accuracy=curacc; // Store the accuracy for this part of the computation

	return resultbubsum+resulttrisum+resultboxsum;
}


template <class PENT,class BOX,class TRI,class BUB> complex<R> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_IR_checked_ep(const eval_param<R>& ep)
{

	if ( ! settings::rational_settings::s_use_IR_in_ratext){
		return eval_rat_ep(ep);
	} else {
		complex<R> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0), kb_bub(0,0), kb_tri(0,0), kb_box(0,0), ksum, kstep;
		R clr_fac=this->template get_colour_fac<R>();
		
		double curacc(MaxDigits<R>()); // keep the lowest recorded accuracy

		// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
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
		_MESSAGE("");
		_MESSAGE("************************ Rat Ext with Check (ep) ******************************");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
#endif

		for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_bub;			

			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at HP as accuracy is ",(*it)->get_accuracy()+acc_shift);
#endif
				vector<Cmom<RHP> > *new_mom=extend_momenta<R,RHP>(ep);
				eval_param<RHP> epHP(*new_mom);
				complex<RHP> HP_result=((*it)->eval(epHP));
				kstep=clr_fac*to_double(HP_result);
				delete new_mom; //delete the new momentum vector created in the extend function

				//Check if this was accurate enough
				if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",(*it)->get_accuracy()+acc_shift);
#endif
					// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
					vector<Cmom<RVHP> > *new_mom=extend_momenta<R,RVHP>(ep);
					eval_param<RVHP> epVHP(*new_mom);

					complex<RVHP> VHP_result=((*it)->eval(epVHP));
					kstep=clr_fac*to_double(VHP_result);
					delete new_mom; //delete the new momentum vector created in the extend function
					 //If this was to fail then there is not much more we can do
#if _VERBOSE==1
					if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
						_MESSAGE3("WARNING : Failed to get ",settings::rational_settings::s_rat_ext_precision," digits of precision even at VHP accuracy.");
					}
#endif
				}
			}
			
			ksum=resultbubsum+kstep;
			kb_bub=(ksum-resultbubsum)-kstep;
			resultbubsum=ksum;
			
			double thisacc=(*it)->get_accuracy();
			if(thisacc+acc_shift<curacc){
				curacc=thisacc+acc_shift;
			}
#if _VERBOSE==1
			_MESSAGE8("BUB : ",*(*it),"=",real(kstep),imag(kstep)>R(0)?"+I ":"-I ",abs(imag(kstep))," estimated acc=",thisacc+acc_shift);
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all bubbles=",resultbubsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

		for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_tri;
			ksum=resulttrisum+kstep;
			kb_tri=(ksum-resulttrisum)-kstep;
			resulttrisum=ksum;
#if _VERBOSE==1
			_MESSAGE6("TRI : ",*(*it),"=",real(kstep),imag(kstep)>R(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all triangles=",resulttrisum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

		for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_box;
			ksum=resultboxsum+kstep;
			kb_box=(ksum-resultboxsum)-kstep;
			resultboxsum=ksum;	
#if _VERBOSE==1
			_MESSAGE6("BOX : ",*(*it),"=",real(kstep),imag(kstep)>R(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
		}
#if _VERBOSE==1
		_MESSAGE2("      Sum of all boxes=",resultboxsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
		for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
			complex<R> resultpent=clr_fac*((*it)->eval(ep));

			_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>R(0)?"+I ":"-I ",abs(imag(resultpent)));
		}
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
#endif

		_accuracy=curacc; // Store the accuracy for this part of the computation

		return resultbubsum+resulttrisum+resultboxsum;
	}
}

template <class PENT,class BOX,class TRI,class BUB> complex<RHP> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_IR_checked_ep(const eval_param<RHP>& ep)
{
	if ( ! settings::rational_settings::s_use_IR_in_ratext){
		return eval_rat_ep(ep);
	} else {
		complex<RHP> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0), kb_bub(0,0), kb_tri(0,0), kb_box(0,0), ksum, kstep;
		RHP clr_fac=this->template get_colour_fac<RHP>();
		
		double curacc(to_double(MaxDigits<RHP>())); // keep the lowest recorded accuracy

		// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
		//  difference of the exponent of the tree result from 1
		double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
		RHP treeres=abs(get_tree(ep)); // Get the absolute power of the tree
		if(treeres > BH::DeltaZero<RHP>()){ // If the tree is zero then there is no error_shift
			if(treeres>1){treeres=RHP(1)/treeres;} // we turn all exponents into positive powers

			acc_shift=to_double(log(treeres)/log(RHP(10)));
		}
#endif

	#if _VERBOSE==1
		_MESSAGE("");
		_MESSAGE("*********************** Rat Ext with Check HP (ep) *****************************");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
	#endif

		for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_bub;

			if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
#if _VERBOSE==1
				_MESSAGE2("WARNING : Recomputing bubble at VHP as accuracy is ",(*it)->get_accuracy()+acc_shift);
#endif
				// We simply extend to VHP as we started with either an mc or a multi_precision_reader and not a multi_precision_reader_HP
				vector<Cmom<RVHP> > *new_mom=extend_momenta<RHP,RVHP>(ep);
				eval_param<RVHP> epVHP(*new_mom);
				complex<RVHP> VHP_result=((*it)->eval(epVHP));
				kstep=clr_fac*to_HP(VHP_result);
				delete new_mom; //delete the new momentum vector created in the extend function
				//If this was to fail then there is not much more we can do
#if _VERBOSE==1
				if((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
					_MESSAGE3("WARNING : Failed to get ",settings::rational_settings::s_rat_ext_precision," digits of precision even at VHP accuracy.");
				}
#endif
			}

			ksum=resultbubsum+kstep;
			kb_bub=(ksum-resultbubsum)-kstep;
			resultbubsum=ksum;
			
			double thisacc=(*it)->get_accuracy();
			if(thisacc+acc_shift<curacc){
				curacc=thisacc+acc_shift;
			}			
#if _VERBOSE==1
			_MESSAGE8("BUB : ",*(*it),"=",real(kstep),imag(kstep)>RHP(0)?"+I ":"-I ",abs(imag(kstep))," estimated acc=",thisacc+acc_shift);
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all bubbles=",resultbubsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

		for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_tri;
			ksum=resulttrisum+kstep;
			kb_tri=(ksum-resulttrisum)-kstep;
			resulttrisum=ksum;
#if _VERBOSE==1
			_MESSAGE6("TRI : ",*(*it),"=",real(kstep),imag(kstep)>RHP(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all triangles=",resulttrisum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

		for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
			kstep=clr_fac*((*it)->eval(ep))-kb_box;
			ksum=resultboxsum+kstep;
			kb_box=(ksum-resultboxsum)-kstep;
			resultboxsum=ksum;
#if _VERBOSE==1
			_MESSAGE6("BOX : ",*(*it),"=",real(kstep),imag(kstep)>RHP(0)?"+I ":"-I ",abs(imag(kstep)));
#endif
		}
#if _VERBOSE==1
		_MESSAGE2("      Sum of all boxes=",resultboxsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
		for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
			complex<RHP> resultpent=clr_fac*((*it)->eval(ep));

			_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>RHP(0)?"+I ":"-I ",abs(imag(resultpent)));
		}
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
#endif

		_accuracy=curacc; // Store the accuracy for this part of the computation

		return resultbubsum+resulttrisum+resultboxsum;

	}

}


template <class PENT,class BOX,class TRI,class BUB> ratext_part<PENT,BOX,TRI,BUB>::~ratext_part()
{
	for(typename  vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
		delete *it;
	}
	for(typename  vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
		delete *it;
	}
	for(typename  vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
		delete *it;
	}
	for(typename  vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
		delete *it;
	}
    delete d_tree_ptr;
}


#if BH_USE_GMP
template <class PENT,class BOX,class TRI,class BUB> complex<RGMP> ratext_part<PENT,BOX,TRI,BUB>::eval_rat_IR_checked(momentum_configuration<RGMP>& mc, const vector<int>& ind)
{

	if ( ! settings::rational_settings::s_use_IR_in_ratext){
		return eval_rat(mc,ind);
	} else {
		complex<RGMP> resultbubsum(0,0), resulttrisum(0,0), resultboxsum(0,0);
		RGMP clr_fac=this->template get_colour_fac<RGMP>();

		double curacc(to_double(MaxDigits<RGMP>())); // keep the lowest recorded accuracy

		// As the result will be proportional to the tree (times a relatively small correction) we estimate an error correction offset based on the
		//  difference of the exponent of the tree result from 1
		double acc_shift(0);
#if _CONSERVATIVE_ERROR_ESTIMATE==1
		RGMP treeres=abs(get_tree(mc,ind)); // Get the absolute power of the tree
		if(treeres > BH::DeltaZero<RGMP>()){ // If the tree is zero then there is no error_shift
			if(treeres>1){treeres=RGMP(1)/treeres;} // we turn all exponents into positive powers

			acc_shift=to_double(log(treeres)/log(RGMP(10)));
		}
#endif


#if _VERBOSE==1
		_MESSAGE("");
		_MESSAGE("************************ Rat Ext with Check (mc) ******************************");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ BUBBLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_bubbles.size(), " Bubbles");
#endif

		for(typename vector<BUB*>::iterator it=d_bubbles.begin();it!=d_bubbles.end();it++){
			complex<RGMP> resultbub=clr_fac*((*it)->eval(mc,ind));
			while((*it)->get_accuracy()+acc_shift<settings::rational_settings::s_rat_ext_precision){
				BH_DEBUG_MESSAGE2("WARNING : Recomputing bubble at higher precision as accuracy is ",(*it)->get_accuracy()+acc_shift);
					int oldPrecision=RGMP::get_current_precision();
					RGMP::set_precision(2*oldPrecision);
					complex<RGMP> HigherPrecisionResult;
					vector<int> new_ind; for (int jj=1;jj<=ind.size();jj++){new_ind.push_back(jj);}
					momentum_configuration<RGMP> mcHP=mc.extend<RGMP>(ind);
					HigherPrecisionResult=((*it)->eval(mcHP,new_ind));
					resultbub=clr_fac*HigherPrecisionResult;
			}
			double thisacc=(*it)->get_accuracy();
			if(thisacc+acc_shift<curacc){
				curacc=thisacc+acc_shift;
			}
			resultbubsum+=resultbub;
#if _VERBOSE==1
			_MESSAGE8("BUB : ",*(*it),"=",real(resultbub),imag(resultbub)>R(0)?"+I ":"-I ",abs(imag(resultbub))," estimated acc=",thisacc+acc_shift);
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all bubbles=",resultbubsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TRIANGLES ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_triangles.size(), " Triangles");
#endif

		for(typename vector<TRI*>::iterator it=d_triangles.begin();it!=d_triangles.end();it++){
#if _VERBOSE==1
			complex<RGMP> resulttri=clr_fac*((*it)->eval(mc,ind));
			resulttrisum+=resulttri;

			_MESSAGE6("TRI : ",*(*it),"=",real(resulttri),imag(resulttri)>R(0)?"+I ":"-I ",abs(imag(resulttri)));
#else
			resulttrisum+=clr_fac*((*it)->eval(mc,ind));
#endif
		}

#if _VERBOSE==1
		_MESSAGE2("      Sum of all triangles=",resulttrisum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Boxes ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("Eval ",d_boxes.size(), " Boxes");
#endif

		for(typename vector<BOX*>::iterator it=d_boxes.begin();it!=d_boxes.end();it++){
#if _VERBOSE==1
			complex<RGMP> resultbox=clr_fac*((*it)->eval(mc,ind));
			resultboxsum+=resultbox;

			_MESSAGE6("BOX : ",*(*it),"=",real(resultbox),imag(resultbox)>R(0)?"+I ":"-I ",abs(imag(resultbox)));
#else
			resultboxsum+=clr_fac*((*it)->eval(mc,ind));
#endif
		}
#if _VERBOSE==1
		_MESSAGE2("      Sum of all boxes=",resultboxsum);
		_MESSAGE("");
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE3("We also have ",d_pentagons.size()," pentagons.");
		for(typename vector<PENT*>::iterator it=d_pentagons.begin();it!=d_pentagons.end();it++){
			complex<RGMP> resultpent=clr_fac*((*it)->eval(mc,ind));

			_MESSAGE6("PENT : ",*(*it),"=",real(resultpent),imag(resultpent)>R(0)?"+I ":"-I ",abs(imag(resultpent)));
		}
		_MESSAGE("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
		_MESSAGE2("      Final result=",resultbubsum+resulttrisum+resultboxsum);
#endif

		_accuracy=curacc; // Store the accuracy for this part of the computation

		return resultbubsum+resulttrisum+resultboxsum;
	}
}


#endif




}

}
