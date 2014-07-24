/*
 * cut_Darren.hpp
 *
 *  Created on: 1 Jul 2009
 *      Author: daniel
 */

#ifndef CUT_DARREN_HPP_
#define CUT_DARREN_HPP_

#include "BH_utilities.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include <iostream>
#include <cassert>
#include "cut_worker.h"
#include "mode_dependent_typedefs.h"  // for TREE_TYPE and TREE_FACTORY_TYPE

#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif

#define _VERBOSE 0

#define USE_NEW_CUT 1 // Set to 1 if you want to use the version of the bubble/triangle/box where we do not use information fed from the level below, otherwise 0

#define use_specs 0

#define use_corner_strat 1

#include "momenta_evaluator.hpp"

namespace BH {

//This function returns 1 if the three-point amplitude is mostly plus or -1 if the amplitude is mostely minus
int Get3ptType_new(corner_type corner);

namespace cut {

namespace Darren {

/*
 *
 *
 * A function to return the maximum number of digits a particular type has to work with
 *
 *
 */

template <class T> inline T MaxDigitsE();

template<> inline R MaxDigitsE(){return R(1e16);}
template<> inline RHP MaxDigitsE(){return RHP(1e32);}
template<> inline RVHP MaxDigitsE(){return RVHP(1e64);}
#if BH_USE_GMP
template<> inline RGMP MaxDigitsE(){return RGMP(exp10(RGMP::get_current_nbr_digits()));}
#endif

template <class T> inline T MaxDigits();

template<> inline R MaxDigits(){return R(16);}
template<> inline RHP MaxDigits(){return RHP(32);}
template<> inline RVHP MaxDigits(){return RVHP(64);}
#if BH_USE_GMP
template<> inline RGMP MaxDigits(){return RGMP(RGMP::get_current_nbr_digits());}
#endif





template <class T> struct do_delete : public std::unary_function<T*,void> {
        void operator()(T* ptr) { delete ptr;}
};




template <class BASE> Darren_wrapper<BASE>::~Darren_wrapper(){
	for(int i=0;i<d_trees.size();i++){
		delete d_trees[i];
	}
}

template <class BASE> Darren_wrapper<BASE>::Darren_wrapper(const BASE& cd): BASE(cd){
	TREE_FACTORY_TYPE TF;
       for (int i=1;i<=BASE::nc();i++){
    	   d_trees.push_back(TF.new_tree(BASE::get_process(i)));
       }

}

template <class BASE> std::complex<R> Darren_wrapper<BASE>::eval_tree(int n,momentum_configuration<R>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
template <class BASE> std::complex<RHP> Darren_wrapper<BASE>::eval_tree(int n,momentum_configuration<RHP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
template <class BASE> std::complex<RVHP> Darren_wrapper<BASE>::eval_tree(int n,momentum_configuration<RVHP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}

template <class BASE> std::complex<R> Darren_wrapper<BASE>::eval_tree(int n,const eval_param<R>& ep){
        return d_trees[n-1]->eval(ep);
}
template <class BASE> std::complex<RHP> Darren_wrapper<BASE>::eval_tree(int n,const eval_param<RHP>& ep){
        return d_trees[n-1]->eval(ep);
}
template <class BASE> std::complex<RVHP> Darren_wrapper<BASE>::eval_tree(int n,const eval_param<RVHP>& ep){
        return d_trees[n-1]->eval(ep);
}
#if BH_USE_GMP
template <class BASE> std::complex<RGMP> Darren_wrapper<BASE>::eval_tree(int n,momentum_configuration<RGMP>& mc, const std::vector<int>& ind){
        return d_trees[n-1]->eval(mc,ind);
}
template <class BASE> std::complex<RGMP> Darren_wrapper<BASE>::eval_tree(int n,const eval_param<RGMP>& ep){
        return d_trees[n-1]->eval(ep);
}
#endif



template <class cutDbase,class TriangleSpecs> triangle_Darren<cutDbase,TriangleSpecs>::~triangle_Darren()
{
        delete coeffkeep;
        delete polekeep;
        delete denpole;

        delete coeffkeep_HP;
        delete polekeep_HP;
        delete denpole_HP;

        delete coeffkeep_VHP;
        delete polekeep_VHP;
        delete denpole_VHP;

    	delete _ep[0];
    	delete _ep_HP[0];
    	delete _ep_VHP[0];
    	delete _ep[1];
    	delete _ep_HP[1];
    	delete _ep_VHP[1];
    	delete _ep[2];
    	delete _ep_HP[2];
    	delete _ep_VHP[2];

#if BH_USE_GMP
        delete coeffkeep_GMP;
        delete polekeep_GMP;
        delete denpole_GMP;
#endif

}



template <class TriangleType> struct box_subtraction_helper {
	template <class T,int CPOINTS> void prepare_tri_legs(TriangleType& self,coeffparam<T,CPOINTS>& tp2,const vector<int>& ind){
		tp2.tri_legs.push_back(self._k1leg-1);
		tp2.tri_legs.push_back(self._k3leg-1);
		tp2.tri_legs.push_back(self._k2leg-1);
	}
};



template <class cutDbase, class TriangleSpecs>  template <class T> void triangle_Darren<cutDbase,TriangleSpecs>::get_coeffs_fn(momentum_configuration<T>& mc,const  vector<int>& ind, coeffparam<T,TriangleSpecs::CPOINTS>& tp)
{

	//
	// Compute the momenta of the three corners of the box
	//

	BH_DEBUG_MESSAGE2("triangle_Darren::get_coeffs_fn for ",*this);
	coeffparam<T,TriangleSpecs::CPOINTS> tp2;
	tp2.mcID=tp.mcID;

	vector<complex<T> > trees_result_1(TriangleSpecs::TPOINTSTRI);
    vector<complex<T> > trees_result_2(TriangleSpecs::TPOINTSTRI);
    vector<complex<T> > trees_result_3(TriangleSpecs::TPOINTSTRI);


    typename TriangleSpecs::CornerTreeStrategy* ptr= static_cast<typename TriangleSpecs::CornerTreeStrategy*>(this);
    ptr->fill_trees(mc,ind,trees_result_1,trees_result_2,trees_result_3,tp2,*this,_k1leg,_k2leg,_k3leg,_masslessleg_type,_reverse);

    BH_DEBUG(
			cout << "trees_result_1: \n" ;
			copy(trees_result_1.begin(),trees_result_1.end(),ostream_iterator<complex<T> >(cout," "));
			cout << "\n";
			cout << "trees_result_2: \n" ;
			copy(trees_result_2.begin(),trees_result_2.end(),ostream_iterator<complex<T> >(cout," "));
			cout << "\n";
			cout << "trees_result_3: \n" ;
			copy(trees_result_3.begin(),trees_result_3.end(),ostream_iterator<complex<T> >(cout," "));
			cout << "\n";
    )


    //
    // Now compute the coefficients
    //

	//Extract the box contributions, the last element in the returned array is the contribution
	// at t=0 if this is a three-mass triangle
	complex<T> boxsub[TriangleSpecs::TPOINTSTRI+1]={complex<T>(0,0)};

    box_subtraction_helper<triangle_Darren<cutDbase,TriangleSpecs> > Nbsh;
    Nbsh.prepare_tri_legs(*this,tp2,ind);

	//Get all the box subtractions adding each new one to the previous ones contained in coeffsret
	// First empty the storage of these boxes in case this coefficient has been computed previously
	tp2.masslessleg=_masslessleg_type;
	tp2.reverse=_reverse;
	empty_box_store(complex<T>(0,0));
    for(size_t box=1; box <= cutDbase::daughters_nbr(); box++){
    	tp2.box_corner=cutDbase::get_opened_corner(box);
#ifndef NDEBUG
    	typename triangle_daughter_type<cutDbase,TriangleSpecs::CPOINTS,TriangleSpecs::TPOINTSTRI>::type* wt=dynamic_cast<typename triangle_daughter_type<cutDbase,TriangleSpecs::CPOINTS,TriangleSpecs::TPOINTSTRI>::type* >(cutDbase::get_daughter(box));
    	assert(wt);
    	wt->get_sub_terms(mc,ind,boxsub,tp2);
#else
    	static_cast< typename triangle_daughter_type<cutDbase,TriangleSpecs::CPOINTS,TriangleSpecs::TPOINTSTRI>::type*>(cutDbase::get_daughter(box))->get_sub_terms(ep,boxsub,tp2);
#endif
    	// Store these results for later use in the bubble coeff subtractions if we have any poles
    	boxcoeff_add(tp2.boxcoeffs[0],tp2.boxcoeffs[1],tp2.Nboxpoles[0],tp2.Nboxpoles[1],tp2.denfac1);
    }

	// Now we calculate the coefficient from examining the triple cut
	complex<T> amp[TriangleSpecs::TPOINTSTRI];
	for(int icirc=0; icirc<TriangleSpecs::TPOINTSTRI; icirc++){

		amp[icirc]=Amp_safe( trees_result_1[icirc]
		                    *trees_result_2[icirc]
		                    *trees_result_3[icirc]
		                     )
                                -boxsub[icirc];

	}




	//We use this CPOINTS*TPOINTSTRI double array contains all the coefficients we use to multiply the 8 terms we generate above
	// to produce the CPOINTS coefficients.
	const complex<T>* coeff_pts;
	triangle_Darren_eval_points<TriangleSpecs::CPOINTS,TriangleSpecs::TPOINTSTRI>::get_coeff_pts(coeff_pts);

	complex<T> triC0res[TriangleSpecs::CPOINTS]={complex<T>(0,0)};
	//Now produce the coeffs by starting from 0 and adding the above terms with the correct signs
	for(int ifin=0;ifin<TriangleSpecs::CPOINTS;ifin++){
		for(int jfin=0;jfin<TriangleSpecs::TPOINTSTRI;jfin++){
			triC0res[ifin]+=coeff_pts[ifin*TriangleSpecs::TPOINTSTRI+jfin]*amp[jfin];
		}
	}

	// Debugging output of all the computed coeffs
#if _VERBOSE
	// Debugging output of all the computed coeffs
	cout << *this << endl << " TRI: {";
	for(int i=0;i<CPOINTS;i++){
		cout << triC0res[i] << ",";
	}
	cout << "}" << endl;
#endif

	// Finally set the actual C0 coefficient, this will be different
	//  from coeffs[(CPOINTS-1)/2] when we have a three mass triangle where we have to add the
	//  "t=0" contribution
	set_C0coeff((triC0res[(TriangleSpecs::CPOINTS-1)/2]+boxsub[TriangleSpecs::TPOINTSTRI]/T(TriangleSpecs::TPOINTSTRI))/T(2));
	BH_DEBUG_PRINT(C0coeff);

	//Save the vec1 and vec2 and opened corner so we can relate the c's we compute here to those
	// computed in the other cases
	set_tri_param_basis_vectors(tp2.vec1c,tp2.vec2c);
	BH_DEBUG_PRINT(tp2.vec1c);
	BH_DEBUG_PRINT(tp2.vec2c);
	set_tri_param_gamma(tp2.gamma);
	BH_DEBUG_PRINT(tp2.gamma);
	coeffkeep_add(triC0res);

//	_PRINT((triC0res[(TriangleSpecs::CPOINTS-1)/2]+boxsub[TriangleSpecs::TPOINTSTRI]/T(TriangleSpecs::TPOINTSTRI))/T(2));

	//Mark that we have computed these using the original mcID not the sub_mom_conf one from the mc that was passed to this function
	mcID=tp.mcID;
	indID=ind;
}


template <class cutDbase, class TriangleSpecs> template <class T> complex<T> triangle_Darren_plusminus<cutDbase,TriangleSpecs>::get_sub_terms_work(momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,TriangleSpecs::CPOINTS>& tp)
{
	tri_sub_info<T> tsi;

	//Set up the K's, Kflats and gamma for this triangle depending upon which leg was opened
	size_t mass_2, mass_3, leg2, leg3;
	// Here K2t can be massless
	momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	BH_DEBUG_PRINT(tp.tri_corner);
	switch(tp.tri_corner)
	{
	case 1:
		for(size_t kiter=1;kiter<=cutDbase::corner_size(1);kiter++){
			K2sum_mom+=mc.mom(ind[cutDbase::corner_ind(1,kiter)-1]);
		}
		leg2=1;
		leg3=2;
		break;

	case 2:
		for(size_t kiter=1;kiter<=cutDbase::corner_size(3);kiter++){
			K2sum_mom-=mc.mom(ind[cutDbase::corner_ind(3,kiter)-1]);
		}
		leg2=3;
		leg3=2;
		break;

	case 3:
		for(size_t kiter=1;kiter<=cutDbase::corner_size(3);kiter++){
			K2sum_mom+=mc.mom(ind[cutDbase::corner_ind(3,kiter)-1]);
		}
		leg2=3;
		leg3=1;
		break;
	}
	tp.K2=K2sum_mom;
	//Save the masses of these legs
	mass_2=cutDbase::corner_size(leg2);
	mass_3=cutDbase::corner_size(leg3);
  	// For a potentially massless leg we check its corner type and check if it is non zero.
 	//  if it is zero then leg 2 is actually massive. We make mass_2 larger than one so that it will avoid
 	//  the massless leg 2 code. As at least one leg must be massless then leg 3 is massless.
    if(mass_2==1){
    	tp.masslessleg=-Get3ptType_new(cutDbase::c_type(leg2));// In this case we have the opposite sign for the type of solution, as stored by masslessleg, than for a massless K3 leg
    	if(tp.masslessleg==0){
    		mass_2=2;
    	}
    }

	//Now construct the constants we will use and pick gamma with the
	// convention that if S1 or S2 would be 0 then it would not vanish
	// this avoids the need to specialize this statement.
	complex<T> K1K2=tp.K1*tp.K2;

	//Now calculate the poles of the triangles sitting above the bubbles
	//  we have a different solution depending upon which leg was massless
	complex<T> Nysolp[TriangleSpecs::TPOINTSBUB];
	complex<T> Nysolm[TriangleSpecs::TPOINTSBUB];

	// Get the points to evaluate the circle on
	bubble_Darren_eval_points<TriangleSpecs::TPOINTSBUB,TriangleSpecs::YPOINTS>::get_eval_pts(tsi.eval_pts);

	// Also set up a variable to keep track of which leg was massless and whether it is a "+" or a "-" leg,
	// 0 means all legs are massive and +1 is a + leg and -1 for a -leg.
	//Find which leg is massless (it can only be 2 or 3), if it is both then we pick the K2 leg solution
	// to determine the sign of masslessleg, though we could pick either.
	if(mass_2==1){// K2 is massless,
        tp.S2=complex<T>(0,0);
        Cmom<T> K2c(tp.K2);
		tp.masslessleg=-Get3ptType_new(cutDbase::c_type(leg2));// In this case we have the opposite sign for the type of solution, as stored by masslessleg, than for a massless K3 leg
		for(int iNysol = 0; iNysol<TriangleSpecs::TPOINTSBUB; iNysol++){
			Nysolp[iNysol]=T(1)+tsi.eval_pts[iNysol]*(tp.K1flatbc.L()*K2c.L())/(tp.chic.L()*K2c.L());
			Nysolm[iNysol]=-tsi.eval_pts[iNysol]*(tp.chic.Lt()*K2c.Lt())/(tp.K1flatbc.Lt()*K2c.Lt());
		}
	}
	else if(mass_3==1){ // K3 is massless
		tp.S2=tp.K2.square();
		tp.masslessleg=Get3ptType_new(cutDbase::c_type(leg3));// Store the type of solution for l we have i.e. wheather the term multiplying t or 1/t vanished
		Cmom<T> K3(tp.K1-tp.K2);
		for(int iNysol = 0; iNysol<TriangleSpecs::TPOINTSBUB; iNysol++){
			Nysolp[iNysol]=T(1)-tsi.eval_pts[iNysol]*(K3.Lt()*tp.chic.Lt())/(K3.Lt()*tp.K1flatbc.Lt());
			Nysolm[iNysol]=-tsi.eval_pts[iNysol]*(K3.L()*tp.K1flatbc.L())/(tp.chic.L()*K3.L());
		}
	}

    // If the triangle has not been computed compute it now
    if(!this->is_eval(tp.mcID,ind)){
    	this->get_coeffs(mc,ind,tp);
    }

	// Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	this->get_boxes(tsi.triboxcoeff,tsi.triboxpoles,tsi.denfac);

	// Get the coefficients and other factors from the previously computed triangle, if not
	//  this will just compute them
	this->coeffkeep_get(tsi.orig_coeffs);
	this->get_tri_param_gamma(tsi.gamma_old);
	this->get_tri_param_basis_vectors(tsi.v1old,tsi.v2old);

    tsi.Nysolm=Nysolm;
    tsi.Nysolp=Nysolp;

    tsi.reverse=this->_reverse;
	BH_DEBUG(
			for(int iNysol = 0; iNysol<TriangleSpecs::TPOINTSBUB; iNysol++){
				cout << "Nysolp[" << iNysol << "]:" << Nysolp[iNysol] << "\n";
				cout << "Nysolm[" << iNysol << "]:" << Nysolm[iNysol] << "\n";
			}
	)
    return this->get_sub_terms_work_pm(*this,mc,ind,tp,tsi);
}


template <class cutDbase, class TriangleSpecs> template <class T> complex<T> triangle_Darren_3mass<cutDbase,TriangleSpecs>::get_sub_terms_work(momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,TriangleSpecs::CPOINTS>& tp)
{
         tp.masslessleg=0;// Mark that this is a three mass triangle

         //Set up the K's, Kflats and gamma for this triangle depending upon which leg was opened
         // Here K2t is never massless
         momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
         switch(tp.tri_corner)
         {
         case 1:
                 for(size_t kiter=1;kiter<=cutDbase::corner_size(1);kiter++){
                         K2sum_mom+=mc.mom(ind[cutDbase::corner_ind(1,kiter)-1]);
                 }
                 break;

         case 2:
                 for(size_t kiter=1;kiter<=cutDbase::corner_size(3);kiter++){
                         K2sum_mom-=mc.mom(ind[cutDbase::corner_ind(3,kiter)-1]);
                 }
                 break;

          case 3:
                 for(size_t kiter=1;kiter<=cutDbase::corner_size(3);kiter++){
                         K2sum_mom+=mc.mom(ind[cutDbase::corner_ind(3,kiter)-1]);
                 }
                 break;
         }
         tp.K2=K2sum_mom;
         tp.S2=tp.K2.square();

 	    // We only need to compute K1flat and K2flat etc if we need to compute the triangle
 	    //  if it has been already computed then we reuse the previous computation
 	    if(!this->is_eval(tp.mcID,ind)){
 	    	this->get_coeffs(mc,ind,tp);
 	    }

         // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
        tri_sub_info<T> tsi;
 	    this->get_boxes(tsi.triboxcoeff,tsi.triboxpoles,tsi.denfac);

         // Get the coefficients and other factors from the previously computed triangle, if not
         //  this will just compute them
         this->coeffkeep_get(tsi.orig_coeffs);
         this->get_tri_param_gamma(tsi.gamma_old);
         this->get_tri_param_basis_vectors(tsi.v1old,tsi.v2old);

         tsi.reverse=this->_reverse;

         return this->get_sub_terms_work_3m(*this,mc,ind,tp,tsi);



}



template <class BASE, int CPOINTS, int TPOINTSTRI> void box_Darren<BASE,CPOINTS,TPOINTSTRI>::init() {
	mcID=-1;
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	indlst[0].assign(BASE::corner_size(1)+2,0);
	indlst[1].assign(BASE::corner_size(2)+2,0);
	indlst[2].assign(BASE::corner_size(3)+2,0);
	indlst[3].assign(BASE::corner_size(4)+2,0);

	//Set up the eval_params for the corners
	_ep[0]=new eval_param<R>(BASE::corner_size(1)+2);
	_ep_HP[0]=new eval_param<RHP>(BASE::corner_size(1)+2);
	_ep_VHP[0]=new eval_param<RVHP>(BASE::corner_size(1)+2);

	_ep[1]=new eval_param<R>(BASE::corner_size(2)+2);
	_ep_HP[1]=new eval_param<RHP>(BASE::corner_size(2)+2);
	_ep_VHP[1]=new eval_param<RVHP>(BASE::corner_size(2)+2);

	_ep[2]=new eval_param<R>(BASE::corner_size(3)+2);
	_ep_HP[2]=new eval_param<RHP>(BASE::corner_size(3)+2);
	_ep_VHP[2]=new eval_param<RVHP>(BASE::corner_size(3)+2);

	_ep[3]=new eval_param<R>(BASE::corner_size(4)+2);
	_ep_HP[3]=new eval_param<RHP>(BASE::corner_size(4)+2);
	_ep_VHP[3]=new eval_param<RVHP>(BASE::corner_size(4)+2);

#if BH_USE_GMP
	_ep_GMP[0]=new eval_param<RGMP>(BASE::corner_size(1)+2);
	_ep_GMP[1]=new eval_param<RGMP>(BASE::corner_size(2)+2);
	_ep_GMP[2]=new eval_param<RGMP>(BASE::corner_size(3)+2);
	_ep_GMP[3]=new eval_param<RGMP>(BASE::corner_size(4)+2);
#endif

}
template <class BASE, int CPOINTS, int TPOINTSTRI> box_Darren<BASE,CPOINTS,TPOINTSTRI>::~box_Darren() {
	for (int i=0;i<4;i++){
		delete _ep[i];
		delete _ep_HP[i];
		delete _ep_VHP[i];
#if BH_USE_GMP
		delete _ep_GMP[i];
#endif
	}
}

template <class BASE, int CPOINTS, int TPOINTSTRI> void box_Darren<BASE,CPOINTS,TPOINTSTRI>::init(const boxD& cd) {


	int i=1;
	while(i<5){
		if((BASE::corner_size(i)>1||(BASE::corner_size(i)==1&&cd.get_process(i).p(2).mass_label()>0))
				&&(BASE::corner_size(((i+2)%4)+1)==1&&cd.get_process(((i+2)%4)+1).p(2).mass_label()==0)){
			break;
		}
		i++;
	};
	_k1leg=(i-1)%4+1;//If no massless K2 with massive K1 was found then i=4 and so we choose simply 1,4,3,2
	_k2leg=((_k1leg+2)%4)+1;
	_k3leg=((_k1leg+1)%4)+1;
	_k4leg=(_k1leg%4)+1;

	if(BASE::corner_size(_k1leg)==1&&cd.get_process(_k1leg).p(2).mass_label()==0){
		_massless_K1=true;
	}
	else{
		_massless_K1=false;
	}

	_masslessleg_type=0; // Start by assuming that k2 is massive, if not then get the corner type
	if((BASE::corner_size(_k2leg)==1)&&(cd.get_process(_k2leg).p(2).mass_label()==0)){
		_masslessleg_type=-Get3ptType_new(BASE::c_type(_k2leg));
	}

	init();
}





template <class cutDbase,int CPOINTS, int TPOINTSTRI> C box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::eval(mom_conf& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a box before computing the bubble");
        }
        return (this->template get_symmetry_factor<double>())
        *
        D0coeff;
}

template <class cutDbase,int CPOINTS, int TPOINTSTRI> CHP box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::eval(mom_conf_HP& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a box before computing the bubble");
        }
        return this->template get_symmetry_factor<RHP>()*D0coeff_HP;
}

template <class cutDbase,int CPOINTS, int TPOINTSTRI> CVHP box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::eval(mom_conf_VHP& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a box before computing the bubble");
        }
        return this->template get_symmetry_factor<RVHP>()*D0coeff_VHP;
}

#if BH_USE_GMP
template <class cutDbase,int CPOINTS, int TPOINTSTRI> CGMP box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a box before computing the bubble");
        }
    	BH_DEBUG_MESSAGE4("symmetry factor: " ,this->template get_symmetry_factor<RGMP>()," prec:",this->template get_symmetry_factor<RGMP>().get_precision())
    	BH_DEBUG_MESSAGE4("box coeff: " ,D0coeff_GMP," prec:",D0coeff_GMP.real().get_precision())

        return this->template get_symmetry_factor<RGMP>()*D0coeff_GMP;
}
#endif

template <class BASE, class TriangleSpecs> void  triangle_Darren<BASE,TriangleSpecs>::init(){

	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	coeffkeep=new vector<C>;
	polekeep=new vector<C>;
	denpole=new vector<C>;

	coeffkeep_HP=new vector<CHP>;
	polekeep_HP=new vector<CHP>;
	denpole_HP=new vector<CHP>;

	coeffkeep_VHP=new vector<CVHP>;
	polekeep_VHP=new vector<CVHP>;
	denpole_VHP=new vector<CVHP>;

	indlst[0].assign(BASE::corner_size(1)+2,0);
	indlst[1].assign(BASE::corner_size(2)+2,0);
	indlst[2].assign(BASE::corner_size(3)+2,0);
	//Set up the eval_params for the corners
	_ep[0]=new eval_param<R>(BASE::corner_size(1)+2);
	_ep_HP[0]=new eval_param<RHP>(BASE::corner_size(1)+2);
	_ep_VHP[0]=new eval_param<RVHP>(BASE::corner_size(1)+2);

	_ep[1]=new eval_param<R>(BASE::corner_size(2)+2);
	_ep_HP[1]=new eval_param<RHP>(BASE::corner_size(2)+2);
	_ep_VHP[1]=new eval_param<RVHP>(BASE::corner_size(2)+2);

	_ep[2]=new eval_param<R>(BASE::corner_size(3)+2);
	_ep_HP[2]=new eval_param<RHP>(BASE::corner_size(3)+2);
	_ep_VHP[2]=new eval_param<RVHP>(BASE::corner_size(3)+2);

#if BH_USE_GMP
	coeffkeep_GMP=new vector<CGMP>;
	polekeep_GMP=new vector<CGMP>;
	denpole_GMP=new vector<CGMP>;
	_ep_GMP[0]=new eval_param<RGMP>(BASE::corner_size(1)+2);
	_ep_GMP[1]=new eval_param<RGMP>(BASE::corner_size(2)+2);
	_ep_GMP[2]=new eval_param<RGMP>(BASE::corner_size(3)+2);
#endif

	// We can have the situation when there are triangles but no boxes but we
	//  still need the static initalised coefficient points from the box so
	//  we initalise it here by creating an object of the correct type
	box_Darren_eval_points<TriangleSpecs::TPOINTSTRI> init_temp;
}





template <class BASE,class BubbleSpecs> void  bubble_Darren<BASE,BubbleSpecs>::init(){
	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	indlst[0].assign(BASE::corner_size(1)+2,0);
	indlst[1].assign(BASE::corner_size(2)+2,0);

	//Set up the eval_params for the corners
	_ep[0]=new eval_param<R>(BASE::corner_size(1)+2);
	_ep_HP[0]=new eval_param<RHP>(BASE::corner_size(1)+2);
	_ep_VHP[0]=new eval_param<RVHP>(BASE::corner_size(1)+2);

	_ep[1]=new eval_param<R>(BASE::corner_size(2)+2);
	_ep_HP[1]=new eval_param<RHP>(BASE::corner_size(2)+2);
	_ep_VHP[1]=new eval_param<RVHP>(BASE::corner_size(2)+2);

#if BH_USE_GMP
	_ep_GMP[0]=new eval_param<RGMP>(BASE::corner_size(1)+2);
	_ep_GMP[1]=new eval_param<RGMP>(BASE::corner_size(2)+2);
#endif

	// We can have the situation when there are bubbles but no triangles but we
	//  still need the static initalised coefficient points from the triangle so
	//  we initalise it here by creating an object of the correct type
	triangle_Darren_eval_points<BubbleSpecs::CPOINTS,BubbleSpecs::TPOINTSTRI> init_temp;
}

template <class BASE,class BubbleSpecs>  bubble_Darren<BASE,BubbleSpecs>::~bubble_Darren(){
	for (int i=0;i<2;i++){
		delete _ep[i];
		delete _ep_HP[i];
		delete _ep_VHP[i];
#if BH_USE_GMP
		delete _ep_GMP[i];
#endif
	}
}

template <class cutDbase, int CPOINTS, int TPOINTSTRI> template <class T> void box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::get_sub_terms(momentum_configuration<T>& mc, const vector<int>& ind, complex<T>* coeffsret, coeffparam<T,CPOINTS>& tp)
{
        //Depending upon the type of triangle this originally came from we will need to choose a specific
        // momentum parameterisation and hence the solution for the poles of the box.
        // To do this we need to get the correct solution depending upon which l_i (l_0,l_1 and l_2) leg appears circularly
        //  to the left of the l_3 leg. We want to solve (l_i-K_3)^2=0 and as the only difference between the l_i's is in the
        //  alpha_i's we simply have to choose these. To do this we produce a hash using (type-1)*3+corner to uniquely identify each corner

        // We need to insert the opened leg in the correct place
        int orig_ltype=(tp.box_corner-1)%3;
        int ltype;
        //If we are going in the opposite direction we must also swap K4->-K4 as well as l->-l
        complex<T> ratiofac(tp.reverse);//If the cut momenta are circulating in the opposite direction then we need to change the sign of l
        if(tp.reverse>0){// if we have 321 ordering for k1, k2, k3 then we have one ordering for l, l1, l2
        	ltype=(3+orig_ltype-tp.tri_legs[0])%3; // Store the location of the leg before l_4 i.e we have l_4=l_[{tp.box_corner-1}]-K_4
        }
        else{// otherwise we have a 123 ordering for k1, k2, k3 then we have one ordering for l, l1, l2
            ltype=2;
            if(tp.tri_legs[0]==orig_ltype){//leg K1 is opened
            	ltype=1;
            }else if(tp.tri_legs[1]==orig_ltype){//leg K2 is opened
            	ltype=0;
            }//leg K3 is opened which is our default choice
    //        _PRINT((3+orig_ltype+tp.tri_legs[0])%3);
    //        _PRINT(ltype);
        }

    	momentum<complex<T> > K4sum_mom(mc.mom(ind[cutDbase::corner_ind(tp.box_corner,1)-1]));
    	for(size_t kiter=2;kiter<=cutDbase::corner_size(tp.box_corner);kiter++){
    		K4sum_mom+=mc.mom(ind[cutDbase::corner_ind(tp.box_corner,kiter)-1]);
    	}
    	momentum<complex<T> > K4=ratiofac*K4sum_mom;
        complex<T> S4=K4.square();


    	//Set up tp.nboxpart1 so that they it matches the L_4=l_i-K_4 (we do not always have L_4=l-K_4)
    	complex<T> nboxpart1;
        switch(ltype){
    	case 1:// l_4=l_1-K4
    		nboxpart1=S4-T(2)*((tp.alp2*tp.K1flatc.P()+tp.alp1*tp.K2flatc.P()-tp.K1)*K4);
    		break;
    	case 2:// l_4=l_2-K4
    		nboxpart1=S4-T(2)*((tp.alp2*tp.K1flatc.P()+tp.alp1*tp.K2flatc.P()-tp.K2)*K4);
    		break;
    	default:// l_4=l-K4
    		nboxpart1=S4-T(2)*((tp.alp2*tp.K1flatc.P()+tp.alp1*tp.K2flatc.P())*K4);
		}

    	//We now calculate the pole solutions
        tp.denfac1=T(2)*(tp.vec1c.P()*K4);
        complex<T> denfac2=T(2)*(tp.vec2c.P()*K4);

        int polepos, nsols;
        if(tp.masslessleg==0){
    		complex<T> nboxpart2=sqrt(pow(nboxpart1,2)-T(4)*tp.alp1*tp.alp2*tp.denfac1*denfac2);

    		// We always have two solutions in the three mass case
    		tp.Nboxpoles[0]=(nboxpart1+nboxpart2)/(T(2)*tp.denfac1);
    		tp.Nboxpoles[1]=(nboxpart1-nboxpart2)/(T(2)*tp.denfac1);
    		nsols=2;
    		polepos=1;
        }
        else{
        	if(tp.nboxpart1.real()<T(0)){
        		tp.Nboxpoles[0]=complex<T>(0,0);
        		tp.Nboxpoles[1]=nboxpart1/tp.denfac1;
        		polepos=1;
        	}
        	else
        	{
        		// Choose the non vanishing solution as the "first" when we have only a single pole i.e. one and two mass cases
        		tp.Nboxpoles[0]=nboxpart1/tp.denfac1;
        		tp.Nboxpoles[1]=complex<T>(0,0);
        		polepos=0;
        	}
        	nsols=1;
        }

        //
    	// Compute the coefficients if we have not already done so
    	//
        if(!((tp.mcID==mcID)&&(ind==indID))){
    		get_coeffs(mc,ind,tp);
    	}

    	//
    	// Now match the solutions to the correct pole
    	//
    	complex<T> ratio, Tn;
        get_denfac(Tn);
        Cmom<T> K1fs, K2fs;
        momentum<complex<T> > K4s;
        get_Kif_coeffs(K1fs,K2fs,K4s);

        //We first select where the first solution should appear and then use the T(q) test to
        //  find the overall ordering relative to the stored solution
        if(polepos==0){
        	momentum<complex<T> > lcmp=tp.alp2*tp.K1flatc.P()+tp.alp1*tp.K2flatc.P()+tp.Nboxpoles[0]*tp.vec1c.P()+tp.alp1*tp.alp2/tp.Nboxpoles[0]*tp.vec2c.P();
        	ratio=((K2fs.Lt()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.Lt())*(K1fs.L()*K2fs.L())
        			-(K2fs.L()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.L())*(K1fs.Lt()*K2fs.Lt()))/Tn;
        }
        else{
        	momentum<complex<T> > lcmp=tp.alp2*tp.K1flatc.P()+tp.alp1*tp.K2flatc.P()+tp.Nboxpoles[1]*tp.vec1c.P()+tp.alp1*tp.alp2/tp.Nboxpoles[1]*tp.vec2c.P();
        	ratio=-((K2fs.Lt()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.Lt())*(K1fs.L()*K2fs.L())
                    -(K2fs.L()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.L())*(K1fs.Lt()*K2fs.Lt()))/Tn;
        }
        // This should always be real and either 1 or -1 the order of the returned solutions is reversed if it is -1
        if((ratiofac*ratio).real()>T(0)){
        	get_coeff(tp.boxcoeffs[0],tp.boxcoeffs[1]);
        }
        else{
        	get_coeff(tp.boxcoeffs[1],tp.boxcoeffs[0]);
        }

        //
        // Now compute the subtraction terms for the triangle
        //

        // Get the points to evaluate at at the correct precision
        const complex<T>* eval_pts;
        box_Darren_eval_points<TPOINTSTRI>::get_eval_pts(eval_pts);
    	if(nsols==2){
    		for(int icirc=0; icirc<TPOINTSTRI; icirc++){
    			coeffsret[icirc]+=complex<T>(0,-1)/(tp.denfac1*(tp.Nboxpoles[1]-tp.Nboxpoles[0]))*(tp.Nboxpoles[0]*tp.boxcoeffs[0]/(eval_pts[icirc]-tp.Nboxpoles[0])
    					-tp.Nboxpoles[1]*tp.boxcoeffs[1]/(eval_pts[icirc]-tp.Nboxpoles[1]));
    		}
    		coeffsret[TPOINTSTRI]+=complex<T>(0,4)*(tp.boxcoeffs[0]-tp.boxcoeffs[1])/(tp.denfac1*(tp.Nboxpoles[1]-tp.Nboxpoles[0]));
    	}
    	else {// We can only have one solution
    		for(int icirc=0; icirc<TPOINTSTRI; icirc++){
    			coeffsret[icirc]+=complex<T>(0,1)*tp.boxcoeffs[polepos]/(tp.denfac1*(eval_pts[icirc]-tp.Nboxpoles[polepos]));
    		}
    		// One solution means that this cannot be a three-mass triangle and so we do not need the t=0 contribution
    	}

}

template <class cutDbase,class BubbleSpecs> void bubble_Darren<cutDbase,BubbleSpecs>::GenIndicesBub(const vector<int>& ind, vector<int*>& pmoml, vector<int*>& pmmoml)
{
        //Insert the legs for the two cut tree amplitudes
        size_t mm;
        for(mm=1;mm<=cutDbase::corner_size(1);mm++){
                (indlst[0])[mm]=ind[cutDbase::corner_ind(1,mm)-1];
        }
        pmmoml[0]=&((indlst[0])[0]);
        pmoml[1]=&((indlst[0])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(2);mm++){
                (indlst[1])[mm]=ind[cutDbase::corner_ind(2,mm)-1];
        }
        pmmoml[1]=&((indlst[1])[0]);
        pmoml[0]=&((indlst[1])[mm]);
}

template <class cutDbase> void Do_GenIndicesBub(cutDbase& cut,const vector<int>& ind, vector<int*>& pmoml, vector<int*>& pmmoml, std::vector<int> (&indlst)[2])
{
        //Insert the legs for the two cut tree amplitudes
        size_t mm;
        for(mm=1;mm<=cut.corner_size(1);mm++){
                (indlst[0])[mm]=ind[cut.corner_ind(1,mm)-1];
        }
        pmmoml[0]=&((indlst[0])[0]);
        pmoml[1]=&((indlst[0])[mm]);

        for(mm=1;mm<=cut.corner_size(2);mm++){
                (indlst[1])[mm]=ind[cut.corner_ind(2,mm)-1];
        }
        pmmoml[1]=&((indlst[1])[0]);
        pmoml[0]=&((indlst[1])[mm]);
}


template <class cutDbase, class TriangleSpecs> void triangle_Darren<cutDbase,TriangleSpecs>::GenIndicesTri(const vector<int>& ind, vector<int*>& pmoml, vector<int*>& pmmoml, vector<int>& leg)
{
        size_t mm;
        for(mm=1;mm<=cutDbase::corner_size(1);mm++){
                (indlst[0])[mm]=ind[cutDbase::corner_ind(1,mm)-1];
        }
        pmmoml[leg[0]]=&((indlst[0])[0]);
        pmoml[leg[1]]=&((indlst[0])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(2);mm++){
                (indlst[1])[mm]=ind[cutDbase::corner_ind(2,mm)-1];
        }
        pmmoml[leg[1]]=&((indlst[1])[0]);
        pmoml[leg[2]]=&((indlst[1])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(3);mm++){
                (indlst[2])[mm]=ind[cutDbase::corner_ind(3,mm)-1];
        }
        pmmoml[leg[2]]=&((indlst[2])[0]);
        pmoml[leg[0]]=&((indlst[2])[mm]);
}


template <class cutDbase, int CPOINTS, int TPOINTSTRI> void box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::GenIndicesBox(const vector<int>& IND, vector<int*>& pmoml, vector<int*>& pmmoml, vector<int>& leg)
{
        size_t mm;
        for(mm=1;mm<=cutDbase::corner_size(1);mm++){
                (indlst[0])[mm]=IND[cutDbase::corner_ind(1,mm)-1];
        }
        pmmoml[leg[0]]=&((indlst[0])[0]);
        pmoml[leg[1]]=&((indlst[0])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(2);mm++){
                (indlst[1])[mm]=IND[cutDbase::corner_ind(2,mm)-1];
        }
        pmmoml[leg[1]]=&((indlst[1])[0]);
        pmoml[leg[2]]=&((indlst[1])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(3);mm++){
                (indlst[2])[mm]=IND[cutDbase::corner_ind(3,mm)-1];
        }
        pmmoml[leg[2]]=&((indlst[2])[0]);
        pmoml[leg[3]]=&((indlst[2])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(4);mm++){
                (indlst[3])[mm]=IND[cutDbase::corner_ind(4,mm)-1];
        }
        pmmoml[leg[3]]=&((indlst[3])[0]);
        pmoml[leg[0]]=&((indlst[3])[mm]);
}

template <class cutDbase, int CPOINTS, int TPOINTSTRI> template <class T> void box_Darren<cutDbase,CPOINTS,TPOINTSTRI>::get_coeffs_fn(momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp)
{

    	//
    	// Compute the momenta of the four corners of the box
    	//

    	momentum<complex<T> > K1sum_mom(mc.mom(ind[cutDbase::corner_ind(_k1leg,1)-1]));
    	size_t k1size=2;
    	for(;k1size<=cutDbase::corner_size(_k1leg);k1size++){
    		K1sum_mom+=mc.mom(ind[cutDbase::corner_ind(_k1leg,k1size)-1]);
    	}
    	momentum<complex<T> > K2sum_mom(-mc.mom(ind[cutDbase::corner_ind(_k2leg,1)-1]));
    	for(size_t kiter=2;kiter<=cutDbase::corner_size(_k2leg);kiter++){
    		K2sum_mom-=mc.mom(ind[cutDbase::corner_ind(_k2leg,kiter)-1]);
    	}
    	momentum<complex<T> > K4sum_mom(mc.mom(ind[cutDbase::corner_ind(_k4leg,1)-1]));
    	for(size_t kiter=2;kiter<=cutDbase::corner_size(_k4leg);kiter++){
    		K4sum_mom+=mc.mom(ind[cutDbase::corner_ind(_k4leg,kiter)-1]);
    	}
    	complex<T> S1=K1sum_mom.square();
    	complex<T> S2=K2sum_mom.square();
    	complex<T> S4=K4sum_mom.square();

    	complex<T> K1K2=K1sum_mom*K2sum_mom;
        complex<T> vec_alp[6];
        complex<T> alp1(1,0), alp2(1,0), f1(1,0), f2(1,0);
        Cmom<T> vec1c, vec2c;
        Cmom<T> K1flatc, K2flatc;
        complex<T> Nboxpoles[2];
        int polepos,nsols;
        if(_masslessleg_type!=0){
        	//K2 is massless so create spinors for the K2 leg
        	Cmom<T> K2m(K2sum_mom);
        	// For gamma we have
        	complex<T> gamfac, gammap(T(2)*K1K2);
        	if(_massless_K1){//In the four point case we can have a massless K1
        		alp1=complex<T>(0,0);
        		f2=complex<T>(2,0);//This is gives the correct answer base don the way the code is setup below
        		gamfac=complex<T>(1,0);
        		K1flatc=Cmom<T>(K1sum_mom);
        	}
        	else{
        		f2=gammap/S1;
        		gamfac=T(1)/f2;
        		K1flatc=Cmom<T>(K1sum_mom-gamfac*K2sum_mom);
        	}
        	alp2=complex<T>(0,0);

        	if(_masslessleg_type>0){
        		// Rescale the spinors of K2 correctly
           		K2flatc=Cmom<T>(K2m.Lt(),gamfac*K2m.L());

           		//Set up the rescaling factors for l1 and l2
        		vec_alp[0]=alp1;
        		vec_alp[1]=complex<T>(0,0);
           		vec_alp[2]=complex<T>(-1,0);
        		vec_alp[3]=complex<T>(0,0);
        		vec_alp[4]=complex<T>(0,0);
        		vec_alp[5]=(T(1)-f2);

        		vec1c=Cmom<T>(K2flatc.Lt(),K1flatc.L());
        		vec2c=Cmom<T>(K1flatc.Lt(),K2flatc.L());
        	}
        	else{
        		// Rescale the spinors of K2 correctly
           		K2flatc=Cmom<T>(gamfac*K2m.Lt(),K2m.L());

           		//Set up the rescaling factors for l1 and l2
        		vec_alp[0]=complex<T>(0,0);
        		vec_alp[1]=alp1;
        		vec_alp[2]=complex<T>(0,0);
        		vec_alp[3]=complex<T>(-1,0);
        		vec_alp[4]=(T(1)-f2);
        		vec_alp[5]=complex<T>(0,0);

           		vec1c=Cmom<T>(K1flatc.Lt(),K2flatc.L());
           		vec2c=Cmom<T>(K2flatc.Lt(),K1flatc.L());
        	}

        	//Set up nboxpart1 so that they it matches the L_4=l_i-K_4 (we do not always have L_4=l-K_4)
        	complex<T> nboxpart1=S4-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom)*K4sum_mom);// l_4=l_1-K4

        	//We now calculate the pole solutions
            complex<T> denfac1=T(2)*(vec1c.P()*K4sum_mom);
            if(nboxpart1.real()<T(0)){
            	Nboxpoles[0]=complex<T>(0,0);
            	Nboxpoles[1]=nboxpart1/denfac1;
            	polepos=1;
            }
            else
            {
            	// Choose the non vanishing solution as the "first" when we have only a single pole i.e. one and two mass cases
            	Nboxpoles[0]=nboxpart1/denfac1;
            	Nboxpoles[1]=complex<T>(0,0);
            	polepos=0;
            }
            nsols=1;
        }
        else{
        	// Compute the gamma of K1flatp and K2flatp so we can define f1 and f2.
     	   complex<T> gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
     	   complex<T> f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
     	   complex<T> f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
     	   // Now convert gammap into gamma for the rest of the computation
    // 	   gamma=gammap/(f1*f2);
     	   vec_alp[0]=alp1;
     	   vec_alp[1]=alp2;
     	   //Set up the rescaling factors for l1 and l2
     	   vec_alp[3]=(T(1)-(gammap*gammap-S1*S2)/(S2*(gammap-S1)));
     	   vec_alp[2]=(T(1)-(gammap*gammap-S1*S2)/(gammap*(gammap-S2)));
     	   vec_alp[5]=(T(1)-(gammap*gammap-S1*S2)/(gammap*(gammap-S1)));
     	   vec_alp[4]=(T(1)-(gammap*gammap-S1*S2)/(S1*(gammap-S2)));

     	   //Create the unscaled K1flat and K2flat
     	   complex<T> denalpha=T(1)/(pow(gammap,2)-S1*S2);
     	   Cmom<T> K1flatcI(denalpha*gammap*(gammap*K1sum_mom-S1*K2sum_mom));
     	   Cmom<T> K2flatcI(denalpha*gammap*(gammap*K2sum_mom-S2*K1sum_mom));

     	   //Now rescale one of the spinors to increase numerical stability
     	   K1flatc=Cmom<T>(K1flatcI.Lt(),(T(1)/f1)*K1flatcI.L());
     	   K2flatc=Cmom<T>((T(1)/f2)*K2flatcI.Lt(),K2flatcI.L());

     	   // Now set up the vectors for the last two terms in the l basis
     	   vec1c=Cmom<T>(K1flatc.Lt(),K2flatc.L());
     	   vec2c=Cmom<T>(K2flatc.Lt(),K1flatc.L());

     	   //Set up nboxpart1 so that they it matches the L_4=l_i-K_4 (we do not always have L_4=l-K_4)
     	   complex<T> nboxpart1=S4-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom)*K4sum_mom);// l_4=l_1-K4

     	   //We now calculate the pole solutions
     	   complex<T> denfac1=T(2)*(vec1c.P()*K4sum_mom);
     	   complex<T> nboxpart2=sqrt(pow(nboxpart1,2)-T(8)*alp1*alp2*denfac1*(vec2c.P()*K4sum_mom));

     	   // We always have two solutions in the three mass case
     	   Nboxpoles[0]=(nboxpart1+nboxpart2)/(T(2)*denfac1);
     	   Nboxpoles[1]=(nboxpart1-nboxpart2)/(T(2)*denfac1);
     	   polepos=1;
     	   nsols=2;
        }

        // Set up the vectors for the trees processes at the corners
        size_t mm;
        for(mm=1;mm<=cutDbase::corner_size(1);mm++){
                (indlst[0])[mm]=ind[cutDbase::corner_ind(1,mm)-1];
        }
        for(mm=1;mm<=cutDbase::corner_size(2);mm++){
                (indlst[1])[mm]=ind[cutDbase::corner_ind(2,mm)-1];
        }
        for(mm=1;mm<=cutDbase::corner_size(3);mm++){
                (indlst[2])[mm]=ind[cutDbase::corner_ind(3,mm)-1];
        }
        for(mm=1;mm<=cutDbase::corner_size(4);mm++){
                (indlst[3])[mm]=ind[cutDbase::corner_ind(4,mm)-1];
        }

        vector<int>::iterator pmoml[4], pmmoml[4];
        pmoml[1]=(indlst[_k1leg-1]).end()-1;
        pmoml[2]=(indlst[_k4leg-1]).end()-1;
        pmoml[3]=(indlst[_k3leg-1]).end()-1;
        pmoml[0]=(indlst[_k2leg-1]).end()-1;
        pmmoml[0]=(indlst[_k1leg-1]).begin();
        pmmoml[1]=(indlst[_k4leg-1]).begin();
        pmmoml[2]=(indlst[_k3leg-1]).begin();
        pmmoml[3]=(indlst[_k2leg-1]).begin();


        //Sum over the number of poles we had
        momentum<complex<T> > mompt1=alp2*K1flatc.P()+alp1*K2flatc.P();
        //Sum over the number of poles
        int polepos_new=polepos-nsols+1;
        complex<T> boxcoeffs[2];
        for(int npole=0; npole<nsols; npole++){
            *(pmoml[0])=mc.insert(vec_alp[1]*vec2c.Lt()*(T(1)/Nboxpoles[polepos_new+npole])+vec1c.Lt(),vec1c.L()*Nboxpoles[polepos_new+npole]+vec_alp[0]*vec2c.L());
            *(pmoml[1])=mc.insert(vec_alp[2]*vec2c.Lt()*(T(1)/Nboxpoles[polepos_new+npole])+vec1c.Lt(),vec1c.L()*Nboxpoles[polepos_new+npole]+vec_alp[3]*vec2c.L());
            *(pmoml[2])=mc.insert(mc.mom(*(pmoml[1]))-K4sum_mom);
            *(pmoml[3])=mc.insert(vec_alp[4]*vec2c.Lt()*(T(1)/Nboxpoles[polepos_new+npole])+vec1c.Lt(),vec1c.L()*Nboxpoles[polepos_new+npole]+vec_alp[5]*vec2c.L());
            *(pmmoml[0])=mc.insert(-mc[*(pmoml[0])]);
            *(pmmoml[1])=mc.insert(-mc[*(pmoml[1])]);
            *(pmmoml[2])=mc.insert(-mc[*(pmoml[2])]);
            *(pmmoml[3])=mc.insert(-mc[*(pmoml[3])]);

            //Construct the three-particle cut subtraction at the above momenta and compute the amplitude
            //If there is only one solution we want to store this in the correct location
            boxcoeffs[polepos_new+npole]=Amp_safe(this->eval_tree(1,mc,indlst[0])*this->eval_tree(2,mc,indlst[1])*this->eval_tree(3,mc,indlst[2])*this->eval_tree(4,mc,indlst[3]));
        }
        //Store the coefficients
        set_coeff(0,boxcoeffs[0]);
        set_coeff(1,boxcoeffs[1]);
        set_D0coeff((boxcoeffs[0]+boxcoeffs[1])/complex<T>(0,-2));

        //Store information allowing us to decide how to choose how to map t_+ and t_- to d_+ and d_- later on.
        // We do this by borrowing the method hidden in the OPP approach basically we can rewrite the OPP expression
        //  for the subtracted box terms as
        //  (1/2)(d_+(1+T(l'_+)/T(l_+))+d_-(1-T(l'_+)/T(l_+))/(t-t_+)+d_+(1+T(l'_-)/T(l_+))+d_-(1-T(l'_-)/T(l_+))/(t-t_-))
        //  with T(l)=Tr(l,K1f,K2f,K4) and l the original solution of l to find d_+ and d_- and l' the new solution for
        //  l we are trying to map d_+ and d_- to.
        //  The ratio of T(l'_+)/T(l_+) will either +1 or -1 depending upon the order of the momenta entering l
        //  and the sign will correctly choose the mapping of t_i->d_i.
        set_Kif_coeffs(K1flatc,K2flatc,K4sum_mom);
        if(polepos==0){
                // compute T(l_+)
                set_denfac((K2flatc.Lt()*smatrix<T>(K4sum_mom)*smatrix<T>(mc.mom(*(pmoml[0])))*K1flatc.Lt())*(K1flatc.L()*K2flatc.L())
                                -(K2flatc.L()*smatrix<T>(K4sum_mom)*smatrix<T>(mc.mom(*(pmoml[0])))*K1flatc.L())*(K1flatc.Lt()*K2flatc.Lt()));
        }
        else{
                //Flip the sign as this is T(l_-) and we store T(l_+)=-T(l_-)
                set_denfac(-(K2flatc.Lt()*smatrix<T>(K4sum_mom)*smatrix<T>(mc.mom(*(pmoml[0])))*K1flatc.Lt())*(K1flatc.L()*K2flatc.L())
                                +(K2flatc.L()*smatrix<T>(K4sum_mom)*smatrix<T>(mc.mom(*(pmoml[0])))*K1flatc.L())*(K1flatc.Lt()*K2flatc.Lt()));
        }

        //Mark that we have computed these using the original mcID not the sub_mom_conf one from the mc that was passed to this function
        mcID=tp.mcID;
        indID=ind;

}


template <class cutDbase, class TriangleSpecs> C triangle_Darren<cutDbase,TriangleSpecs>::eval(mom_conf& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a triangle before computing the bubble");
        }
        return this->template get_symmetry_factor<R>()*C0coeff;
}

template <class cutDbase, class TriangleSpecs> CHP triangle_Darren<cutDbase,TriangleSpecs>::eval(mom_conf_HP& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a triangle before computing the bubble");
        }
        return this->template get_symmetry_factor<RHP>()*C0coeff_HP;
}

template <class cutDbase, class TriangleSpecs> CVHP triangle_Darren<cutDbase,TriangleSpecs>::eval(mom_conf_VHP& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a triangle before computing the bubble");
        }
        return this->template get_symmetry_factor<RVHP>()*C0coeff_VHP;
}

#if BH_USE_GMP
template <class cutDbase, class TriangleSpecs> CGMP triangle_Darren<cutDbase,TriangleSpecs>::eval(momentum_configuration<RGMP>& mc,const vector<int>& ind)
{
        // If we have not previously computed the box compute it
        if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING2(*this," : Attempting to compute a triangle before computing the bubble");
        }
        return this->template get_symmetry_factor<RGMP>()*C0coeff_GMP;
}

#endif




template <class cutDbase, class BubbleSpecs> template <class T> complex<T> bubble_Darren<cutDbase,BubbleSpecs>::get_coeffs(momentum_configuration<T>& mc, const vector<int>& ind)
{
	//We will store all information that needs to be passed between the bubble, triangle and box in the coeffparam structure
	coeffparam<T,BubbleSpecs::CPOINTS> tp;

	//Mark that we have computed these and save the actual mcID so we can use this in the triangle and boxes later
	tp.mcID=mc.get_ID();
	mcID=tp.mcID;
	indID=ind;

	// Create a sub momentum configuration so that at the end of the computation we can dump all the
	//  newly created momenta that are no longer needed
	sub_momentum_configuration<T> sub_mc(mc);

    vector<complex<T> > trees_result_1(BubbleSpecs::YPOINTS*BubbleSpecs::TPOINTSBUB);
    vector<complex<T> > trees_result_2(BubbleSpecs::YPOINTS*BubbleSpecs::TPOINTSBUB);


    typename BubbleSpecs::CornerTreeStrategy* ptr= static_cast<typename BubbleSpecs::CornerTreeStrategy*>(this);
    ptr->fill_trees(sub_mc,ind,trees_result_1,trees_result_2,tp,*this);

//    copy(trees_result_1.begin(),trees_result_1.end(),ostream_iterator<complex<T> >(cout,"\n"));
//    copy(trees_result_2.begin(),trees_result_2.end(),ostream_iterator<complex<T> >(cout,"\n"));


    const complex<T>* eval_pts;
    MomentaEvaluatorType::get_eval_pts(eval_pts);


    complex<T> ampl[BubbleSpecs::YPOINTS]={complex<T>(0,0)}, ampl_err[BubbleSpecs::YPOINTS]={complex<T>(0,0)}, temp[BubbleSpecs::YPOINTS];

    for(int j=0;j<BubbleSpecs::YPOINTS;j++){
            for(int i=0;i<BubbleSpecs::TPOINTSBUB;i++){
                int current_index=j*BubbleSpecs::TPOINTSBUB+i;
                    temp[j]=trees_result_1[current_index]*trees_result_2[current_index];
					BH_DEBUG_PRINT(temp[j]);
                    ampl[j]+=temp[j];
           			ampl_err[j]+=temp[j]*eval_pts[i];
        }
    }

	//Sum over all the triangles getting their subtraction pieces
	complex<T> trisuby(0,0), amp_error(0,0);
	RVHP tri_sub_acc_comp(to_double(MaxDigitsE<T>()));
	for(size_t tri=1; tri <= cutDbase::daughters_nbr(); tri++){
		tp.tri_corner=cutDbase::get_opened_corner(tri);
#ifndef NDEBUG
		typename bubble_daughter_type<cutDbase,BubbleSpecs>::type* wt= dynamic_cast<typename bubble_daughter_type<cutDbase,BubbleSpecs>::type* >(cutDbase::get_daughter(tri));
		assert(wt);
#else
		typename bubble_daughter_type<cutDbase,CPOINTS,TPOINTSBUB,TPOINTSTRI,YPOINTS>::type* wt= static_cast<typename bubble_daughter_type<cutDbase,CPOINTS,TPOINTSBUB,TPOINTSTRI,YPOINTS>::type* >(cutDbase::get_daughter(tri));
#endif
		trisuby+=wt->get_sub_terms(sub_mc,ind,tp);
		BH_DEBUG_PRINT(trisuby);
		amp_error+=tp.amp_err;
		if(tp.tri_sub_acc<tri_sub_acc_comp){
			tri_sub_acc_comp=tp.tri_sub_acc;
		}
	}



	RVHP accuracy;
//	typename BubbleSpecs::BubbleCombiner* BCptr= static_cast<typename BubbleSpecs::BubbleCombiner*>(this);
	complex<T> result= BubbleSpecs::BubbleCombiner::combine(ampl,
			ampl_err,
			trisuby,
			amp_error,
			tri_sub_acc_comp,
			accuracy);

	this->set_accuracy(to_double(log(accuracy))/log(10.));

	set_B0coeff(result);

	//We have already calculated this so return the value
	return result;
}


template <class cutDbase, class BubbleSpecs> C bubble_Darren<cutDbase,BubbleSpecs>::eval(mom_conf& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<R>()*get_coeffs(mc,ind);
}

template <class cutDbase, class BubbleSpecs> CHP bubble_Darren<cutDbase,BubbleSpecs>::eval(mom_conf_HP& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<RHP>()*get_coeffs(mc,ind);
}

template <class cutDbase, class BubbleSpecs> CVHP bubble_Darren<cutDbase,BubbleSpecs>::eval(mom_conf_VHP& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<RVHP>()*get_coeffs(mc,ind);
}
#if BH_USE_GMP
template <class cutDbase, class BubbleSpecs> std::complex<RGMP> bubble_Darren<cutDbase,BubbleSpecs>::eval(momentum_configuration<RGMP>& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<RGMP>()*get_coeffs(mc,ind);
}
#endif


/*
 *
 *
 * Specialized versions
 *
 *
 *
 */

#if ugly_specialization

template <class BASE> void  bubble_Darren<
	BASE,
	Normal_Bubble_Specification<BASE> >::init(){
	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	indlst[0].assign(BASE::corner_size(1)+2,0);
	indlst[1].assign(BASE::corner_size(2)+2,0);

	//Set up the eval_params for the corners
	_ep[0]=new eval_param<R>(BASE::corner_size(1)+2);
	_ep_HP[0]=new eval_param<RHP>(BASE::corner_size(1)+2);
	_ep_VHP[0]=new eval_param<RVHP>(BASE::corner_size(1)+2);

	_ep[1]=new eval_param<R>(BASE::corner_size(2)+2);
	_ep_HP[1]=new eval_param<RHP>(BASE::corner_size(2)+2);
	_ep_VHP[1]=new eval_param<RVHP>(BASE::corner_size(2)+2);

	// We can have the situation when there are bubbles but no triangles but we
	//  still need the static initalised coefficient points from the triangle so
	//  we initalise it here by creating an object of the correct type
	triangle_Darren_eval_points<7,TRIPOINTS_CD> init_temp;

}

template <class cutDbase> void bubble_Darren<cutDbase,Normal_Bubble_Specification<cutDbase> >::GenIndicesBub(const vector<int>& ind, vector<int*>& pmoml, vector<int*>& pmmoml)
{
        //Insert the legs for the two cut tree amplitudes
        size_t mm;
        for(mm=1;mm<=cutDbase::corner_size(1);mm++){
                (indlst[0])[mm]=ind[cutDbase::corner_ind(1,mm)-1];
        }
        pmmoml[0]=&((indlst[0])[0]);
        pmoml[1]=&((indlst[0])[mm]);

        for(mm=1;mm<=cutDbase::corner_size(2);mm++){
                (indlst[1])[mm]=ind[cutDbase::corner_ind(2,mm)-1];
        }
        pmmoml[1]=&((indlst[1])[0]);
        pmoml[0]=&((indlst[1])[mm]);
}

template <class cutDbase> C bubble_Darren<
	cutDbase,
	Normal_Bubble_Specification<cutDbase>
 >::eval(mom_conf& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<R>()*get_coeffs(mc,ind);
}

template <class cutDbase> CHP bubble_Darren<
	cutDbase,
	Normal_Bubble_Specification<cutDbase>
>::eval(mom_conf_HP& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<RHP>()*get_coeffs(mc,ind);
}

template <class cutDbase> CVHP bubble_Darren<cutDbase,
Normal_Bubble_Specification<cutDbase> >::eval(mom_conf_VHP& mc, const vector<int>& ind){
        return this->template get_symmetry_factor<RVHP>()*get_coeffs(mc,ind);
}

template <class cutDbase> template <class T> complex<T> bubble_Darren<cutDbase,
Normal_Bubble_Specification<cutDbase> >::get_coeffs(momentum_configuration<T>& mc, const vector<int>& ind)
{
	typedef Normal_Bubble_Specification<cutDbase> BubbleSpecs;
	//We will store all information that needs to be passed between the bubble, triangle and box in the coeffparam structure
        coeffparam<T,7> tp;

        //Mark that we have computed these and save the actual mcID so we can use this in the triangle and boxes later
        tp.mcID=mc.get_ID();
        mcID=tp.mcID;
        indID=ind;

        // Create a sub momentum configuration so that at the end of the computation we can dump all the
        //  newly created momenta that are no longer needed
        sub_momentum_configuration<T> sub_mc(mc);
#if use_old
        // Construct K1 by summing over all the legs of c(1) using ind to get their location in the momconf
        momentum<complex<T> > K1sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
        for(size_t k1iter=1;k1iter<=cutDbase::corner_size(1);k1iter++){
                K1sum_mom+=sub_mc.mom(ind[cutDbase::corner_ind(1,k1iter)-1]);
        }
        tp.K1=K1sum_mom;
        tp.S1=tp.K1.square();

        //Construct chi so that all the coefficients of l^{\mu} are of order 1
        // this is achieved by scaling the vector (1,0,1,0) (we want to avoid a chi proportional to one of the external legs)
        // so that gamma=S1
        //	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,-1));
        //	tp.chic=Cmom<T>(tp.S1/(T(2)*(chi_init*tp.K1))*chi_init);
        momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,-1));

    	// Now rescale to the actual chi
    	tp.chic=Cmom<T>(tp.S1/(T(2)*(chi_init*tp.K1))*chi_init);
        tp.gammab=tp.S1;
        tp.K1flatbc=Cmom<T>(tp.K1-tp.chic.P());

#endif


#if use_older == 0
        vector<complex<T> > trees_result_1(2*BUBPOINTS_CD);
        vector<complex<T> > trees_result_2(2*BUBPOINTS_CD);

        this->template fill_trees(sub_mc,ind,trees_result_1,trees_result_2,tp,*this);

        copy(trees_result_1.begin(),trees_result_1.end(),ostream_iterator<complex<T> >(cout,"\n"));
        copy(trees_result_2.begin(),trees_result_2.end(),ostream_iterator<complex<T> >(cout,"\n"));

        const complex<T>* eval_pts;
        MomentaEvaluatorType::get_eval_pts(eval_pts);


#else
        // ORIGINAL VERSION
        //Now set up the three points we evaluate around the circle
        const complex<T>* eval_pts;
        bubble_Darren_eval_points<BUBPOINTS_CD,2>::get_eval_pts(eval_pts);


    	const complex<T>* ypoint;
    	bubble_Darren_eval_points<BUBPOINTS_CD,2>::get_yeval_pts(ypoint);

    	MomentaEvaluatorType::MomentaGridType momenta[2];
    	MomentaEvaluatorType::MomentaGridType momentam[2];

    	this->template compute_momenta<T,7,bubble_Darren<cutDbase,7,BUBPOINTS_CD,TRIPOINTS_CD,2> >(sub_mc,ind,momenta,momentam, tp,*this);


    	//We check the error of our computation by looking at the 1/t component of the bubble numerator structure.
    	// Due to our parameterization this should always have a coefficient of zero.
        vector<int*> pmoml(2),pmmoml(2);
        GenIndicesBub(ind, pmoml, pmmoml);

        vector<complex<T> > trees_result_1(2*BUBPOINTS_CD);
        vector<complex<T> > trees_result_2(2*BUBPOINTS_CD);

        for(int j=0;j<2;j++){
                for(int i=0;i<BUBPOINTS_CD;i++){
#if use_old
                    //Construct the cut-momenta
                	*(pmoml[0])=sub_mc.insert((ypoint[j]/eval_pts[i])*tp.K1flatbc.Lt()+tp.chic.Lt(),eval_pts[i]*tp.K1flatbc.L()+(T(1)-ypoint[j])*tp.chic.L());
                        *(pmoml[1])=sub_mc.insert((ypoint[j]-T(1))*tp.K1flatbc.Lt()+eval_pts[i]*tp.chic.Lt(),tp.K1flatbc.L()-(ypoint[j]/eval_pts[i])*tp.chic.L());
                        *(pmmoml[0])=sub_mc.insert(-sub_mc.Lt(*(pmoml[0])),sub_mc.L(*(pmoml[0])));
                        *(pmmoml[1])=sub_mc.insert(-sub_mc.Lt(*(pmoml[1])),sub_mc.L(*(pmoml[1])));
                        //Construct the two-particle cut at the momenta above
#else
            			//Construct the cut-momenta
            			*(pmoml[0])=momenta[j][i][0];
            			*(pmoml[1])=momenta[j][i][1];
            			*(pmmoml[0])=momentam[j][i][0];
            			*(pmmoml[1])=momentam[j][i][1];
            			//Construct the two-particle cut at the momenta above
#endif
            			int current_index=j*BUBPOINTS_CD+i;
                        trees_result_1[current_index]=eval_tree(1,sub_mc,indlst[0]);
                        trees_result_2[current_index]=eval_tree(2,sub_mc,indlst[1]);
            }
        }
#endif



        complex<T> ampl[2]={complex<T>(0,0)}, ampl_err[2]={complex<T>(0,0)}, temp[2];


        for(int j=0;j<2;j++){
                for(int i=0;i<BUBPOINTS_CD;i++){
                    int current_index=j*BUBPOINTS_CD+i;
                        temp[j]=trees_result_1[current_index]*trees_result_2[current_index];
                        ampl[j]+=temp[j];
               			ampl_err[j]+=temp[j]*eval_pts[i];
            }
        }


//NEW VERSION
//        // At most we can have up to 2 powers of alpha2 or alpha3 so we need 3 values of each naively
//        //  first test if setting alpha2 and alpha3 to zero works
//       vector<int*> pmoml(2),pmmoml(2);
//       GenIndicesBub(ind, pmoml, pmmoml);
//
//       //Generate the eval points for the two alaphas, they are symetric so we use the same points for both
//       complex<T> eval_alpha_pts[3];
//       for(int i=0;i<3;i++){
//    	   eval_alpha_pts[i]=exp(complex<T>(0,2)*pi<T>()*T(i)/T(3));
//       }
//       complex<T> eval_alpha_pts_matrix[3*3];
//       for(int j=0;j<3;j++){
//    	   for(int i=0;i<3;i++){
//    		   eval_alpha_pts_matrix[i+j*3]=pow(eval_alpha_pts[i],-j)/T(3);
//    	   }
//       }
//
//       momentum<complex<T> > n1=complex<T>(0,1)*(tp.K1flatbc.P()-tp.chic.P())/sqrt(tp.S1);
//       Cmom<T> nA(tp.K1flatbc.L(),tp.chic.Lt());
//       Cmom<T> nB(tp.chic.L(),tp.K1flatbc.Lt());
//       momentum<complex<T> > n2=complex<T>(0,1)*(nA.P()+nB.P())/sqrt(tp.gammab);
//       momentum<complex<T> > n3=(nA.P()-nB.P())/sqrt(tp.gammab);
//       complex<T> amp(0,0);
//       for(int j=0;j<3;j++){
//    	   for(int i=0;i<3;i++){
//    	       *(pmoml[0])=sub_mc.insert(tp.K1/complex<T>(2,0)
//							   +sqrt(-tp.S1/T(4)-eval_alpha_pts[i]*eval_alpha_pts[i]-eval_alpha_pts[j]*eval_alpha_pts[j])*n1
//							   +eval_alpha_pts[i]*n2+eval_alpha_pts[j]*n3);
//    	       *(pmoml[1])=sub_mc.insert(sub_mc.mom(*(pmoml[0]))-tp.K1);
//    	       *(pmmoml[0])=sub_mc.insert(-sub_mc.Lt(*(pmoml[0])),sub_mc.L(*(pmoml[0])));
//    	       *(pmmoml[1])=sub_mc.insert(-sub_mc.Lt(*(pmoml[1])),sub_mc.L(*(pmoml[1])));
//    	       _MESSAGE4("l0^2=",sub_mc[*(pmoml[0])].P().square(),", l0^2=",sub_mc[*(pmoml[1])].P().square());
//    	       //Construct the two-particle cut at the momenta above
//    	       amp+=eval_tree(1,sub_mc,indlst[0])*eval_tree(2,sub_mc,indlst[1]);
//    	   }
//       }
//
////       *(pmoml[0])=sub_mc.insert((-tp.K1flatbc.P()+complex<T>(3)*tp.chic.P())/complex<T>(2));
////       *(pmoml[1])=sub_mc.insert(sub_mc.mom(*(pmoml[0]))-tp.K1);
////       *(pmmoml[0])=sub_mc.insert(-sub_mc.Lt(*(pmoml[0])),sub_mc.L(*(pmoml[0])));
////       *(pmmoml[1])=sub_mc.insert(-sub_mc.Lt(*(pmoml[1])),sub_mc.L(*(pmoml[1])));
////       //Construct the two-particle cut at the momenta above
////       complex<T> amp0=eval_tree(1,sub_mc,indlst[0])*eval_tree(2,sub_mc,indlst[1]);

        //Sum over all the triangles getting their subtraction pieces
        complex<T> trisuby(0,0), amp_error(0,0);
        for(size_t tri=1; tri <= cutDbase::daughters_nbr(); tri++){
                tp.tri_corner=cutDbase::get_opened_corner(tri);
#ifndef NDEBUG
                typename bubble_daughter_type<cutDbase,Normal_Bubble_Specification<cutDbase> >::type* wt= dynamic_cast<typename bubble_daughter_type<cutDbase,Normal_Bubble_Specification<cutDbase> >::type* >(cutDbase::get_daughter(tri));
                assert(wt);
#else
                worker_triangleDarren* wt= static_cast<worker_triangleDarren*>(wcd);
#endif
                trisuby+=wt->get_sub_terms(sub_mc,ind,tp);
        		amp_error+=tp.amp_err;
        }



#if dont_use_combiner == 0
	double accuracy;
	typename BubbleSpecs::BubbleCombiner* BCptr= static_cast<typename BubbleSpecs::BubbleCombiner*>(this);
	complex<T> result= Normal_Bubble_Specification<cutDbase>::BubbleCombiner::combine(ampl,
			ampl_err,
			trisuby,
			amp_error,
			accuracy);

	this->set_accuracy(accuracy);

#else
//ORIGINAL VERSION
    	set_accuracy(-to_double(log(abs(amp_error-ampl_err[0]))/log(10.)));

        //Calculate the bubble coefficient including all subtractions
#if BUBPOINTS_CD==4

    	_PRINT(trisuby);
    	_PRINT(ampl[0]);
    	_PRINT(ampl[1]);
    	_PRINT((ampl[0]+T(3)*ampl[1]));
    	complex<T> result=(complex<T>(0,-1)*(ampl[0]+T(3)*ampl[1])+trisuby)/T(16);
#else
        complex<T> result=(complex<T>(0,-2)*(ampl[0]/T(3)+ampl[1])+trisuby)/T(8);
#endif /*BUBPOINTS==4*/

// NEW VERSION
//    	set_accuracy(to_double(T(16.)));// Temporarily turn off the error check
//
//    	_MESSAGE4("BUB :",*this,"=",amp);
//        //Calculate the bubble coefficient including all subtractions
//        complex<T> result=amp;
#endif
        set_B0coeff(result);

//    	_MESSAGE3(*this,"==",(complex<T>(0,-1)*(ampl[0]+T(3)*ampl[1])+trisuby)/T(16));
////    	_MESSAGE6("BUB ycoeff [0]:",ycoeffs[0],", ycoeff [1]:",ycoeffs[1],", ycoeff [2]:",ycoeffs[2]);
//    	_MESSAGE6("2 BUB ycoeffs [0] : ",ampl[0]," ycoeffs y=2/3 : ",ampl[1]," trisub y:",trisuby);
//    	_MESSAGE4("2 BUB err BUB:",ampl_err[0]," TRI:",amp_error);


        //We have already calculated this so return the value
        return result;
}

#endif

template <class TriangleType,int CPOINTS,int TPOINTSBUB> template <class T> complex<T> Normal_Triangle_Subtraction<TriangleType,CPOINTS,TPOINTSBUB>::get_sub_terms_work_3m(const TriangleType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi)
{


         complex<T> K1K2=tp.K1*tp.K2;


         // We construct the "gamma-" contributing terms by using the symmetry of the momentum param t->a01 a02/t
         //  to relate the gamma+ solution 1/t^3, 1/t^2 and 1/t terms to t^3, t^2 and t.
         momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
         momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

         complex<T> chiK2Klflat=T(2)*(chik1b*tp.K2);//mc.spab(tp.chi, tp.K2t, tp.K1flatb);
         complex<T> VPK1flatK2=tp.K1flatbc.P()*tp.K2;
         complex<T> VPchiK2=tp.chic.P()*tp.K2;
         complex<T> Delta3=(K1K2*K1K2-tp.S1*tp.S2);

         complex<T> ampfull[2]={complex<T>(0,0)}, invcircposp, invcircposm;
         complex<T> Nysolp, Nysolm;
         complex<T> reddenfac, ressqrt, resa, resb;
         complex<T> facover=T(-2)/tsi.gamma_old;
         complex<T> ampp, ampm, ypfac;
         complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermpbox, subtermmbox, subtermptri, subtermmtri;

         // Get the points to evaluate the circle on
         const complex<T>* eval_pts;
         bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::get_eval_pts(eval_pts);

         complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
         complex<T> fac2[4]={tsi.v1old.P()*tp.K1flatbc.P(),tsi.v1old.P()*tp.chic.P(),tsi.v1old.P()*k1bchi,tsi.v1old.P()*chik1b};

         // Reset the precision error gathering variable
         tp.amp_err=complex<T>(0,0);
         complex<T> temp[2];

    	 // Calculate the subtraction pieces
    	 for(int icirc=0; icirc<BUBPOINTS_STD; icirc++){
    		 // First compute the solutions for the poles in y
    		 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-tp.S2)*chiK2Klflat/Delta3+T(1)/T(4)*pow(chiK2Klflat,2)/Delta3;
    		 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat*T(1)/T(2);
    		 resb=sqrt(Delta3)*sqrt(ressqrt);

    		 Nysolp=(resa+resb)/chiK2Klflat;
    		 Nysolm=(resa-resb)/chiK2Klflat;


    		 // Now construct the t in the bubble parameterisation
    		 invcircposp=Nysolp/eval_pts[icirc]*(T(1)-Nysolp);
    		 invcircposm=Nysolm/eval_pts[icirc]*(T(1)-Nysolm);
    		 bubtsp=T(tsi.reverse)*facover*(Nysolp*fac1[0]
    			              +(T(1)-Nysolp)*fac1[1]
    						  +eval_pts[icirc]*fac1[2]
    						  +invcircposp*fac1[3]);

    		 bubtisp=T(tsi.reverse)*facover*(Nysolp*fac2[0]
    						  +(T(1)-Nysolp)*fac2[1]
    						  +eval_pts[icirc]*fac2[2]
    						  +invcircposp*fac2[3]);

    		 bubtsm=T(tsi.reverse)*facover*(Nysolm*fac1[0]
    						  +(T(1)-Nysolm)*fac1[1]
    						  +eval_pts[icirc]*fac1[2]
    						  +invcircposm*fac1[3]);

    		 bubtism=T(tsi.reverse)*facover*(Nysolm*fac2[0]
    						  +(T(1)-Nysolm)*fac2[1]
    						  +eval_pts[icirc]*fac2[2]
    						  +invcircposm*fac2[3]);

    		 // Reconstruct the subtracted numerator
    		 subtermptri=complex<T>(0,-1)*(pow(bubtsp,3)*tsi.orig_coeffs[0]+pow(bubtsp,2)*tsi.orig_coeffs[1]+bubtsp*tsi.orig_coeffs[2]+tsi.orig_coeffs[3]
    				            +bubtisp*tsi.orig_coeffs[4]+pow(bubtisp,2)*tsi.orig_coeffs[5]+pow(bubtisp,3)*tsi.orig_coeffs[6]);
    		 subtermmtri=complex<T>(0,-1)*(pow(bubtsm,3)*tsi.orig_coeffs[0]+pow(bubtsm,2)*tsi.orig_coeffs[1]+bubtsm*tsi.orig_coeffs[2]+tsi.orig_coeffs[3]
    				            +bubtism*tsi.orig_coeffs[4]+pow(bubtism,2)*tsi.orig_coeffs[5]+pow(bubtism,3)*tsi.orig_coeffs[6]);

    		 //Now add the box contributions to the triple cut subtraction term
    		 complex<T> subtermpbox(0,0),subtermmbox(0,0);
    		 for(size_t ibox=0;ibox<tsi.denfac->size();ibox++){
    			 reddenfac=T(1)/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
    			 subtermpbox-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsp-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsp-(*tsi.triboxpoles)[2*ibox+1]))*reddenfac;
    			 subtermmbox-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsm-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsm-(*tsi.triboxpoles)[2*ibox+1]))*reddenfac;
    		 }

    		 // Now construct the actual subtraction term
    		 ypfac=eval_pts[icirc]/(T(2)*resb);
    		 temp[0]=-(T(1)/Nysolp*(subtermpbox+subtermptri)-T(1)/Nysolm*(subtermmbox+subtermmtri))*ypfac;
			 BH_DEBUG_PRINT(temp[0]);
    		 temp[1]=(T(3)/(T(2)/T(3)-Nysolp)*(subtermpbox+subtermptri)-T(3)/(T(2)/T(3)-Nysolm)*(subtermmbox+subtermmtri))*ypfac;
			 BH_DEBUG_PRINT(temp[1]);

    		 ampfull[0]+=temp[0];
			 BH_DEBUG_PRINT(ampfull[0]);
			 ampfull[1]+=temp[1];
			 BH_DEBUG_PRINT(ampfull[1]);
    		 tp.amp_err+=temp[0]*eval_pts[icirc]; //We use the y=0 component in the computation of the error from the 1/t coeff
    	 }

         //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
         // triangle triple cut using the analytic formula below.
         complex<T> subdenfac=T(1)/(tsi.gamma_old*chiK2Klflat);
         complex<T> subfacb=T(2)*(tsi.v1old.P()*chik1b)*subdenfac;
         complex<T> subfac=T(2)*(tsi.v2old.P()*chik1b)*subdenfac;
         complex<T> subfac2=-tp.S2+K1K2;
         complex<T> subdelta=(pow(K1K2,2)-tp.S1*tp.S2)/T(3);
         complex<T> pow3=(pow(subfac2,2)+subdelta);

         complex<T> subres=(tsi.orig_coeffs[0]*pow(subfac,3)*pow3+tsi.orig_coeffs[1]*pow(subfac,2)*subfac2+tsi.orig_coeffs[2]*subfac
                                        +tsi.orig_coeffs[4]*subfacb+tsi.orig_coeffs[5]*pow(subfacb,2)*subfac2+tsi.orig_coeffs[6]*pow(subfacb,3)*pow3);

    	 // Because errors show up when we multiply a large number against
    	 //   what should be zero but is 1e-16 at best due to the number of digits
    	 //   then the error on this result will be the difference of the exponent of the
    	 //   two quantities
    	 // To judge weather we are in this situation we apply this error if computation
    	 //   when the dynamic range of the computation is over 16, which should correspond
    	 //   to situations of this type.
         if(real(subfac*conj(subfac))<real(subfacb*conj(subfacb))){
        	 T comp(abs(pow3*pow(subfacb,3)/tsi.orig_coeffs[6]));
        	 if(comp>MaxDigitsE<T>()){
        		 tp.tri_sub_acc=convertToRVHP(T(1e2)/(abs(pow3*pow(subfacb,3)*tsi.orig_coeffs[6])));
        	 }
        	 else{
        		 tp.tri_sub_acc=convertToRVHP(MaxDigitsE<T>());
        	 }
//        	 _MESSAGE4("For :",subres," accuracy ",to_double(log(tp.tri_sub_acc))/log(10.0));
//        	 tp.tri_sub_acc=to_double(MaxDigits<T>()-log(abs(pow3*pow(subfacb,3)))/log(T(10.0)))-2.0;
         }
         else{
           	 T comp(abs(pow3*pow(subfac,3)/tsi.orig_coeffs[0]));
           	 if(comp>MaxDigitsE<T>()){
           		 tp.tri_sub_acc=convertToRVHP(T(1e2)/(abs(pow3*pow(subfac,3)*tsi.orig_coeffs[0])));
           	 }
           	 else{
           		 tp.tri_sub_acc=convertToRVHP(MaxDigitsE<T>());
           	 }
//           	 _MESSAGE4("For :",subres," accuracy ",to_double(log(tp.tri_sub_acc))/log(10.0));
//        	 tp.tri_sub_acc=to_double(MaxDigits<T>()-log(abs(pow3*pow(subfac,3)))/log(T(10.0)))-2.0;
         }

//    	 _MESSAGE7(*this," ycoeffs[0] : ",ampfull[0]," ycoeffs[2/3] : ",ampfull[1]," subres=",subres);
//    	 _MESSAGE2("Error : ",tp.amp_err);

         //Now return the complete subtraction result
         return complex<T>(0,1)*(ampfull[0]+ampfull[1])+T(16)*subres;
}


template <class cutDType,int CPOINTS,int TPOINTSBUB> template <class T> complex<T> Normal_Triangle_Subtraction<cutDType,CPOINTS,TPOINTSBUB>::get_sub_terms_work_pm(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi)
{
    complex<T> K1K2=tp.K1*tp.K2;

    momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
    momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

    //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
    // triangle triple cut using the analytic formula below.
    complex<T> facover=T(-2)/tsi.gamma_old;
    complex<T> ampp, ampm, ypfac;
    complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermtri;
    complex<T> chiK2Klflat=T(2)*(chik1b*tp.K2);

    // Reset the precision error gathering variable
    tp.amp_err=complex<T>(0,0);
    complex<T> temp[2];
    // Contains the value of the triangle pole subtraction for the bubble
    complex<T> ampfull[2]={complex<T>(0,0)};
    complex<T> subres(0,0);

	 // Now reconstruct the triple-cut subtraction terms
	 if(tp.masslessleg<0){
		 // This is a plus vertex triangle, so we only have 1...4 of triC0res as non-zero
		 //  we therefore need only the <K2flat||K1flat> pieces of the three mass case

		 // Construct the triangle pole subtraction
		 complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
		 for(int icirc=0; icirc<BUBPOINTS_STD; icirc++){
			 // Where vec2=<K2flat-||K1flat-> so that to get t we multiply t<K1flat-||K2flat->* vec2
			 bubtsp=T(tsi.reverse)*facover*(tsi.Nysolp[icirc]*fac1[0]
			              +(T(1)-tsi.Nysolp[icirc])*fac1[1]
		 		    	  +tsi.eval_pts[icirc]*fac1[2]
		 		    	  +tsi.Nysolp[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolp[icirc])*fac1[3]);
			 subtermtri=(pow(bubtsp,3)*tsi.orig_coeffs[0]+pow(bubtsp,2)*tsi.orig_coeffs[1]+bubtsp*tsi.orig_coeffs[2]+tsi.orig_coeffs[3])/complex<T>(0,1);

			 //Now add the box contributions to the triple cut
			 complex<T> subtermbox(0,0);
			 for(size_t ibox=0;ibox<tsi.denfac->size();ibox++){
				 subtermbox-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsp-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsp-(*tsi.triboxpoles)[2*ibox+1]))/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
			 }

			 // Now construct the actual subtraction term
			 ypfac=tsi.eval_pts[icirc]/(chiK2Klflat*(tsi.Nysolm[icirc]-tsi.Nysolp[icirc]));
			 temp[0]=T(1)/tsi.Nysolp[icirc]*(subtermtri+subtermbox)*ypfac;
			 temp[1]=-T(3)/(T(2)/T(3)-tsi.Nysolp[icirc])*(subtermtri+subtermbox)*ypfac;
			 ampfull[0]+=temp[0];
			 ampfull[1]+=temp[1];
			 tp.amp_err+=temp[0]*tsi.eval_pts[icirc]; //We use the y=0 component in the computation of the error from the 1/t coeff
		 }

		 //Now subtract off the triple-cut contributions
		 // Here subfac=vec1*<K1flatb||chi>
		 complex<T> subfac=T(2)*(tsi.v2old.P()*chik1b)/(tsi.gamma_old*chiK2Klflat);
		 complex<T> subfac2=-tp.S2+K1K2;
		 subres=(tsi.orig_coeffs[0]*pow(subfac,3)*(pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))+tsi.orig_coeffs[1]*pow(subfac,2)*subfac2+tsi.orig_coeffs[2]*subfac);

      	 T comp(abs((pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))*pow(subfac,3)/tsi.orig_coeffs[0]));
      	 if(comp>MaxDigitsE<T>()){
      		 tp.tri_sub_acc=convertToRVHP(T(1e2)/(abs((pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))*pow(subfac,3)*tsi.orig_coeffs[0])));
      	 }
      	 else{
      		 tp.tri_sub_acc=convertToRVHP(MaxDigitsE<T>());
      	 }
//    		 tp.tri_sub_acc=to_double(MaxDigits<T>()-log(abs(pow(subfac,3)*(pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))))/log(T(10)))-4.0;
	 }
	 else{
		 // This is a minus vertex triangle and so we get only 0...4 of triC0res filled
		 //  but the solution is the "inverse" of the plus case so we must multiply
		 //  using factors with K1flat<->K2flat which is the same as t<->a01*a02/t

		 // Construct the triangle pole subtraction
		 complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
		 for(int icirc=0; icirc<TPOINTSBUB; icirc++){
			 // Where vec2=<K1flat-||K2flat-> so that to get t we multiply t<K2flat-||K1flat->*vec2
			 bubtism=T(tsi.reverse)*facover*(tsi.Nysolm[icirc]*fac1[0]
			              +(T(1)-tsi.Nysolm[icirc])*fac1[1]
			   			  +tsi.eval_pts[icirc]*fac1[2]
	 	    			  +tsi.Nysolm[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolm[icirc])*fac1[3]);
			 BH_DEBUG_PRINT(bubtism);

			 subtermtri=complex<T>(0,-1)*(pow(bubtism,3)*tsi.orig_coeffs[0]+pow(bubtism,2)*tsi.orig_coeffs[1]+bubtism*tsi.orig_coeffs[2]+tsi.orig_coeffs[3]);
			 BH_DEBUG_PRINT(tsi.orig_coeffs[0]);
			 BH_DEBUG_PRINT(tsi.orig_coeffs[1]);
			 BH_DEBUG_PRINT(tsi.orig_coeffs[2]);
			 BH_DEBUG_PRINT(tsi.orig_coeffs[3]);
			 BH_DEBUG_PRINT(subtermtri);

			 // Subtract off the box terms from the bubble
			 complex<T> subtermbox(0,0);
			 for(size_t ibox=0;ibox<tsi.denfac->size();ibox++){
				 subtermbox-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtism-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtism-(*tsi.triboxpoles)[2*ibox+1]))/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
			 }
			 BH_DEBUG_PRINT(subtermbox);

			 // Now construct the actual subtraction term
			 ypfac=tsi.eval_pts[icirc]/(chiK2Klflat*(tsi.Nysolm[icirc]-tsi.Nysolp[icirc]));
			 temp[0]=-T(1)/tsi.Nysolm[icirc]*(subtermtri+subtermbox)*ypfac;
			 temp[1]=T(3)/(T(2)/T(3)-tsi.Nysolm[icirc])*(subtermtri+subtermbox)*ypfac;
			 BH_DEBUG_PRINT(temp[0]);
			 BH_DEBUG_PRINT(temp[1]);
			 ampfull[0]+=temp[0];
			 BH_DEBUG_PRINT(ampfull[0]);
			 ampfull[1]+=temp[1];
			 BH_DEBUG_PRINT(ampfull[1]);
			 tp.amp_err+=temp[0]*tsi.eval_pts[icirc]; //We use the y=0 component in the computation of the error from the 1/t coeff
		 }

		 //Now subtract off the triple-cut contributions
		 // Here subfacb=vec1*<K1flatb||chi>
		 complex<T> subfacb=T(2)*(tsi.v2old.P()*chik1b)/(tsi.gamma_old*chiK2Klflat);
		 complex<T> subfac2=-tp.S2+K1K2;
		 subres=(tsi.orig_coeffs[2]*subfacb+tsi.orig_coeffs[1]*pow(subfacb,2)*subfac2+tsi.orig_coeffs[0]*pow(subfacb,3)*(pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3)));

      	 T comp(abs((pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))*pow(subfacb,3)/tsi.orig_coeffs[0]));
      	 if(comp>MaxDigitsE<T>()){
      		 tp.tri_sub_acc=convertToRVHP(T(1e2)/(abs((pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))*pow(subfacb,3)*tsi.orig_coeffs[0])));
      	 }
      	 else{
      		 tp.tri_sub_acc=convertToRVHP(MaxDigitsE<T>());
      	 }
//    		 tp.tri_sub_acc=to_double(MaxDigits<T>()-log(abs(pow(subfacb,3)*(pow(subfac2,2)+(tp.S1*tp.S2-pow(K1K2,2))/T(-3))))/log(T(10)))-4.0;
	 }
//    	 _MESSAGE7(*this," ycoeffs[0] : ",ampfull[0]," ycoeffs[2/3] : ",ampfull[1]," subres=",subres);
//       	 _MESSAGE2("Error : ",tp.amp_err);

    //Now return the complete subtraction result
    return complex<T>(0,1)*(ampfull[0]+ampfull[1])+T(16)*subres;

}

template <class cutDbase, int CPOINTS,int TPOINTSBUB,int YPOINTS> template <class T> complex<T> General_Triangle_Subtraction<cutDbase,CPOINTS,TPOINTSBUB,YPOINTS>::get_sub_terms_work_3m(const cutDbase& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi)
{
		complex<T> K1K2=tp.K1*tp.K2;

          // We construct the "gamma-" contributing terms by using the symmetry of the momentum param t->a01 a02/t
         //  to relate the gamma+ solution 1/t^3, 1/t^2 and 1/t terms to t^3, t^2 and t.
         momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
         momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

         complex<T> chiK2Klflat=T(2)*(chik1b*tp.K2);//mc.spab(tp.chi, tp.K2t, tp.K1flatb);
         complex<T> VPK1flatK2=tp.K1flatbc.P()*tp.K2;
         complex<T> VPchiK2=tp.chic.P()*tp.K2;
         complex<T> Delta3=(K1K2*K1K2-tp.S1*tp.S2);

         // Get the points to evaluate the circle on
         const complex<T>* eval_pts;
         bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_eval_pts(eval_pts);

         complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
         complex<T> fac2[4]={tsi.v1old.P()*tp.K1flatbc.P(),tsi.v1old.P()*tp.chic.P(),tsi.v1old.P()*k1bchi,tsi.v1old.P()*chik1b};

    	 complex<T> part, temp[YPOINTS]={complex<T>(0,0)}, temp_err[YPOINTS]={complex<T>(0,0)};

    	 const complex<T>* eval_ypts;
    	 bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_yeval_pts(eval_ypts);

    	 complex<T> temp_subres[YPOINTS]={complex<T>(0,0)};
         complex<T> invcircposp, invcircposm;
         complex<T> Nysolp, Nysolm;
         complex<T> reddenfac, ressqrt, resa, resb;
         complex<T> facover=T(-2)/tsi.gamma_old;
         complex<T> ampp, ampm, ypfac;
         complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermptri, subtermmtri, bubts, bubtis, subterm, invcircpos;
    	 // Calculate the subtraction pieces
    	 for(int icirc=0; icirc<TPOINTSBUB; icirc++){
    		 // First compute the solutions for the poles in y
    		 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-tp.S2)*chiK2Klflat/Delta3+T(1)/T(4)*pow(chiK2Klflat,2)/Delta3;
    		 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat*T(1)/T(2);
    		 resb=sqrt(Delta3)*sqrt(ressqrt);

    		 Nysolp=(resa+resb)/chiK2Klflat;
    		 Nysolm=(resa-resb)/chiK2Klflat;

    		 // Now construct the t in the bubble parameterisation
    		 bubtsp=T(tsi.reverse)*facover*(Nysolp*fac1[0]
    				              +(T(1)-Nysolp)*fac1[1]
    							  +eval_pts[icirc]*fac1[2]
    							  +Nysolp/eval_pts[icirc]*(T(1)-Nysolp)*fac1[3]);

    		 bubtsm=T(tsi.reverse)*facover*(Nysolm*fac1[0]
    							  +(T(1)-Nysolm)*fac1[1]
    							  +eval_pts[icirc]*fac1[2]
    							  +Nysolm/eval_pts[icirc]*(T(1)-Nysolm)*fac1[3]);

//    		 /************************************************************************/
//       		 bubtisp=T(tsi.reverse)*facover*(Nysolp*fac2[0]
//        				              +(T(1)-Nysolp)*fac2[1]
//        							  +eval_pts[icirc]*fac2[2]
//        							  +Nysolp/eval_pts[icirc]*(T(1)-Nysolp)*fac2[3]);
//
//       		 bubtism=T(tsi.reverse)*facover*(Nysolm*fac2[0]
//        							  +(T(1)-Nysolm)*fac2[1]
//        							  +eval_pts[icirc]*fac2[2]
//        							  +Nysolm/eval_pts[icirc]*(T(1)-Nysolm)*fac2[3]);
//
//			 subtermmtri=tsi.orig_coeffs[(CPOINTS-1)/2];
//			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				 subtermmtri+=pow(bubtsm,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
//			 }
//			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				 subtermmtri+=pow(bubtism,icoeffs+1)*tsi.orig_coeffs[(CPOINTS-1)/2+icoeffs+1];
//			 }
//			 subtermmtri*=complex<T>(0,-1);
//
//			 subtermptri=tsi.orig_coeffs[(CPOINTS-1)/2];
//			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				 subtermptri+=pow(bubtsp,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
//			 }
//			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				 subtermptri+=pow(bubtisp,icoeffs+1)*tsi.orig_coeffs[(CPOINTS-1)/2+icoeffs+1];
//			 }
//			 subtermptri*=complex<T>(0,-1);
//    		 /************************************************************************/

    		 //Now add the box contributions to the triple cut subtraction term
			 reddenfac=T(1)/((*tsi.denfac)[0]*((*tsi.triboxpoles)[1]-(*tsi.triboxpoles)[0]));
			 complex<T> subtermp(0,0), subtermm(0,0);
    		 for(size_t ibox=0;ibox<tsi.denfac->size();ibox++){
    			 reddenfac=T(1)/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
    			 subtermp-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsp-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsp-(*tsi.triboxpoles)[2*ibox+1]))*reddenfac;
    			 subtermm-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsm-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsm-(*tsi.triboxpoles)[2*ibox+1]))*reddenfac;
    		 }

    		 ypfac=eval_pts[icirc]/(T(2)*resb);
    		 for(int iy=0;iy<YPOINTS;iy++){
    			 // Now construct the t in the bubble parameterisation
    			 invcircpos=eval_ypts[iy]/eval_pts[icirc]*(T(1)-eval_ypts[iy]);
    			 bubts=T(tsi.reverse)*facover*(eval_ypts[iy]*fac1[0]
    			              +(T(1)-eval_ypts[iy])*fac1[1]
    						  +eval_pts[icirc]*fac1[2]
    						  +invcircpos*fac1[3]);

    			 bubtis=T(tsi.reverse)*facover*(eval_ypts[iy]*fac2[0]
    						  +(T(1)-eval_ypts[iy])*fac2[1]
    						  +eval_pts[icirc]*fac2[2]
    						  +invcircpos*fac2[3]);

    			 // Reconstruct the subtracted numerator
    			 subterm=tsi.orig_coeffs[(CPOINTS-1)/2];
    			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
    				 subterm+=pow(bubts,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
    			 }
    			 for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
    				 subterm+=pow(bubtis,icoeffs+1)*tsi.orig_coeffs[(CPOINTS-1)/2+icoeffs+1];
    			 }
    			 subterm*=complex<T>(0,-1);

    			 // Now construct the actual subtraction term we add and subtract the triangle poles and so we end up just adding the box
    			 //  pole terms at the triangle triple cut to the points around the circle we are using to compute the contribution of the triangle
    			 //  to the bubble at infinity (we do not include the box terms at the box pole in y or those from the eval around the points for
    			 //  the extra contribution also all these box terms (including the ones we use) cancel each other (as there is no boundary))
    			 momentum<complex<T> > l2=eval_ypts[iy]*tp.K1flatbc.P()+(T(1)-eval_ypts[iy])*tp.chic.P()+eval_pts[icirc]*k1bchi+eval_ypts[iy]/eval_pts[icirc]*(T(1)-eval_ypts[iy])*chik1b-tp.K2;
    			 part=subterm/l2.square()+(subtermp/(eval_ypts[iy]-Nysolp)-subtermm/(eval_ypts[iy]-Nysolm))*ypfac;
    			 temp[iy]+=part;
//    			 temp_err[iy]+=part*eval_pts[icirc]; // We will only take the error from the y=0 component
//    			 temp[iy]+=((subtermp+subtermptri)/(eval_ypts[iy]-Nysolp)-(subtermm+subtermmtri)/(eval_ypts[iy]-Nysolm))*ypfac;
//    			 temp_subres[iy]+=subterm/l2.square()-((subtermptri)/(eval_ypts[iy]-Nysolp)-(subtermmtri)/(eval_ypts[iy]-Nysolm))*ypfac;
    			 temp_err[iy]+=((subtermp+subtermptri)/(eval_ypts[iy]-Nysolp)-(subtermm+subtermmtri)/(eval_ypts[iy]-Nysolm))*ypfac*eval_pts[icirc];
    		 }
    	 }

    	 // Extract the coefficients of y
    	 complex<T> ampfull(0,0);
    	 const complex<T>* ymatrixpoint;
    	 bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_y_matrix_eval_pts(ymatrixpoint);
    	 for(int i=0;i<YPOINTS;i++){
    		 for(int j=0;j<YPOINTS;j++){
    			 ampfull+=ymatrixpoint[j+i*(YPOINTS)]*temp[j];
    		 }
    	 }

    	 //For the error we only need one component, the y=0 coefficient
         tp.amp_err=complex<T>(0,0); // Reset the precision error gathering variable
    	 for(int inum=0;inum<YPOINTS;inum++){
    		 tp.amp_err+=ymatrixpoint[inum]*temp_err[inum];
    	 }

//    	 _MESSAGE6("3m : ",*this," ampfull : ",ampfull,", Error : ",tp.amp_err);

    	 //Now return the complete subtraction result
    	 return ampfull;


}

template <class cutDType,int CPOINTS,int TPOINTSBUB,int YPOINTS> template <class T> complex<T> General_Triangle_Subtraction<cutDType,CPOINTS,TPOINTSBUB,YPOINTS>::get_sub_terms_work_pm(const cutDType& cd,momentum_configuration<T>& mc, const vector<int>& ind, coeffparam<T,CPOINTS>& tp,tri_sub_info<T>& tsi)
{
	BH_DEBUG_MESSAGE("HERE");
	momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
	momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

	//The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
	// triangle triple cut using the analytic formula below.
	complex<T> facover=T(-2)/tsi.gamma_old;
	complex<T> ampp, ampm, ypfac;
	complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermptri, subtermmtri, subterm;
	complex<T> chiK2Klflat=T(2)*(chik1b*tp.K2);

	complex<T> temp[YPOINTS]={complex<T>(0,0)}, temp_err[YPOINTS]={complex<T>(0,0)};
	// Contains the value of the triangle pole subtraction for the bubble
	complex<T> part;

	complex<T> temp_subres[YPOINTS]={complex<T>(0,0)};

	const complex<T>* eval_ypts;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_yeval_pts(eval_ypts);

	// Now reconstruct the triple-cut subtraction terms
	if(tp.masslessleg<0){
		// This is a plus vertex triangle, so we only have 1...4 of triC0res as non-zero
		//  we therefore need only the <K2flat||K1flat> pieces of the three mass case

		// Construct the triangle pole subtraction
		complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
		for(int icirc=0; icirc<TPOINTSBUB; icirc++){
			// Where vec2=<K2flat-||K1flat-> so that to get t we multiply t<K1flat-||K2flat->* vec2
			bubtsp=T(tsi.reverse)*facover*(tsi.Nysolp[icirc]*fac1[0]
    			              +(T(1)-tsi.Nysolp[icirc])*fac1[1]
    		 		    	  +tsi.eval_pts[icirc]*fac1[2]
    		 		    	  +tsi.Nysolp[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolp[icirc])*fac1[3]);

//			/*****************************************************/
//			subtermptri=tsi.orig_coeffs[(CPOINTS-1)/2];
//			for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				subtermptri+=pow(bubtsp,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
//			}
//			subtermptri*=complex<T>(0,-1);
//			/*****************************************************/
			//We dont include the minus solution normally as this corresponds to a vanishing subtraction term
			// but in computing around t=0 we include it so we need to subtract it out
			bubtsm=T(tsi.reverse)*facover*(tsi.Nysolm[icirc]*fac1[0]
    			              +(T(1)-tsi.Nysolm[icirc])*fac1[1]
    		 		    	  +tsi.eval_pts[icirc]*fac1[2]
    		 		    	  +tsi.Nysolm[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolm[icirc])*fac1[3]);

			subtermmtri=tsi.orig_coeffs[(CPOINTS-1)/2];
			for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
				subtermmtri+=pow(bubtsm,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
			}
			subtermmtri*=complex<T>(0,-1);

			//Now add the box contributions to the triple cut
			complex<T> subtermp(0,0);
			for(size_t ibox=1;ibox<tsi.denfac->size();ibox++){
				subtermp-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtsp-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtsp-(*tsi.triboxpoles)[2*ibox+1]))/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
			}

			// Now construct the actual subtraction term
			ypfac=tsi.eval_pts[icirc]/(chiK2Klflat*(tsi.Nysolm[icirc]-tsi.Nysolp[icirc]));
			for(int iy=0;iy<YPOINTS;iy++){
				// Where vec2=<K2flat-||K1flat-> so that to get t we multiply t<K1flat-||K2flat->* vec2
				bubtsp=T(tsi.reverse)*facover*(eval_ypts[iy]*fac1[0]
    			              +(T(1)-eval_ypts[iy])*fac1[1]
    		 		    	  +tsi.eval_pts[icirc]*fac1[2]
    		 		    	  +eval_ypts[iy]/tsi.eval_pts[icirc]*(T(1)-eval_ypts[iy])*fac1[3]);

				subterm=tsi.orig_coeffs[(CPOINTS-1)/2];
				for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
					subterm+=pow(bubtsp,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
				}
				subterm*=complex<T>(0,-1);

				// Now construct the actual subtraction term
				momentum<complex<T> > l2=eval_ypts[iy]*tp.K1flatbc.P()+(T(1)-eval_ypts[iy])*tp.chic.P()+tsi.eval_pts[icirc]*k1bchi+eval_ypts[iy]/tsi.eval_pts[icirc]*(T(1)-eval_ypts[iy])*chik1b-tp.K2;
				part=subterm/l2.square()-(subtermp/(eval_ypts[iy]-tsi.Nysolp[icirc])+subtermmtri/(eval_ypts[iy]-tsi.Nysolm[icirc]))*ypfac;
    			temp[iy]+=part;
//    			temp_err[iy]+=part*eval_pts[icirc];
//				temp[iy]+=-ypfac*(subtermptri+subtermp)/(eval_ypts[iy]-Nysolp[icirc]);
//				temp_subres[iy]+=subterm/l2.square()+((subtermptri)/(eval_ypts[iy]-Nysolp[icirc])-subtermmtri/(eval_ypts[iy]-Nysolm[icirc]))*ypfac;
				temp_err[iy]-=tsi.eval_pts[icirc]*ypfac*(subtermptri+subtermp)/(eval_ypts[iy]-tsi.Nysolp[icirc]);
			}
		}
	}
	else{
		// This is a minus vertex triangle and so we get only 0...4 of triC0res filled
		//  but the solution is the "inverse" of the plus case so we must multiply
		//  using factors with K1flat<->K2flat which is the same as t<->a01*a02/t

		// Construct the triangle pole subtraction
		complex<T> fac1[4]={tsi.v2old.P()*tp.K1flatbc.P(),tsi.v2old.P()*tp.chic.P(),tsi.v2old.P()*k1bchi,tsi.v2old.P()*chik1b};
		for(int icirc=0; icirc<TPOINTSBUB; icirc++){
			// Where vec2=<K1flat-||K2flat-> so that to get t we multiply t<K2flat-||K1flat->*vec2
			bubtism=T(tsi.reverse)*facover*(tsi.Nysolm[icirc]*fac1[0]
    			              +(T(1)-tsi.Nysolm[icirc])*fac1[1]
    			   			  +tsi.eval_pts[icirc]*fac1[2]
    	 	    			  +tsi.Nysolm[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolm[icirc])*fac1[3]);

//			/*****************************************************/
//			subtermmtri=tsi.orig_coeffs[(CPOINTS-1)/2];
//			for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
//				subtermmtri+=pow(bubtism,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
//			}
//			subtermmtri*=complex<T>(0,-1);
//			/*****************************************************/
			//We dont include the minus solution normally as this corresponds to a vanishing subtraction term
			// but in computing around t=0 we include it so we need to subtract it out
			bubtisp=T(tsi.reverse)*facover*(tsi.Nysolp[icirc]*fac1[0]
    			              +(T(1)-tsi.Nysolp[icirc])*fac1[1]
    			   			  +tsi.eval_pts[icirc]*fac1[2]
    	 	    			  +tsi.Nysolp[icirc]/tsi.eval_pts[icirc]*(T(1)-tsi.Nysolp[icirc])*fac1[3]);

			subtermptri=tsi.orig_coeffs[(CPOINTS-1)/2];
			for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
				subtermptri+=pow(bubtisp,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
			}
			subtermptri*=complex<T>(0,-1);

			// Subtract off the box terms from the bubble
			complex<T> subtermm(0,0);
			for(size_t ibox=0;ibox<tsi.denfac->size();ibox++){
				subtermm-=((*tsi.triboxpoles)[2*ibox]*(*tsi.triboxcoeff)[2*ibox]/(bubtism-(*tsi.triboxpoles)[2*ibox])-(*tsi.triboxpoles)[2*ibox+1]*(*tsi.triboxcoeff)[2*ibox+1]/(bubtism-(*tsi.triboxpoles)[2*ibox+1]))/((*tsi.denfac)[ibox]*((*tsi.triboxpoles)[2*ibox+1]-(*tsi.triboxpoles)[2*ibox]));
			}

			// Now construct the actual subtraction term
			ypfac=tsi.eval_pts[icirc]/(chiK2Klflat*(tsi.Nysolm[icirc]-tsi.Nysolp[icirc]));
			for(int iy=0;iy<YPOINTS;iy++){
				bubtism=T(tsi.reverse)*facover*(eval_ypts[iy]*fac1[0]
    				              +(T(1)-eval_ypts[iy])*fac1[1]
    				   			  +tsi.eval_pts[icirc]*fac1[2]
    		 	    			  +eval_ypts[iy]/tsi.eval_pts[icirc]*(T(1)-eval_ypts[iy])*fac1[3]);
				subterm=tsi.orig_coeffs[(CPOINTS-1)/2];
				for(int icoeffs=0;icoeffs<(CPOINTS-1)/2;icoeffs++){
					subterm+=pow(bubtism,(CPOINTS-1)/2-icoeffs)*tsi.orig_coeffs[icoeffs];
				}
				subterm*=complex<T>(0,-1);

				// Now construct the actual subtraction term
				momentum<complex<T> > l2=eval_ypts[iy]*tp.K1flatbc.P()+(T(1)-eval_ypts[iy])*tp.chic.P()+tsi.eval_pts[icirc]*k1bchi+eval_ypts[iy]/tsi.eval_pts[icirc]*(T(1)-eval_ypts[iy])*chik1b-tp.K2;
    			part=subterm/l2.square()+(subtermptri/(eval_ypts[iy]-tsi.Nysolp[icirc])+subtermm/(eval_ypts[iy]-tsi.Nysolm[icirc]))*ypfac;
    			temp[iy]+=part;
//    			temp_err[iy]+=part*eval_pts[icirc];
//				temp[iy]+=ypfac*(subtermmtri+subtermm)/(eval_ypts[iy]-Nysolm[icirc]);
//				temp_subres[iy]+=subterm/l2.square()+(subtermptri/(eval_ypts[iy]-Nysolp[icirc])-subtermmtri/(eval_ypts[iy]-Nysolm[icirc]))*ypfac;
				temp_err[iy]+=tsi.eval_pts[icirc]*ypfac*(subtermmtri+subtermm)/(eval_ypts[iy]-tsi.Nysolm[icirc]);
			}
		}
	}

	// Extract the coefficients of y
	complex<T> ampfull(0,0);
	const complex<T>* ymatrixpoint;
	bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_y_matrix_eval_pts(ymatrixpoint);
	for(int i=0;i<YPOINTS;i++){
		for(int j=0;j<YPOINTS;j++){
			ampfull+=ymatrixpoint[j+i*(YPOINTS)]*temp[j];
			BH_DEBUG_PRINT(ampfull);
		}
	}

	//For the error we only need one component
	tp.amp_err=complex<T>(0,0); // Reset the precision error gathering variable
	for(int inum=0;inum<YPOINTS;inum++){
		tp.amp_err+=ymatrixpoint[inum]*temp_err[inum];
	}

//	_MESSAGE7(tp.masslessleg," : ",*this," ampfull : ",ampfull,", Error : ",tp.amp_err);

	//Now return the complete subtraction result
	return ampfull;
}



template <class cutDbase,class TriangleSpecs> void triangle_Darren<cutDbase,TriangleSpecs>::change_klegs(int k1,int k2,int k3){
	_k1leg=k1;
	_k2leg=k2;
	_k3leg=k3;

	// Find out which way around the loop the momentum circulates
	int facdec=_k1leg-_k2leg;
	if(facdec==-1||facdec==2){// if we have 321 ordering for k1, k2, k3 then we have one ordering for l, l1, l2
		_reverse=-1;
	}
	else{
		_reverse=1;
	}
}

}
}
}



#endif /* CUT_DARREN_HPP_ */
