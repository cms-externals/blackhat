/*
 * rational_eval.cpp
 *
 *  Created on: Apr 7, 2009
 *      Author: dforde
 */

#include "BH_utilities.h"
#include "rational.h"
#include "BH_A0.h"
#include "cut_part_factory.h"
#ifndef BH_PUBLIC
#include "partial_order.h"
#include "ordering.h"
#endif




#define  _VERBOSE 0

#if _VERBOSE
//Definitions needed for debugging purposes
#include "known_rational.h"
#include "rec_rational.h"
#include "inf_rat.h"
#endif


using namespace std;

namespace BH {

/*
 *
 *
 * Code for Rec_Pair with eval_param
 *
 *
 *
 */

template <class T> complex<T>  Rec_Pair::Rec_Pair_eval(const eval_param<T>& ep){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_eval ",get_part().get_code()," -----------------------------------");
#endif
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(get_l_eval_param<T>());
	eval_param<T>& epr(get_r_eval_param<T>());
	
	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom(ep.p(c1[0].ind())->P());
	epl.set(0,ep.p(c1[0].ind()));
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		size_t ind=c1[psiter].ind();
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}
	
	// Create the shifted momenta
	momentum<complex<T> > ij_mom;
	Spinor_to_momentum(ep.p(j)->Lt(),ep.p(i)->L(),ij_mom);
	const complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));
	
	// Now compute the on-shell propagator momentum
	ij_mom.mult_by(z);
	ij_mom.add_to(Psum_mom);
	get_phat<T>().set_to_U(ij_mom);
	get_mphat<T>().set_to(get_phat<T>().L(),-get_phat<T>().Lt());

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		epr.set(kr,ep.p(c2[kr].ind()));
	}

	// Insert the shifted legs
	get_i_leg<T>().set_to(ep.p(i)->L(),ep.p(i)->Lt()-z*ep.p(j)->Lt());
	get_j_leg<T>().set_to(ep.p(j)->L()+z*ep.p(i)->L(),ep.p(j)->Lt());
	epr.set(shifted_ind_i,&get_i_leg<T>());  // The shifted legs are massless
	epl.set(shifted_ind_j,&get_j_leg<T>());

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempresl=Amp_safe(left()->eval(epl));
	complex<T> tempresr=Amp_safe(right()->eval(epr));
	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(fact()->eval(epr))/Psqr);
	_MESSAGE7(tempres*tempresl*tempresr," : From ",tempresl," * (",T(1.)/Psqr,") * ",tempresr);
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	
	return tempres*tempresl*tempresr;
#else
	complex<T> A1=left()->eval(epl);
	complex<T> A2=right()->eval(epr);
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2*(fact()->eval(epr)))/Psqr);
#endif

}

template <class T> complex<T>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(const eval_param<T>& ep){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_massive_prop_massless_shift ",get_part().get_code()," -----------------------------------");
#endif
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(get_l_eval_param<T>());
	eval_param<T>& epr(get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom(ep.p(c1[0].ind())->P());
	epl.set(0,ep.p(c1[0].ind()));
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		size_t ind=c1[psiter].ind();
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}

	// Create the shifted momenta
	momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),ep.p(i)->L());
	const complex<T> Psqr=Psum_mom.square()-eval_param<T>::mass2(_mass_leg);
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));
    
	// Now compute the on-shell propagator momentum
	Cmom<T> PHat(Psum_mom+z*ij_mom,_mt_massive);
	Cmom<T> mPHat(-PHat.P(),_mt_massive);
    
	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&PHat);
	epl.set(maxl-1,&mPHat);
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		epr.set(kr,ep.p(c2[kr].ind()));
	}

	// Insert the shifted legs
	Cmom<T> shifti(ep.p(i)->L(),ep.p(i)->Lt()-z*ep.p(j)->Lt());
	Cmom<T> shiftj(ep.p(j)->L()+z*ep.p(i)->L(),ep.p(j)->Lt());
	epr.set(shifted_ind_i,&shifti); // The shifted legs are massless
	epl.set(shifted_ind_j,&shiftj);

	//Pass down reference vector and masses
	epr.set_ref(ep.ref());
	epl.set_ref(ep.ref());

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempresl=Amp_safe(left()->eval(epl));
	complex<T> tempresr=Amp_safe(right()->eval(epr));
	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(fact()->eval(epr))/Psqr);
	_MESSAGE7(tempres*tempresl*tempresr," : From ",tempresl," * (",T(1.)/Psqr,") * ",tempresr);
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");

	return tempres*tempresl*tempresr;
#else
	complex<T> A1=left()->eval(epl);
	complex<T> A2=right()->eval(epr);
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2*(fact()->eval(epr)))/Psqr);
#endif
}

template <class T> complex<T>  Rec_Pair_massive_unshifted::Rec_Pair_eval(const eval_param<T>& ep){
#if _VERBOSE
	_MESSAGE3("-------------------------Begin Rec_Pair_massive_unshifted ",get_part().get_code()," -----------------------------------");
#endif
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(get_l_eval_param<T>());
	eval_param<T>& epr(get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom(ep.p(c1[0].ind())->P());
	epl.set(0,ep.p(c1[0].ind()));
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		size_t ind(c1[psiter].ind());
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}

	// Create the shifted momenta
	momentum<complex<T> > ij_mom=PfLLt(ep.p(j)->Lt(),ep.p(i)->L());
	const complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));

	// Now compute the on-shell propagator momentum
	Cmom<T> PHat(Psum_mom+z*ij_mom);
	Cmom<T> mPHat(PHat.L(),-PHat.Lt());
    
	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&PHat); // This is a massless propagator
	epl.set(maxl-1,&mPHat);
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		epr.set(kr,ep.p(c2[kr].ind()));
	}

	// Insert the shifted legs
	Cmom<T> shifti(ep.p(i)->L(),ep.p(i)->Lt()-z*ep.p(j)->Lt());
	Cmom<T> shiftj(ep.p(j)->L()+z*ep.p(i)->L(),ep.p(j)->Lt());
	epr.set(shifted_ind_i,&shifti); // The shifted legs are massless
	epl.set(shifted_ind_j,&shiftj);

	//Pass down reference vector and masses
	epr.set_ref(ep.ref());
	epl.set_ref(ep.ref());

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};
	
	complex<T> tempresl=Amp_safe(left()->eval(epl));
	complex<T> tempresr=Amp_safe(right()->eval(epr));
	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(fact()->eval(epr))/Psqr);
	_MESSAGE7(tempres*tempresl*tempresr," : From ",tempresl," * (",T(1.)/Psqr,") * ",tempresr);
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	
	return tempres*tempresl*tempresr;
#else
	complex<T> A1=left()->eval(epl);
	complex<T> A2=right()->eval(epr);
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2*(fact()->eval(epr)))/Psqr);
#endif
}

template <class T> complex<T> Rec_Pair_massive::Rec_Pair_eval(const eval_param<T>& ep){
#if _VERBOSE
	_MESSAGE3("--------------------------Begin Rec_Pair_massive ",get_part().get_code()," ----------------------------------");
#endif
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(get_l_eval_param<T>());
	eval_param<T>& epr(get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom(ep.p(c1[0].ind())->P());
	epl.set(0,ep.p(c1[0].ind()));
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		size_t ind(c1[psiter].ind());
		Psum_mom.add_to(ep.p(ind)->P());
        epl.set(psiter,ep.p(ind));
	}

	// Create the shifted momenta
	const complex<T> Psqr=Psum_mom.square();
	Cmom<T> shifti, shiftj, Phat;
    const Cmom<T> *ref_i, *ref_j;
	get_shifted_ij(ep,shifti,shiftj,Phat,Psum_mom,Psqr,ref_i,ref_j);
	Cmom<T> mPHat(Phat.L(),-Phat.Lt());

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&Phat); // This is a massless propagator
	epl.set(maxl-1,&mPHat);
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		epr.set(kr,ep.p(c2[kr].ind()));
	}

	// Insert the shifted legs
	epr.set(shifted_ind_i,&shifti);
	epl.set(shifted_ind_j,&shiftj);

	//Pass down reference vector and masses
//  Chenged to the default values to make the recursion work, not sure what the problem is here.
//	epr.set_ref(ref_i);
//	epl.set_ref(ref_j);
	epr.set_ref(ep.ref());
	epl.set_ref(ep.ref());

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempresl=Amp_safe(left()->eval(epl));
	complex<T> tempresr=Amp_safe(right()->eval(epr));
	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(fact()->eval(epr))/Psqr);
	_MESSAGE7(tempres*tempresl*tempresr," : From ",tempresl," * (",T(1.)/Psqr,") * ",tempresr);
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	
	return tempres*tempresl*tempresr;
#else
	complex<T> A1=left()->eval(epl);
	complex<T> A2=right()->eval(epr);
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2*(fact()->eval(epr)))/Psqr);
#endif
}

template <class T> complex<T> Rec_Pair_massive_prop::Rec_Pair_eval(const eval_param<T>& ep){
#if _VERBOSE
	_MESSAGE3("--------------------------Begin Rec_Pair_massive_prop ",get_part().get_code()," ----------------------------------");
#endif
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(get_l_eval_param<T>());
	eval_param<T>& epr(get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	const std::vector<plabel>& c1=get_part().c(1);
	momentum<complex<T> > Psum_mom(ep.p(c1[0].ind())->P());
	epl.set(0,ep.p(c1[0].ind()));
	for(size_t psiter=1;psiter<maxl-1;psiter++){
		size_t ind(c1[psiter].ind());
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}

	// Create the shifted momenta
	const complex<T> Psqr=Psum_mom.square()-eval_param<T>::mass2(_mass_leg);
	Cmom<T> shifti, shiftj, Phat;
    const Cmom<T> *ref_i, *ref_j;
	get_shifted_ij(ep,shifti,shiftj,Phat,Psum_mom,Psqr,ref_i,ref_j);
	Cmom<T> mPHat=-Phat;

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&Phat); // This is a massless propagator
	epl.set(maxl-1,&mPHat);
	const std::vector<plabel>& c2=get_part().c(2);
	for(size_t kr=1;kr<maxr;kr++){
		epr.set(kr,ep.p(c2[kr].ind()));
	}

	// Insert the shifted legs
	epr.set(shifted_ind_i,&shifti);
	epl.set(shifted_ind_j,&shiftj);

	//  Chenged to the default values to make the recursion work, not sure what the problem is here.
	//	epr.set_ref(ref_i);
	//	epl.set_ref(ref_j);
	epr.set_ref(ep.ref());
	epl.set_ref(ep.ref());

	//Evaluate the recursive terms and check if the returned value is nan, if so then
	// the contribution should be zero as we this comes from a 0/0 situation
#if _VERBOSE
	char llabel[]="T", rlabel[]="T";
	if(typeid(*(left()))==typeid(Known_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Unknown_Rec_Rational)){llabel[0]='R';};
	if(typeid(*(left()))==typeid(Known_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(left()))==typeid(Unknown_Rec_InfRational)){llabel[0]='I';};
	if(typeid(*(right()))==typeid(Known_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Unknown_Rec_Rational)){rlabel[0]='R';};
	if(typeid(*(right()))==typeid(Known_Rec_InfRational)){rlabel[0]='I';};
	if(typeid(*(right()))==typeid(Unknown_Rec_InfRational)){rlabel[0]='I';};

	complex<T> tempres=Amp_safe(complex<T>(0.0,-1.0)*(left()->eval(epl)*right()->eval(epr)*(fact()->eval(epr)))/Psqr);
	_MESSAGE7(tempres," : From ",left()->eval(epl)," * (",T(1.)/Psqr,") * ",right()->eval(epr));
	_MESSAGE7(get_part().get_code()," : ",llabel,get_part().c(1),"*",rlabel,get_part().c(2));
	_MESSAGE3("--------------------------End Eval ",get_part().get_code()," ----------------------------------");
	return tempres;
#else
	complex<T> A1=left()->eval(epl);
	complex<T> A2=right()->eval(epr);
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2*(fact()->eval(epr)))/Psqr);
#endif
}

//Explicit Instantiation
template complex<R>  Rec_Pair_massive_unshifted::Rec_Pair_eval(const eval_param<R>& ep);
template complex<RHP>  Rec_Pair_massive_unshifted::Rec_Pair_eval(const eval_param<RHP>& ep);
template complex<RVHP>  Rec_Pair_massive_unshifted::Rec_Pair_eval(const eval_param<RVHP>& ep);

template complex<R>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(const eval_param<R>& ep);
template complex<RHP>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(const eval_param<RHP>& ep);
template complex<RVHP>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(const eval_param<RVHP>& ep);

template complex<R>  Rec_Pair::Rec_Pair_eval(const eval_param<R>& ep);
template complex<RHP>  Rec_Pair::Rec_Pair_eval(const eval_param<RHP>& ep);
template complex<RVHP>  Rec_Pair::Rec_Pair_eval(const eval_param<RVHP>& ep);

template complex<R>  Rec_Pair_massive::Rec_Pair_eval(const eval_param<R>& ep);
template complex<RHP>  Rec_Pair_massive::Rec_Pair_eval(const eval_param<RHP>& ep);
template complex<RVHP>  Rec_Pair_massive::Rec_Pair_eval(const eval_param<RVHP>& ep);

template complex<R>  Rec_Pair_massive_prop::Rec_Pair_eval(const eval_param<R>& ep);
template complex<RHP>  Rec_Pair_massive_prop::Rec_Pair_eval(const eval_param<RHP>& ep);
template complex<RVHP>  Rec_Pair_massive_prop::Rec_Pair_eval(const eval_param<RVHP>& ep);


#if BH_USE_GMP

template complex<RGMP>  Rec_Pair_massive_unshifted::Rec_Pair_eval(const eval_param<RGMP>& ep);
template complex<RGMP>  Rec_Pair_massive_prop_massless_shift::Rec_Pair_eval(const eval_param<RGMP>& ep);
template complex<RGMP>  Rec_Pair::Rec_Pair_eval(const eval_param<RGMP>& ep);
template complex<RGMP>  Rec_Pair_massive::Rec_Pair_eval(const eval_param<RGMP>& ep);
template complex<RGMP>  Rec_Pair_massive_prop::Rec_Pair_eval(const eval_param<RGMP>& ep);

#endif

}



