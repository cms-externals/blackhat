/*
 * pentagon_ratext.hpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef PENTAGON_RATEXT_HPP_
#define PENTAGON_RATEXT_HPP_

#include "ratext/pentagon_ratext.h"
#include "ratext/rat_ext.h"
#include "ratext/box_ratext.h"
#include <cassert>

#define _USE_ORIGINAL 1

using std::complex;

namespace BH {

namespace ratext {

template <class BASE, class RatPentSpecs> void  pentagon_Rat<BASE,RatPentSpecs>::init()
{
	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	indlst[0].assign(BASE::corner_size(1)+4,0);
	indlst[1].assign(BASE::corner_size(2)+4,0);
	indlst[2].assign(BASE::corner_size(3)+4,0);
	indlst[3].assign(BASE::corner_size(4)+4,0);
	indlst[4].assign(BASE::corner_size(5)+4,0);

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

	_ep[4]=new eval_param<R>(BASE::corner_size(5)+2);
	_ep_HP[4]=new eval_param<RHP>(BASE::corner_size(5)+2);
	_ep_VHP[4]=new eval_param<RVHP>(BASE::corner_size(5)+2);

#if BH_USE_GMP
	_ep_GMP[0]=new eval_param<RGMP>(BASE::corner_size(1)+2);
	_ep_GMP[1]=new eval_param<RGMP>(BASE::corner_size(2)+2);
	_ep_GMP[2]=new eval_param<RGMP>(BASE::corner_size(3)+2);
	_ep_GMP[3]=new eval_param<RGMP>(BASE::corner_size(4)+2);
	_ep_GMP[4]=new eval_param<RGMP>(BASE::corner_size(5)+2);
#endif
}

template <class BASE, class RatPentSpecs> void pentagon_Rat<BASE,RatPentSpecs>::add_mass(const cutD& pcd)
{
	add_mass(pcd.l(1).mass_label(),
			pcd.l(2).mass_label(),
			pcd.l(3).mass_label(),
			pcd.l(4).mass_label(),
			pcd.l(5).mass_label());

}

template <class BASE, class RatPentSpecs> void pentagon_Rat<BASE,RatPentSpecs>::add_mass(int m1, int m2, int m3 , int m4,int m5)
{
	if(find(_cut_mass.begin(),_cut_mass.end(),m1)==_cut_mass.end()){_cut_mass.push_back(m1);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m2)==_cut_mass.end()){_cut_mass.push_back(m2);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m3)==_cut_mass.end()){_cut_mass.push_back(m3);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m4)==_cut_mass.end()){_cut_mass.push_back(m4);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m5)==_cut_mass.end()){_cut_mass.push_back(m5);};

	//Set up the masses for the legs
	_leg_masses->set(0,m1);
	_leg_masses->set(1,m2);
	_leg_masses->set(2,m3);
	_leg_masses->set(3,m4);
	_leg_masses->set(4,m5);
}


template <class BASE, class RatPentSpecs> pentagon_Rat<BASE,RatPentSpecs>::~pentagon_Rat()
{
	delete _ep[0];
	delete _ep_HP[0];
	delete _ep_VHP[0];
	delete _ep[1];
	delete _ep_HP[1];
	delete _ep_VHP[1];
	delete _ep[2];
	delete _ep_HP[2];
	delete _ep_VHP[2];
	delete _ep[3];
	delete _ep_HP[3];
	delete _ep_VHP[3];
	delete _ep[4];
	delete _ep_HP[4];
	delete _ep_VHP[4];
#if BH_USE_GMP
	delete _ep_GMP[0];
	delete _ep_GMP[1];
	delete _ep_GMP[2];
	delete _ep_GMP[3];
	delete _ep_GMP[4];
#endif

	delete _leg_masses;
}


template <class BASE, class RatPentSpecs> C pentagon_Rat<BASE,RatPentSpecs>::eval(const eval_param<R>& ep)
{
	return this->template get_symmetry_factor<R>()*E0coeff;
}

template <class BASE, class RatPentSpecs> CHP pentagon_Rat<BASE,RatPentSpecs>::eval(const eval_param<RHP>& ep)
{
	return this->template get_symmetry_factor<RHP>()*E0coeff_HP;
}
template <class BASE, class RatPentSpecs> CVHP pentagon_Rat<BASE,RatPentSpecs>::eval(const eval_param<RVHP>& ep)
{
	return this->template get_symmetry_factor<RVHP>()*E0coeff_VHP;
}

template <class BASE, class RatPentSpecs> C pentagon_Rat<BASE,RatPentSpecs>::eval(mom_conf& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the pentagon compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("pentagon_Rat<BASE> : Attempting to compute a pentagon before computing the bubble");
	}
	return this->template get_symmetry_factor<R>()*E0coeff;
}

template <class BASE, class RatPentSpecs> CHP pentagon_Rat<BASE,RatPentSpecs>::eval(mom_conf_HP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the pentagon compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("pentagon_Rat<BASE> : Attempting to compute a pentagon before computing the bubble");
	}
	return this->template get_symmetry_factor<RHP>()*E0coeff_HP;
}

template <class BASE, class RatPentSpecs> CVHP pentagon_Rat<BASE,RatPentSpecs>::eval(mom_conf_VHP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the pentagon compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("pentagon_Rat<BASE> : Attempting to compute a pentagon before computing the bubble");
	}
	return this->template get_symmetry_factor<RVHP>()*E0coeff_VHP;
}

#if BH_USE_GMP
template <class BASE, class RatPentSpecs> CGMP pentagon_Rat<BASE,RatPentSpecs>::eval(const eval_param<RGMP>& ep)
{
	return this->template get_symmetry_factor<RGMP>()*E0coeff_GMP;
}


template <class BASE, class RatPentSpecs> CGMP pentagon_Rat<BASE,RatPentSpecs>::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the pentagon compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("pentagon_Rat<BASE> : Attempting to compute a pentagon before computing the bubble");
	}
	return this->template get_symmetry_factor<RGMP>()*E0coeff_GMP;
}

#endif

template <class BASE, class RatPentSpecs> template <class T> void pentagon_Rat<BASE,RatPentSpecs>::get_sub_terms(momentum_configuration<T>& mc, const std::vector<int>& ind, pent_param<T,RatPentSpecs::MUPOINTS>& pp)
{
	momentum<complex<T> > K5sum_mom(mc.mom(ind[BASE::corner_ind(pp.pent_corner,1)-1]));
	for(size_t kiter=2;kiter<=BASE::corner_size(pp.pent_corner);kiter++){
		K5sum_mom+=mc.mom(ind[BASE::corner_ind(pp.pent_corner,kiter)-1]);
	}
	complex<T> S5=K5sum_mom.square();

	complex<T> npentpart;
	switch(pp.pent_corner){
	case 1://l_5=l-K_5
	case 5:
		npentpart=S5-T(2)*((pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P())*K5sum_mom);
		break;
	case 2://l_5=l_1-K_5=l-K_1-K_5
		npentpart=S5+T(2)*((pp.K1-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;
	case 3://l_5=l_4-K_5=l_1-K_4-K_5=l-K_1-K_4-K_5
		npentpart=S5+T(2)*((pp.K1+pp.K4-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;//l_4=l_1-K4
	case 4://l_5=l_2-K_5=l-K_2-K_5
		npentpart=S5+T(2)*((pp.K2-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;
	}

	// Construct the subtraction terms

	// Now compute the pentagon coeffs if needed
	if(!((pp.mcID==mcID)&&(ind==indID))){
		get_coeffs(mc,ind,pp.mcID,pp.vertex_ref);
	}

	//Set up the m0 masses we evaluate around
	const complex<T>* m0mass;
	(static_cast<typename RatPentSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatPentSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix);

	complex<T> m0sol0, sum[RatPentSpecs::MUPOINTS-2]={complex<T>(0,0)};
	complex<T> vec1cK5(pp.vec1c*K5sum_mom), vec2cK5(pp.vec2c*K5sum_mom);
	for(int m0pole=0;m0pole<RatPentSpecs::MUPOINTS;m0pole++){
		// the mu^2 subtraction terms, we evaluate at the value of t for the box
		m0sol0=pp.gamma*(T(2)*pow(pp.boxpoles[2*m0pole],2)*vec1cK5-npentpart*pp.boxpoles[2*m0pole]+T(2)*pp.alp1*pp.alp2*vec2cK5)/(T(2)*vec2cK5);
		for(int imupoints=0;imupoints<RatPentSpecs::MUPOINTS-2;imupoints++){
			sum[imupoints]+=m0mass_matrix[m0pole+(imupoints+2)*RatPentSpecs::MUPOINTS]*pp.boxpoles[2*m0pole]/(m0mass[m0pole]-m0sol0);
		}

		m0sol0=pp.gamma*(T(2)*pow(pp.boxpoles[1+2*m0pole],2)*vec1cK5-npentpart*pp.boxpoles[1+2*m0pole]+T(2)*pp.alp1*pp.alp2*vec2cK5)/(T(2)*vec2cK5);
		for(int imupoints=0;imupoints<RatPentSpecs::MUPOINTS-2;imupoints++){
			sum[imupoints]+=m0mass_matrix[m0pole+(imupoints+2)*RatPentSpecs::MUPOINTS]*pp.boxpoles[1+2*m0pole]/(m0mass[m0pole]-m0sol0);
		}
	}
	complex<T> pentcoeff_mass;
	get_E0coeff(pentcoeff_mass);
	for(int imupoints=0;imupoints<RatPentSpecs::MUPOINTS-2;imupoints++){
		pp.pent_sub_ret[imupoints]-=sum[imupoints]*pp.gamma*pentcoeff_mass/(complex<T>(0,2)*vec2cK5);
	}

	// Estimate the accuracy of the computation and if it is worse than the current value store it
	double pent_acc=to_double(MaxDigits<T>()-(log(abs(pentcoeff_mass)))/log(T(10)));
	if(pp.accuracy>pent_acc){
		pp.accuracy=pent_acc;
	}
}

template <class BASE, class RatPentSpecs> template <class T> void pentagon_Rat<BASE,RatPentSpecs>::get_sub_terms(const eval_param<T>& ep, pent_param<T,RatPentSpecs::MUPOINTS>& pp)
{
	momentum<complex<T> > K5sum_mom(ep.p(BASE::corner_ind(pp.pent_corner,1)-1)->P());
	for(size_t kiter=2;kiter<=BASE::corner_size(pp.pent_corner);kiter++){
		K5sum_mom+=ep.p(BASE::corner_ind(pp.pent_corner,kiter)-1)->P();
	}
	complex<T> S5=K5sum_mom.square();

	complex<T> npentpart;
	switch(pp.pent_corner){
	case 1://l_5=l-K_5
	case 5:
		npentpart=S5-T(2)*((pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P())*K5sum_mom);
		break;
	case 2://l_5=l_1-K_5=l-K_1-K_5
		npentpart=S5+T(2)*((pp.K1-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;
	case 3://l_5=l_4-K_5=l_1-K_4-K_5=l-K_1-K_4-K_5
		npentpart=S5+T(2)*((pp.K1+pp.K4-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;//l_4=l_1-K4
	case 4://l_5=l_2-K_5=l-K_2-K_5
		npentpart=S5+T(2)*((pp.K2-(pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()))*K5sum_mom);
		break;
	}

	// Construct the subtraction terms

	// Now compute the pentagon coeffs if needed
	if(!(ep.get_ID()==epID)){
		get_coeffs(ep);
	}

	//Set up the m0 masses we evaluate around
	const complex<T>* m0mass;
	(static_cast<typename RatPentSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatPentSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix);


	complex<T> pentcoeff_mass;
	get_E0coeff(pentcoeff_mass);



	complex<T> m0sol0, sum[RatPentSpecs::MUPOINTS-2]={complex<T>(0,0)};
	complex<T> vec1cK5(pp.vec1c*K5sum_mom), vec2cK5(pp.vec2c*K5sum_mom); // v1/v2 propto scale denfac propto scale boxpoles proptto 1/scale
	complex<T> fac1=pp.gamma/(T(2)*vec2cK5);
	complex<T> fac2=T(2)*pp.alp1*pp.alp2*vec2cK5;

    
	complex<T> fac3=-pp.gamma*pentcoeff_mass/(complex<T>(0,2)*vec2cK5);
    for(int m0pole=0;m0pole<RatPentSpecs::MUPOINTS;m0pole++){
		// the mu^2 subtraction terms, we evaluate at the value of t for the box
		m0sol0=fac1*(pp.boxpoles[2*m0pole]*(T(2)*pp.boxpoles[2*m0pole]*vec1cK5-npentpart)+fac2);
	for(int imupoints=0;imupoints<RatPentSpecs::MUPOINTS-2;imupoints++){
            pp.pent_sub_ret[imupoints]+=fac3*m0mass_matrix[m0pole+(imupoints+2)*RatPentSpecs::MUPOINTS]*pp.boxpoles[2*m0pole]/(m0mass[m0pole]-m0sol0);
	}

		m0sol0=fac1*(pp.boxpoles[1+2*m0pole]*(T(2)*pp.boxpoles[1+2*m0pole]*vec1cK5-npentpart)+fac2);
        for(int imupoints=0;imupoints<RatPentSpecs::MUPOINTS-2;imupoints++){
			pp.pent_sub_ret[imupoints]+=fac3*m0mass_matrix[m0pole+(imupoints+2)*RatPentSpecs::MUPOINTS]*pp.boxpoles[1+2*m0pole]/(m0mass[m0pole]-m0sol0);
		}
	}
    
	
	
	// Estimate the accuracy of the computation and if it is worse than the current value store it
	double pent_acc=to_double(MaxDigits<T>()-(log(abs(pentcoeff_mass)))/log(T(10)));
	if(pp.accuracy>pent_acc){
		pp.accuracy=pent_acc;
	}
}


template <class BASE, class RatPentSpecs> template <class T> void pentagon_Rat<BASE,RatPentSpecs>::get_coeffs(momentum_configuration<T>& mc, const std::vector<int>& ind, const long int original_mcID, const size_t vertex_ref)
{
	momentum<complex<T> > K1sum_mom(mc.mom(ind[BASE::corner_ind(1,1)-1]));
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		K1sum_mom+=mc.mom(ind[BASE::corner_ind(1,k1size+1)-1]);
	}
	momentum<complex<T> > K2sum_mom(-mc.mom(ind[BASE::corner_ind(5,1)-1]));
	size_t k2size=1;
	for(;k2size<BASE::corner_size(5);k2size++){
		K2sum_mom-=mc.mom(ind[BASE::corner_ind(5,k2size+1)-1]);
	}
	momentum<complex<T> > K4sum_mom(mc.mom(ind[BASE::corner_ind(2,1)-1]));
	for(size_t k4size=1;k4size<BASE::corner_size(2);k4size++){
		K4sum_mom+=mc.mom(ind[BASE::corner_ind(2,k4size+1)-1]);
	}
	momentum<complex<T> > K5sum_mom(mc.mom(ind[BASE::corner_ind(3,1)-1]));
	for(size_t k5size=1;k5size<BASE::corner_size(3);k5size++){
		K5sum_mom+=mc.mom(ind[BASE::corner_ind(3,k5size+1)-1]);
	}
	complex<T> S1=K1sum_mom.square();
	complex<T> S2=K2sum_mom.square();
	complex<T> S4=K4sum_mom.square();
	complex<T> S5=K5sum_mom.square();

	complex<T> K1K2=K1sum_mom*K2sum_mom;

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	complex<T> alp1(1,0), alp2(1,0);
	complex<T> f1, f2;
	if(!_k1massive){
		alp1=complex<T>(0,0);
		f2=complex<T>(1,0);
	}
	else{
		f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
	}
	if(!_k2massive){
		alp2=complex<T>(0,0);
		f1=complex<T>(1,0);
	}
	else{
		f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	complex<T> gamma=gammap/(f1*f2);
	complex<T> gamfac=gamma/(pow(gammap,2)-S1*S2);
	Cmom<T> K1flatc(f2*gamfac*(gammap*K1sum_mom-S1*K2sum_mom));
	Cmom<T> K2flatc(f1*gamfac*(gammap*K2sum_mom-S2*K1sum_mom));
	momentum<complex<T> > vec1c_mass=PfLLt(K1flatc.Lt(),K2flatc.L());
	momentum<complex<T> > vec2c_mass=PfLLt(K2flatc.Lt(),K1flatc.L());

	// Compute the l_5^2-mu^2=0 solution in terms of mu^2, we have l_5=l_4-K_5 and l_4=l_1-K_4, with l_1=l-K_1 and l_2=l+K_2
	complex<T> v1K4=vec1c_mass*K4sum_mom;
	complex<T> v2K4=vec2c_mass*K4sum_mom;
	complex<T> v1K5=vec1c_mass*K5sum_mom;
	complex<T> v2K5=vec2c_mass*K5sum_mom;

	complex<T> a1=S5-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom-K4sum_mom)*K5sum_mom);
	complex<T> a2=S4-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom)*K4sum_mom);
	complex<T> tsol_at_m0sol=T(0.5)*(a1*v2K4-a2*v2K5)/(v1K5*v2K4-v1K4*v2K5);
	complex<T> m0sol=gamma*(pow(tsol_at_m0sol,2)*v1K5-T(0.5)*a1*tsol_at_m0sol+alp1*alp2*v2K5)/v2K5;

	// Always put the branch cut in the same place for the sqrt(mass) computation
	if(abs(imag(m0sol))<DeltaZero<T>()){
		m0sol=complex<T>(real(m0sol),0);
	}

	// Set up the vectors for the trees processes at the corners
	for(size_t mm=1;mm<=BASE::corner_size(1);mm++){
		(indlst[0])[mm]=ind[BASE::corner_ind(1,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(2);mm++){
		(indlst[1])[mm]=ind[BASE::corner_ind(2,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(3);mm++){
		(indlst[2])[mm]=ind[BASE::corner_ind(3,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(4);mm++){
		(indlst[3])[mm]=ind[BASE::corner_ind(4,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(5);mm++){
		(indlst[4])[mm]=ind[BASE::corner_ind(5,mm)-1];
	}

	// We get the wrong sign on mass from sqrt(mass) in some cases as we are
	//  not consistently on one side of the sqrt branch cut so the code below fixes this.
	//  Ideally we would find a way of moving off m_0^2 off the negative real axis, not sure how to
	//  do this though as it is a physical pole.
	*(indlst[4].end()-3)=mc.insert(alp2*K1flatc.P()+alp1*K2flatc.P()+tsol_at_m0sol*vec1c_mass+(alp1*alp2-m0sol/gamma)/tsol_at_m0sol*vec2c_mass,_mt_massive);
//	*(indlst[4].end()-3)=mc.insert(alp2*K1flatc.P()+alp1*K2flatc.P()
//												+T(0.5)*((a1*v2K4-a2*v2K5)/(v1K5*v2K4-v1K4*v2K5))*vec1c_mass
//												+T(0.5)*((a2*v1K5-a1*v1K4)/(v1K5*v2K4-v1K4*v2K5))*vec2c_mass,_mt_massive);
	*(indlst[0].end()-3)=mc.insert(mc.mom(*(indlst[4].end()-3))-K1sum_mom,_mt_massive);
	*(indlst[1].end()-3)=mc.insert(mc.mom(*(indlst[0].end()-3))-K4sum_mom,_mt_massive);
	*(indlst[2].end()-3)=mc.insert(mc.mom(*(indlst[1].end()-3))-K5sum_mom,_mt_massive);
	*(indlst[3].end()-3)=mc.insert(mc.mom(*(indlst[4].end()-3))-K2sum_mom,_mt_massive);

	*(indlst[0].begin())=mc.insert(-mc.mom(*(indlst[4].end()-3)),_mt_massive);
	*(indlst[1].begin())=mc.insert(-mc.mom(*(indlst[0].end()-3)),_mt_massive);
	*(indlst[2].begin())=mc.insert(-mc.mom(*(indlst[1].end()-3)),_mt_massive);
	*(indlst[3].begin())=mc.insert(-mc.mom(*(indlst[2].end()-3)),_mt_massive);
	*(indlst[4].begin())=mc.insert(-mc.mom(*(indlst[3].end()-3)),_mt_massive);

//	size_t vertex_ref=mc.insert(momentum<complex<T> >(sqrt(complex<T>(-2,2)),complex<T>(0,1),complex<T>(1,1),complex<T>(0,1)),_mt_unknown);
	*(indlst[0].end()-2)=vertex_ref;
	*(indlst[1].end()-2)=vertex_ref;
	*(indlst[2].end()-2)=vertex_ref;
	*(indlst[3].end()-2)=vertex_ref;
	*(indlst[4].end()-2)=vertex_ref;

	size_t ref_mass=mc.insert(momentum<complex<T> >(m0sol,sqrt(m0sol),complex<T>(0,0),complex<T>(0,0)),_mt_massive);
	*(indlst[0].end()-1)=ref_mass;
	*(indlst[1].end()-1)=ref_mass;
	*(indlst[2].end()-1)=ref_mass;
	*(indlst[3].end()-1)=ref_mass;
	*(indlst[4].end()-1)=ref_mass;

	// Compute the penta-cut
	complex<T> pentcoeff_mass(0,0);
	//Construct the two-particle cut at the momenta above
	for(int i=0;i<BASE::decendant_nbr();i++){
		pentcoeff_mass+=BASE::eval_tree(i,0,mc,indlst[0])*BASE::eval_tree(i,1,mc,indlst[1])*BASE::eval_tree(i,2,mc,indlst[2])*BASE::eval_tree(i,3,mc,indlst[3])*BASE::eval_tree(i,4,mc,indlst[4]);
	}

	//Store the coefficient
	set_E0coeff(pentcoeff_mass);

	//Mark that we have computed these, we use the mcID passed from the bubble as this will
	// have been computed using a sub_mom_conf
	mcID=original_mcID;
	indID=ind;
}

template <class BASE, class RatPentSpecs> template <class T> void pentagon_Rat<BASE,RatPentSpecs>::get_coeffs(const eval_param<T>& ep)
{
#if _USE_ORIGINAL==1
	momentum<complex<T> > K1sum_mom(ep.p(BASE::corner_ind(1,1)-1)->P());
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		K1sum_mom+=ep.p(BASE::corner_ind(1,k1size+1)-1)->P();
	}
	momentum<complex<T> > K2sum_mom(-ep.p(BASE::corner_ind(5,1)-1)->P());
	size_t k2size=1;
	for(;k2size<BASE::corner_size(5);k2size++){
		K2sum_mom-=ep.p(BASE::corner_ind(5,k2size+1)-1)->P();
	}
	momentum<complex<T> > K4sum_mom(ep.p(BASE::corner_ind(2,1)-1)->P());
	for(size_t k4size=1;k4size<BASE::corner_size(2);k4size++){
		K4sum_mom+=ep.p(BASE::corner_ind(2,k4size+1)-1)->P();
	}
	momentum<complex<T> > K5sum_mom(ep.p(BASE::corner_ind(3,1)-1)->P());
	for(size_t k5size=1;k5size<BASE::corner_size(3);k5size++){
		K5sum_mom+=ep.p(BASE::corner_ind(3,k5size+1)-1)->P();
	}
	complex<T> S1=K1sum_mom.square();
	complex<T> S2=K2sum_mom.square();
	complex<T> S4=K4sum_mom.square();
	complex<T> S5=K5sum_mom.square();

	complex<T> K1K2=K1sum_mom*K2sum_mom;

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	complex<T> alp1(1,0), alp2(1,0);
	complex<T> f1, f2;
	if(!_k1massive){
		alp1=complex<T>(0,0);
		f2=complex<T>(1,0);
	}
	else{
		f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
	}
	if(!_k2massive){
		alp2=complex<T>(0,0);
		f1=complex<T>(1,0);
	}
	else{
		f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	complex<T> gamma=gammap/(f1*f2);
	complex<T> gamfac=gamma/(pow(gammap,2)-S1*S2);
	Cmom<T> K1flatc(f2*gamfac*(gammap*K1sum_mom-S1*K2sum_mom));
	Cmom<T> K2flatc(f1*gamfac*(gammap*K2sum_mom-S2*K1sum_mom));
	momentum<complex<T> > vec1c_mass=PfLLt(K1flatc.Lt(),K2flatc.L());
	momentum<complex<T> > vec2c_mass=PfLLt(K2flatc.Lt(),K1flatc.L());

	// Compute the l_5^2-mu^2=0 solution in terms of mu^2, we have l_5=l_4-K_5 and l_4=l_1-K_4, with l_1=l-K_1 and l_2=l+K_2
	complex<T> v1K4=vec1c_mass*K4sum_mom;
	complex<T> v2K4=vec2c_mass*K4sum_mom;
	complex<T> v1K5=vec1c_mass*K5sum_mom;
	complex<T> v2K5=vec2c_mass*K5sum_mom;

	complex<T> a1=S5-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom-K4sum_mom)*K5sum_mom);
	complex<T> a2=S4-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom)*K4sum_mom);
	complex<T> tsol_at_m0sol=T(0.5)*(a1*v2K4-a2*v2K5)/(v1K5*v2K4-v1K4*v2K5);
	complex<T> m0sol=gamma*(pow(tsol_at_m0sol,2)*v1K5-T(0.5)*a1*tsol_at_m0sol+alp1*alp2*v2K5)/v2K5;

	// Always put the branch cut in the same place for the sqrt(mass) computation
	if(abs(imag(m0sol))<DeltaZero<T>()){
		m0sol=complex<T>(real(m0sol),0);
	}

	eval_param<T>** epc;
	get_ep(epc);

	for(int mm=1;mm<=BASE::corner_size(1);mm++){
		epc[0]->set(mm,ep.p(BASE::corner_ind(1,mm)-1));
	}
	for(int mm=1;mm<=BASE::corner_size(2);mm++){
		epc[1]->set(mm,ep.p(BASE::corner_ind(2,mm)-1));
	}
	for(int mm=1;mm<=BASE::corner_size(3);mm++){
		epc[2]->set(mm,ep.p(BASE::corner_ind(3,mm)-1));
	}
	for(int mm=1;mm<=BASE::corner_size(4);mm++){
		epc[3]->set(mm,ep.p(BASE::corner_ind(4,mm)-1));
	}
	for(int mm=1;mm<=BASE::corner_size(5);mm++){
		epc[4]->set(mm,ep.p(BASE::corner_ind(5,mm)-1));
	}

	for(int icutmass=0;icutmass<_cut_mass.size();icutmass++){
		eval_param<T>::set_dynamic2(_cut_mass[icutmass], m0sol);
	}

	Cmom<T> l=Cmom<T>(alp2*K1flatc.P()+alp1*K2flatc.P()+tsol_at_m0sol*vec1c_mass+(alp1*alp2-m0sol/gamma)/tsol_at_m0sol*vec2c_mass,_mt_massive);
	Cmom<T> l1=Cmom<T>(l.P()-K1sum_mom,_mt_massive);
	Cmom<T> l2=Cmom<T>(l1.P()-K4sum_mom,_mt_massive);
	Cmom<T> l4=Cmom<T>(l2.P()-K5sum_mom,_mt_massive);
	Cmom<T> l5=Cmom<T>(l.P()-K2sum_mom,_mt_massive);
	Cmom<T> ml=Cmom<T>(-l.P(),_mt_massive);
	Cmom<T> ml1=Cmom<T>(-l1.P(),_mt_massive);
	Cmom<T> ml2=Cmom<T>(-l2.P(),_mt_massive);
	Cmom<T> ml4=Cmom<T>(-l4.P(),_mt_massive);
	Cmom<T> ml5=Cmom<T>(-l5.P(),_mt_massive);

	// We get the wrong sign on mass from sqrt(mass) in some cases as we are
	//  not consistently on one side of the sqrt branch cut so the code below fixes this.
	//  Ideally we would find a way of moving off m_0^2 off the negative real axis, not sure how to
	//  do this though as it is a physical pole.
	*(epc[4]->back())=&l;
	*(epc[0]->back())=&l1;
	*(epc[1]->back())=&l2;
	*(epc[2]->back())=&l4;
	*(epc[3]->back())=&l5;
	*(epc[0]->begin())=&ml;
	*(epc[1]->begin())=&ml1;
	*(epc[2]->begin())=&ml2;
	*(epc[3]->begin())=&ml4;
	*(epc[4]->begin())=&ml5;

	// Compute the penta-cut
	complex<T> pentcoeff_mass(0,0);
	//Construct the two-particle cut at the momenta above
	for(int i=0;i<BASE::decendant_nbr();i++){
		pentcoeff_mass+=BASE::eval_tree(i,0,*(epc[0]))*BASE::eval_tree(i,1,*(epc[1]))*BASE::eval_tree(i,2,*(epc[2]))*BASE::eval_tree(i,3,*(epc[3]))*BASE::eval_tree(i,4,*(epc[4]));
	}

	//Store the coefficient
	set_E0coeff(pentcoeff_mass);

	//Mark that we have computed these, we set the epID
	epID=ep.get_ID();
#else
	// Compute the trees for each corner
	vector<complex<T> > tr[5]={vector<complex<T> >(BASE::decendant_nbr())};
    typename RatPentSpecs::CornerTreeStrategy* ptr= static_cast<typename RatPentSpecs::CornerTreeStrategy*>(this);
    ptr->fill_trees(ep,tr,_cut_mass,*this);

    // Compute the coefficient
    complex<T> result(0,0);
	for(int i=0;i<BASE::decendant_nbr();i++){
		result+=(tr[i])[0]*(tr[i])[1]*(tr[i])[2]*(tr[i])[3]*(tr[i])[4];
	}

	//Store the coefficient
	set_E0coeff(result/T(2));

	//Mark that we have computed these, we set the epID
	epID=ep.get_ID();

#endif
}



template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> template <class T,class cutDbase> void pentagon_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::compute_momenta(const eval_param<T>& ep,Cmom<T> (&momenta)[5],Cmom<T> (&momenta_m)[5], complex<T>& pm0sol,const cutDbase& cb){
	momentum<complex<T> > K1sum_mom(ep.p(cb.corner_ind(1,1)-1)->P());
	size_t k1size=1;
	for(;k1size<cb.corner_size(1);k1size++){
		K1sum_mom+=ep.p(cb.corner_ind(1,k1size+1)-1)->P();
	}
	momentum<complex<T> > K2sum_mom(-ep.p(cb.corner_ind(5,1)-1)->P());
	size_t k2size=1;
	for(;k2size<cb.corner_size(5);k2size++){
		K2sum_mom-=ep.p(cb.corner_ind(5,k2size+1)-1)->P();
	}
	momentum<complex<T> > K4sum_mom(ep.p(cb.corner_ind(2,1)-1)->P());
	for(size_t k4size=1;k4size<cb.corner_size(2);k4size++){
		K4sum_mom+=ep.p(cb.corner_ind(2,k4size+1)-1)->P();
	}
	momentum<complex<T> > K5sum_mom(ep.p(cb.corner_ind(3,1)-1)->P());
	for(size_t k5size=1;k5size<cb.corner_size(3);k5size++){
		K5sum_mom+=ep.p(cb.corner_ind(3,k5size+1)-1)->P();
	}
	complex<T> S1=K1sum_mom.square();
	complex<T> S2=K2sum_mom.square();
	complex<T> S4=K4sum_mom.square();
	complex<T> S5=K5sum_mom.square();

	complex<T> K1K2=K1sum_mom*K2sum_mom;

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	complex<T> gammap;
	if((k1size==1)||(k2size==1)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	complex<T> alp1(1,0), alp2(1,0);
	complex<T> f1, f2;
	if(k1size==1){
		alp1=complex<T>(0,0);
		f2=complex<T>(1,0);
	}
	else{
		f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
	}
	if(k2size==1){
		alp2=complex<T>(0,0);
		f1=complex<T>(1,0);
	}
	else{
		f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	complex<T> gamma=gammap/(f1*f2);
	complex<T> gamfac=gamma/(pow(gammap,2)-S1*S2);
	Cmom<T> K1flatc(f2*gamfac*(gammap*K1sum_mom-S1*K2sum_mom));
	Cmom<T> K2flatc(f1*gamfac*(gammap*K2sum_mom-S2*K1sum_mom));
	momentum<complex<T> > vec1c_mass=PfLLt(K1flatc.Lt(),K2flatc.L());
	momentum<complex<T> > vec2c_mass=PfLLt(K2flatc.Lt(),K1flatc.L());

	// Compute the l_5^2-mu^2=0 solution in terms of mu^2, we have l_5=l_4-K_5 and l_4=l_1-K_4, with l_1=l-K_1 and l_2=l+K_2
	complex<T> v1K4=vec1c_mass*K4sum_mom;
	complex<T> v2K4=vec2c_mass*K4sum_mom;
	complex<T> v1K5=vec1c_mass*K5sum_mom;
	complex<T> v2K5=vec2c_mass*K5sum_mom;

	complex<T> a1=S5-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom-K4sum_mom)*K5sum_mom);
	complex<T> a2=S4-T(2)*((alp2*K1flatc.P()+alp1*K2flatc.P()-K1sum_mom)*K4sum_mom);
	complex<T> tsol_at_m0sol=T(0.5)*(a1*v2K4-a2*v2K5)/(v1K5*v2K4-v1K4*v2K5);
	complex<T> m0sol=gamma*(pow(tsol_at_m0sol,2)*v1K5-T(0.5)*a1*tsol_at_m0sol+alp1*alp2*v2K5)/v2K5;

	// Always put the branch cut in the same place for the sqrt(mass) computation
	if(abs(imag(m0sol))<DeltaZero<T>()){
		pm0sol=complex<T>(real(m0sol),0);
	}

	momenta[0]=Cmom<T>(alp2*K1flatc.P()+alp1*K2flatc.P()+tsol_at_m0sol*vec1c_mass+(alp1*alp2-pm0sol/gamma)/tsol_at_m0sol*vec2c_mass,_mt_massive);
	momenta[1]=Cmom<T>(momenta[0].P()-K1sum_mom,_mt_massive);
	momenta[2]=Cmom<T>(momenta[1].P()-K4sum_mom,_mt_massive);
	momenta[3]=Cmom<T>(momenta[2].P()-K5sum_mom,_mt_massive);
	momenta[4]=Cmom<T>(momenta[0].P()-K2sum_mom,_mt_massive);
	momenta_m[0]=Cmom<T>(-momenta[0].P(),_mt_massive);
	momenta_m[1]=Cmom<T>(-momenta[1].P(),_mt_massive);
	momenta_m[2]=Cmom<T>(-momenta[2].P(),_mt_massive);
	momenta_m[3]=Cmom<T>(-momenta[3].P(),_mt_massive);
	momenta_m[4]=Cmom<T>(-momenta[4].P(),_mt_massive);
}



template <class MomentumEvaluator,class CutType> template <class T> void Normal_RatPent_Corner_Tree_Strategy<MomentumEvaluator,CutType>::fill_trees(
		const eval_param<T>& ep,
		std::vector<std::complex<T> > (&trees_result)[5],
		std::vector<size_t>& cut_mass,
		CutType& self)
{
	Cmom<T> momenta[5];
 	Cmom<T> momenta_m[5];
 	complex<T> m0sol;
	this->template compute_momenta<T,CutType>(ep,momenta,momenta_m,m0sol,self);

	eval_param<T>* epc[5];
	epc[0]=new eval_param<T>(self.corner_size(1)+2);
	epc[1]=new eval_param<T>(self.corner_size(2)+2);
	epc[2]=new eval_param<T>(self.corner_size(3)+2);
	epc[3]=new eval_param<T>(self.corner_size(4)+2);
	epc[4]=new eval_param<T>(self.corner_size(5)+2);
//	eval_param<T>** epc;
//	self.get_ep(epc);

	for(int mm=1;mm<=self.corner_size(1);mm++){
		epc[0]->set(mm,ep.p(self.corner_ind(1,mm)-1));
	}
	for(int mm=1;mm<=self.corner_size(2);mm++){
		epc[1]->set(mm,ep.p(self.corner_ind(2,mm)-1));
	}
	for(int mm=1;mm<=self.corner_size(3);mm++){
		epc[2]->set(mm,ep.p(self.corner_ind(3,mm)-1));
	}
	for(int mm=1;mm<=self.corner_size(4);mm++){
		epc[3]->set(mm,ep.p(self.corner_ind(4,mm)-1));
	}
	for(int mm=1;mm<=self.corner_size(5);mm++){
		epc[4]->set(mm,ep.p(self.corner_ind(5,mm)-1));
	}

	for(int icutmass=0;icutmass<cut_mass.size();icutmass++){
		eval_param<T>::set_dynamic2(cut_mass[icutmass], m0sol);
	}

	// We get the wrong sign on mass from sqrt(mass) in some cases as we are
	//  not consistently on one side of the sqrt branch cut so the code below fixes this.
	//  Ideally we would find a way of moving off m_0^2 off the negative real axis, not sure how to
	//  do this though as it is a physical pole.
	*(epc[4]->back())=&momenta[0];
	*(epc[0]->back())=&momenta[1];
	*(epc[1]->back())=&momenta[2];
	*(epc[2]->back())=&momenta[3];
	*(epc[3]->back())=&momenta[4];
	*(epc[0]->begin())=&momenta_m[0];
	*(epc[1]->begin())=&momenta_m[1];
	*(epc[2]->begin())=&momenta_m[2];
	*(epc[3]->begin())=&momenta_m[3];
	*(epc[4]->begin())=&momenta_m[4];


	//Compute all the possible trees
	for(int n=0;n<5;n++){
		for(int i=0;i<self.decendant_nbr();i++){
			(trees_result[i])[n]=self.eval_tree(i,n,*(epc[n]));
		}
	}

	delete epc[0];
	delete epc[1];
	delete epc[2];
	delete epc[3];
	delete epc[4];
}




}
}
#endif /* PENTAGON_RATEXT_HPP_ */
