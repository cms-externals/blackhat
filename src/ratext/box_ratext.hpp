/*
 * box_ratext.hpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef BOX_RATEXT_HPP_
#define BOX_RATEXT_HPP_


#include "ratext/box_ratext.h"
#include "ratext/rat_ext.h"
#include "ratext/triangle_ratext.h"
#include <cassert>

using std::complex;

namespace BH {

#define _USE_ORIGINAL 1

namespace ratext {

template <class BASE, class RatBoxSpecs> C box_Rat<BASE,RatBoxSpecs>::eval(const eval_param<R>& ep)
{
	return this->template get_symmetry_factor<R>()*D0coeff;
}

template <class BASE, class RatBoxSpecs> CHP box_Rat<BASE,RatBoxSpecs>::eval(const eval_param<RHP>& ep)
{
	return this->template get_symmetry_factor<RHP>()*D0coeff_HP;
}
template <class BASE, class RatBoxSpecs> CVHP box_Rat<BASE,RatBoxSpecs>::eval(const eval_param<RVHP>& ep)
{
	return this->template get_symmetry_factor<RVHP>()*D0coeff_VHP;
}

template <class BASE, class RatBoxSpecs> C box_Rat<BASE,RatBoxSpecs>::eval(mom_conf& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("box_Rat<BASE,RatBoxSpecs> : Attempting to compute a box before computing the bubble");
	}
	return this->template get_symmetry_factor<R>()*D0coeff;
}

template <class BASE, class RatBoxSpecs> CHP box_Rat<BASE,RatBoxSpecs>::eval(mom_conf_HP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("box_Rat<BASE,RatBoxSpecs> : Attempting to compute a box before computing the bubble");
	}
	return this->template get_symmetry_factor<RHP>()*D0coeff_HP;
}

template <class BASE, class RatBoxSpecs> CVHP box_Rat<BASE,RatBoxSpecs>::eval(mom_conf_VHP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("box_Rat<BASE,RatBoxSpecs> : Attempting to compute a box before computing the bubble");
	}
	return this->template get_symmetry_factor<RVHP>()*D0coeff_VHP;
}


#if BH_USE_GMP

template <class BASE, class RatBoxSpecs> CGMP box_Rat<BASE,RatBoxSpecs>::eval(const eval_param<RGMP>& ep)
{
	return this->template get_symmetry_factor<RGMP>()*D0coeff_GMP;
}
template <class BASE, class RatBoxSpecs> CGMP box_Rat<BASE,RatBoxSpecs>::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)&&(ind!=indID)){
//		_WARNING("box_Rat<BASE,RatBoxSpecs> : Attempting to compute a box before computing the bubble");
	}
	return this->template get_symmetry_factor<RGMP>()*D0coeff_GMP;
}
#endif

template <class BASE, class RatBoxSpecs> template <class T> void box_Rat<BASE,RatBoxSpecs>::get_sub_terms(const eval_param<T>& ep, box_param<T,RatBoxSpecs::MUTRIPOINTS,RatBoxSpecs::CPOINTS>& bp)
{
	momentum<complex<T> > K4sum_mom(ep.p(BASE::corner_ind(bp.box_corner,1)-1)->P());
	for(size_t kiter=2;kiter<=BASE::corner_size(bp.box_corner);kiter++){
		K4sum_mom+=ep.p(BASE::corner_ind(bp.box_corner,kiter)-1)->P();
	}
	
	//Set up the alphas so that they match the L_4=l_i-K_4 (we do not always have L_4=l-K_4)
	complex<T> nboxpart1;
	switch(bp.box_corner){
	case 1:	//l_4=l-K-4
	case 4:
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P())*K4sum_mom);
		break;
	case 2: //l_4=l_1-K_4=l-K_1-K_4
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()-bp.K1)*K4sum_mom);
		break;
	case 3: //l_4=l_2-K_4=l-K_2-K_4
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()-bp.K2)*K4sum_mom);
		break;
	}

	//We now find the number of pole solutions as well as calculating them
	bp.denfac1=T(2)*(bp.vec1c*K4sum_mom);

	//Set up the m0 masses and t values so sum over all RATPOINT terms
	const complex<T>* m0mass;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_eval_pts(m0mass);

	complex<T> nboxpart2;
	for(int m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		nboxpart2=sqrt(pow(nboxpart1,2)-T(8)*(bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)*bp.denfac1*(bp.vec2c*K4sum_mom));
		//There will always be two poles when we have a massive leg
		bp.boxpoles[2*m0pole]=(nboxpart1+nboxpart2)/(T(2)*bp.denfac1);
		bp.boxpoles[1+2*m0pole]=(nboxpart1-nboxpart2)/(T(2)*bp.denfac1);
	}

	// Get the points to evaluate at for the correct precision
	const complex<T>* eval_pts;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->TriPoints.get_t_eval_pts(eval_pts);

	// Compute the coefficients if we have not already done so
	if(!(ep.get_ID()==epID)){
		get_coeffs(ep,bp.accuracy);
	}

	// Now match the solutions to the correct pole
	complex<T> ratio, Tn;
	Cmom<T> K1fs, K2fs;
	momentum<complex<T> > K4s;
	//We first select where the first solution should appear and then use the T(q) test to
	//  find the overall ordering relative to the stored solution
	// This is done separately for each m0pole term as the ratio can vary as we go around the circle
	for(size_t m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		get_denfac(m0pole,Tn);
		get_Kif_coeffs(m0pole,K1fs,K2fs,K4s);
		momentum<complex<T> > lcmp=bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()+bp.boxpoles[1+2*m0pole]*bp.vec1c+(bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)/bp.boxpoles[1+2*m0pole]*bp.vec2c;
		ratio=-((K2fs.Lt()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.Lt())*(K1fs.L()*K2fs.L())
				-(K2fs.L()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.L())*(K1fs.Lt()*K2fs.Lt()))/Tn;
		// ratio should always be real and either 1 or -1, the order of the returned solutions is reversed if it is -1
		if(ratio.real()>T(0)){
			get_coeff(m0pole,bp.boxcoeffs[2*m0pole],bp.boxcoeffs[1+2*m0pole]);
		}
		else{
			get_coeff(m0pole,bp.boxcoeffs[1+2*m0pole],bp.boxcoeffs[2*m0pole]);
		}
	}

	//Now construct the t and mu^2 subtraction terms, here we subtract from the triangle so use RATPOINTS_SUB terms
	for(int m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		for(int icirc=0; icirc<RatBoxSpecs::CPOINTS; icirc++){
			//CHANGE : Replace this with a ((d_+ + d_-)+(d_-+d_+)/something)/l_4^2 so that we subtract off any potential double poles correctly
			bp.tri_sub_ret[icirc+(RatBoxSpecs::CPOINTS+1)*m0pole]+=complex<T>(0,-1)*(bp.boxpoles[2*m0pole]*bp.boxcoeffs[2*m0pole]/(eval_pts[icirc]-bp.boxpoles[2*m0pole])
					-bp.boxpoles[1+2*m0pole]*bp.boxcoeffs[1+2*m0pole]/(eval_pts[icirc]-bp.boxpoles[1+2*m0pole]))
					/(bp.denfac1*(bp.boxpoles[1+2*m0pole]-bp.boxpoles[2*m0pole]));
		}
		bp.tri_sub_ret[RatBoxSpecs::CPOINTS+(RatBoxSpecs::CPOINTS+1)*m0pole]+=(bp.boxcoeffs[2*m0pole]-bp.boxcoeffs[1+2*m0pole])/(complex<T>(0,-2)*bp.denfac1*(bp.boxpoles[1+2*m0pole]-bp.boxpoles[2*m0pole]));
	}
}

template <class BASE, class RatBoxSpecs> template <class T> void box_Rat<BASE,RatBoxSpecs>::get_sub_terms(momentum_configuration<T>& mc, const std::vector<int>& ind, box_param<T,RatBoxSpecs::MUTRIPOINTS,RatBoxSpecs::CPOINTS>& bp)
{
	momentum<complex<T> > K4sum_mom(mc.mom(ind[BASE::corner_ind(bp.box_corner,1)-1]));
	for(size_t kiter=2;kiter<=BASE::corner_size(bp.box_corner);kiter++){
		K4sum_mom+=mc.mom(ind[BASE::corner_ind(bp.box_corner,kiter)-1]);
	}

	//Set up the alphas so that they match the L_4=l_i-K_4 (we do not always have L_4=l-K_4)
	complex<T> nboxpart1;
	switch(bp.box_corner){
	case 1:	//l_4=l-K-4
	case 4:
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P())*K4sum_mom);
		break;
	case 2: //l_4=l_1-K_4=l-K_1-K_4
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()-bp.K1)*K4sum_mom);
		break;
	case 3: //l_4=l_2-K_4=l-K_2-K_4
		nboxpart1=K4sum_mom.square()-T(2)*((bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()-bp.K2)*K4sum_mom);
		break;
	}

	//We now find the number of pole solutions as well as calculating them
	bp.denfac1=T(2)*(bp.vec1c*K4sum_mom);

	//Set up the m0 masses and t values so sum over all RATPOINT terms
	const complex<T>* m0mass;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_eval_pts(m0mass);

	complex<T> nboxpart2;
	for(int m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		nboxpart2=sqrt(pow(nboxpart1,2)-T(8)*(bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)*bp.denfac1*(bp.vec2c*K4sum_mom));
		//There will always be two poles when we have a massive leg
		bp.boxpoles[2*m0pole]=(nboxpart1+nboxpart2)/(T(2)*bp.denfac1);
		bp.boxpoles[1+2*m0pole]=(nboxpart1-nboxpart2)/(T(2)*bp.denfac1);
	}

	// Get the points to evaluate at for the correct precision
	const complex<T>* eval_pts;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->TriPoints.get_t_eval_pts(eval_pts);

	// Compute the coefficients if we have not already done so
	if(!((bp.mcID==mcID)&&(ind==indID))){
		get_coeffs(mc,ind,bp.accuracy,bp.mcID,bp.vertex_ref);
	}

	// Now match the solutions to the correct pole
	complex<T> ratio, Tn;
	Cmom<T> K1fs, K2fs;
	momentum<complex<T> > K4s;
	//We first select where the first solution should appear and then use the T(q) test to
	//  find the overall ordering relative to the stored solution
	// This is done separately for each m0pole term as the ratio can vary as we go around the circle
	for(size_t m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		get_denfac(m0pole,Tn);
		get_Kif_coeffs(m0pole,K1fs,K2fs,K4s);
		momentum<complex<T> > lcmp=bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P()+bp.boxpoles[1+2*m0pole]*bp.vec1c+(bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)/bp.boxpoles[1+2*m0pole]*bp.vec2c;
		ratio=-((K2fs.Lt()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.Lt())*(K1fs.L()*K2fs.L())
				-(K2fs.L()*smatrix<T>(K4s)*smatrix<T>(lcmp)*K1fs.L())*(K1fs.Lt()*K2fs.Lt()))/Tn;
		// ratio should always be real and either 1 or -1, the order of the returned solutions is reversed if it is -1
		if(ratio.real()>T(0)){
			get_coeff(m0pole,bp.boxcoeffs[2*m0pole],bp.boxcoeffs[1+2*m0pole]);
		}
		else{
			get_coeff(m0pole,bp.boxcoeffs[1+2*m0pole],bp.boxcoeffs[2*m0pole]);
		}
	}

	//Now construct the t and mu^2 subtraction terms, here we subtract from the triangle so use RATPOINTS_SUB terms
	for(int m0pole=0;m0pole<RatBoxSpecs::MUTRIPOINTS;m0pole++){
		for(int icirc=0; icirc<RatBoxSpecs::CPOINTS; icirc++){
			//CHANGE : Replace this with a ((d_+ + d_-)+(d_-+d_+)/something)/l_4^2 so that we subtract off any potential double poles correctly
			bp.tri_sub_ret[icirc+(RatBoxSpecs::CPOINTS+1)*m0pole]+=complex<T>(0,-1)*(bp.boxpoles[2*m0pole]*bp.boxcoeffs[2*m0pole]/(eval_pts[icirc]-bp.boxpoles[2*m0pole])
					-bp.boxpoles[1+2*m0pole]*bp.boxcoeffs[1+2*m0pole]/(eval_pts[icirc]-bp.boxpoles[1+2*m0pole]))
					/(bp.denfac1*(bp.boxpoles[1+2*m0pole]-bp.boxpoles[2*m0pole]));
		}
		bp.tri_sub_ret[RatBoxSpecs::CPOINTS+(RatBoxSpecs::CPOINTS+1)*m0pole]+=(bp.boxcoeffs[2*m0pole]-bp.boxcoeffs[1+2*m0pole])/(complex<T>(0,-2)*bp.denfac1*(bp.boxpoles[1+2*m0pole]-bp.boxpoles[2*m0pole]));
	}
}


template <class BASE, class RatBoxSpecs> template <class T> void box_Rat<BASE,RatBoxSpecs>::get_coeffs_fn(const eval_param<T>& ep, double& accuracy)
{
	pent_param<T,RatBoxSpecs::MUPOINTS> pp;
	pp.accuracy=accuracy;

	pp.K1=ep.p(BASE::corner_ind(1,1)-1)->P();
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		pp.K1+=ep.p(BASE::corner_ind(1,k1size+1)-1)->P();
	}
	pp.K2=-ep.p(BASE::corner_ind(4,1)-1)->P();
	size_t k2size=1;
	for(;k2size<BASE::corner_size(4);k2size++){
		pp.K2-=ep.p(BASE::corner_ind(4,k2size+1)-1)->P();
	}
	pp.K4=ep.p(BASE::corner_ind(2,1)-1)->P();
	for(size_t k4size=1;k4size<BASE::corner_size(2);k4size++){
		pp.K4+=ep.p(BASE::corner_ind(2,k4size+1)-1)->P();
	}
	complex<T> S1=pp.K1.square();
	complex<T> S2=pp.K2.square();

	complex<T> K1K2=pp.K1*pp.K2;

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	complex<T> f2;
	if(!_k1massive){
		pp.alp1=complex<T>(0,0);
		f2=complex<T>(1,0);
	}
	else{
		pp.alp1=complex<T>(1,0);
		f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
	}
	complex<T> f1;
	if(!_k2massive){
		pp.alp2=complex<T>(0,0);
		f1=complex<T>(1,0);
	}
	else{
		pp.alp2=complex<T>(1,0);
		f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	pp.gamma=gammap/(f1*f2);
	complex<T> gamfac=pp.gamma/(pow(gammap,2)-S1*S2);
	pp.K1flatc=Cmom<T>(f2*gamfac*(gammap*pp.K1-S1*pp.K2));
	pp.K2flatc=Cmom<T>(f1*gamfac*(gammap*pp.K2-S2*pp.K1));
	pp.vec1c=PfLLt(pp.K1flatc.Lt(),pp.K2flatc.L());
	pp.vec2c=PfLLt(pp.K2flatc.Lt(),pp.K1flatc.L());

	//We now find the number of pole solutions as well as calculating them
	complex<T> denfac1=T(2)*(pp.vec1c*pp.K4);
	complex<T> denfac2=T(2)*(pp.vec2c*pp.K4);
	complex<T> nboxpart1=pp.K4.square()-T(2)*((pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()-pp.K1)*pp.K4);

	const complex<T>* m0mass;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_eval_pts(m0mass);

	complex<T> nboxpart2;
	for(int m0pole=0;m0pole<RatBoxSpecs::MUPOINTS;m0pole++){
		nboxpart2=sqrt(pow(nboxpart1,2)-T(4)*(pp.alp1*pp.alp2-m0mass[m0pole]/pp.gamma)*denfac1*denfac2);
		//There will always be two poles when we have a massive leg
		pp.boxpoles[2*m0pole]=(nboxpart1+nboxpart2)/(T(2)*denfac1);
		pp.boxpoles[1+2*m0pole]=(nboxpart1-nboxpart2)/(T(2)*denfac1);
	}

	//Clear out the previous computation
	for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS*(RatBoxSpecs::MUPOINTS-2);imupoints++){
		pp.pent_sub_ret[imupoints]=complex<T>(0,0);
	}
	//Get all the pentagon subtractions
	for(size_t pent=1; pent <= BASE::daughters_nbr(); pent++){
		pp.pent_corner=BASE::get_opened_corner(pent);
#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(pent));
		assert(dau);
#else
		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(pent));
#endif
		dau->get_sub_terms(ep,pp);
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

	Cmom<T> l[4],ml[4];
	*(epc[0]->back())=&l[1];
	*(epc[1]->back())=&l[2];
	*(epc[2]->back())=&l[3];
	*(epc[3]->back())=&l[0];
	*(epc[0]->begin())=&ml[0];
	*(epc[1]->begin())=&ml[1];
	*(epc[2]->begin())=&ml[2];
	*(epc[3]->begin())=&ml[3];

	//Sum over the number of poles we had
	momentum<complex<T> > mompt1=pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P();

	complex<T> test[RatBoxSpecs::MUPOINTS-2]={complex<T>(0,0)};
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_matrix_eval_pts(m0mass_matrix);
	for(size_t m0pole=0;m0pole<RatBoxSpecs::MUPOINTS;m0pole++){
		// Clear any previously stored results
		complex<T> boxcoeffs_mass[2]={complex<T>(0,0)};
		for(int icutmass=0;icutmass<_cut_mass.size();icutmass++){
			eval_param<T>::set_dynamic2(_cut_mass[icutmass], m0mass[m0pole]);
		}
		for(size_t npole=0; npole<2; npole++){
			l[0]=Cmom<T>(mompt1+pp.boxpoles[npole+2*m0pole]*pp.vec1c+(pp.alp1*pp.alp2-m0mass[m0pole]/pp.gamma)/pp.boxpoles[npole+2*m0pole]*pp.vec2c,_mt_massive);
			l[1]=Cmom<T>(l[0].P()-pp.K1,_mt_massive);
			l[2]=Cmom<T>(l[1].P()-pp.K4,_mt_massive);
			l[3]=Cmom<T>(l[0].P()-pp.K2,_mt_massive);
			ml[0]=Cmom<T>(-l[0].P(),_mt_massive);
			ml[1]=Cmom<T>(-l[1].P(),_mt_massive);
			ml[2]=Cmom<T>(-l[2].P(),_mt_massive);
			ml[3]=Cmom<T>(-l[3].P(),_mt_massive);

			//Construct the three-particle cut subtraction at the above momenta and compute the amplitude
			//If there is only one solution we want to store this in the correct location
			//Construct the two-particle cut at the momenta above
			complex<T> temp_res(0,0);
			for(int i=0;i<BASE::decendant_nbr();i++){
				temp_res+=BASE::eval_tree(i,0,*(epc[0]))*BASE::eval_tree(i,1,*(epc[1]))*BASE::eval_tree(i,2,*(epc[2]))*BASE::eval_tree(i,3,*(epc[3]));
			}
			// We only use the first two m0pole results in the triangle subtraction
			boxcoeffs_mass[npole]+=temp_res;
			
			for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS-2;imupoints++){
				test[imupoints]+=m0mass_matrix[m0pole+(imupoints+2)*RatBoxSpecs::MUPOINTS]*temp_res;
			}
		}

		//Store information allowing us to decide how to choose how to map t_+ and t_- to d_+ and d_- later on.
		// We do this by borrowing the method hidden in the OPP approach basically we can rewrite the OPP expresion
		//  for the subtracted box terms as
		//  (1/2)((d_+(1+T(l'_+)/T(l_+))+d_-(1-T(l'_+)/T(l_+)))/(t-t_+)+(d_+(1+T(l'_-)/T(l_+))+d_-(1-T(l'_-)/T(l_+)))/(t-t_-))
		//  with T(l)=Tr(l,K1f,K2f,K4) and l the original solution of l to find d_+ and d_- and l' the new solution for
		//  l we are trying to map d_+ and d_- to.
		//  The ratio of T(l'_+)/T(l_+) will be either +1 or -1 depending upon the order of the momenta entering l
		//  and the sign will correctly choose the mapping of t_i->d_i.
		set_Kif_coeffs(m0pole,pp.K1flatc,pp.K2flatc,pp.K4);
		//Flip the sign as this is T(l_-) and we store T(l_+)=-T(l_-)
		// Here l[0] corresponds to the result at m0point[RatBoxSpecs::MUPOINTS-1]
		set_denfac(m0pole,-(pp.K2flatc.Lt()*smatrix<T>(pp.K4)*smatrix<T>(l[0].P())*pp.K1flatc.Lt())*(pp.K1flatc.L()*pp.K2flatc.L())
				+(pp.K2flatc.L()*smatrix<T>(pp.K4)*smatrix<T>(l[0].P())*pp.K1flatc.L())*(pp.K1flatc.Lt()*pp.K2flatc.Lt()));

		//Store the coefficients
		set_coeff(0,m0pole,boxcoeffs_mass[0]);
		set_coeff(1,m0pole,boxcoeffs_mass[1]);
	}
	complex<T> fullboxcoeff(0,0);

	for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS-2;imupoints++){
		fullboxcoeff+=(test[imupoints]+pp.pent_sub_ret[imupoints])*(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1);
	}
	set_D0coeff(fullboxcoeff);

	//Mark that we have computed these, we use the mcID passed from the bubble as this will
	// have been computed using a sub_mom_conf
	epID=ep.get_ID();
	accuracy=pp.accuracy;
}

template <class BASE, class RatBoxSpecs> template <class T> void box_Rat<BASE,RatBoxSpecs>::get_coeffs_fn(momentum_configuration<T>& mc, const std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref)
{
	pent_param<T,RatBoxSpecs::MUPOINTS> pp;
	pp.accuracy=accuracy;
	pp.mcID=original_mcID;
	pp.vertex_ref=vertex_ref;

	pp.K1=mc.mom(ind[BASE::corner_ind(1,1)-1]);
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		pp.K1+=mc.mom(ind[BASE::corner_ind(1,k1size+1)-1]);
	}
	pp.K2=-mc.mom(ind[BASE::corner_ind(4,1)-1]);
	size_t k2size=1;
	for(;k2size<BASE::corner_size(4);k2size++){
		pp.K2-=mc.mom(ind[BASE::corner_ind(4,k2size+1)-1]);
	}
	pp.K4=mc.mom(ind[BASE::corner_ind(2,1)-1]);
	for(size_t k4size=1;k4size<BASE::corner_size(2);k4size++){
		pp.K4+=mc.mom(ind[BASE::corner_ind(2,k4size+1)-1]);
	}
	complex<T> S1=pp.K1.square();
	complex<T> S2=pp.K2.square();

	complex<T> K1K2=pp.K1*pp.K2;

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	complex<T> f2;
	if(!_k1massive){
		pp.alp1=complex<T>(0,0);
		f2=complex<T>(1,0);
	}
	else{
		pp.alp1=complex<T>(1,0);
		f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
	}
	complex<T> f1;
	if(!_k2massive){
		pp.alp2=complex<T>(0,0);
		f1=complex<T>(1,0);
	}
	else{
		pp.alp2=complex<T>(1,0);
		f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	pp.gamma=gammap/(f1*f2);
	complex<T> gamfac=pp.gamma/(pow(gammap,2)-S1*S2);
	pp.K1flatc=Cmom<T>(f2*gamfac*(gammap*pp.K1-S1*pp.K2));
	pp.K2flatc=Cmom<T>(f1*gamfac*(gammap*pp.K2-S2*pp.K1));
	pp.vec1c=PfLLt(pp.K1flatc.Lt(),pp.K2flatc.L());
	pp.vec2c=PfLLt(pp.K2flatc.Lt(),pp.K1flatc.L());

	//We now find the number of pole solutions as well as calculating them
	complex<T> denfac1=T(2)*(pp.vec1c*pp.K4);
	complex<T> denfac2=T(2)*(pp.vec2c*pp.K4);
	complex<T> nboxpart1=pp.K4.square()-T(2)*((pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P()-pp.K1)*pp.K4);

	const complex<T>* m0mass;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_eval_pts(m0mass);

	complex<T> nboxpart2;
	for(int m0pole=0;m0pole<RatBoxSpecs::MUPOINTS;m0pole++){
		nboxpart2=sqrt(pow(nboxpart1,2)-T(4)*(pp.alp1*pp.alp2-m0mass[m0pole]/pp.gamma)*denfac1*denfac2);
		//There will always be two poles when we have a massive leg
		pp.boxpoles[2*m0pole]=(nboxpart1+nboxpart2)/(T(2)*denfac1);
		pp.boxpoles[1+2*m0pole]=(nboxpart1-nboxpart2)/(T(2)*denfac1);
	}

	//Clear out the previous computation
	for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS-2;imupoints++){
		pp.pent_sub_ret[imupoints]=complex<T>(0,0);
	}
	//Get all the pentagon subtractions
	for(size_t pent=1; pent <= BASE::daughters_nbr(); pent++){
		pp.pent_corner=BASE::get_opened_corner(pent);

		typedef typename daughter_info<box_Rat<BASE,RatBoxSpecs> >::daughter_type* dau_type_p ;
#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(pent));
		assert(dau);
#else
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(pent));
#endif
		dau->get_sub_terms(mc,ind,pp);
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

	*(indlst[0].end()-2)=vertex_ref;
	*(indlst[1].end()-2)=vertex_ref;
	*(indlst[2].end()-2)=vertex_ref;
	*(indlst[3].end()-2)=vertex_ref;

	//Sum over the number of poles we had
	momentum<complex<T> > mompt1=pp.alp2*pp.K1flatc.P()+pp.alp1*pp.K2flatc.P();
	complex<T> test[2*RatBoxSpecs::MUPOINTS]={complex<T>(0,0)};
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_mu_matrix_eval_pts(m0mass_matrix);
	for(size_t m0pole=0;m0pole<RatBoxSpecs::MUPOINTS;m0pole++){
		// Clear any previously stored results
		complex<T> boxcoeffs_mass[2]={complex<T>(0,0)};
		size_t ref_mass=mc.insert(momentum<complex<T> >(m0mass[m0pole],sqrt(m0mass[m0pole]),complex<T>(0,0),complex<T>(0,0)),_mt_massive);

		*(indlst[0].end()-1)=ref_mass;
		*(indlst[1].end()-1)=ref_mass;
		*(indlst[2].end()-1)=ref_mass;
		*(indlst[3].end()-1)=ref_mass;

		for(size_t npole=0; npole<2; npole++){
			*(indlst[3].end()-3)=mc.insert(mompt1+pp.boxpoles[npole+2*m0pole]*pp.vec1c+(pp.alp1*pp.alp2-m0mass[m0pole]/pp.gamma)/pp.boxpoles[npole+2*m0pole]*pp.vec2c,_mt_massive);
			*(indlst[0].end()-3)=mc.insert(mc.mom(*(indlst[3].end()-3))-pp.K1,_mt_massive);
			*(indlst[1].end()-3)=mc.insert(mc.mom(*(indlst[0].end()-3))-pp.K4,_mt_massive);
			*(indlst[2].end()-3)=mc.insert(mc.mom(*(indlst[3].end()-3))-pp.K2,_mt_massive);
			*(indlst[0].begin())=mc.insert(-mc.mom(*(indlst[3].end()-3)),_mt_massive);
			*(indlst[1].begin())=mc.insert(-mc.mom(*(indlst[0].end()-3)),_mt_massive);
			*(indlst[2].begin())=mc.insert(-mc.mom(*(indlst[1].end()-3)),_mt_massive);
			*(indlst[3].begin())=mc.insert(-mc.mom(*(indlst[2].end()-3)),_mt_massive);

			//Construct the three-particle cut subtraction at the above momenta and compute the amplitude
			//If there is only one solution we want to store this in the correct location
			complex<T> temp_res(0,0);
			//Construct the two-particle cut at the momenta above
			for(int i=0;i<BASE::decendant_nbr();i++){
				temp_res+=BASE::eval_tree(i,0,mc,indlst[0])*BASE::eval_tree(i,1,mc,indlst[1])*BASE::eval_tree(i,2,mc,indlst[2])*BASE::eval_tree(i,3,mc,indlst[3]);
			}
			boxcoeffs_mass[npole]+=temp_res; // We only use the first two m0pole results in the triangle subtraction
			for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS-2;imupoints++){
				test[npole+imupoints*2]+=m0mass_matrix[m0pole+(imupoints+2)*RatBoxSpecs::MUPOINTS]*temp_res;
			}
		}

		//Store information allowing us to decide how to choose how to map t_+ and t_- to d_+ and d_- later on.
		// We do this by borrowing the method hidden in the OPP approach basically we can rewrite the OPP expresion
		//  for the subtracted box terms as
		//  (1/2)((d_+(1+T(l'_+)/T(l_+))+d_-(1-T(l'_+)/T(l_+)))/(t-t_+)+(d_+(1+T(l'_-)/T(l_+))+d_-(1-T(l'_-)/T(l_+)))/(t-t_-))
		//  with T(l)=Tr(l,K1f,K2f,K4) and l the original solution of l to find d_+ and d_- and l' the new solution for
		//  l we are trying to map d_+ and d_- to.
		//  The ratio of T(l'_+)/T(l_+) will be either +1 or -1 depending upon the order of the momenta entering l
		//  and the sign will correctly choose the mapping of t_i->d_i.
		set_Kif_coeffs(m0pole,pp.K1flatc,pp.K2flatc,pp.K4);
		//Flip the sign as this is T(l_-) and we store T(l_+)=-T(l_-)
		// Here *(indlst[3].end()-3) points to l of m0point[RATPOINTS-1]
		set_denfac(m0pole,-(pp.K2flatc.Lt()*smatrix<T>(pp.K4)*smatrix<T>(mc.mom(*(indlst[3].end()-3)))*pp.K1flatc.Lt())*(pp.K1flatc.L()*pp.K2flatc.L())
				+(pp.K2flatc.L()*smatrix<T>(pp.K4)*smatrix<T>(mc.mom(*(indlst[3].end()-3)))*pp.K1flatc.L())*(pp.K1flatc.Lt()*pp.K2flatc.Lt()));

		//Store the coefficients
		set_coeff(0,m0pole,boxcoeffs_mass[0]);
		set_coeff(1,m0pole,boxcoeffs_mass[1]);
	}
	complex<T> fullboxcoeff(0,0);
	for(int imupoints=0;imupoints<RatBoxSpecs::MUPOINTS-2;imupoints++){
		fullboxcoeff+=(test[0+2*imupoints]+test[1+2*imupoints]+pp.pent_sub_ret[imupoints])*(static_cast<typename RatBoxSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1);
	}
	set_D0coeff(fullboxcoeff);

	//Mark that we have computed these, we use the mcID passed from the bubble as this will
	// have been computed using a sub_mom_conf
	mcID=original_mcID;
	indID=ind;

	accuracy=pp.accuracy;
}


template <class BASE, class RatBoxSpecs> void box_Rat<BASE,RatBoxSpecs>::init()
{
	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	indlst[0].assign(BASE::corner_size(1)+4,0);
	indlst[1].assign(BASE::corner_size(2)+4,0);
	indlst[2].assign(BASE::corner_size(3)+4,0);
	indlst[3].assign(BASE::corner_size(4)+4,0);

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

template <class BASE, class RatBoxSpecs> void box_Rat<BASE,RatBoxSpecs>::add_mass(const cutD& pcd)
{
	add_mass(pcd.l(1).mass_label(),pcd.l(2).mass_label(),pcd.l(3).mass_label(),pcd.l(4).mass_label());

}

template <class BASE, class RatBoxSpecs> void box_Rat<BASE,RatBoxSpecs>::add_mass(int m1,int m2,int m3,int m4)
{
	if(find(_cut_mass.begin(),_cut_mass.end(),m1)==_cut_mass.end()){_cut_mass.push_back(m1);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m2)==_cut_mass.end()){_cut_mass.push_back(m2);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m3)==_cut_mass.end()){_cut_mass.push_back(m3);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m4)==_cut_mass.end()){_cut_mass.push_back(m4);};

	//Set up the masses for the legs
	_leg_masses->set(0,m1);
	_leg_masses->set(1,m2);
	_leg_masses->set(2,m3);
	_leg_masses->set(3,m4);
}


template <class BASE, class RatBoxSpecs> box_Rat<BASE,RatBoxSpecs>::~box_Rat()
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
#if BH_USE_GMP

	delete _ep_GMP[0];
	delete _ep_GMP[1];
	delete _ep_GMP[2];
	delete _ep_GMP[3];

#endif

	delete _leg_masses;
}




}
}


#endif /* BOX_RATEXT_HPP_ */
