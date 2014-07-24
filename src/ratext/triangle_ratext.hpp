/*
 * triangle_ratext.hpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef TRIANGLE_RATEXT_HPP_
#define TRIANGLE_RATEXT_HPP_

#include "ratext/triangle_ratext.h"
#include "ratext/rat_ext.h"
#include "ratext/bubble_ratext.h"
#include <cassert>

namespace BH {

namespace ratext {

template <class BASE,  class RatTriSpecs> void triangle_Rat<BASE,RatTriSpecs>::init()
{
	mcID=-1;// We've not computed anything yet so set this to a value that mc.get_ID() will not return
	epID=~1; // We set this to the largest unsigned long int (using a unary not operator) so that we do not accidentally think we have evaluated anything

	//Set up the inds for the cut legs
	size_t mm;
	indlst[0].assign(BASE::corner_size(1)+4,0);
	indlst[1].assign(BASE::corner_size(2)+4,0);
	indlst[2].assign(BASE::corner_size(3)+4,0);

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
	_ep_GMP[0]=new eval_param<RGMP>(BASE::corner_size(1)+2);
	_ep_GMP[1]=new eval_param<RGMP>(BASE::corner_size(2)+2);
	_ep_GMP[2]=new eval_param<RGMP>(BASE::corner_size(3)+2);
#endif

}

template <class BASE,  class RatTriSpecs> void triangle_Rat<BASE,RatTriSpecs>::final_initialisation(){
	// Add another base cut base cut
	for(int nmbr=0;nmbr<BASE::daughters_nbr();nmbr++){
		for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
			coeffkeep[m0pole].push_back(C(0,0));
			coeffkeep[m0pole].push_back(C(0,0));
			polekeep[m0pole].push_back(C(0,0));
			polekeep[m0pole].push_back(C(0,0));
			coeffkeep_HP[m0pole].push_back(CHP(0,0));
			coeffkeep_HP[m0pole].push_back(CHP(0,0));
			polekeep_HP[m0pole].push_back(CHP(0,0));
			polekeep_HP[m0pole].push_back(CHP(0,0));
			coeffkeep_VHP[m0pole].push_back(CVHP(0,0));
			coeffkeep_VHP[m0pole].push_back(CVHP(0,0));
			polekeep_VHP[m0pole].push_back(CVHP(0,0));
			polekeep_VHP[m0pole].push_back(CVHP(0,0));
#if BH_USE_GMP
			coeffkeep_GMP[m0pole].push_back(CGMP(0,0));
			coeffkeep_GMP[m0pole].push_back(CGMP(0,0));
			polekeep_GMP[m0pole].push_back(CGMP(0,0));
			polekeep_GMP[m0pole].push_back(CGMP(0,0));
#endif
		}
		denpole.push_back(C(0,0));
		denpole_HP.push_back(CHP(0,0));
		denpole_VHP.push_back(CVHP(0,0));
#if BH_USE_GMP
		denpole_GMP.push_back(CGMP(0,0));
#endif
	}
}

template <class BASE,  class RatTriSpecs> void triangle_Rat<BASE,RatTriSpecs>::add_mass(const cutD& pcd)
{
	add_mass(pcd.l(1).mass_label(),pcd.l(2).mass_label(),pcd.l(3).mass_label());
}

template <class BASE,  class RatTriSpecs> void triangle_Rat<BASE,RatTriSpecs>::add_mass(int m1,int m2,int m3)
{
	if(find(_cut_mass.begin(),_cut_mass.end(),m1)==_cut_mass.end()){_cut_mass.push_back(m1);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m2)==_cut_mass.end()){_cut_mass.push_back(m2);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m3)==_cut_mass.end()){_cut_mass.push_back(m3);};

	//Set up the masses for the legs
	_leg_masses->set(0,m1);
	_leg_masses->set(1,m2);
	_leg_masses->set(2,m3);
}


template <class BASE,  class RatTriSpecs> triangle_Rat<BASE,RatTriSpecs>::~triangle_Rat()
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
#if BH_USE_GMP
	delete _ep_GMP[0];
	delete _ep_GMP[1];
	delete _ep_GMP[2];
#endif
	 delete _leg_masses;
}



template <class BASE,  class RatTriSpecs> C triangle_Rat<BASE,RatTriSpecs>::eval(mom_conf& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("basetriangleRat_comp : Attempting to compute a triangle before computing the bubble");
	}
	return this->template get_symmetry_factor<R>()*C0coeff;
}

template <class BASE,  class RatTriSpecs> CHP triangle_Rat<BASE,RatTriSpecs>::eval(mom_conf_HP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("basetriangleRat_comp : Attempting to compute a triangle before computing the bubble");
	}
	return this->template get_symmetry_factor<RHP>()*C0coeff_HP;
}

template <class BASE,  class RatTriSpecs> CVHP triangle_Rat<BASE,RatTriSpecs>::eval(mom_conf_VHP& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("basetriangleRat_comp : Attempting to compute a triangle before computing the bubble");
	}
	return this->template get_symmetry_factor<RVHP>()*C0coeff_VHP;
}


template <class BASE,  class RatTriSpecs> C triangle_Rat<BASE,RatTriSpecs>::eval(const eval_param<R>& ep)
{
//	// If we have not previously computed the box compute it
//	if(ep.get_ID()!=mcID){
//		coeffparam_mass<R> tp;
//		get_coeffs(mc,ind,tp);
//	}
	return this->template get_symmetry_factor<R>()*C0coeff;
}

template <class BASE,  class RatTriSpecs> CHP triangle_Rat<BASE,RatTriSpecs>::eval(const eval_param<RHP>& ep)
{
	// If we have not previously computed the box compute it
//	if(ep.get_ID()!=mcID){
//		coeffparam_mass<RHP> tp;
//		get_coeffs(mc,ind,tp);
//	}
	return this->template get_symmetry_factor<RHP>()*C0coeff_HP;
}

template <class BASE,  class RatTriSpecs> CVHP triangle_Rat<BASE,RatTriSpecs>::eval(const eval_param<RVHP>& ep)
{
	// If we have not previously computed the box compute it
//	if(ep.get_ID()!=mcID){
//		coeffparam_mass<RVHP> tp;
//		get_coeffs(mc,ind,tp);
//	}
	return this->template get_symmetry_factor<RVHP>()*C0coeff_VHP;
}


#if BH_USE_GMP
template <class BASE,  class RatTriSpecs> CGMP triangle_Rat<BASE,RatTriSpecs>::eval(momentum_configuration<RGMP>& mc,const std::vector<int>& ind)
{
	// If we have not previously computed the box compute it
	if((mc.get_ID()!=mcID)||(ind!=indID)){
//		_WARNING("basetriangleRat_comp : Attempting to compute a triangle before computing the bubble");
	}
	return this->template get_symmetry_factor<RGMP>()*C0coeff_GMP;
}

template <class BASE,  class RatTriSpecs> CGMP triangle_Rat<BASE,RatTriSpecs>::eval(const eval_param<RGMP>& ep)
{
	return this->template get_symmetry_factor<RGMP>()*C0coeff_GMP;
}

#endif

template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat_plusminus<BASE,RatTriSpecs>::get_sub_terms_work(momentum_configuration<T>& mc, const std::vector<int>& ind, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)
{
	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	 if(tp.tri_corner==2){
		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
			 K2sum_mom-=mc.mom(ind[BASE::corner_ind(3,kiter+1)-1]);
		 }
	 }
	 else
	 {
		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
			 K2sum_mom+=mc.mom(ind[BASE::corner_ind(tp.tri_corner,kiter+1)-1]);
		 }
	 }

	 //Now construct the constants we will use and pick gamma with the
	 // convention that if S1 or S2 would be 0 then it would not vanish
	 // this avoids the need to specialize this statement.
	 complex<T> S2=K2sum_mom.square();
	 complex<T> K1K2=tp.K1*K2sum_mom;

	 if(!this->is_eval(tp.mcID,ind)){
		 this->get_coeffs(mc,ind,tp.accuracy,tp.mcID,tp.vertex_ref);
	 }

	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
	 // triangle triple cut using the analytic formula below.
	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);//mc.spab(tp.chi, K2sum_momt, tp.K1flatb);
	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);

	 complex<T> ampp, ampm, ypfac;
	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermptri, subtermmtri;
	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;

	 // Contains the value of the triangle pole subtraction for the bubble
	 complex<T> ampfull[RatTriSpecs::YPOINTS*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> ampfull_mass[RatTriSpecs::YPOINTS]={complex<T>(0,0)};
	 complex<T> subres_1a(0,0),subres_2a(0,0),subres_1b(0,0),subres_2b(0,0);

	 // Get the points to evaluate m0 around
	 const complex<T>* m0mass;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	 const complex<T>* m0mass_matrix_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);

	 // Get the points to evaluate the circle on
	 const complex<T>* eval_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
	 const complex<T>* ypoint;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_eval_pts(ypoint);

	 // Get the coefficients and other factors from the previously computed triangle, if not
	 //  this will just compute them
	 complex<T> gamma_old,gammap_old,*orig_coeffs;
	 this->coeffkeep_get(orig_coeffs);
	 this->get_tri_param_gamma(gamma_old,gammap_old);
	 momentum<complex<T> > v1old, v2old;
	 this->get_tri_param_basis_vectors(v1old,v2old);

	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
	 this->get_boxes(triboxcoeff,triboxpoles,denfac);

	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};

	 complex<T> temp[RatTriSpecs::YPOINTS*RatTriSpecs::MUBUBPOINTS]={complex<T>(0,0)};
	 complex<T> temp_mass[RatTriSpecs::YPOINTS*RatTriSpecs::MUBUBPOINTS]={complex<T>(0,0)};
	 complex<T> invcircpos,bubts,bubtis,subterm;

	 for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
		 int indexbase=m0pole*RatTriSpecs::CPOINTS;

		 // We are subtracting the mu^2(t+c+1/t) contribution and we only want the c piece
		 for(int icirc=0; icirc<RatTriSpecs::CBUBPOINTS; icirc++){
			 //Compute the y+/- poles
			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat/T(2);
			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
			 resb=sqrt(Delta3)*sqrt(ressqrt);

			 Nysolp=(resa+resb)/chiK2Klflat;
			 Nysolm=(resa-resb)/chiK2Klflat;

			 // Construct the new t coefficient
			 bubtsp=-T(2)*(Nysolp*fac1[0]
			              +(T(1)-Nysolp)*fac1[1]
		 		    	  +eval_pts[icirc]*fac1[2]
		 		    	  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 bubtsm=-T(2)*(Nysolm*fac1[0]
				          +(T(1)-Nysolm)*fac1[1]
				          +eval_pts[icirc]*fac1[2]
				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

//			 /************************************************************************/
//			 bubtism=T(2)*(Nysolm*fac2[0]
//				 		  +(T(1)-Nysolm)*fac2[1]
//				 		  +eval_pts[icirc]*fac2[2]
//				 	 	  +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/m0mass[m0pole];
//
//			 bubtisp=T(2)*(Nysolp*fac2[0]
//				 		  +(T(1)-Nysolp)*fac2[1]
//				 		  +eval_pts[icirc]*fac2[2]
//				 		  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/m0mass[m0pole];
//
//			 subtermmtri=orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2];
//			 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//				 subtermmtri+=pow(bubtsm,(RatTriSpecs::CPOINTS-1)/2-icoeffs)*orig_coeffs[indexbase+icoeffs];
//			 }
//			 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//				 subtermmtri+=pow(bubtism,icoeffs+1)*orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2+icoeffs+1];
//			 }
//			 subtermmtri*=complex<T>(0,-1);
//
//			 subtermptri=orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2];
//			 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//				 subtermptri+=pow(bubtsp,(RatTriSpecs::CPOINTS-1)/2-icoeffs)*orig_coeffs[indexbase+icoeffs];
//			 }
//			 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//				 subtermptri+=pow(bubtisp,icoeffs+1)*orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2+icoeffs+1];
//			 }
//			 subtermptri*=complex<T>(0,-1);
////			 /************************************************************************/
			 //Now add the box contributions to the triple cut
			 complex<T> subtermp(0,0), subtermm(0,0);
			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
			 }

			 // Now construct the actual subtraction term
			 ypfac=eval_pts[icirc]/(T(2)*resb);
			 for(int iy=0;iy<RatTriSpecs::YPOINTS;iy++){
				 // Now construct the t in the bubble parameterisation
				 invcircpos=(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)/eval_pts[icirc];
				 bubts=-T(2)*(ypoint[iy]*fac1[0]
				              +(T(1)-ypoint[iy])*fac1[1]
							  +eval_pts[icirc]*fac1[2]
							  +invcircpos*fac1[3])/gamma_old;

				 bubtis=T(2)*(ypoint[iy]*fac2[0]
							  +(T(1)-ypoint[iy])*fac2[1]
							  +eval_pts[icirc]*fac2[2]
							  +invcircpos*fac2[3])/m0mass[m0pole];

				 // Reconstruct the subtracted numerator
				 subterm=orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2];
				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
					 subterm+=pow(bubts,(RatTriSpecs::CPOINTS-1)/2-icoeffs)*orig_coeffs[indexbase+icoeffs];
				 }
				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
					 subterm+=pow(bubtis,icoeffs+1)*orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2+icoeffs+1];
				 }
				 subterm*=complex<T>(0,-1);

				 // Now construct the actual subtraction term we add and subtract the triangle poles and so we end up just adding the box
				 //  pole terms at the triangle triple cut to the points around the circle we are using to compute the contribution of the triangle
				 //  to the bubble at infinity (we do not include the box terms at the box pole in y or those from the eval around the points for
				 //  the extra contribution also all these box terms (including the ones we use) cancel each other (as there is no boundary))
				 momentum<complex<T> > l2=ypoint[iy]*tp.K1flatbc.P()+(T(1)-ypoint[iy])*tp.chic.P()+eval_pts[icirc]*k1bchi+(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)*chik1b/eval_pts[icirc]-K2sum_mom;
				 complex<T> part=subterm/(l2.square()-m0mass[m0pole])+(subtermp/(ypoint[iy]-Nysolp)-subtermm/(ypoint[iy]-Nysolm))*ypfac;

				 complex<T> amp_at=((subtermp+subtermptri)/(ypoint[iy]-Nysolp)-(subtermm+subtermmtri)/(ypoint[iy]-Nysolm))*ypfac;
//				 ampfull_mass[iy]+=amp_at;
				 temp_mass[iy]+=part;
				 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
//					 ampfull[iy+imupoints*RatTriSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*amp_at;
					 temp[iy+imupoints*RatTriSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*part;
				 }
			 }
		 }

//		 subres_1a+=orig_coeffs[indexbase];
//		 subres_1b+=orig_coeffs[indexbase+6]/pow(m0mass[m0pole],3);
//		 subres_2b+=orig_coeffs[indexbase+4]/pow(m0mass[m0pole],2);
//		 subres_2a+=orig_coeffs[indexbase+2]/m0mass[m0pole];
	 }
//	 complex<T> subfacb=-(v1old*chik1b)*T(2)/chiK2Klflat;
//	 complex<T> subfac=(v2old*chik1b)*T(2)/(gamma_old*chiK2Klflat);
//	 complex<T> submufac=-(pow(K1K2-tp.S1*S2/tp.S1,2)-(tp.S1*S2-pow(K1K2,2))/T(3))/tp.S1;
//	 complex<T> subres=((subres_1a*pow(subfac,3)+subres_1b*pow(subfacb,3))*submufac+subres_2a*subfac+subres_2b*subfacb);

	 const complex<T> *y_matrix_point;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_matrix_eval_pts(y_matrix_point);
//	 complex<T> mass_add[RatTriSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
	 complex<T> full_mass_add[RatTriSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
	 for(int ext_int=0;ext_int<RatTriSpecs::MUBUBPOINTS-1;ext_int++){
		 for(int ext=0;ext<RatTriSpecs::YPOINTS;ext++){
			 // Compute the y^2 contribution to the result
//		 	 mass_add[ext_int]+=y_matrix_point[ext+(2+ext_int)*RatTriSpecs::YPOINTS]*ampfull_mass[ext];
			 full_mass_add[ext_int]+=y_matrix_point[ext+(2+ext_int)*RatTriSpecs::YPOINTS]*temp_mass[ext];
		 }
	 }
//	 mass_add[0]/=(-T(18)*tp.S1);
	 full_mass_add[0]*=(T(1)/(-T(3)*tp.S1))/T(6);

	 //To estimate the computation error return the extra components
	 for(int ret=0;ret<RatTriSpecs::YPOINTS;ret++){
//		 tp.bub_acc_pieces[ret]+=ampfull[ret];
		 tp.bub_acc_pieces[ret]+=temp[ret];
	 }

	 //Now return the complete subtraction result
	 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
//		 tp.tri_sub_pieces[imupoints]+=(ampfull[imupoints*RatTriSpecs::YPOINTS]/T(RatTriSpecs::YPOINTS)+mass_add[0]+subres/(complex<T>(0,1)*T(RatTriSpecs::MUBUBPOINTS)));
		 tp.tri_sub_pieces[imupoints]+=temp[imupoints*RatTriSpecs::YPOINTS]/T(RatTriSpecs::YPOINTS)+full_mass_add[0];
	 }
}



template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat_3mass<BASE,RatTriSpecs>::get_sub_terms_work(momentum_configuration<T>& mc, const std::vector<int>& ind, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)
{
	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	 if(tp.tri_corner==2){
		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
			 K2sum_mom-=mc.mom(ind[BASE::corner_ind(3,kiter+1)-1]);
		 }
	 }
	 else
	 {
		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
			 K2sum_mom+=mc.mom(ind[BASE::corner_ind(tp.tri_corner,kiter+1)-1]);
		 }
	 }

	 //Now construct the constants we will use and pick gamma with the
	 // convention that if S1 or S2 would be 0 then it would not vanish
	 // this avoids the need to specialise this statement.
	 complex<T> S2=K2sum_mom.square();
	 complex<T> K1K2=tp.K1*K2sum_mom;

	 //For gamma we need to choose this such that alpha1 and alpha2 are not small and
	 // ideally of ~1. We can do this by checking that we choose the correct sign in gamma
	 // when S2~0 and similarly when S3~0 that gamma gamma!=S1 or S2 depending on the solution

//	 // Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
//	 tp.gammap=(K1K2+sqrt(pow(K1K2,2)-tp.S1*tp.S2));

	 if(!this->is_eval(tp.mcID,ind)){
		 this->get_coeffs(mc,ind,tp.accuracy,tp.mcID,tp.vertex_ref);
	 }

	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
	 // triangle triple cut using the analytic formula below.
	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);//mc.spab(tp.chi, K2sum_momt, tp.K1flatb);
	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);

	 complex<T> ampp, ampm, ypfac;
	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermptri, subtermmtri;
	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;

	 // Contains the value of the triangle pole subtraction for the bubble
	 complex<T> ampfull[RatTriSpecs::YPOINTS*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> ampfull_mass[RatTriSpecs::YPOINTS]={complex<T>(0,0)};
	 complex<T> subres_1a(0,0),subres_1b(0,0),subres_2a(0,0),subres_2b(0,0);

	 // Get the points to evaluate m0 around
	 const complex<T>* m0mass;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	 const complex<T>* m0mass_matrix_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);

	 // Get the points to evaluate the circle on
	 const complex<T>* eval_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
	 const complex<T>* ypoint;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_eval_pts(ypoint);

	 // Get the coefficients and other factors from the previously computed triangle, if not
	 //  this will just compute them
	 complex<T> gamma_old,gammap_old,*orig_coeffs;
	 this->coeffkeep_get(orig_coeffs);
	 this->get_tri_param_gamma(gamma_old,gammap_old);
	 momentum<complex<T> > v1old, v2old;
	 this->get_tri_param_basis_vectors(v1old,v2old);

	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
	 this->get_boxes(triboxcoeff,triboxpoles,denfac);

	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};

	 complex<T> temp[RatTriSpecs::YPOINTS*RatTriSpecs::MUBUBPOINTS]={complex<T>(0,0)};
	 complex<T> temp_mass[RatTriSpecs::YPOINTS*RatTriSpecs::MUBUBPOINTS]={complex<T>(0,0)};
	 complex<T> invcircpos,bubts,bubtis,subterm;

	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	 for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
		 int indexbase=m0pole*RatTriSpecs::CPOINTS;

		 // We are subtracting the mu^2(t+c+1/t) contribution and we only want the c piece
		 for(int icirc=0; icirc<RatTriSpecs::CBUBPOINTS; icirc++){
			 //Compute the y+/- poles
			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat*T(0.5);
			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
			 resb=sqrt(Delta3)*sqrt(ressqrt);

			 Nysolp=(resa+resb)/chiK2Klflat;
			 Nysolm=(resa-resb)/chiK2Klflat;

			 // Construct the new t coefficient
			 bubtsp=-T(2)*(Nysolp*fac1[0]
			                 +(T(1)-Nysolp)*fac1[1]
	 		 		    	 +eval_pts[icirc]*fac1[2]
	 		 		    	 +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 bubtsm=-T(2)*(Nysolm*fac1[0]
	 				          +(T(1)-Nysolm)*fac1[1]
	 				          +eval_pts[icirc]*fac1[2]
	 				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 complex<T> subtermp(0,0), subtermm(0,0);

			 //Now add the box contributions to the triple cut
			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
				 complex<T> reddenfac=T(1)/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))*reddenfac;
				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))*reddenfac;
			 }

			 // Now construct the actual subtraction term
			 ypfac=eval_pts[icirc]/(T(2)*resb);
			 for(int iy=0;iy<RatTriSpecs::YPOINTS;iy++){
//				 // Now construct the t in the bubble parameterisation
//				 invcircpos=((ypoint[iy]-a)*(T(1)-(ypoint[iy]-a))-m0mass[m0pole]/tp.gammab)/eval_pts[icirc];
//				 bubts=-T(2)*((ypoint[iy]-a)*fac1[0]
//				              +(T(1)-(ypoint[iy]-a))*fac1[1]
//							  +eval_pts[icirc]*fac1[2]
//							  +invcircpos*fac1[3])/gamma_old;
//
//				 bubtis=T(2)*((ypoint[iy]-a)*fac2[0]
//							  +(T(1)-(ypoint[iy]-a))*fac2[1]
//							  +eval_pts[icirc]*fac2[2]
//							  +invcircpos*fac2[3])/m0mass[m0pole];
//
//				 // Reconstruct the subtracted numerator
//				 subterm=orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2];
//				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//					 subterm+=pow(bubts,(RatTriSpecs::CPOINTS-1)/2-icoeffs)*orig_coeffs[indexbase+icoeffs];
//				 }
//				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
//					 subterm+=pow(bubtis,icoeffs+1)*orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2+icoeffs+1];
//				 }
//				 subterm*=complex<T>(0,-1);
//
//				 // Now construct the actual subtraction term we add and subtract the triangle poles and so we end up just adding the box
//				 //  pole terms at the triangle triple cut to the points around the circle we are using to compute the contribution of the triangle
//				 //  to the bubble at infinity (we do not include the box terms at the box pole in y or those from the eval around the points for
//				 //  the extra contribution also all these box terms (including the ones we use) cancel each other (as there is no boundary))
//				 momentum<complex<T> > l2=(ypoint[iy]-a)*tp.K1flatbc.P()+(T(1)-(ypoint[iy]-a))*tp.chic.P()+eval_pts[icirc]*k1bchi+((ypoint[iy]-a)*(T(1)-(ypoint[iy]-a))-m0mass[m0pole]/tp.gammab)*chik1b/eval_pts[icirc]-K2sum_mom;
//				 complex<T> part=subterm/(l2.square()-m0mass[m0pole])+(subtermp/(ypoint[iy]-Nysolp)-subtermm/(ypoint[iy]-Nysolm))*ypfac;
//
				 complex<T> amp_at=((subtermp+subtermptri)/(ypoint[iy]-Nysolp)-(subtermm+subtermmtri)/(ypoint[iy]-Nysolm))*ypfac;
				 ampfull_mass[iy]+=amp_at;
//				 temp_mass[iy]+=part;
				 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
					 ampfull[iy+imupoints*RatTriSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*amp_at;
//					 temp[iy+imupoints*RatTriSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*part;
				 }
			 }

//			 ypfac=eval_pts[icirc]/(T(2)*resb);
//			 for(int yloop=0;yloop<RatTriSpecs::YPOINTS;yloop++){
//				 complex<T> amp_at=(subtermp/(ypoint[yloop]-Nysolp)-subtermm/(ypoint[yloop]-Nysolm))*ypfac;
//				 ampfull_mass[yloop]+=amp_at;
//				 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
//					 ampfull[yloop+imupoints*RatTriSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*amp_at;
//				 }
//			 }
		 }

		 //Now subtract off the triple-cut contributions
		 subres_1a+=orig_coeffs[indexbase];
		 subres_1b+=orig_coeffs[indexbase+6]/pow((gamma_old-m0mass[m0pole]),3);
		 subres_2b+=orig_coeffs[indexbase+4]/(m0mass[m0pole]*(gamma_old-m0mass[m0pole]));
		 subres_2a+=orig_coeffs[indexbase+2]/m0mass[m0pole];
	 }

	 complex<T> subfacb=v1old*chik1b*T(2)/chiK2Klflat;
	 complex<T> subfac=v2old*chik1b*T(2)/(gamma_old*chiK2Klflat);
	 complex<T> gamma_fac=(K1K2+sqrt(pow(K1K2,2)-tp.S1*S2));//This is gammap for this triangle subtraction we need it for the line below
	 complex<T> submufac=-(pow(K1K2-tp.S1*S2/gamma_fac,2)-(tp.S1*S2-pow(K1K2,2))/T(3))/tp.S1;
	 complex<T> subres=((subres_1a*pow(subfac,3)+subres_1b*pow(subfacb,3))*submufac+subres_2a*subfac+subres_2b*subfacb);

	 //Extract the different components we need
	 const complex<T> *y_matrix_point;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_matrix_eval_pts(y_matrix_point);
	 complex<T> mass_add[RatTriSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
//	 complex<T> full_mass_add[RatTriSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
	 for(int ext_int=0;ext_int<RatTriSpecs::MUBUBPOINTS-1;ext_int++){
		 for(int ext=0;ext<RatTriSpecs::YPOINTS;ext++){
			 // Compute the y^2 contribution to the result
		 	 mass_add[ext_int]+=y_matrix_point[ext+(2+ext_int)*RatTriSpecs::YPOINTS]*ampfull_mass[ext];
//			 full_mass_add[ext_int]+=y_matrix_point[ext+(2+ext_int)*RatTriSpecs::YPOINTS]*temp_mass[ext];
		 }
	 }
	 mass_add[0]/=(-T(18)*tp.S1);
//	 full_mass_add[0]*=(T(1)/(-T(3)*tp.S1))/T(6);

	 //To estimate the computation error return the extra components
	 for(int ret=0;ret<RatTriSpecs::YPOINTS;ret++){
		 tp.bub_acc_pieces[ret]+=ampfull[ret];
//		 tp.bub_acc_pieces[ret]+=temp[ret];
	 }

	 //Now return the complete subtraction result
	 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
		 tp.tri_sub_pieces[imupoints]+=(ampfull[imupoints*RatTriSpecs::YPOINTS]/T(RatTriSpecs::YPOINTS)+mass_add[0]+subres/(complex<T>(0,1)*T(RatTriSpecs::MUBUBPOINTS)));
//		 tp.tri_sub_pieces[imupoints]+=temp[imupoints*RatTriSpecs::YPOINTS]/T(RatTriSpecs::YPOINTS)+full_mass_add[0];
	 }
}


template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat_plusminus<BASE,RatTriSpecs>::get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)
{
	//_MESSAGE4("Sub +/- bub:",BASE::corner_size(1),BASE::corner_size(2),BASE::corner_size(3));
	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	 if(tp.tri_corner==2){
		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
			 K2sum_mom-=ep.p(BASE::corner_ind(3,kiter+1)-1)->P();
		 }
	 }
	 else
	 {
		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
			 K2sum_mom+=ep.p(BASE::corner_ind(tp.tri_corner,kiter+1)-1)->P();
		 }
	 }

	 //Now construct the constants we will use and pick gamma with the
	 // convention that if S1 or S2 would be 0 then it would not vanish
	 // this avoids the need to specialise this statement.
	 complex<T> S2=K2sum_mom.square();
	 complex<T> K1K2=tp.K1*K2sum_mom;

	 //For gamma we need to choose this such that alpha1 and alpha2 are not small and
	 // ideally of ~1. We can do this by checking that we choose the correct sign in gamma
	 // when S2~0 and similarly when S3~0 that gamma gamma!=S1 or S2 depending on the solution

	 if(!this->is_eval(ep.get_ID())){
		 this->get_coeffs(ep,tp.accuracy);
	 }

	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
	 // triangle triple cut using the analytic formula below.
	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);
	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);

	 complex<T> ampp, ampm, ypfac;
	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermptri, subtermmtri;
	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;

	 // Get the points to evaluate m0 around
	 const complex<T>* m0mass;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	 const complex<T>* m0mass_matrix_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);

	 // Get the points to evaluate the circle on
	 const complex<T>* eval_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
	 const complex<T>* ypoint;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_eval_pts(ypoint);

	 // Get the coefficients and other factors from the previously computed triangle, if not
	 //  this will just compute them
	 complex<T> gamma_old,gammap_old,*orig_coeffs;
	 this->coeffkeep_get(orig_coeffs);
	 this->get_tri_param_gamma(gamma_old,gammap_old);
	 momentum<complex<T> > v1old, v2old;
	 this->get_tri_param_basis_vectors(v1old,v2old);

	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
	 this->get_boxes(triboxcoeff,triboxpoles,denfac);

	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};

	 complex<T> temp[(RatTriSpecs::YPOINTS-2)*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> temp_mass[(RatTriSpecs::YPOINTS-2)*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> invcircpos,bubts,bubtis,subterm,massform;
		 
	 const complex<T> *y_matrix_point;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_matrix_eval_pts(y_matrix_point);

	 for(int m0pole=0;m0pole<RatTriSpecs::MUBUBPOINTS;m0pole++){
		 int indexbase=m0pole*RatTriSpecs::CPOINTS;

		 // We are subtracting the mu^2(t+c+1/t) contribution and we only want the c piece
		 for(int icirc=0; icirc<RatTriSpecs::CBUBPOINTS; icirc++){
			 //Compute the y+/- poles
			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat/T(2);
			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
			 resb=sqrt(Delta3*ressqrt);

			 Nysolp=(resa+resb)/chiK2Klflat;
			 Nysolm=(resa-resb)/chiK2Klflat;

			 
			 complex<T> bubtsp=-T(2)*(Nysolp*fac1[0]
			              +(T(1)-Nysolp)*fac1[1]
		 		    	  +eval_pts[icirc]*fac1[2]
		 		    	  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 complex<T> bubtsm=-T(2)*(Nysolm*fac1[0]
				          +(T(1)-Nysolm)*fac1[1]
				          +eval_pts[icirc]*fac1[2]
				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 //Now add the box contributions to the triple cut
			 complex<T> subtermp(0,0), subtermm(0,0);
			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
			 }

			 // Now construct the actual subtraction term
			 ypfac=eval_pts[icirc]/(T(2)*resb);
			 for(int iy=0;iy<RatTriSpecs::YPOINTS;iy++){
				 invcircpos=(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)/eval_pts[icirc];

				 complex<T> bubts=-T(2)*(ypoint[iy]*fac1[0]
				              +(T(1)-ypoint[iy])*fac1[1]
							  +eval_pts[icirc]*fac1[2]
							  +invcircpos*fac1[3])/gamma_old;

				 complex<T> bubtis=T(2)*(ypoint[iy]*fac2[0]
							  +(T(1)-ypoint[iy])*fac2[1]
							  +eval_pts[icirc]*fac2[2]
							  +invcircpos*fac2[3])/m0mass[m0pole];

				 complex<T> subterm_pos(0,0),subterm_neg(0,0);
				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
					 subterm_pos=bubts*(subterm_pos+orig_coeffs[indexbase+icoeffs]);
					 subterm_neg=bubtis*(subterm_neg+orig_coeffs[indexbase+RatTriSpecs::CPOINTS-1-icoeffs]);
				 }
				 subterm=complex<T>(0,-1)*(subterm_pos+subterm_neg+orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2]);

				 // Now construct the actual subtraction term we add and subtract the triangle poles and so we end up just adding the box
				 //  pole terms at the triangle triple cut to the points around the circle we are using to compute the contribution of the triangle
				 //  to the bubble at infinity (we do not include the box terms at the box pole in y or those from the eval around the points for
				 //  the extra contribution also all these box terms (including the ones we use) cancel each other (as there is no boundary))
				 momentum<complex<T> > l2=ypoint[iy]*tp.K1flatbc.P()+(T(1)-ypoint[iy])*tp.chic.P()+eval_pts[icirc]*k1bchi+(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)*chik1b/eval_pts[icirc]-K2sum_mom;
				 complex<T> part=subterm/(l2.square()-m0mass[m0pole])+(subtermp/(ypoint[iy]-Nysolp)-subtermm/(ypoint[iy]-Nysolm))*ypfac;

				 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
					 // For the piece independent of m^2 we can simply add all the terms together as the conversion matrix for m is just 1
					 //  we need to get all terms from y^2 upwards and m^2 y^2 upwards etc.
					 // TODO: fix this for higher mu powers. Need to rescale by number of mu points I think
					 for(int ylvl=0;ylvl<RatTriSpecs::YPOINTS-2;ylvl++){
						 temp_mass[imupoints+ylvl*(RatTriSpecs::MUBUBPOINTS-1)]+=m0mass_matrix_pts[m0pole+imupoints*RatTriSpecs::MUBUBPOINTS]*y_matrix_point[iy+(2+ylvl)*RatTriSpecs::YPOINTS]*part;
					 }

					 massform=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*part;
//					 // We want the y^2 and y component of these terms
//					 temp_mass[imupoints]+=y_matrix_point[iy+(2+imupoints)*RatTriSpecs::YPOINTS]*massform*(T(1)/T(3)*(T(1)-m0mass[m0pole]/tp.S1)+a+a*a)/T(6);
//					 temp_mass2[imupoints]+=y_matrix_point[iy+(1+imupoints)*RatTriSpecs::YPOINTS]*massform*(T(1)/T(2)+a);
					 // Here we just want the first element in the y array and not y^2 etc so we do not multiply by y_matrix_point
//					 temp[iy+imupoints*RatTriSpecs::YPOINTS]+=massform;
					 // We want to sum over all the m^2(y^0+y^1) terms 
					 for(int ylp=0;ylp<RatTriSpecs::YPOINTS-2;ylp++){
						 temp[imupoints+ylp*(RatTriSpecs::MUBUBPOINTS-1)]+=y_matrix_point[iy+ylp*RatTriSpecs::YPOINTS]*massform/T(ylp+1);
					 }
//					 //Compute the m^2 y contribution to be used in the error analysis
//					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(1+imupoints)*RatTriSpecs::YPOINTS]*massform;
//					 //For the error analysis look at the y^3 component
//					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(3+imupoints)*RatTriSpecs::YPOINTS]*part;
//					 // Compute the m^2 y^2 contribution to be used in the error analysis
					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(2+imupoints)*RatTriSpecs::YPOINTS]*massform;
				 }
			 }
		 }
	 }

	 //Now return the complete subtraction result
	 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
		 for(int iypts=0;iypts<RatTriSpecs::YPOINTS-2;iypts++){
			 tp.tri_sub_pieces[imupoints]+=temp[imupoints+iypts*(RatTriSpecs::MUBUBPOINTS-1)]+temp_mass[imupoints+iypts*(RatTriSpecs::MUBUBPOINTS-1)]*(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_extra_fac(imupoints,tp.S1);
		 }
	 }
}

template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat_3mass<BASE,RatTriSpecs>::get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,RatTriSpecs::MUBUBPOINTS,RatTriSpecs::YPOINTS>& tp)
{
	//_MESSAGE4("Sub 3mass bub:",BASE::corner_size(1),BASE::corner_size(2),BASE::corner_size(3));
	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
	 if(tp.tri_corner==2){
		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
			 K2sum_mom-=ep.p(BASE::corner_ind(3,kiter+1)-1)->P();
		 }
	 }
	 else
	 {
		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
			 K2sum_mom+=ep.p(BASE::corner_ind(tp.tri_corner,kiter+1)-1)->P();
		 }
	 }

	 if(!this->is_eval(ep.get_ID())){
		 this->get_coeffs(ep,tp.accuracy);
	 }

	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
	 // triangle triple cut using the analytic formula below.
	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());

     complex<T> S2=K2sum_mom.square();
     complex<T> K1K2=tp.K1*K2sum_mom;
	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);
	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);

	 complex<T> ampp, ampm, ypfac;
	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermp, subtermm;
	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;

	 // Get the points to evaluate m0 around
	 const complex<T>* m0mass;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	 const complex<T>* m0mass_matrix_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);

	 // Get the points to evaluate the circle on
	 const complex<T>* eval_pts;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
	 const complex<T>* ypoint;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_eval_pts(ypoint);

	 // Get the coefficients and other factors from the previously computed triangle, if not
	 //  this will just compute them
	 complex<T> gamma_old,gammap_old,*orig_coeffs;
	 this->coeffkeep_get(orig_coeffs);
	 this->get_tri_param_gamma(gamma_old,gammap_old);
	 momentum<complex<T> > v1old, v2old;
	 this->get_tri_param_basis_vectors(v1old,v2old);

	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
	 this->get_boxes(triboxcoeff,triboxpoles,denfac);

	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};

	 complex<T> temp[(RatTriSpecs::YPOINTS-2)*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> temp_mass[(RatTriSpecs::YPOINTS-2)*(RatTriSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	 complex<T> invcircpos,bubts,bubtis,subterm,massform;

	 const complex<T> *y_matrix_point;
	 (static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_y_matrix_eval_pts(y_matrix_point);

	 for(int m0pole=0;m0pole<RatTriSpecs::MUBUBPOINTS;m0pole++){
		 int indexbase=m0pole*RatTriSpecs::CPOINTS;

		 // We are subtracting the mu^2(t+c+1/t) contribution and we only want the c piece
		 for(int icirc=0; icirc<RatTriSpecs::CBUBPOINTS; icirc++){

			 //Compute the y+/- poles
			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat*T(0.5);
			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
			 resb=sqrt(Delta3*ressqrt);

			 Nysolp=(resa+resb)/chiK2Klflat;
			 Nysolm=(resa-resb)/chiK2Klflat;

			 // Construct the new t coefficient
			 bubtsp=-T(2)*(Nysolp*fac1[0]
			                 +(T(1)-Nysolp)*fac1[1]
	 		 		    	 +eval_pts[icirc]*fac1[2]
	 		 		    	 +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 bubtsm=-T(2)*(Nysolm*fac1[0]
	 				          +(T(1)-Nysolm)*fac1[1]
	 				          +eval_pts[icirc]*fac1[2]
	 				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;

			 //Now add the box contributions to the triple cut
			 complex<T> subtermp(0,0), subtermm(0,0);
			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
				 complex<T> reddenfac=T(1)/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))*reddenfac;
				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))*reddenfac;
			 }

			 ypfac=eval_pts[icirc]/(T(2)*resb);
			 for(int iy=0;iy<RatTriSpecs::YPOINTS;iy++){
				 // Now construct the t in the bubble parameterisation
				 invcircpos=(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)/eval_pts[icirc];
				 
				 bubts=-T(2)*(ypoint[iy]*fac1[0]
				              +(T(1)-ypoint[iy])*fac1[1]
							  +eval_pts[icirc]*fac1[2]
							  +invcircpos*fac1[3])/gamma_old;

				 bubtis=-T(2)*(ypoint[iy]*fac2[0]
							  +(T(1)-ypoint[iy])*fac2[1]
							  +eval_pts[icirc]*fac2[2]
							  +invcircpos*fac2[3])/(gamma_old-m0mass[m0pole]);

				 complex<T> subterm_pos(0,0),subterm_neg(0,0);
				 for(int icoeffs=0;icoeffs<(RatTriSpecs::CPOINTS-1)/2;icoeffs++){
					 subterm_pos=bubts*(subterm_pos+orig_coeffs[indexbase+icoeffs]);
					 subterm_neg=bubtis*(subterm_neg+orig_coeffs[indexbase+RatTriSpecs::CPOINTS-1-icoeffs]);
				 }
				 subterm=complex<T>(0,-1)*(subterm_pos+subterm_neg+orig_coeffs[indexbase+(RatTriSpecs::CPOINTS-1)/2]);

				 // Now construct the actual subtraction term we add and subtract the triangle poles and so we end up just adding the box
				 //  pole terms at the triangle triple cut to the points around the circle we are using to compute the contribution of the triangle
				 //  to the bubble at infinity (we do not include the box terms at the box pole in y or those from the eval around the points for
				 //  the extra contribution also all these box terms (including the ones we use) cancel each other (as there is no boundary))
				 momentum<complex<T> > l2=ypoint[iy]*tp.K1flatbc.P()+(T(1)-ypoint[iy])*tp.chic.P()+eval_pts[icirc]*k1bchi+(ypoint[iy]*(T(1)-ypoint[iy])-m0mass[m0pole]/tp.gammab)*chik1b/eval_pts[icirc]-K2sum_mom;
				 complex<T> part=subterm/(l2.square()-m0mass[m0pole])+(subtermp/(ypoint[iy]-Nysolp)-subtermm/(ypoint[iy]-Nysolm))*ypfac;

				 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
					 // For the piece independent of m^2 we can simply add all the terms together as the conversion matrix is just 1
					 for(int ylvl=0;ylvl<RatTriSpecs::YPOINTS-2;ylvl++){
						 temp_mass[imupoints+ylvl*(RatTriSpecs::MUBUBPOINTS-1)]+=m0mass_matrix_pts[m0pole+imupoints*RatTriSpecs::MUBUBPOINTS]*y_matrix_point[iy+(2+ylvl)*RatTriSpecs::YPOINTS]*part;
					 }

					 massform=m0mass_matrix_pts[m0pole+(imupoints+1)*RatTriSpecs::MUPOINTS]*part;
					 // We want to sum over all the m^2(y^0+y^1+y^2+...) terms 
					 for(int ylp=0;ylp<RatTriSpecs::YPOINTS-2;ylp++){
						 temp[imupoints+ylp*(RatTriSpecs::MUBUBPOINTS-1)]+=y_matrix_point[iy+ylp*RatTriSpecs::YPOINTS]*massform/T(ylp+1);
					 }
					 // Compute the m^2 y contribution to be used in the error analysis
//					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(1+imupoints)*RatTriSpecs::YPOINTS]*massform;
//					 //For the error analysis look at the y^3 component
//					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(3+imupoints)*RatTriSpecs::YPOINTS]*part;

//					 // Compute the m^2 y^2 contribution to be used in the error analysis
					 tp.bub_acc_pieces[imupoints]+=y_matrix_point[iy+(2+imupoints)*RatTriSpecs::YPOINTS]*massform;
				 }
			 }
		 }
	 }

	 //Now return the complete subtraction result
	 for(int imupoints=0;imupoints<RatTriSpecs::MUBUBPOINTS-1;imupoints++){
		 for(int iypts=0;iypts<RatTriSpecs::YPOINTS-2;iypts++){
			tp.tri_sub_pieces[imupoints]+=temp[imupoints+iypts*(RatTriSpecs::MUBUBPOINTS-1)]+temp_mass[imupoints+iypts*(RatTriSpecs::MUBUBPOINTS-1)]*(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BubPoints.get_extra_fac(imupoints,tp.S1);
		 }
	 }
}



template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat<BASE,RatTriSpecs>::get_coeffs_fn(momentum_configuration<T>& mc,const  std::vector<int>& ind, double& accuracy, const long int original_mcID, const size_t vertex_ref)
{
	box_param<T,RatTriSpecs::MUTRIPOINTS,RatTriSpecs::CPOINTS> bp;
	bp.accuracy=accuracy;
	bp.mcID=original_mcID;
	bp.vertex_ref=vertex_ref;

	bp.K1=mc.mom(ind[BASE::corner_ind(1,1)-1]);
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		bp.K1+=mc.mom(ind[BASE::corner_ind(1,k1size+1)-1]);
	}
	bp.K2=-mc.mom(ind[BASE::corner_ind(3,1)-1]);
	size_t k2size=1;
	for(;k2size<BASE::corner_size(3);k2size++){
		bp.K2-=mc.mom(ind[BASE::corner_ind(3,k2size+1)-1]);
	}
	complex<T> S1=bp.K1.square();
	complex<T> S2=bp.K2.square();

	complex<T> K1K2=bp.K1*bp.K2;

	//Get the number of legs in the third leg
	size_t k3size=BASE::corner_size(2);

	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	// if K3^2=0 then (K1-K2)^2=0 => S1+S2=2*K1K2 => K1K2^2-S1*S2=(1/4)(S1-S2)^2
	//  => 2*tp.gammap=S1+S2+/-(S1-S2) we can have either tp.gammap=S1 or S2, we choose
	//  the solution so that tp.gammap equals whichever of S1 or S2 is not zero,
	//  unless both are positive in which case we choose the S2 solution.
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}


	//Set up the scaling factors for the momenta
	//We need to be careful if K_1 and K_3 are massless then gamap=S2 and not S1
	complex<T> f1, f2;
	if(!_k3massive){
		if(!_k1massive){// gamma=S2
			bp.alp1=complex<T>(0,0);
			bp.alp2=complex<T>(1,0);
			f2=complex<T>(1,0);
			f1=complex<T>(1,0);
		}
		else{ //gamma=S1 or gamma=S1
			gammap=S1;
			bp.alp1=complex<T>(1,0);
			bp.alp2=complex<T>(0,0);
			f2=complex<T>(1,0);
			f1=complex<T>(1,0);
		}
	}
	else{
		if(!_k1massive){
			bp.alp1=complex<T>(0,0);
			f2=complex<T>(1,0);
		}
		else{
			bp.alp1=complex<T>(1,0);
			f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
		}
		if(!_k2massive){
			bp.alp2=complex<T>(0,0);
			f1=complex<T>(1,0);
		}
		else{
			bp.alp2=complex<T>(1,0);
			f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
		}
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	bp.gamma=gammap/(f1*f2);
	complex<T> gamfac=bp.gamma/(pow(gammap,2)-S1*S2);
	bp.K1flatc=Cmom<T>(f2*gamfac*(gammap*bp.K1-S1*bp.K2));
	bp.K2flatc=Cmom<T>(f1*gamfac*(gammap*bp.K2-S2*bp.K1));
	bp.vec1c=PfLLt(bp.K1flatc.Lt(),bp.K2flatc.L());
	bp.vec2c=PfLLt(bp.K2flatc.Lt(),bp.K1flatc.L());

	//Extract the box contributions, the last element in the returned array is the contribution
	// at t=0 if this is a three-mass triangle
	complex<T> zero(0,0);
	for(int ipole=0;ipole<=(RatTriSpecs::CPOINTS+1)*RatTriSpecs::MUTRIPOINTS;ipole++){
		bp.tri_sub_ret[ipole]=zero;
	}
	//Get all the box subtractions adding each new one to the previous ones contained in coeffsret_mass
	for(size_t box=1; box <= BASE::daughters_nbr(); box++){
		bp.box_corner=BASE::get_opened_corner(box);
#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(box));
		assert(dau);
#else
		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(box));
#endif
		dau->get_sub_terms(mc,ind,bp);

		// Store these results for later use in the bubble coeff subtractions if we have any poles
		boxcoeff_add(box-1,bp.boxcoeffs,bp.boxpoles,bp.denfac1);
	}

	// Now we calculate the coefficient from examining the triple cut
	for(size_t mm=1;mm<=BASE::corner_size(1);mm++){
		(indlst[0])[mm]=ind[BASE::corner_ind(1,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(2);mm++){
		(indlst[1])[mm]=ind[BASE::corner_ind(2,mm)-1];
	}

	for(size_t mm=1;mm<=BASE::corner_size(3);mm++){
		(indlst[2])[mm]=ind[BASE::corner_ind(3,mm)-1];
	}

	// Get the points to evaluate at at the correct precision
	const complex<T>* eval_pts;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_t_eval_pts(eval_pts);

	//We use this RatTriSpecs::CPOINTS*RatTriSpecs::CPOINTS double array contains all the coefficients we use to multiply the 7 terms we generate above
	// to produce the RatTriSpecs::CPOINTS coefficients.
	const complex<T>* coeff_pts;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_t_matrix_eval_pts(coeff_pts);

	// Get the points to evaluate m0 around
	const complex<T>* m0mass;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix);

	*(indlst[0].end()-2)=vertex_ref;
	*(indlst[1].end()-2)=vertex_ref;
	*(indlst[2].end()-2)=vertex_ref;

	momentum<complex<T> > mompt1(bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P());
	complex<T> c0res[RatTriSpecs::MUTRIPOINTS-1]={complex<T>(0,0)}, amp[RatTriSpecs::CPOINTS];
	complex<T> triC0res_mass[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
		momentum<complex<T> > mompt3((bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)*bp.vec2c);
		size_t ref_mass=mc.insert(momentum<complex<T> >(m0mass[m0pole],sqrt(m0mass[m0pole]),complex<T>(0,0),complex<T>(0,0)),_mt_massive);
		for(int icirc=0; icirc<RatTriSpecs::CPOINTS; icirc++){
			//Construct the cut-momenta simply by recomputing the triple cut at each pole
			//  this is the slowest method as it does not reuse any information
			*(indlst[2].end()-3)=mc.insert(mompt1+eval_pts[icirc]*bp.vec1c+mompt3/eval_pts[icirc],_mt_massive);
			*(indlst[0].end()-3)=mc.insert(mc.mom(*(indlst[2].end()-3))-bp.K1,_mt_massive);
			*(indlst[1].end()-3)=mc.insert(mc.mom(*(indlst[2].end()-3))-bp.K2,_mt_massive);

			*(indlst[0].begin())=mc.insert(-mc.mom(*(indlst[2].end()-3)),_mt_massive);
			*(indlst[1].begin())=mc.insert(-mc.mom(*(indlst[0].end()-3)),_mt_massive);
			*(indlst[2].begin())=mc.insert(-mc.mom(*(indlst[1].end()-3)),_mt_massive);

			*(indlst[0].end()-1)=ref_mass;
			*(indlst[1].end()-1)=ref_mass;
			*(indlst[2].end()-1)=ref_mass;

			amp[icirc]=complex<T>(0,0);
			//Construct the two-particle cut at the momenta above
			for(int i=0;i<BASE::decendant_nbr();i++){
				amp[icirc]+=BASE::eval_tree(i,0,mc,indlst[0])*BASE::eval_tree(i,1,mc,indlst[1])*BASE::eval_tree(i,2,mc,indlst[2]);
			}
			amp[icirc]-=bp.tri_sub_ret[icirc+(RatTriSpecs::CPOINTS+1)*m0pole];
		}

		//Now produce the coeffs by starting from 0 and adding the above terms with the correct signs
		for(int ifin=0;ifin<7;ifin++){
			int indexC0=ifin+7*m0pole;
			triC0res_mass[indexC0]=zero;
			for(int jfin=0;jfin<RatTriSpecs::CPOINTS;jfin++){
				triC0res_mass[indexC0]+=coeff_pts[ifin*RatTriSpecs::CPOINTS+jfin]*amp[jfin];
			}
		}
		//We compute C0 using 2c0bare-2(t_sub terms)-(t->0_terms)
		for(int iMupoints=0;iMupoints<RatTriSpecs::MUTRIPOINTS-1;iMupoints++){
			c0res[iMupoints]+=m0mass_matrix[m0pole+(iMupoints+1)*RatTriSpecs::MUPOINTS]*(triC0res_mass[3+RatTriSpecs::CPOINTS*m0pole]+bp.tri_sub_ret[RatTriSpecs::CPOINTS+(RatTriSpecs::CPOINTS+1)*m0pole]);
		}
	}

	// Finally set the actual C0 coefficient, this will be different
	//  from coeffs[3] when we have a three mass triangle where we have to add the
	//  "t=0" contribution
	complex<T> fullC0res(0,0);
	momentum<complex<T> > K3=-(bp.K1+bp.K2);
	complex<T> S3=K3.square();
	for(int imupoints=0;imupoints<RatTriSpecs::MUTRIPOINTS-1;imupoints++){
		fullC0res+=c0res[imupoints]*(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1,S2,S3);
	}
	set_C0coeff(fullC0res);

	//Save the vec1 and vec2 and opened corner so we can relate the c's we compute here to those
	// computed in the other cases
	set_tri_param_basis_vectors(bp.vec1c,bp.vec2c);
	set_tri_param_gamma(bp.gamma,gammap);
	coeffkeep_add(triC0res_mass);

	//Mark that we have computed these, we use the mcID passed from the bubble as this will
	// have been computed using a sub_mom_conf
	mcID=original_mcID;
	indID=ind;
	accuracy=bp.accuracy;
}


template <class BASE,  class RatTriSpecs> template <class T> void triangle_Rat<BASE,RatTriSpecs>::get_coeffs_fn(const eval_param<T>& ep, double& accuracy)
{
	box_param<T,RatTriSpecs::MUTRIPOINTS,RatTriSpecs::CPOINTS> bp;
	bp.accuracy=accuracy;

	bp.K1=ep.p(BASE::corner_ind(1,1)-1)->P();
	size_t k1size=1;
	for(;k1size<BASE::corner_size(1);k1size++){
		bp.K1+=ep.p(BASE::corner_ind(1,k1size+1)-1)->P();
	}
	bp.K2=-ep.p(BASE::corner_ind(3,1)-1)->P();
	size_t k2size=1;
	for(;k2size<BASE::corner_size(3);k2size++){
		bp.K2-=ep.p(BASE::corner_ind(3,k2size+1)-1)->P();
	}
		
	complex<T> S1=bp.K1.square();
	complex<T> S2=bp.K2.square();

	complex<T> K1K2=bp.K1*bp.K2;


	//Compute gammap, if we have a massless leg we need to avoid the "zero" solution
	// if K3^2=0 then (K1-K2)^2=0 => S1+S2=2*K1K2 => K1K2^2-S1*S2=(1/4)(S1-S2)^2
	//  => 2*tp.gammap=S1+S2+/-(S1-S2) we can have either tp.gammap=S1 or S2, we choose
	//  the solution so that tp.gammap equals whichever of S1 or S2 is not zero,
	//  unless both are positive in which case we choose the S2 solution.
	complex<T> gammap;
	if(!(_k1massive&&_k2massive)){
		gammap=T(2)*K1K2;
	}else{
		gammap=(K1K2+sqrt(pow(K1K2,2)-S1*S2));
	}

	//Set up the scaling factors for the momenta
	//We need to be careful if K_1 and K_3 are massless then gammap=S2 and not S1
	complex<T> f1, f2;
	if(!_k3massive){
		if(!_k1massive){// gamma=S2
			bp.alp1=complex<T>(0,0);
			bp.alp2=complex<T>(1,0);
			f2=complex<T>(1,0);
			f1=complex<T>(1,0);
		}
		else{ //gamma=S1
			gammap=S1;
			bp.alp1=complex<T>(1,0);
			bp.alp2=complex<T>(0,0);
			f2=complex<T>(1,0);
			f1=complex<T>(1,0);
		}
	}
	else{
		if(!_k1massive){
			bp.alp1=complex<T>(0,0);
			f2=complex<T>(1,0);
		}
		else{
			bp.alp1=complex<T>(1,0);
			f2=(S1*S2-pow(gammap,2))/(S1*(S2-gammap));
		}
		if(!_k2massive){
			bp.alp2=complex<T>(0,0);
			f1=complex<T>(1,0);
		}
		else{
			bp.alp2=complex<T>(1,0);
			f1=(S1*S2-pow(gammap,2))/(S2*(S1-gammap));
		}
	}

	// Compute the gammap of K1flatp and K2flatp so we can define f1 and f2.
	bp.gamma=gammap/(f1*f2);
	complex<T> gamfac=bp.gamma/(pow(gammap,2)-S1*S2);
	bp.K1flatc=Cmom<T>(f2*gamfac*(gammap*bp.K1-S1*bp.K2));
	bp.K2flatc=Cmom<T>(f1*gamfac*(gammap*bp.K2-S2*bp.K1));

//	complex<T> denalpha=T(1)/(pow(gammap,2)-S1*S2);
//	Cmom<T> K1flatcI(denalpha*gammap*(gammap*bp.K1-S1*bp.K2));
//	Cmom<T> K2flatcI(denalpha*gammap*(gammap*bp.K2-S2*bp.K1));
//
//	bp.K1flatc=Cmom<T>(K1flatcI.Lt(),((T(1)/f1))*K1flatcI.L());
//	bp.K2flatc=Cmom<T>(((T(1)/f2))*K2flatcI.Lt(),K2flatcI.L());

	bp.vec1c=PfLLt(bp.K1flatc.Lt(),bp.K2flatc.L());
	bp.vec2c=PfLLt(bp.K2flatc.Lt(),bp.K1flatc.L());

	//Extract the box contributions, the last element in the returned array is the contribution
	// at t=0 if this is a three-mass triangle
	complex<T> zero(0,0);
	for(int ipole=0;ipole<=(RatTriSpecs::CPOINTS+1)*RatTriSpecs::MUTRIPOINTS;ipole++){
		bp.tri_sub_ret[ipole]=zero;
	}
	
	//Get all the box subtractions adding each new one to the previous ones contained in coeffsret_mass
	for(size_t box=1; box <= BASE::daughters_nbr(); box++){
		bp.box_corner=BASE::get_opened_corner(box);

#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(box));
		assert(dau);
#else
		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(box));
#endif
		dau->get_sub_terms(ep,bp);
//		_MESSAGE4("      Box:",*dau," corner:",bp.box_corner);

		// Store these results for later use in the bubble coeff subtractions if we have any poles
		boxcoeff_add(box-1,bp.boxcoeffs,bp.boxpoles,bp.denfac1);
	}

	// Get the points to evaluate at at the correct precision
	const complex<T>* eval_pts;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_t_eval_pts(eval_pts);

	//We use this RatTriSpecs::CPOINTS*RatTriSpecs::CPOINTS double array contains all the coefficients we use to multiply the 7 terms we generate above
	// to produce the RatTriSpecs::CPOINTS coefficients.
	const complex<T>* coeff_pts;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_t_matrix_eval_pts(coeff_pts);

	// Get the points to evaluate m0 around
	const complex<T>* m0mass;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix;
	(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix);

	// Using the eval_param, for now we have to put the incoming momenta into the left or right side
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

	Cmom<T> l[3],ml[3];
	*(epc[0]->back())=&l[1];
	*(epc[1]->back())=&l[2];
	*(epc[2]->back())=&l[0];
	*(epc[0]->begin())=&ml[0];
	*(epc[1]->begin())=&ml[1];
	*(epc[2]->begin())=&ml[2];

	complex<T> triC0res_mass[RatTriSpecs::CPOINTS*RatTriSpecs::MUTRIPOINTS];
	momentum<complex<T> > mompt1(bp.alp2*bp.K1flatc.P()+bp.alp1*bp.K2flatc.P());
	complex<T> c0res[RatTriSpecs::MUTRIPOINTS-1]={complex<T>(0,0)}, amp[RatTriSpecs::CPOINTS];
	for(int m0pole=0;m0pole<RatTriSpecs::MUTRIPOINTS;m0pole++){
		momentum<complex<T> > mompt3((bp.alp1*bp.alp2-m0mass[m0pole]/bp.gamma)*bp.vec2c);
		for(int icutmass=0;icutmass<_cut_mass.size();icutmass++){
			eval_param<T>::set_dynamic2(_cut_mass[icutmass], m0mass[m0pole]);
		}

		for(int icirc=0; icirc<RatTriSpecs::CPOINTS; icirc++){
			l[0]=Cmom<T>(mompt1+eval_pts[icirc]*bp.vec1c+mompt3/eval_pts[icirc],_mt_massive);
			l[1]=Cmom<T>(l[0].P()-bp.K1,_mt_massive);
			l[2]=Cmom<T>(l[0].P()-bp.K2,_mt_massive);
			ml[0]=Cmom<T>(-l[0].P(),_mt_massive);
			ml[1]=Cmom<T>(-l[1].P(),_mt_massive);
			ml[2]=Cmom<T>(-l[2].P(),_mt_massive);

			amp[icirc]=-bp.tri_sub_ret[icirc+(RatTriSpecs::CPOINTS+1)*m0pole];
			//Construct the two-particle cut at the momenta above
			for(int i=0;i<BASE::decendant_nbr();i++){
				amp[icirc]+=BASE::eval_tree(i,0,*(epc[0]))*BASE::eval_tree(i,1,*(epc[1]))*BASE::eval_tree(i,2,*(epc[2]));
			}
		}

		//Now produce the coeffs by starting from 0 and adding the above terms with the correct signs
		for(int ifin=0;ifin<RatTriSpecs::CPOINTS;ifin++){
			int indexC0=ifin+RatTriSpecs::CPOINTS*m0pole;
			triC0res_mass[indexC0]=coeff_pts[ifin*RatTriSpecs::CPOINTS]*amp[0];
			for(int jfin=1;jfin<RatTriSpecs::CPOINTS;jfin++){
				triC0res_mass[indexC0]+=coeff_pts[ifin*RatTriSpecs::CPOINTS+jfin]*amp[jfin];
			}
		}

		//We compute C0 using 2c0bare-2(t_sub terms)-(t->0_terms)
		for(int iMupoints=0;iMupoints<RatTriSpecs::MUTRIPOINTS-1;iMupoints++){
			c0res[iMupoints]+=m0mass_matrix[m0pole+(iMupoints+1)*RatTriSpecs::MUPOINTS]*(triC0res_mass[(RatTriSpecs::CPOINTS-1)/2+RatTriSpecs::CPOINTS*m0pole]+bp.tri_sub_ret[RatTriSpecs::CPOINTS+(RatTriSpecs::CPOINTS+1)*m0pole]);
		}
	}

	// Finally set the actual C0 coefficient, this will be different
	//  from coeffs[3] when we have a three mass triangle where we have to add the
	//  "t=0" contribution
	complex<T> fullC0res(0,0);
	for(int imupoints=0;imupoints<RatTriSpecs::MUTRIPOINTS-1;imupoints++){
		// We compute S3 from K3^2=(-K1+K2)^2
//		_MESSAGE8("   Tri bit:",imupoints," gives ",c0res[imupoints]," * ",(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1,S2,(-bp.K1+bp.K2).square()),"=",c0res[imupoints]*(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1,S2,(-bp.K1+bp.K2).square()));
		fullC0res+=c0res[imupoints]*(static_cast<typename RatTriSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,S1,S2,(-bp.K1+bp.K2).square());
	}

	set_C0coeff(fullC0res);

	//Save the vec1 and vec2 and opened corner so we can relate the c's we compute here to those
	// computed in the other cases
	set_tri_param_basis_vectors(bp.vec1c,bp.vec2c);
	set_tri_param_gamma(bp.gamma,gammap);
	coeffkeep_add(triC0res_mass);

	//Mark that we have computed these, we use the mcID passed from the bubble as this will
	// have been computed using a sub_mom_conf
	epID=ep.get_ID();
	accuracy=bp.accuracy;
}


///*
// *
// *
// *
// *
// * Specialised versions of the triangle amplitudes
// *
// *
// *
// */
//
//
//template <class BASE> template <class T> complex<T> triangle_Rat_plusminus<BASE,typename Normal_RatTri_Specification<BASE> >::get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,3>& tp)
//{
//	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
//	 if(tp.tri_corner==2){
//		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
//			 K2sum_mom-=ep.p(BASE::corner_ind(3,kiter+1)-1)->P();
//		 }
//	 }
//	 else
//	 {
//		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
//			 K2sum_mom+=ep.p(BASE::corner_ind(tp.tri_corner,kiter+1)-1)->P();
//		 }
//	 }
//
//	 //Now construct the constants we will use and pick gamma with the
//	 // convention that if S1 or S2 would be 0 then it would not vanish
//	 // this avoids the need to specialise this statement.
//	 complex<T> S2=K2sum_mom.square();
//	 complex<T> K1K2=tp.K1*K2sum_mom;
//
//	 if(!is_eval(ep.get_ID())){
//		 get_coeffs(ep,tp.accuracy);
//	 }
//
//	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
//	 // triangle triple cut using the analytic formula below.
//	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
//	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());
//
//	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);//mc.spab(tp.chi, K2sum_momt, tp.K1flatb);
//	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
//	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
//	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);
//
//	 complex<T> ampp, ampm, ypfac;
//	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermp, subtermm;
//	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;
//
//	 // Contains the value of the triangle pole subtraction for the bubble
//	 complex<T> ampfull[3]={complex<T>(0,0),complex<T>(0,0),complex<T>(0,0)};
//	 complex<T> ampfull_mass[3]={complex<T>(0,0),complex<T>(0,0),complex<T>(0,0)};
//	 complex<T> amp_at[3];
//	 complex<T> subres_1a(0,0),subres_2a(0,0),subres_1b(0,0),subres_2b(0,0);
//
//	 // Get the points to evaluate m0 around
//	 const complex<T>* m0mass;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
//	 const complex<T>* m0mass_matrix_pts;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);
//
//	 // Get the points to evaluate the circle on
//	 const complex<T>* eval_pts;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
//
//	 //Set up the points to evaluate y around
//	 const complex<T> ypoint[3]={complex<T>(1,0),
//	 								complex<T>(-1,sqrt(T(3)))/T(2),
//	 								complex<T>(-1,-sqrt(T(3)))/T(2)};
//
//	 // Get the coefficients and other factors from the previously computed triangle, if not
//	 //  this will just compute them
//	 complex<T> gamma_old,gammap_old,*orig_coeffs;
//	 coeffkeep_get(orig_coeffs);
//	 get_tri_param_gamma(gamma_old,gammap_old);
//	 momentum<complex<T> > v1old, v2old;
//	 get_tri_param_basis_vectors(v1old,v2old);
//
//	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
//	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};
//
//	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
//	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
//	 get_boxes(triboxcoeff,triboxpoles,denfac);
//	 for(int m0pole=0;m0pole<Normal_RatTri_Specification::MUTRIPOINTS;m0pole++){
//		 int indexbase=m0pole*Normal_RatTri_Specification::CPOINTS;
//
//		 for(int icirc=0; icirc<Normal_RatTri_Specification::CBUBPOINTS; icirc++){
//			 //Compute the y+/- poles
//			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat/T(2);
//			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
//			 resb=sqrt(Delta3)*sqrt(ressqrt);
//			 Nysolp=(resa+resb)/chiK2Klflat;
//			 Nysolm=(resa-resb)/chiK2Klflat;
//			 // Construct the new t coefficient
//			 bubtsp=-T(2)*(Nysolp*fac1[0]
//			              +(T(1)-Nysolp)*fac1[1]
//		 		    	  +eval_pts[icirc]*fac1[2]
//		 		    	  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;
//
//			 bubtsm=-T(2)*(Nysolm*fac1[0]
//				          +(T(1)-Nysolm)*fac1[1]
//				          +eval_pts[icirc]*fac1[2]
//				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;
//
//			 bubtism=T(2)*(Nysolm*fac2[0]
//				 		  +(T(1)-Nysolm)*fac2[1]
//				 		  +eval_pts[icirc]*fac2[2]
//				 	 	  +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/m0mass[m0pole];
//
//			 bubtisp=T(2)*(Nysolp*fac2[0]
//				 		  +(T(1)-Nysolp)*fac2[1]
//				 		  +eval_pts[icirc]*fac2[2]
//				 		  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/m0mass[m0pole];
//
//			 subtermp=-(pow(bubtsp,3)*orig_coeffs[indexbase]+pow(bubtsp,2)*orig_coeffs[indexbase+1]+bubtsp*orig_coeffs[indexbase+2]+orig_coeffs[indexbase+3]
//				 		     +bubtisp*orig_coeffs[indexbase+4]+pow(bubtisp,2)*orig_coeffs[indexbase+5]+pow(bubtisp,3)*orig_coeffs[indexbase+6])*complex<T>(0,1)/*complex<T>(0,0.125)*/;
//			 subtermm=-(pow(bubtsm,3)*orig_coeffs[indexbase]+pow(bubtsm,2)*orig_coeffs[indexbase+1]+bubtsm*orig_coeffs[indexbase+2]+orig_coeffs[indexbase+3]
//				 			 +bubtism*orig_coeffs[indexbase+4]+pow(bubtism,2)*orig_coeffs[indexbase+5]+pow(bubtism,3)*orig_coeffs[indexbase+6])*complex<T>(0,1)/*complex<T>(0,0.125)*/;
//
//			 //Now add the box contributions to the triple cut
//			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
//				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
//				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
//			 }
//
//			 ampfull[0]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[0];
//			 ampfull_mass[0]+=amp_at[0];
//			 ampfull_mass[1]+=amp_at[1];
//			 ampfull_mass[2]+=amp_at[2];
//
//			 //Compute to test the accuracy
//			 ampfull[1]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[1];
//			 ampfull[2]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[2];
//		 }
//
//		 subres_1a+=orig_coeffs[indexbase];
//		 subres_1b+=orig_coeffs[indexbase+6]/pow(m0mass[m0pole],3);
//		 subres_2b+=orig_coeffs[indexbase+4]/pow(m0mass[m0pole],2);
//		 subres_2a+=orig_coeffs[indexbase+2]/m0mass[m0pole];
//	 }
//	 complex<T> subfacb=-(v1old*chik1b)*T(2)/chiK2Klflat;
//	 complex<T> subfac=(v2old*chik1b)*T(2)/(gamma_old*chiK2Klflat);
//	 complex<T> submufac=-(pow(K1K2-tp.S1*S2/tp.S1,2)-(tp.S1*S2-pow(K1K2,2))/T(3))/tp.S1;
//	 complex<T> subres=((subres_1a*pow(subfac,3)+subres_1b*pow(subfacb,3))*submufac+subres_2a*subfac+subres_2b*subfacb);
//
//	 //Extract the different components we need
//	 complex<T> mass_add=-(ypoint[1]*ampfull_mass[1]+ypoint[2]*ampfull_mass[2]+ampfull_mass[0])/(T(9)*tp.S1);
//
//	 //To estimate the computation error return the extra components
//   tp.bub_acc_pieces[0]+=ampfull[0];
//   tp.bub_acc_pieces[1]+=ampfull[1];
//   tp.bub_acc_pieces[2]+=ampfull[2];
//
//	 //Now return the complete subtraction result
//#if BUBPOINTS_EXT==4
//	 return complex<T>(0,1)*(ampfull[0]+mass_add)+subres/T(2);
//#else
//	 return ((ampfull[0]+mass_add)/T(Normal_RatTri_Specification::YPOINTS)+complex<T>(0,-1)*subres);
//#endif /*BUBPOINTS_EXT==4*/
//}
//
//template <class BASE> template <class T> complex<T> triangle_Rat_3mass<BASE,typename Normal_RatTri_Specification<BASE> >::get_sub_terms_work(const eval_param<T>& ep, triangle_param<T,3>& tp)
//{
//	 momentum<complex<T> > K2sum_mom(complex<T>(0,0),complex<T>(0,0),complex<T>(0,0),complex<T>(0,0));
//	 if(tp.tri_corner==2){
//		 for(size_t kiter=0;kiter<BASE::corner_size(3);kiter++){
//			 K2sum_mom-=ep.p(BASE::corner_ind(3,kiter+1)-1)->P();
//		 }
//	 }
//	 else
//	 {
//		 for(size_t kiter=0;kiter<BASE::corner_size(tp.tri_corner);kiter++){
//			 K2sum_mom+=ep.p(BASE::corner_ind(tp.tri_corner,kiter+1)-1)->P();
//		 }
//	 }
//
//	 //Now construct the constants we will use and pick gamma with the
//	 // convention that if S1 or S2 would be 0 then it would not vanish
//	 // this avoids the need to specialise this statement.
//	 complex<T> S2=K2sum_mom.square();
//	 complex<T> K1K2=tp.K1*K2sum_mom;
//
//	 //For gamma we need to choose this such that alpha1 and alpha2 are not small and
//	 // ideally of ~1. We can do this by checking that we choose the correct sign in gamma
//	 // when S2~0 and similarly when S3~0 that gamma gamma!=S1 or S2 depending on the solution
//
//	 if(!is_eval(ep.get_ID())){
//		 get_coeffs(ep,tp.accuracy);
//	 }
//
//	 //The triple-cut contributions to the bubble can be computed from the 7 coefficients of the
//	 // triangle triple cut using the analytic formula below.
//	 momentum<complex<T> > chik1b=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());
//	 momentum<complex<T> > k1bchi=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());
//
//	 complex<T> chiK2Klflat=T(2)*(chik1b*K2sum_mom);//mc.spab(tp.chi, K2sum_momt, tp.K1flatb);
//	 complex<T> VPK1flatK2=tp.K1flatbc.P()*K2sum_mom;
//	 complex<T> VPchiK2=tp.chic.P()*K2sum_mom;
//	 complex<T> Delta3=(K1K2*K1K2-tp.S1*S2);
//
//	 complex<T> ampp, ampm, ypfac;
//	 complex<T> bubtsp, bubtsm, bubtisp, bubtism, subtermp, subtermm;
//	 complex<T> Nysolp, Nysolm, ressqrt, resa, resb;
//
//	 // Contains the value of the triangle pole subtraction for the bubble
//	 complex<T> ampfull[Normal_RatTri_Specification::YPOINTS]={complex<T>(0,0)};
//	 complex<T> ampfull_mass[Normal_RatTri_Specification::YPOINTS]={complex<T>(0,0)};
//	 complex<T> amp_at[Normal_RatTri_Specification::YPOINTS];
//	 complex<T> subres_1a(0,0),subres_1b(0,0),subres_2a(0,0),subres_2b(0,0);
//
//	 // Get the points to evaluate m0 around
//	 const complex<T>* m0mass;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
//	 const complex<T>* m0mass_matrix_pts;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);
//
//	 // Get the points to evaluate the circle on
//	 const complex<T>* eval_pts;
//	 (static_cast<typename Normal_RatTri_Specification::CornerTreeStrategy*>(this))->BubPoints.get_t_eval_pts(eval_pts);
//
//	 // Get the coefficients and other factors from the previously computed triangle, if not
//	 //  this will just compute them
//	 complex<T> gamma_old,gammap_old,*orig_coeffs;
//	 coeffkeep_get(orig_coeffs);
//	 get_tri_param_gamma(gamma_old,gammap_old);
//	 momentum<complex<T> > v1old, v2old;
//	 get_tri_param_basis_vectors(v1old,v2old);
//
//	 complex<T> fac1[4]={v2old*tp.K1flatbc.P(),v2old*tp.chic.P(),v2old*k1bchi,v2old*chik1b};
//	 complex<T> fac2[4]={v1old*tp.K1flatbc.P(),v1old*tp.chic.P(),v1old*k1bchi,v1old*chik1b};
//
//
//	 // Get the points to evaluate y around
//	 const complex<T> ypoint[3]={complex<T>(1,0),
//								 complex<T>(-1,sqrt(T(3)))/T(2),
//								 complex<T>(-1,-sqrt(T(3)))/T(2)};
//	 // Get the box terms for this triangle, we must have already called the get_coeffs function for this to work
//	 std::vector<complex<T> > *triboxcoeff, *triboxpoles, *denfac;
//	 get_boxes(triboxcoeff,triboxpoles,denfac);
//	 for(int m0pole=0;m0pole<Normal_RatTri_Specification::MUTRIPOINTS;m0pole++){
//		 int indexbase=m0pole*Normal_RatTri_Specification::CPOINTS;
//
//		 for(int icirc=0; icirc<Normal_RatTri_Specification::CBUBPOINTS; icirc++){
//			 //Compute the y+/- poles
//			 resa=(VPK1flatK2-VPchiK2)*eval_pts[icirc]+chiK2Klflat*T(0.5);
//			 ressqrt=pow(eval_pts[icirc],2)+eval_pts[icirc]*(K1K2-S2)*chiK2Klflat/Delta3+(T(1)/T(4)-m0mass[m0pole]/tp.gammab)*pow(chiK2Klflat,2)/Delta3;
//			 resb=sqrt(Delta3)*sqrt(ressqrt);
//
//			 Nysolp=(resa+resb)/chiK2Klflat;
//			 Nysolm=(resa-resb)/chiK2Klflat;
//
//			 // Construct the new t coefficient
//			 bubtsp=-T(2)*(Nysolp*fac1[0]
//			                 +(T(1)-Nysolp)*fac1[1]
//	 		 		    	 +eval_pts[icirc]*fac1[2]
//	 		 		    	 +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;
//
//			 bubtsm=-T(2)*(Nysolm*fac1[0]
//	 				          +(T(1)-Nysolm)*fac1[1]
//	 				          +eval_pts[icirc]*fac1[2]
//	 				          +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac1[3]/eval_pts[icirc])/gamma_old;
//
//			 bubtism=-T(2)*(Nysolm*fac2[0]
//	 				 		  +(T(1)-Nysolm)*fac2[1]
//	 				 		  +eval_pts[icirc]*fac2[2]
//	 				 	 	  +(Nysolm*(T(1)-Nysolm)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/(gamma_old-m0mass[m0pole]);
//
//			 bubtisp=-T(2)*(Nysolp*fac2[0]
//	 				 		  +(T(1)-Nysolp)*fac2[1]
//	 				 		  +eval_pts[icirc]*fac2[2]
//	 				 		  +(Nysolp*(T(1)-Nysolp)-m0mass[m0pole]/tp.gammab)*fac2[3]/eval_pts[icirc])/(gamma_old-m0mass[m0pole]);
//
//			 subtermp=-(pow(bubtsp,3)*orig_coeffs[indexbase]+pow(bubtsp,2)*orig_coeffs[indexbase+1]+bubtsp*orig_coeffs[indexbase+2]+orig_coeffs[indexbase+3]
//	 				 		     +bubtisp*orig_coeffs[indexbase+4]+pow(bubtisp,2)*orig_coeffs[indexbase+5]+pow(bubtisp,3)*orig_coeffs[indexbase+6])*complex<T>(0,1)/*complex<T>(0,0.125)*/;
//			 subtermm=-(pow(bubtsm,3)*orig_coeffs[indexbase]+pow(bubtsm,2)*orig_coeffs[indexbase+1]+bubtsm*orig_coeffs[indexbase+2]+orig_coeffs[indexbase+3]
//	 				 			 +bubtism*orig_coeffs[indexbase+4]+pow(bubtism,2)*orig_coeffs[indexbase+5]+pow(bubtism,3)*orig_coeffs[indexbase+6])*complex<T>(0,1)/*complex<T>(0,0.125)*/;
//			 //Now add the box contributions to the triple cut
//			 for(size_t ibox=0;ibox<BASE::daughters_nbr();ibox++){
//				 subtermp-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsp-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsp-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
//				 subtermm-=((triboxpoles[m0pole])[2*ibox]*(triboxcoeff[m0pole])[2*ibox]/(bubtsm-(triboxpoles[m0pole])[2*ibox])-(triboxpoles[m0pole])[2*ibox+1]*(triboxcoeff[m0pole])[2*ibox+1]/(bubtsm-(triboxpoles[m0pole])[2*ibox+1]))/((*denfac)[ibox]*((triboxpoles[m0pole])[2*ibox+1]-(triboxpoles[m0pole])[2*ibox]));
//			 }
//
//			 // Now construct the actual subtraction term
//			 ypfac=eval_pts[icirc]/(T(2)*resb);
//			 amp_at[0]=(subtermp/(ypoint[0]-Nysolp)-subtermm/(ypoint[0]-Nysolm))*ypfac;
//			 amp_at[1]=(subtermp/(ypoint[1]-Nysolp)-subtermm/(ypoint[1]-Nysolm))*ypfac;
//			 amp_at[2]=(subtermp/(ypoint[2]-Nysolp)-subtermm/(ypoint[2]-Nysolm))*ypfac;
//			 ampfull[0]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[0];
//			 ampfull_mass[0]+=amp_at[0];
//			 ampfull_mass[1]+=amp_at[1];
//			 ampfull_mass[2]+=amp_at[2];
//
//			 //Compute to test the accuracy
//			 ampfull[1]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[1];
//			 ampfull[2]+=m0mass_matrix_pts[m0pole+RatTriSpecs::MUPOINTS]*amp_at[2];
//		 }
//
//		 //Now subtract off the triple-cut contributions
//		 subres_1a+=orig_coeffs[indexbase];
//		 subres_1b+=orig_coeffs[indexbase+6]/pow((gamma_old-m0mass[m0pole]),3);
//		 subres_2b+=orig_coeffs[indexbase+4]/(m0mass[m0pole]*(gamma_old-m0mass[m0pole]));
//		 subres_2a+=orig_coeffs[indexbase+2]/m0mass[m0pole];
//	 }
//	 complex<T> subfacb=v1old*chik1b*T(2)/chiK2Klflat;
//	 complex<T> subfac=v2old*chik1b*T(2)/(gamma_old*chiK2Klflat);
//	 complex<T> gamma_fac=(K1K2+sqrt(pow(K1K2,2)-tp.S1*S2));//This is gammap for this triangle subtraction we need it for the line below
//	 complex<T> submufac=-(pow(K1K2-tp.S1*S2/gamma_fac,2)-(tp.S1*S2-pow(K1K2,2))/T(3))/tp.S1;
//	 complex<T> subres=((subres_1a*pow(subfac,3)+subres_1b*pow(subfacb,3))*submufac+subres_2a*subfac+subres_2b*subfacb);
//
//	 //Extract the different components we need
//	 complex<T> mass_add=-(ypoint[1]*ampfull_mass[1]+ypoint[2]*ampfull_mass[2]+ampfull_mass[0])/(T(9)*tp.S1);
//
//	 //To estimate the computation error return the extra components
//   tp.bub_acc_pieces[0]+=ampfull[0];
//   tp.bub_acc_pieces[1]+=ampfull[1];
//   tp.bub_acc_pieces[2]+=ampfull[2];
//
//	 //Now return the complete subtraction result
//#if BUBPOINTS_EXT==4
//	 return complex<T>(0,1)*(ampfull[0]+mass_add)+subres/T(2);
//#else
//	 return ((ampfull[0]+mass_add)/T(Normal_RatTri_Specification::YPOINTS)+complex<T>(0,-1)*subres);
//#endif /*BUBPOINTS_EXT==4*/
//}


}

}

#endif /* TRIANGLE_RATEXT_HPP_ */
