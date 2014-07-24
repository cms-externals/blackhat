/*
 * bubble_ratext.cpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#ifndef BUBBLE_RATEXT_HPP_
#define BUBBLE_RATEXT_HPP_

#include "ratext/bubble_ratext.h"
#include "ratext/rat_ext.h"
#include <cassert>
#include "BH_debug.h"

#define _OLD_PARAM 1

namespace BH {

namespace ratext {

template <class BASE, class RatBubSpecs> C bubble_Rat<BASE, RatBubSpecs>::eval(const eval_param<R>& ep){
	return get_coeffs(ep);
}

template <class BASE, class RatBubSpecs> CHP bubble_Rat<BASE, RatBubSpecs>::eval(const eval_param<RHP>& ep){
	return get_coeffs(ep);
}

template <class BASE, class RatBubSpecs> CVHP bubble_Rat<BASE, RatBubSpecs>::eval(const eval_param<RVHP>& ep){
	return get_coeffs(ep);
}


template <class BASE, class RatBubSpecs> C bubble_Rat<BASE, RatBubSpecs>::eval(mom_conf& mc, const std::vector<int>& ind){
	return get_coeffs(mc,ind);
}

template <class BASE, class RatBubSpecs> CHP bubble_Rat<BASE, RatBubSpecs>::eval(mom_conf_HP& mc, const std::vector<int>& ind){
	return get_coeffs(mc,ind);
}

template <class BASE, class RatBubSpecs> CVHP bubble_Rat<BASE, RatBubSpecs>::eval(mom_conf_VHP& mc, const std::vector<int>& ind){
	return get_coeffs(mc,ind);
}

#if BH_USE_GMP
template <class BASE, class RatBubSpecs> CGMP bubble_Rat<BASE, RatBubSpecs>::eval(const eval_param<RGMP>& ep){
	return get_coeffs(ep);
}
template <class BASE, class RatBubSpecs> CGMP bubble_Rat<BASE, RatBubSpecs>::eval(momentum_configuration<RGMP>& mc, const std::vector<int>& ind){
	return get_coeffs(mc,ind);
}
#endif

template <class BASE, class RatBubSpecs> void bubble_Rat<BASE, RatBubSpecs>::init()
{
	//Set up the ind's for the cut legs
	indlst[0].assign(BASE::corner_size(1)+4,0);
	indlst[1].assign(BASE::corner_size(2)+4,0);

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

}

template <class BASE, class RatBubSpecs> void bubble_Rat<BASE, RatBubSpecs>::add_mass(const cutD& pcd)
{
	add_mass(pcd.l(1).mass_label(),pcd.l(2).mass_label());
}

template <class BASE, class RatBubSpecs> void bubble_Rat<BASE, RatBubSpecs>::add_mass(int m1,int m2)
{
	if(find(_cut_mass.begin(),_cut_mass.end(),m1)==_cut_mass.end()){_cut_mass.push_back(m1);};
	if(find(_cut_mass.begin(),_cut_mass.end(),m2)==_cut_mass.end()){_cut_mass.push_back(m2);};

	//Set up the masses for the legs (we reset these each time this is uneccassary fix this later)
	_leg_masses->set(0,m1);
	_leg_masses->set(1,m2);
}


template <class BASE, class RatBubSpecs> bubble_Rat<BASE, RatBubSpecs>::~bubble_Rat()
{
	delete _ep[0];
	delete _ep_HP[0];
	delete _ep_VHP[0];
	delete _ep[1];
	delete _ep_HP[1];
	delete _ep_VHP[1];
#if BH_USE_GMP
	delete _ep_GMP[0];
	delete _ep_GMP[1];
#endif
	 delete _leg_masses;
}


template <class BASE, class RatBubSpecs> template <class T> complex<T> bubble_Rat<BASE, RatBubSpecs>::get_coeffs(const eval_param<T>& ep)
{
	//We will store all information that needs to be passed between the bubble, triangle and box in the coeffparam_mass structure
	triangle_param<T,RatBubSpecs::MUBUBPOINTS,RatBubSpecs::YPOINTS> tp;

	//Set/Reset the accuracy parameter to the maximum number of digits we expect to get
	tp.accuracy=to_double(MaxDigits<T>());

	//Mark that we have computed this bubble
	epID=ep.get_ID();

	// Construct K1 by summing over all the legs of c(1) using ind to get their location in the momconf
	momentum<complex<T> > K1sum_mom(ep.p(BASE::corner_ind(1,1)-1)->P());
	for(size_t k1iter=2;k1iter<=BASE::corner_size(1);k1iter++){
		K1sum_mom+=ep.p(BASE::corner_ind(1,k1iter)-1)->P();
	}
	tp.K1=K1sum_mom;
	tp.S1=tp.K1.square();

	//Construct chi so that all the coefficients of l^{\mu} are of order 1
	// this is achieved by scaling the vector (1,0,1,0) (we want to avoid a chi proportional to one of the external legs)
	// so that gamma=S1
//	tp.chi=sub_mc.insert(momentum<complex<T> >(sqrt(complex<T>(2,0)),complex<T>(2,0),complex<T>(0,1),complex<T>(0,1)),_mt_unknown);
//	tp.chi=sub_mc.insert(momentum<complex<T> >(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,1)),_mt_unknown);
	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,1));

	// Now rescale to the actual chi
	tp.chic=Cmom<T>(tp.S1/(T(2)*(chi_init*tp.K1))*chi_init);
	tp.gammab=tp.S1;
	tp.K1flatbc=Cmom<T>(tp.K1-tp.chic.P());
	momentum<complex<T> > vec_bub1=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());
	momentum<complex<T> > vec_bub2=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());

	// Using the eval_param, for now we have to put the incoming momenta into the left or right side
	eval_param<T>** epc;
	get_ep(epc);

	for(int mm=1;mm<=BASE::corner_size(1);mm++){
		epc[0]->set(mm,ep.p(BASE::corner_ind(1,mm)-1));
	}
	for(int mm=1;mm<=BASE::corner_size(2);mm++){
		epc[1]->set(mm,ep.p(BASE::corner_ind(2,mm)-1));
	}
	Cmom<T> l[2], ml[2];
	*(epc[0]->back())=&l[1];
	*(epc[1]->back())=&l[0];
	*(epc[0]->begin())=&ml[0];
	*(epc[1]->begin())=&ml[1];

	// Get the points to evaluate the circle on
	const complex<T>* eval_pts;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_t_eval_pts(eval_pts);
	//Set up the m0 masses
	const complex<T>* m0mass;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix_pts;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);
	// Get the points to evaluate in y
	const complex<T> *ypoint;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_y_eval_pts(ypoint);

//#if _OLD_PARAM==0
//	 momentum<complex<T> > n1=complex<T>(0,1)*(tp.K1flatbc.P()-tp.chic.P())/sqrt(tp.S1);
//	 Cmom<T> nA(tp.K1flatbc.L(),tp.chic.Lt());
//	 Cmom<T> nB(tp.chic.L(),tp.K1flatbc.Lt());
//	 momentum<complex<T> > n2=complex<T>(0,1)*(nA.P()+nB.P())/sqrt(tp.gammab);
//	 momentum<complex<T> > n3=(nA.P()-nB.P())/sqrt(tp.gammab);
//
//	complex<T> amplm[RatBubSpecs::YPOINTS*(RatBubSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
//	complex<T> amplp[RatBubSpecs::YPOINTS*(RatBubSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
//	complex<T> amptemp;
//	for(int m0pole=0;m0pole<RatBubSpecs::MUBUBPOINTS;m0pole++){
//		for(int icutmass=0;icutmass<_cut_mass.size();icutmass++){
//			eval_param<T>::set_dynamic2(_cut_mass[icutmass], m0mass[m0pole]);
//		}
//		for(int j=0;j<RatBubSpecs::YPOINTS;j++){ // We will want different y components
//			for(int i=0;i<RatBubSpecs::CBUBPOINTS;i++){ // We will only ever want the t^0 component
//				complex<T> sqrsol=sqrt(m0mass[m0pole]-tp.S1/T(4)-eval_pts[i]*eval_pts[i]-ypoint[j]*ypoint[j]);
//				l[0]=Cmom<T>(tp.K1+sqrsol*n1+eval_pts[i]*n2+ypoint[j]*n3,_mt_massive);
//				l[1]=Cmom<T>(l[0].P()-tp.K1,_mt_massive);
//				ml[0]=Cmom<T>(-l[0].P(),_mt_massive);
//				ml[1]=Cmom<T>(-l[1].P(),_mt_massive);
//				
////				l[0]=Cmom<T>(tp.K1-sqrsol*n1+eval_pts[i]*n2+ypoint[j]*n3,_mt_massive);
////				l[1]=Cmom<T>(l[0].P()-tp.K1,_mt_massive);
////				ml[0]=Cmom<T>(-l[0].P(),_mt_massive);
////				ml[1]=Cmom<T>(-l[1].P(),_mt_massive);				
//
//				//Construct the two-particle cut at the momenta above, we must have at least one bubble
//				amptemp=BASE::eval_tree(0,0,*(epc[0]))*BASE::eval_tree(0,1,*(epc[1]));
//				//Construct the two-particle cut at the momenta above
//				for(int ii=1;ii<BASE::decendant_nbr();ii++){
//					amptemp+=BASE::eval_tree(ii,0,*(epc[0]))*BASE::eval_tree(ii,1,*(epc[1]));
//				}
//				for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
//					amplp[j+imupoints*RatBubSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatBubSpecs::MUPOINTS]*amptemp;
//				}
//			}
//		}
//	}
//
//	//Sum over all the triangles getting their subtraction pieces
//	complex<T> piece[RatBubSpecs::YPOINTS]={complex<T>(0,0)};
//	for(int imupoints=0;imupoints<RatBubSpecs::YPOINTS;imupoints++){
//		tp.tri_sub_pieces_p[imupoints]=complex<T>(0,0);
//		tp.tri_sub_pieces_m[imupoints]=complex<T>(0,0);
//		tp.tri_sub_pieces_0p[imupoints]=complex<T>(0,0);
//		tp.tri_sub_pieces_0m[imupoints]=complex<T>(0,0);
//	}
//	for(size_t tri=1; tri <= BASE::daughters_nbr(); tri++){
//		tp.tri_corner=BASE::get_opened_corner(tri);
//#ifndef NDEBUG
//		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(tri));
//		assert(dau);
//#else
//		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(tri));
//#endif
//		dau->get_sub_terms(ep,tp);
//	}
//
//	//Extract the different components we need
//	const complex<T> *y_matrix_point;
//	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_y_matrix_eval_pts(y_matrix_point);
//	complex<T> err_bubp(0,0), err_bubm(0,0), err_tri(0,0), err_trip(0,0), err_trim(0,0), err_tri0p(0,0), err_tri0m(0,0), resp[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)},
//			resm[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)},
//			res_subp[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)},
//			res_subm[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)},
//			res_sub0p[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)},
//			res_sub0m[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
//	for(int ext=0;ext<RatBubSpecs::YPOINTS;ext++){
//		// Compute the m^2 contribution which is our result
//		for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
//			resp[imupoints]+=amplp[ext+imupoints*RatBubSpecs::YPOINTS];
//			resm[imupoints]+=amplm[ext+imupoints*RatBubSpecs::YPOINTS];
//			res_subp[imupoints]+=tp.tri_sub_pieces_p[ext+imupoints*RatBubSpecs::YPOINTS];
//			res_subm[imupoints]+=tp.tri_sub_pieces_m[ext+imupoints*RatBubSpecs::YPOINTS];
//			res_sub0p[imupoints]+=tp.tri_sub_pieces_0p[ext+imupoints*RatBubSpecs::YPOINTS];
//			res_sub0m[imupoints]+=tp.tri_sub_pieces_0m[ext+imupoints*RatBubSpecs::YPOINTS];
//		}
//		// Compute the m^2 y contribution to be used in the error analysis
//		err_bubp+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*amplp[ext];
//		err_bubm+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*amplm[ext];
//		err_trip+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*tp.tri_sub_pieces_p[ext];
//		err_trim+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*tp.tri_sub_pieces_m[ext];
//		err_tri0p+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*tp.tri_sub_pieces_0p[ext];
//		err_tri0m+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*tp.tri_sub_pieces_0m[ext];
//	}
//
//	// Check the cancellation of the m^2 y for our error analysis
//	double bubble_acc=-to_double(log(abs(err_bubp-err_trip-err_tri0p-err_tri0m))/log(10.));
//	// If the accuracy of the bubble pieces (and consequentially the triangle and the boxes pieces) is less than
//	//  the accuracy of the pentagons (which only effect the box result and do not filter down to the triangle and bubble)
//	//  then that is the minimum accuracy
//	if(bubble_acc<tp.accuracy){
//		set_accuracy(bubble_acc);
//	}
//	else{
//		set_accuracy(tp.accuracy);
//	}
//	_PRINT(err_bubp);
//	_PRINT(err_bubm);
//	_PRINT(err_trip);
//	_PRINT(err_trim);
//	_PRINT(err_tri0p);
//	_PRINT(err_tri0m);
//
//	//Calculate the bubble coefficient including all subtractions
//	complex<T> bubfullresult(0,0);
//	for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
//		_PRINT(resp[imupoints]/T(RatBubSpecs::YPOINTS));
//		_PRINT(resm[imupoints]/T(RatBubSpecs::YPOINTS));
//		_PRINT(res_subp[imupoints]/T(RatBubSpecs::YPOINTS));
//		_PRINT(res_subm[imupoints]/T(RatBubSpecs::YPOINTS));
//		_PRINT(res_sub0p[imupoints]/T(RatBubSpecs::YPOINTS));
//		_PRINT(res_sub0m[imupoints]/T(RatBubSpecs::YPOINTS));
//		bubfullresult+=((resp[imupoints]+resm[imupoints]+res_subp[imupoints]+res_subm[imupoints]+res_sub0p[imupoints]+res_sub0m[imupoints])/T(RatBubSpecs::YPOINTS))*(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,tp.S1);
//	}
//	set_B0coeff(bubfullresult);
//
//	//We have already calculated this so return the value
//	return bubfullresult;
//
//#else
	//Extract the different components we need
	const complex<T> *y_matrix_point;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_y_matrix_eval_pts(y_matrix_point);

	complex<T> ampl[(RatBubSpecs::YPOINTS-2)*(RatBubSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)}, amplm[(RatBubSpecs::YPOINTS-2)*(RatBubSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	complex<T> amptemp,massform,err_bub[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)};
	for(int m0pole=0;m0pole<RatBubSpecs::MUBUBPOINTS;m0pole++){
		for(int icutmass=0;icutmass<_cut_mass.size();icutmass++){
			eval_param<T>::set_dynamic2(_cut_mass[icutmass], m0mass[m0pole]);
		}
		
		for(int i=0;i<RatBubSpecs::CBUBPOINTS;i++){ // We will only ever want the t^0 component
			for(int j=0;j<RatBubSpecs::YPOINTS;j++){ // We will want different y components
				l[0]=Cmom<T>(ypoint[j]*tp.K1flatbc.P()+(T(1)-ypoint[j])*tp.chic.P()+eval_pts[i]*vec_bub1
							+(ypoint[j]*(T(1)-ypoint[j])-m0mass[m0pole]/tp.gammab)/eval_pts[i]*vec_bub2,_mt_massive);
				l[1]=Cmom<T>(l[0].P()-tp.K1,_mt_massive);
				ml[0]=Cmom<T>(-l[0].P(),_mt_massive);
				ml[1]=Cmom<T>(-l[1].P(),_mt_massive);

				//Construct the two-particle cut at the momenta above, we must have at least one bubble

				amptemp=BASE::eval_tree(0,0,*(epc[0]))*BASE::eval_tree(0,1,*(epc[1]));
				//Construct the two-particle cut at the momenta above
				for(int iamp=1;iamp<BASE::decendant_nbr();iamp++){
					amptemp+=BASE::eval_tree(iamp,0,*(epc[0]))*BASE::eval_tree(iamp,1,*(epc[1]));
				}
				for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
					// For the piece independent of m^2 we can simply add all the terms together as the conversion matrix for m is just 1
					//  we need to get all terms from y^2 upwards
					// TODO: fix this for higher mu powers. Need to rescale by number of mu points I think
					for(int ylvl=0;ylvl<RatBubSpecs::YPOINTS-2;ylvl++){
						ampl[imupoints+ylvl*(RatBubSpecs::MUBUBPOINTS-1)]+=m0mass_matrix_pts[m0pole+imupoints*RatBubSpecs::MUBUBPOINTS]*y_matrix_point[j+(2+ylvl)*RatBubSpecs::YPOINTS]*amptemp;
					}

					// Extract the m^2 piece of the result
					massform=m0mass_matrix_pts[m0pole+(imupoints+1)*RatBubSpecs::MUPOINTS]*amptemp;
					// We want to sum over all the m^2(y^0+y^1+y^2+...) terms 
					for(int ylp=0;ylp<RatBubSpecs::YPOINTS-2;ylp++){
						amplm[imupoints+ylp*(RatBubSpecs::MUBUBPOINTS-1)]+=y_matrix_point[j+ylp*RatBubSpecs::YPOINTS]*massform/T(ylp+1);
					}
					// Compute the m^2 y contribution to be used in the error analysis
//					err_bub[imupoints]+=y_matrix_point[j+(1+imupoints)*RatBubSpecs::YPOINTS]*massform;
//					//For the error analysis look at the y^3 component
//					err_bub[imupoints]+=y_matrix_point[j+(3+imupoints)*RatBubSpecs::YPOINTS]*amptemp;
//					// Compute the m^2 y^2 contribution to be used in the error analysis
					err_bub[imupoints]+=y_matrix_point[j+(2+imupoints)*RatBubSpecs::YPOINTS]*massform;
				}
			}
		}
	}

	//Sum over all the triangles getting their subtraction pieces
	for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
		tp.bub_acc_pieces[imupoints]=complex<T>(0,0);
		tp.tri_sub_pieces[imupoints]=complex<T>(0,0);
	}

//	_MESSAGE2("Bub:",*this);
	for(size_t tri=1; tri <= BASE::daughters_nbr(); tri++){
		tp.tri_corner=BASE::get_opened_corner(tri);
#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(tri));
		assert(dau);
#else
		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(tri));
#endif
		dau->get_sub_terms(ep,tp);
	}

	// Check the cancellation of the m^2 y for our error analysis
	double bubble_acc=-to_double(log(abs(err_bub[0]-tp.bub_acc_pieces[0])))/log(10.);
	// If the accuracy of the bubble pieces (and consequentially the triangle and the boxes pieces) is less than
	//  the accuracy of the pentagons (which only effect the box result and do not filter down to the triangle and bubble)
	//  then that is the minimum accuracy
	if(bubble_acc<tp.accuracy){
		set_accuracy(bubble_acc);
	}
	else{
		set_accuracy(tp.accuracy);
	}

	//Calculate the bubble coefficient including all subtractions
	complex<T> bubfullresult(0,0);
	for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
		complex<T> sumtmp=-tp.tri_sub_pieces[imupoints];
		for(int iypts=0;iypts<RatBubSpecs::YPOINTS-2;iypts++){
			sumtmp+=amplm[imupoints+iypts*(RatBubSpecs::MUBUBPOINTS-1)]+ampl[imupoints+iypts*(RatBubSpecs::MUBUBPOINTS-1)]*(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_extra_fac(iypts,tp.S1);
		}
		
		// Multiply the coefficient by the rational integral for the bubble
		bubfullresult+=sumtmp*(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,tp.S1)/T(RatBubSpecs::CBUBPOINTS);
	}

	set_B0coeff(bubfullresult);

	//We have already calculated this so return the value
	return bubfullresult;
//#endif
}


template <class BASE, class RatBubSpecs> template <class T> complex<T> bubble_Rat<BASE, RatBubSpecs>::get_coeffs(momentum_configuration<T>& mc, const std::vector<int>& ind)
{
	//We will store all information that needs to be passed between the bubble and triangle in the structure
	triangle_param<T,RatBubSpecs::MUBUBPOINTS,RatBubSpecs::YPOINTS> tp;

	//Set/Reset the accuracy parameter to the maximum number of digits we expect to get
	tp.accuracy=to_double(MaxDigits<T>());

	// Create a sub momentum configuration so that at the end of the computation we can dump all the
	//  newly created momenta that are no longer needed
	sub_momentum_configuration<T> sub_mc(mc);

	//Mark that we have computed this bubble
	tp.mcID=mc.get_ID();
	mcID=tp.mcID;
	indID=ind;

	// Construct K1 by summing over all the legs of c(1) using ind to get their location in the momconf
	momentum<complex<T> > K1sum_mom(sub_mc.mom(ind[BASE::corner_ind(1,1)-1]));
	for(size_t k1iter=2;k1iter<=BASE::corner_size(1);k1iter++){
		K1sum_mom+=sub_mc.mom(ind[BASE::corner_ind(1,k1iter)-1]);
	}
	tp.K1=K1sum_mom;
	tp.S1=tp.K1.square();

	//Construct chi so that all the coefficients of l^{\mu} are of order 1
	// this is achieved by scaling the vector (1,0,1,0) (we want to avoid a chi proportional to one of the external legs)
	// so that gamma=S1
//	tp.chi=sub_mc.insert(momentum<complex<T> >(sqrt(complex<T>(2,0)),complex<T>(2,0),complex<T>(0,1),complex<T>(0,1)),_mt_unknown);
//	tp.chi=sub_mc.insert(momentum<complex<T> >(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,1)),_mt_unknown);
	momentum<complex<T> > chi_init(sqrt(complex<T>(-2,2)),complex<T>(1,1),complex<T>(0,1),complex<T>(0,1));

	// Now rescale to the actual chi
	tp.chic=Cmom<T>(tp.S1/(T(2)*(chi_init*tp.K1))*chi_init);
	tp.gammab=tp.S1;
	tp.K1flatbc=Cmom<T>(tp.K1-tp.chic.P());
	momentum<complex<T> > vec_bub1=PfLLt(tp.chic.Lt(),tp.K1flatbc.L());
	momentum<complex<T> > vec_bub2=PfLLt(tp.K1flatbc.Lt(),tp.chic.L());

	// Get the points to evaluate the circle on
	const complex<T>* eval_pts;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_t_eval_pts(eval_pts);
	//Set up the m0 masses
	const complex<T>* m0mass;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_eval_pts(m0mass);
	const complex<T>* m0mass_matrix_pts;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->BoxPoints.get_mu_matrix_eval_pts(m0mass_matrix_pts);
	// Get the points to evaluate in y
	const complex<T> *ypoint;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_y_eval_pts(ypoint);

	// Set up the reference momenta for the vertices
//	tp.vertex_ref=vec_bub2;
	tp.vertex_ref=sub_mc.insert(momentum<complex<T> >(sqrt(complex<T>(-2,2)),complex<T>(0,1),complex<T>(1,1),complex<T>(0,1)),_mt_unknown);

	//Insert the legs for the two cut tree amplitudes
	for(size_t mm=1;mm<=BASE::corner_size(1);mm++){
		(indlst[0])[mm]=ind[BASE::corner_ind(1,mm)-1];
	}
	for(size_t mm=1;mm<=BASE::corner_size(2);mm++){
		(indlst[1])[mm]=ind[BASE::corner_ind(2,mm)-1];
	}
	*(indlst[0].end()-2)=tp.vertex_ref;
	*(indlst[1].end()-2)=tp.vertex_ref;

	complex<T> ampl[RatBubSpecs::YPOINTS]={complex<T>(0,0)}, amplm[RatBubSpecs::YPOINTS*(RatBubSpecs::MUBUBPOINTS-1)]={complex<T>(0,0)};
	complex<T> amptemp;
	for(int m0pole=0;m0pole<RatBubSpecs::MUBUBPOINTS;m0pole++){
		size_t ref_mass=sub_mc.insert(momentum<complex<T> >(m0mass[m0pole],sqrt(m0mass[m0pole]),complex<T>(0,0),complex<T>(0,0)),_mt_massive);
		for(int j=0;j<RatBubSpecs::YPOINTS;j++){ // We will want different y components
			for(int i=0;i<RatBubSpecs::CBUBPOINTS;i++){ // We will only ever want the t^0 component
				//Construct the cut-momenta
				*(indlst[1].end()-3)=sub_mc.insert(ypoint[j]*tp.K1flatbc.P()
								                      +(T(1)-ypoint[j])*tp.chic.P()
								                      +eval_pts[i]*vec_bub1
								                      +(ypoint[j]*(T(1)-ypoint[j])-m0mass[m0pole]/tp.gammab)/eval_pts[i]*vec_bub2,_mt_massive);
				*(indlst[0].end()-3)=sub_mc.insert(sub_mc.mom(*((indlst[1].end())-3))-tp.K1,_mt_massive);
				*(indlst[0].begin())=sub_mc.insert(-sub_mc.mom(*(indlst[1].end()-3)),_mt_massive);
				*(indlst[1].begin())=sub_mc.insert(-sub_mc.mom(*(indlst[0].end()-3)),_mt_massive);

				*(indlst[0].end()-1)=ref_mass;
				*(indlst[1].end()-1)=ref_mass;

				amptemp=BASE::eval_tree(0,0,sub_mc,indlst[0])*BASE::eval_tree(0,1,sub_mc,indlst[1]);
				//Construct the two-particle cut at the momenta above
				for(int ii=1;ii<BASE::decendant_nbr();ii++){
					amptemp+=BASE::eval_tree(ii,0,sub_mc,indlst[0])*BASE::eval_tree(ii,1,sub_mc,indlst[1]);
				}

				ampl[j]+=amptemp;//The inversion matrix for the mu independent terms is always just a sum of 1's so we can skip multiplying by it
				for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
					amplm[j+imupoints*RatBubSpecs::YPOINTS]+=m0mass_matrix_pts[m0pole+(imupoints+1)*RatBubSpecs::MUPOINTS]*amptemp;
				}
			}
		}
	}

	//Sum over all the triangles getting their subtraction pieces
	complex<T> piece[RatBubSpecs::YPOINTS]={complex<T>(0,0)};
	for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
		tp.tri_sub_pieces[imupoints]=complex<T>(0,0);
	}
	for(size_t tri=1; tri <= BASE::daughters_nbr(); tri++){
		tp.tri_corner=BASE::get_opened_corner(tri);
#ifndef NDEBUG
		dau_type_p dau=dynamic_cast<dau_type_p>(BASE::get_daughter(tri));
		assert(dau);
#else
		dau_type_p dau=static_cast<dau_type_p>(BASE::get_daughter(tri));
#endif
		dau->get_sub_terms(sub_mc,ind,tp);
	}

	//Extract the different components we need
	const complex<T> *y_matrix_point;
	(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_y_matrix_eval_pts(y_matrix_point);
	complex<T> mass_add[RatBubSpecs::MUBUBPOINTS-1]={complex<T>(0,0)}, err_bub(0,0), err_tri(0,0);
	for(int ext=0;ext<RatBubSpecs::YPOINTS;ext++){
		// Compute the y^2 contribution to the result
		mass_add[0]+=y_matrix_point[ext+2*RatBubSpecs::YPOINTS]*ampl[ext];

		// Compute the m^2 y contribution to be used in the error analysis
		err_bub+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*amplm[ext];
		err_tri+=y_matrix_point[ext+RatBubSpecs::YPOINTS]*tp.bub_acc_pieces[ext];
	}
	mass_add[0]/=(-T(18)*tp.S1);

	// Check the cancellation of the m^2 y for our error analysis
	double bubble_acc=-to_double(log(abs(err_bub-err_tri)))/log(10.);
	// If the accuracy of the bubble pieces (and consequentially the triangle and the boxes pieces) is less than
	//  the accuracy of the pentagons (which only effect the box result and do not filter down to the triangle and bubble)
	//  then that is the minimum accuracy
	if(bubble_acc<tp.accuracy){
		set_accuracy(bubble_acc);
	}
	else{
		set_accuracy(tp.accuracy);
	}

	//Calculate the bubble coefficient including all subtractions
	complex<T> bubfullresult(0,0);
	BH_DEBUG_PRINT(typeid(T).name());
	BH_DEBUG_PRINT(RatBubSpecs::MUBUBPOINTS);
	for(int imupoints=0;imupoints<RatBubSpecs::MUBUBPOINTS-1;imupoints++){
//		_MESSAGE8("Bub at:",imupoints,"=",(amplm[imupoints*RatBubSpecs::YPOINTS]+mass_add[imupoints])/T(RatBubSpecs::YPOINTS),", ",-tp.tri_sub_pieces[imupoints]," * ",((static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,tp.S1)));
		BH_DEBUG_PRINT(bubfullresult);
		bubfullresult+=(amplm[imupoints*RatBubSpecs::YPOINTS]/T(RatBubSpecs::YPOINTS)+mass_add[imupoints]-tp.tri_sub_pieces[imupoints])*(static_cast<typename RatBubSpecs::CornerTreeStrategy*>(this))->get_rat_integral(imupoints,tp.S1);
	}
	set_B0coeff(bubfullresult);

	//We have already calculated this so return the value
	return bubfullresult;
}




}

}
#endif /* BUBBLE_RATEXT_HPP_ */
