/*
 * R2q3g_r.cpp
 *
 *  Created on: 07-Nov-2008
 *      Author: daniel
 */
#include "R_2q3g_r_eval.h"
#include "eval_param.h" 


#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

#define _VERBOSE 0


using namespace std;

namespace BH  {


template <class T> complex<T> Rqgqgg_mpppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqgg : -++++ R");
#endif
	return
	- (complex<T>(0,-1)/complex<T>(2,0)*SPA(1,3)*(SPA(1,2)*SPA(1,3)*SPB(3,2) +
		       SPA(1,2)*SPA(1,4)*SPB(4,2) + SPA(1,3)*SPA(1,4)*SPB(4,3)))/
		   (SPA(1,2)*SPA(1,5)*SPA(2,3)*SPA(3,4)*SPA(4,5))
	;
}

template <class T> complex<T> Rqgqgg_pmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqgg : +---- R");
#endif
	return
	-(complex<T>(0,-1)/complex<T>(2,0)*SPB(1,3)*(SPB(1,2)*SPB(1,3)*SPA(3,2) +
		       SPB(1,2)*SPB(1,4)*SPA(4,2) + SPB(1,3)*SPB(1,4)*SPA(4,3)))/
		   (SPB(1,2)*SPB(1,5)*SPB(2,3)*SPB(3,4)*SPB(4,5))
	;
}

template <class T> complex<T> Rqgqgg_mmpmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
// obtained by using charge parity, parity, and (3.2) of hep-ph/9409393

	#if _VERBOSE
	_MESSAGE("R2q3g : --+-- R sl");
#endif
return
+ (complex<T>(0,-1)/complex<T>(2,0)*SPB(3,1)*(SPA(1,2)*SPB(3,1)*SPB(3,2) +
       SPA(1,5)*SPB(3,1)*SPB(5,3) + SPA(2,5)*SPB(3,2)*SPB(5,3)))/
   (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,1)*SPB(5,4))
	;
	}

template <class T> complex<T> Rqgqgg_ppmpp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
// obtained by using charge parity, parity, and (3.2) of hep-ph/9409393

	#if _VERBOSE
	_MESSAGE("R2q3g : --+-- R sl");
#endif
return
 (complex<T>(0,-1)/complex<T>(2,0)*SPA(3,1)*(SPB(1,2)*SPA(3,1)*SPA(3,2) +
       SPB(1,5)*SPA(3,1)*SPA(5,3) + SPB(2,5)*SPA(3,2)*SPA(5,3)))/
   (SPA(2,1)*SPA(3,2)*SPA(4,3)*SPA(5,1)*SPA(5,4))
	;
}


template <class T> complex<T> R2q3g0(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : -----");
#endif
	return complex<T>(0,0);
}



template <class T> complex<T> R2q3g_R_1005(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqgg_mpppp_R(ep,mpc);}   // qbm p qp p p
template <class T> complex<T> R2q3g_R_18(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqgg_pmmmm_R(ep,mpc);}   // qbp m qm m m
template <class T> complex<T> R2q3g_R_33  (const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqgg_mmpmm_R(ep,mpc);}   // qbm p qp p p
template <class T> complex<T> R2q3g_R_990(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqgg_ppmpp_R(ep, mpc);}


template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R2q3g_mpppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : -++++ R");
#endif
	return -complex<T>(0,1)*(-(ep.spa(i1,i2)*ep.spa(i1,i3)*ep.spb(i3,i2) +
		       ep.spa(i1,i4)*ep.spa(i1,i5)*
		        ep.spb(i5,i4))/
		    (complex<T>(2,0)*ep.spa(i1,i5)*ep.spa(i2,i3)*
		      ep.spa(i3,i4)*ep.spa(i4,i5))


		     );
}


template <int i1, int i2, int i3, int i4, int i5, class T> complex<T> R2q3g_pmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : +---- R");
#endif
	// one minus from inverting all helicities, and one minus because spa(i,j)-> spb(i,j) instead of spb(j,i) and odd number of spa/b
	return -complex<T>(0,+1)*( -(ep.spb(i1,i2)*ep.spb(i1,i3)*ep.spa(i3,i2) +
		       ep.spb(i1,i4)*ep.spb(i1,i5)*
		        ep.spa(i5,i4))/
		    (complex<T>(2,0)*ep.spb(i1,i5)*ep.spb(i2,i3)*
		      ep.spb(i3,i4)*ep.spb(i4,i5))


		 );
}




template <class T> complex<T> R2q3g_mpmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : -+--- R");
#endif
	return

	-complex<T>(0,1)*( -
	     (-(SPA(3,4)*SPB(3,2)*SPB(4,2)) -
	        SPA(1,5)*SPB(2,1)*SPB(5,2))/
	      (complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,1)*
	        SPB(5,4)))

		     ;
}



template <class T> complex<T> R2q3g_pmppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : -+--- R");
#endif
	return

	+ complex<T>(0,1)*( -
	     (-(SPB(3,4)*SPA(3,2)*SPA(4,2)) -
	        SPB(1,5)*SPA(2,1)*SPA(5,2))/
	      (complex<T>(2,0)*SPA(3,2)*SPA(4,3)*SPA(5,1)*
	        SPA(5,4)))

		     ;
}

template <class T> complex<T> R2q3g_pmppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q3g : -+--- R");
#endif
	return

	+ complex<T>(0,1)*(((SPB(4,5)*pow(SPA(4,2),2)*
	           SPA(5,2))/
	         (SPA(2,1)*SPA(3,2)*SPA(4,3)*
	           pow(SPA(5,4),2)) -
	        (SPB(1,3)*SPB(1,5))/
	         (SPB(1,2)*SPA(4,3)*SPA(5,4)) +
	        (SPB(3,4)*SPA(3,2)*SPA(4,1)*
	           SPA(4,2))/
	         (SPA(2,1)*pow(SPA(4,3),2)*
	           SPA(5,1)*SPA(5,4)))/complex<T>(3,0) )

		     ;
}





template <class T> complex<T> R2q3g_R_1017(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q3g_mpppp_R<0,1,2,3,4>(ep,mpc);}   // qm qp p p p
template <class T> complex<T> R2q3g_R_6   (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q3g_pmmmm_R<0,1,2,3,4>(ep,mpc);}   // qp qm m m m
template <class T> complex<T> R2q3g_R_9(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q3g_mpmmm_R(ep,mpc);}   // qm qp m m m
template <class T> complex<T> R2q3g_R_1014(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q3g_pmppp_R(ep,mpc);}   // qp qm p p p



#define _CASE_R2q3g_R_Ptr_eval(K) case K : return &R2q3g_R_ ## K

template <class T> complex<T> (*R2q3g_SLC_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {


	switch (hc) {

	_CASE_R2q3g_R_Ptr_eval(1017);
	_CASE_R2q3g_R_Ptr_eval(6);
	_CASE_R2q3g_R_Ptr_eval(9);
	_CASE_R2q3g_R_Ptr_eval(1014);
	_CASE_R2q3g_R_Ptr_eval(990);
	_CASE_R2q3g_R_Ptr_eval(33);
	_CASE_R2q3g_R_Ptr_eval(1005);
	_CASE_R2q3g_R_Ptr_eval(18);

	default: return 0;
	}
}



 template complex<R> (*R2q3g_SLC_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R2q3g_SLC_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R2q3g_SLC_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R2q3g_SLC_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif





}

