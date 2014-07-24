/*
 * R_2q4g.cpp
 *
 *  Created on: Jul 26, 2008
 *      Author: daniel
 */

#include "R_2q4g_eval.h"
#include "eval_param.h"

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)


using namespace std;
namespace BH  {

template <class T> complex<T> Rqgqggg_mppppp_sl(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : -+++++ L");
#endif
	return

	(complex<T>(0,1)/complex<T>(6,0)*
		     (complex<T>(2,0)*ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
		        pow(ep.spa(4,5),2)*ep.spb(4,3)*
		        pow(ep.spa(0,3)*ep.spb(5,3) +
		          ep.spa(0,4)*ep.spb(5,4),2) -
		       ep.s(3,4)*ep.s(0,1,2)*
		        (complex<T>(3,0)*pow(ep.spa(0,2),2)*ep.spa(3,4)*
		           ep.spa(4,5)*
		           (ep.spa(0,3)*ep.spb(3,2) +
		             ep.spa(0,4)*ep.spb(4,2)) +
		          complex<T>(3,0)*ep.spa(0,2)*ep.spa(3,4)*
		           ep.spa(4,5)*
		           (ep.spa(0,3)*ep.spa(0,4)*
		              ep.spb(4,3) -
		             ep.spa(0,1)*ep.spa(0,5)*
		              ep.spb(5,1)) +
		          complex<T>(2,0)*ep.spa(0,4)*
		           (ep.spa(0,3)*ep.spa(2,3)*
		              ep.spa(4,5)*
		              (ep.spa(0,4)*ep.spb(4,3) +
		                ep.spa(0,5)*ep.spb(5,3)) +
		             ep.spa(0,5)*
		              (ep.spa(0,5)*ep.spa(2,4)*
		                 ep.spa(3,4) +
		                ep.spa(0,4)*ep.spa(2,3)*
		                 ep.spa(4,5))*ep.spb(5,4)))))/
		   (ep.s(3,4)*ep.s(0,1,2)*ep.spa(0,1)*
		     ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
		     pow(ep.spa(3,4),2)*pow(ep.spa(4,5),2));
}

template <class T> complex<T> Rqgqggg_mppppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : -+++++ R");
#endif
	return
	(complex<T>(0,1)/complex<T>(2,0)*SPA(1,3)*(SPA(1,2)*SPA(1,3)*SPB(3,2) +
	       SPA(1,2)*SPA(1,4)*SPB(4,2) + SPA(1,3)*SPA(1,4)*SPB(4,3) +
	       SPA(1,2)*SPA(1,5)*SPB(5,2) + SPA(1,3)*SPA(1,5)*SPB(5,3) +
	       SPA(1,4)*SPA(1,5)*SPB(5,4)))/
	   (SPA(1,2)*SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6));

}

template <class T> complex<T> Rqgqggg_mppppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : -+++++ nf");
#endif
	return
	(complex<T>(0,1)/complex<T>(3,0)*
	     ((SPA(1,5)*pow(SPA(1,6),2)*SPA(3,5)*SPB(6,5))/SPA(5,6) +
	       (SPA(1,5)*SPA(3,4)*(SPA(1,4)*SPA(1,5)*SPB(5,4) +
	            SPA(1,4)*SPA(1,6)*SPB(6,4) + SPA(1,5)*SPA(1,6)*SPB(6,5)))/
	        SPA(4,5) - (SPA(3,4)*SPA(5,6)*
	          pow(SPA(1,4)*SPA(1,6)*SPB(6,4) + SPA(1,5)*SPA(1,6)*SPB(6,5),
	           2)*(SPA(1,6)*SPA(3,4)*SPB(6,4) + SPA(1,6)*SPA(3,5)*SPB(6,5))*
	          (-(SPA(1,2)*SPA(1,6)*pow(SPA(4,5),2)*SPB(4,2)*SPB(5,4)*
	               SPB(6,5)) - SPA(1,3)*SPA(1,6)*pow(SPA(4,5),2)*SPB(4,3)*
	             SPB(5,4)*SPB(6,5)))/
	        (S(4,5)*SPA(1,6)*SPA(4,5)*
	          (-(SPA(1,2)*SPA(4,5)*SPB(4,2)) - SPA(1,3)*SPA(4,5)*SPB(4,3))*
	          (-(SPA(1,2)*SPA(4,6)*SPB(4,2)) - SPA(1,3)*SPA(4,6)*SPB(4,3) -
	            SPA(1,2)*SPA(5,6)*SPB(5,2) - SPA(1,3)*SPA(5,6)*SPB(5,3))*
	          SPB(6,5)*(-(SPA(1,6)*SPA(3,4)*SPB(6,4)) -
	            SPA(1,6)*SPA(3,5)*SPB(6,5)))))/
	   (SPA(1,2)*SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6))	;

}

template <class T> complex<T> Rqgqggg_pmmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : +----- R");
#endif
	return

	(complex<T>(0,1)/complex<T>(-2,0)*SPB(3,1)*(SPA(2,3)*SPB(2,1)*SPB(3,1) +
	       SPA(2,4)*SPB(2,1)*SPB(4,1) + SPA(3,4)*SPB(3,1)*SPB(4,1) +
	       SPA(2,5)*SPB(2,1)*SPB(5,1) + SPA(3,5)*SPB(3,1)*SPB(5,1) +
	       SPA(4,5)*SPB(4,1)*SPB(5,1)))/
	   (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,1)*SPB(6,5))
	;

}

template <class T> complex<T> Rqgqggg_pmmmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : +----- nf");
#endif

	return (complex<T>(0,-1)/complex<T>(3,0)*
		     ((SPB(4,3)*SPB(5,1)*(SPA(4,5)*SPB(4,1)*SPB(5,1) +
		            SPA(4,6)*SPB(4,1)*SPB(6,1) + SPA(5,6)*SPB(5,1)*SPB(6,1)))/
		        SPB(5,4) + (SPA(5,6)*SPB(5,1)*SPB(5,3)*pow(SPB(6,1),2))/
		        SPB(6,5) - (SPB(4,3)*pow(SPA(4,6)*SPB(4,1)*SPB(6,1) +
		            SPA(5,6)*SPB(5,1)*SPB(6,1),2)*
		          (SPA(4,6)*SPB(4,3)*SPB(6,1) + SPA(5,6)*SPB(5,3)*SPB(6,1))*
		          (-(SPA(2,4)*SPA(4,5)*SPA(5,6)*SPB(2,1)*pow(SPB(5,4),2)*
		               SPB(6,1)) - SPA(3,4)*SPA(4,5)*SPA(5,6)*SPB(3,1)*
		             pow(SPB(5,4),2)*SPB(6,1))*SPB(6,5))/
		        (S(4,5)*SPA(5,6)*SPB(5,4)*
		          (-(SPA(2,4)*SPB(2,1)*SPB(5,4)) - SPA(3,4)*SPB(3,1)*SPB(5,4))*
		          SPB(6,1)*(-(SPA(4,6)*SPB(4,3)*SPB(6,1)) -
		            SPA(5,6)*SPB(5,3)*SPB(6,1))*
		          (-(SPA(2,4)*SPB(2,1)*SPB(6,4)) - SPA(3,4)*SPB(3,1)*SPB(6,4) -
		            SPA(2,5)*SPB(2,1)*SPB(6,5) - SPA(3,5)*SPB(3,1)*SPB(6,5)))))/
		   (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,1)*SPB(6,5))

;

}

template <class T> complex<T> Rqgqggg_pmmmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
return 	-Rqgqggg_pmmmmm_R(ep,mpc)- Rqgqggg_pmmmmm_nf(ep,mpc);
}

template <class T> complex<T> Rqgqggg_mmpmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : --+--- R");
#endif
	return
	(complex<T>(0,1)/complex<T>(2,0)*SPB(3,1)*(-(SPA(1,2)*SPB(3,1)*SPB(3,2)) -
	       SPA(1,5)*SPB(3,1)*SPB(5,3) - SPA(2,5)*SPB(3,2)*SPB(5,3) -
	       SPA(1,6)*SPB(3,1)*SPB(6,3) - SPA(2,6)*SPB(3,2)*SPB(6,3) -
	       SPA(5,6)*SPB(5,3)*SPB(6,3)))/
	   (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,1)*SPB(6,5))

	;

}

template <class T> complex<T> Rqgqggg_mmpmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : --+--- nf");
#endif
	return
	(complex<T>(0,-1)/complex<T>(3,0)*
	     ((SPA(4,5)*pow(SPB(4,3),2)*SPB(5,1)*SPB(5,3))/SPB(5,4) -
	       (SPB(5,3)*SPB(6,1)*(-(SPA(4,5)*SPB(4,3)*SPB(5,3)) -
	            SPA(4,6)*SPB(4,3)*SPB(6,3) - SPA(5,6)*SPB(5,3)*SPB(6,3)))/
	        SPB(6,5) + (SPB(5,4)*SPB(6,1)*
	          (-(SPA(4,5)*SPB(4,3)*SPB(5,1)) - SPA(4,6)*SPB(4,3)*SPB(6,1))*
	          pow(-(SPA(4,5)*SPB(4,3)*SPB(5,3)) -
	            SPA(4,6)*SPB(4,3)*SPB(6,3),2)*
	          (SPA(1,6)*SPA(4,5)*SPA(5,6)*SPB(3,1)*SPB(4,3)*
	             pow(SPB(6,5),2) +
	            SPA(2,6)*SPA(4,5)*SPA(5,6)*SPB(3,2)*SPB(4,3)*
	             pow(SPB(6,5),2)))/
	        (S(5,6)*SPA(4,5)*SPB(4,3)*
	          (SPA(4,5)*SPB(4,3)*SPB(5,1) + SPA(4,6)*SPB(4,3)*SPB(6,1))*
	          (-(SPA(1,5)*SPB(3,1)*SPB(5,4)) - SPA(2,5)*SPB(3,2)*SPB(5,4) -
	            SPA(1,6)*SPB(3,1)*SPB(6,4) - SPA(2,6)*SPB(3,2)*SPB(6,4))*
	          SPB(6,5)*(-(SPA(1,6)*SPB(3,1)*SPB(6,5)) -
	            SPA(2,6)*SPB(3,2)*SPB(6,5)))))/
	   (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,1)*SPB(6,5))

	;

}

template <class T> complex<T> Rqgqggg_mmpmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
return 	-Rqgqggg_mmpmmm_R(ep,mpc)- Rqgqggg_mmpmmm_nf(ep,mpc);
}



template <class T> complex<T> Rqgqggg_ppmppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : ++-+++ R");
#endif
	return

	(complex<T>(0,1)/complex<T>(-2,0)*SPA(1,3)*(-(SPA(1,3)*SPA(2,3)*SPB(2,1)) -
	       SPA(1,3)*SPA(3,5)*SPB(5,1) - SPA(2,3)*SPA(3,5)*SPB(5,2) -
	       SPA(1,3)*SPA(3,6)*SPB(6,1) - SPA(2,3)*SPA(3,6)*SPB(6,2) -
	       SPA(3,5)*SPA(3,6)*SPB(6,5)))/
	   (SPA(1,2)*SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6))
	;

}
template <class T> complex<T> Rqgqggg_ppmppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("Rqgqggg : ++-+++ nf");
#endif
	return

	(complex<T>(0,1)/complex<T>(3,0)*
	     ((SPA(1,5)*pow(SPA(3,4),2)*SPA(3,5)*SPB(5,4))/SPA(4,5) -
	       (SPA(1,6)*SPA(3,5)*(-(SPA(3,4)*SPA(3,5)*SPB(5,4)) -
	            SPA(3,4)*SPA(3,6)*SPB(6,4) - SPA(3,5)*SPA(3,6)*SPB(6,5)))/
	        SPA(5,6) + (SPA(1,6)*SPA(4,5)*
	          (-(SPA(1,5)*SPA(3,4)*SPB(5,4)) - SPA(1,6)*SPA(3,4)*SPB(6,4))*
	          pow(-(SPA(3,4)*SPA(3,5)*SPB(5,4)) -
	            SPA(3,4)*SPA(3,6)*SPB(6,4),2)*
	          (SPA(1,3)*SPA(3,4)*pow(SPA(5,6),2)*SPB(5,4)*SPB(6,1)*
	             SPB(6,5) + SPA(2,3)*SPA(3,4)*pow(SPA(5,6),2)*SPB(5,4)*
	             SPB(6,2)*SPB(6,5)))/
	        (S(5,6)*SPA(3,4)*SPA(5,6)*SPB(5,4)*
	          (-(SPA(1,3)*SPA(4,5)*SPB(5,1)) - SPA(2,3)*SPA(4,5)*SPB(5,2) -
	            SPA(1,3)*SPA(4,6)*SPB(6,1) - SPA(2,3)*SPA(4,6)*SPB(6,2))*
	          (-(SPA(1,3)*SPA(5,6)*SPB(6,1)) - SPA(2,3)*SPA(5,6)*SPB(6,2))*
	          (SPA(1,5)*SPA(3,4)*SPB(5,4) + SPA(1,6)*SPA(3,4)*SPB(6,4)))))/
	   (SPA(1,2)*SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6))
	;

}

template <class T> complex<T> Rqgqggg_ppmppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
return 	-Rqgqggg_ppmppp_R(ep,mpc)- Rqgqggg_ppmppp_nf(ep,mpc);
}



template <class T> complex<T> R2q4g_mppppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ R");
#endif
	return
	-(complex<T>(0,-1)/complex<T>(2,0)*(SPA(1,2)*SPA(1,3)*SPB(3,2) +
	       SPA(1,2)*SPA(1,4)*SPB(4,2) + SPA(1,3)*SPA(1,4)*SPB(4,3) +
	       SPA(1,2)*SPA(1,5)*SPB(5,2) + SPA(1,3)*SPA(1,5)*SPB(5,3) +
	       SPA(1,4)*SPA(1,5)*SPB(5,4)))/
	   (SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6)) ;
}

template <class T> complex<T> R2q4g_mppppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ nf");
#endif
	return - complex<T>(0,1)/complex<T>(3,0)*
	   (-(pow(SPA(1,3)*SPB(3,2) + SPA(1,4)*SPB(4,2),2)/
	        (SPA(1,6)*pow(SPA(3,4),2)*SPA(5,6)*
	          (-(SPA(3,5)*SPB(3,2)) - SPA(4,5)*SPB(4,2)))) -
	     (SPA(1,4)*(SPA(1,2)*SPA(1,3)*SPB(3,2) +
	          SPA(1,2)*SPA(1,4)*SPB(4,2) + SPA(1,3)*SPA(1,4)*SPB(4,3)))/
	      (SPA(1,2)*SPA(1,6)*pow(SPA(3,4),2)*SPA(4,5)*SPA(5,6)) -
	     (SPA(1,5)*SPA(1,6)*SPA(2,5)*SPB(6,5))/
	      (SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*pow(SPA(5,6),2)) -
	     (SPA(1,5)*SPA(2,4)*(SPA(1,4)*SPA(1,5)*SPB(5,4) +
	          SPA(1,4)*SPA(1,6)*SPB(6,4) + SPA(1,5)*SPA(1,6)*SPB(6,5)))/
	      (SPA(1,2)*SPA(1,6)*SPA(2,3)*SPA(3,4)*pow(SPA(4,5),2)*SPA(5,6)) +
	     (pow(-(SPA(1,4)*SPB(6,4)) - SPA(1,5)*SPB(6,5),2)*
	        (-(SPA(2,4)*SPB(6,4)) - SPA(2,5)*SPB(6,5)))/
	      (SS(4,5,6)*SPA(1,2)*SPA(2,3)*pow(SPA(4,5),2)*
	        (-(SPA(3,4)*SPB(6,4)) - SPA(3,5)*SPB(6,5))) +
	     (pow(SPB(6,2),2)*(pow(SPA(3,4),2)*SPB(3,2)*SPB(4,3)*SPB(6,4) +
	          SPA(3,4)*SPA(3,5)*SPB(3,2)*SPB(5,3)*SPB(6,4) +
	          SPA(3,4)*SPA(4,5)*SPB(4,2)*SPB(5,3)*SPB(6,4) +
	          SPA(3,4)*SPA(3,5)*SPB(3,2)*SPB(4,3)*SPB(6,5) +
	          pow(SPA(3,5),2)*SPB(3,2)*SPB(5,3)*SPB(6,5) +
	          SPA(3,5)*SPA(4,5)*SPB(4,2)*SPB(5,3)*SPB(6,5) +
	          SPA(3,5)*SPA(4,5)*SPB(3,2)*SPB(5,4)*SPB(6,5) +
	          pow(SPA(4,5),2)*SPB(4,2)*SPB(5,4)*SPB(6,5)))/
	      (SS(3,4,5)*SPA(3,4)*SPA(4,5)*SPB(2,1)*
	        (-(SPA(3,5)*SPB(3,2)) - SPA(4,5)*SPB(4,2))*
	        (-(SPA(3,4)*SPB(6,4)) - SPA(3,5)*SPB(6,5))))
		     ;
}

template <class T> complex<T> R2q4g_mppppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ L");
#endif
	return - R2q4g_mppppp_nf(ep,mpc) - R2q4g_mppppp_R(ep,mpc)
		     ;
}



template <class T> complex<T> R2q4g_pmmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ R");
#endif
	return
	-(complex<T>(0,-1)/complex<T>(2,0)*(SPB(1,2)*SPB(1,3)*SPA(3,2) +
	       SPB(1,2)*SPB(1,4)*SPA(4,2) + SPB(1,3)*SPB(1,4)*SPA(4,3) +
	       SPB(1,2)*SPB(1,5)*SPA(5,2) + SPB(1,3)*SPB(1,5)*SPA(5,3) +
	       SPB(1,4)*SPB(1,5)*SPA(5,4)))/
	   (SPB(1,6)*SPB(2,3)*SPB(3,4)*SPB(4,5)*SPB(5,6)) ;
}

template <class T> complex<T> R2q4g_pmmmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ nf");
#endif
	return - complex<T>(0,1)/complex<T>(3,0)*
	   (-(pow(SPB(1,3)*SPA(3,2) + SPB(1,4)*SPA(4,2),2)/
	        (SPB(1,6)*pow(SPB(3,4),2)*SPB(5,6)*
	          (-(SPB(3,5)*SPA(3,2)) - SPB(4,5)*SPA(4,2)))) -
	     (SPB(1,4)*(SPB(1,2)*SPB(1,3)*SPA(3,2) +
	          SPB(1,2)*SPB(1,4)*SPA(4,2) + SPB(1,3)*SPB(1,4)*SPA(4,3)))/
	      (SPB(1,2)*SPB(1,6)*pow(SPB(3,4),2)*SPB(4,5)*SPB(5,6)) -
	     (SPB(1,5)*SPB(1,6)*SPB(2,5)*SPA(6,5))/
	      (SPB(1,2)*SPB(2,3)*SPB(3,4)*SPB(4,5)*pow(SPB(5,6),2)) -
	     (SPB(1,5)*SPB(2,4)*(SPB(1,4)*SPB(1,5)*SPA(5,4) +
	          SPB(1,4)*SPB(1,6)*SPA(6,4) + SPB(1,5)*SPB(1,6)*SPA(6,5)))/
	      (SPB(1,2)*SPB(1,6)*SPB(2,3)*SPB(3,4)*pow(SPB(4,5),2)*SPB(5,6)) +
	     (pow(-(SPB(1,4)*SPA(6,4)) - SPB(1,5)*SPA(6,5),2)*
	        (-(SPB(2,4)*SPA(6,4)) - SPB(2,5)*SPA(6,5)))/
	      (SS(4,5,6)*SPB(1,2)*SPB(2,3)*pow(SPB(4,5),2)*
	        (-(SPB(3,4)*SPA(6,4)) - SPB(3,5)*SPA(6,5))) +
	     (pow(SPA(6,2),2)*(pow(SPB(3,4),2)*SPA(3,2)*SPA(4,3)*SPA(6,4) +
	          SPB(3,4)*SPB(3,5)*SPA(3,2)*SPA(5,3)*SPA(6,4) +
	          SPB(3,4)*SPB(4,5)*SPA(4,2)*SPA(5,3)*SPA(6,4) +
	          SPB(3,4)*SPB(3,5)*SPA(3,2)*SPA(4,3)*SPA(6,5) +
	          pow(SPB(3,5),2)*SPA(3,2)*SPA(5,3)*SPA(6,5) +
	          SPB(3,5)*SPB(4,5)*SPA(4,2)*SPA(5,3)*SPA(6,5) +
	          SPB(3,5)*SPB(4,5)*SPA(3,2)*SPA(5,4)*SPA(6,5) +
	          pow(SPB(4,5),2)*SPA(4,2)*SPA(5,4)*SPA(6,5)))/
	      (SS(3,4,5)*SPB(3,4)*SPB(4,5)*SPA(2,1)*
	        (-(SPB(3,5)*SPA(3,2)) - SPB(4,5)*SPA(4,2))*
	        (-(SPB(3,4)*SPA(6,4)) - SPB(3,5)*SPA(6,5))))
		     ;
}

template <class T> complex<T> R2q4g_pmmmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -++++ L");
#endif
	return - R2q4g_pmmmmm_nf(ep,mpc) - R2q4g_pmmmmm_R(ep,mpc)
		     ;
}

template <class T> complex<T> R2q4g_mpmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -+---- R");
#endif
	return
	(complex<T>(0,1)/complex<T>(2,0)*(-(SPA(1,4)*SPB(2,1)*SPB(4,2)) -
	        SPA(1,5)*SPB(2,1)*SPB(5,2) - SPA(4,5)*SPB(4,2)*SPB(5,2) -
	        SPA(1,6)*SPB(2,1)*SPB(6,2) - SPA(4,6)*SPB(4,2)*SPB(6,2) -
	        SPA(5,6)*SPB(5,2)*SPB(6,2)))/
	    (SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,1)*SPB(6,5)) 	;
}

template <class T> complex<T> R2q4g_mpmmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -+---- nf");
#endif
	return -
	   complex<T>(0,1)/complex<T>(3,0)*
	    ((SPB(5,2)*(-(SPA(1,5)*SPB(2,1)*SPB(5,2)) -
	           SPA(1,6)*SPB(2,1)*SPB(6,2) - SPA(5,6)*SPB(5,2)*SPB(6,2)))/
	       (SPB(2,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*pow(SPB(6,5),2)) +
	      pow(SPA(1,5)*SPB(5,2) + SPA(1,6)*SPB(6,2),2)/
	       (SPB(3,2)*SPB(4,3)*(SPA(1,5)*SPB(5,4) + SPA(1,6)*SPB(6,4))*
	         pow(SPB(6,5),2)) +
	      (SPB(4,2)*SPB(5,1)*(-(SPA(3,4)*SPB(3,2)*SPB(4,2)) -
	           SPA(3,5)*SPB(3,2)*SPB(5,2) - SPA(4,5)*SPB(4,2)*SPB(5,2)))/
	       (SPB(2,1)*SPB(3,2)*SPB(4,3)*pow(SPB(5,4),2)*SPB(6,1)*SPB(6,5))\
	       - (SPA(3,4)*SPB(3,2)*SPB(4,1)*SPB(4,2))/
	       (SPB(2,1)*pow(SPB(4,3),2)*SPB(5,4)*SPB(6,1)*SPB(6,5)) -
	      ((SPA(3,4)*SPB(4,1) + SPA(3,5)*SPB(5,1))*
	         pow(SPA(3,4)*SPB(4,2) + SPA(3,5)*SPB(5,2),2))/
	       (SS(3,4,5)*SPB(2,1)*pow(SPB(5,4),2)*SPB(6,1)*
	         (-(SPA(3,4)*SPB(6,4)) - SPA(3,5)*SPB(6,5))) -
	      (pow(SPA(1,3),2)*(SPA(1,5)*SPA(3,4)*SPA(4,5)*
	            pow(SPB(5,4),2) +
	           SPA(1,6)*SPA(3,4)*SPA(4,5)*SPB(5,4)*SPB(6,4) +
	           SPA(1,5)*SPA(3,4)*SPA(4,6)*SPB(5,4)*SPB(6,4) +
	           SPA(1,6)*SPA(3,4)*SPA(4,6)*pow(SPB(6,4),2) +
	           SPA(1,5)*SPA(3,5)*SPA(4,6)*SPB(5,4)*SPB(6,5) +
	           SPA(1,6)*SPA(3,5)*SPA(4,6)*SPB(6,4)*SPB(6,5) +
	           SPA(1,6)*SPA(3,4)*SPA(5,6)*SPB(6,4)*SPB(6,5) +
	           SPA(1,6)*SPA(3,5)*SPA(5,6)*pow(SPB(6,5),2)))/
	       (SS(4,5,6)*SPA(1,2)*SPB(5,4)*
	         (SPA(1,5)*SPB(5,4) + SPA(1,6)*SPB(6,4))*SPB(6,5)*
	         (-(SPA(3,4)*SPB(6,4)) - SPA(3,5)*SPB(6,5))))
		     ;
}

template <class T> complex<T> R2q4g_mpmmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : -+---- L");
#endif
	return - R2q4g_mpmmmm_nf(ep,mpc) - R2q4g_mpmmmm_R(ep,mpc)
		     ;
}


template <class T> complex<T> R2q4g_pmpppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : +-++++ R");
#endif
	return
	(complex<T>(0,1)/complex<T>(2,0)*(-(SPB(1,4)*SPA(2,1)*SPA(4,2)) -
	        SPB(1,5)*SPA(2,1)*SPA(5,2) - SPB(4,5)*SPA(4,2)*SPA(5,2) -
	        SPB(1,6)*SPA(2,1)*SPA(6,2) - SPB(4,6)*SPA(4,2)*SPA(6,2) -
	        SPB(5,6)*SPA(5,2)*SPA(6,2)))/
	    (SPA(3,2)*SPA(4,3)*SPA(5,4)*SPA(6,1)*SPA(6,5)) 	;
}

template <class T> complex<T> R2q4g_pmpppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : +-++++ nf");
#endif
	return -
	   complex<T>(0,1)/complex<T>(3,0)*
	    ((SPA(5,2)*(-(SPB(1,5)*SPA(2,1)*SPA(5,2)) -
	           SPB(1,6)*SPA(2,1)*SPA(6,2) - SPB(5,6)*SPA(5,2)*SPA(6,2)))/
	       (SPA(2,1)*SPA(3,2)*SPA(4,3)*SPA(5,4)*pow(SPA(6,5),2)) +
	      pow(SPB(1,5)*SPA(5,2) + SPB(1,6)*SPA(6,2),2)/
	       (SPA(3,2)*SPA(4,3)*(SPB(1,5)*SPA(5,4) + SPB(1,6)*SPA(6,4))*
	         pow(SPA(6,5),2)) +
	      (SPA(4,2)*SPA(5,1)*(-(SPB(3,4)*SPA(3,2)*SPA(4,2)) -
	           SPB(3,5)*SPA(3,2)*SPA(5,2) - SPB(4,5)*SPA(4,2)*SPA(5,2)))/
	       (SPA(2,1)*SPA(3,2)*SPA(4,3)*pow(SPA(5,4),2)*SPA(6,1)*SPA(6,5))\
	       - (SPB(3,4)*SPA(3,2)*SPA(4,1)*SPA(4,2))/
	       (SPA(2,1)*pow(SPA(4,3),2)*SPA(5,4)*SPA(6,1)*SPA(6,5)) -
	      ((SPB(3,4)*SPA(4,1) + SPB(3,5)*SPA(5,1))*
	         pow(SPB(3,4)*SPA(4,2) + SPB(3,5)*SPA(5,2),2))/
	       (SS(3,4,5)*SPA(2,1)*pow(SPA(5,4),2)*SPA(6,1)*
	         (-(SPB(3,4)*SPA(6,4)) - SPB(3,5)*SPA(6,5))) -
	      (pow(SPB(1,3),2)*(SPB(1,5)*SPB(3,4)*SPB(4,5)*
	            pow(SPA(5,4),2) +
	           SPB(1,6)*SPB(3,4)*SPB(4,5)*SPA(5,4)*SPA(6,4) +
	           SPB(1,5)*SPB(3,4)*SPB(4,6)*SPA(5,4)*SPA(6,4) +
	           SPB(1,6)*SPB(3,4)*SPB(4,6)*pow(SPA(6,4),2) +
	           SPB(1,5)*SPB(3,5)*SPB(4,6)*SPA(5,4)*SPA(6,5) +
	           SPB(1,6)*SPB(3,5)*SPB(4,6)*SPA(6,4)*SPA(6,5) +
	           SPB(1,6)*SPB(3,4)*SPB(5,6)*SPA(6,4)*SPA(6,5) +
	           SPB(1,6)*SPB(3,5)*SPB(5,6)*pow(SPA(6,5),2)))/
	       (SS(4,5,6)*SPB(1,2)*SPA(5,4)*
	         (SPB(1,5)*SPA(5,4) + SPB(1,6)*SPA(6,4))*SPA(6,5)*
	         (-(SPB(3,4)*SPA(6,4)) - SPB(3,5)*SPA(6,5))))
		     ;
}


template <class T> complex<T> R2q4g_pmpppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R2q4g qqgggg : +-++++ L");
#endif
	return - R2q4g_pmpppp_nf(ep,mpc) - R2q4g_pmpppp_R(ep,mpc)
		     ;
}


template <class T> complex<T> Rqgqggg_mmpppp_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{

return

-complex<T>(0,1) * (
		(pow(SS(1,2,3),2)*pow(SPA(2,4),2)*SPB(4,3))/
    (complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6)*
      (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))*SPB(3,2)*
      (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))) +
   (SS(1,2,3)*SPB(3,1)*(SPA(1,4)*SPB(3,1) + SPA(2,4)*SPB(3,2))*SPB(4,3))/
    (complex<T>(2,0)*SPA(4,5)*SPA(5,6)*SPB(2,1)*
      (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))*SPB(3,2)*
      (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))) +
   (SPB(3,1)*(SPA(1,5)*SPB(3,1) + SPA(2,5)*SPB(3,2))*SPB(6,5))/
    (complex<T>(2,0)*SPA(4,5)*SPA(5,6)*SPB(2,1)*
      (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))*SPB(3,2)) +
   ((SPA(2,3)*SPB(3,1) + SPA(2,4)*SPB(4,1))*
      (-(SPA(2,3)*SPB(6,3)) - SPA(2,4)*SPB(6,4))*SPB(6,5))/
    (complex<T>(2,0)*SS(2,3,4)*SPA(3,4)*SPA(5,6)*
      (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))*SPB(6,1)) -
   (pow(SPA(2,5),2)*pow(SPB(6,5),2))/
    (complex<T>(2,0)*SS(2,3,4)*SPA(3,4)*SPA(4,5)*SPA(5,6)*SPB(6,1)) -
   (pow(SPA(1,2),2)*pow(SPA(2,4),2)*pow(SPB(5,4),2))/
    (complex<T>(2,0)*SPA(1,6)*SPA(2,3)*SPA(3,4)*SPA(4,5)*
      (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))*
      (SPA(1,2)*SPB(5,1) + SPA(2,6)*SPB(6,5))) +
   (pow(SPA(1,2),2)*SPB(5,3)*SPB(5,4)*
      (SPA(1,2)*SPB(3,1) + SPA(2,6)*SPB(6,3)))/
    (complex<T>(2,0)*SS(3,4,5)*SPA(1,6)*SPA(4,5)*
      (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))*
      (SPA(1,2)*SPB(5,1) + SPA(2,6)*SPB(6,5))) -
   (pow(SPA(1,2),2)*(pow(-(SPA(1,5)*SPA(2,4)*SPB(5,4)) -
           SPA(1,6)*SPA(2,4)*SPB(6,4),2)/
         (S(2,3)*(S(2,3) - SS(1,5,6))*pow(SPA(1,2),2)) +
        (pow(SPA(1,6),2)*pow(SPA(2,5),2)*pow(SPB(6,5),2))/
         (SS(1,5,6)*(-S(1,6) + SS(1,5,6))*pow(SPA(1,2),2))))/
    (complex<T>(2,0)*SPA(1,6)*SPA(3,4)*SPA(4,5)*SPA(5,6))

    )

    ;
}

template <class T> complex<T> Rqgqggg_mmpppp_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{

return - complex<T>(0,1) * (
		((pow(SPA(1,2),2)*SPB(4,3)*SPB(5,3))/
		      (SS(3,4,5)*SPA(1,6)*SPA(4,5)*
		        (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))) -
		     ((SPA(1,4)*SPB(3,1) + SPA(2,4)*SPB(3,2))*
		        pow(SPA(1,5)*SPB(3,1) + SPA(2,5)*SPB(3,2),2)*SPB(5,4))/
		      (SS(1,2,3)*pow(SPA(4,5),2)*SPA(5,6)*SPB(2,1)*SPB(3,2)*
		        (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))) +
		     (SPB(4,3)*SPB(6,3))/(SPA(4,5)*SPA(5,6)*SPB(2,1)*SPB(3,2)) -
		     ((-(SPA(2,3)*SPB(5,3)) - SPA(2,4)*SPB(5,4))*
		        (-(SPA(2,3)*SPB(6,3)) - SPA(2,4)*SPB(6,4)))/
		      (SS(2,3,4)*SPA(3,4)*SPA(5,6)*
		        (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))) -
		     ((-(SPA(2,5)*SPB(2,1)) - SPA(3,5)*SPB(3,1))*
		        (SPA(1,5)*SPB(3,1) + SPA(2,5)*SPB(3,2))*
		        (SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))*SPB(6,5))/
		      (SS(1,2,3)*SPA(4,5)*pow(SPA(5,6),2)*SPB(2,1)*
		        (-(SPA(2,4)*SPB(2,1)) - SPA(3,4)*SPB(3,1))*SPB(3,2)))/complex<T>(3,0)
)
;

}

template <class T> complex<T> Rqgqggg_mmpppp_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
return 	-Rqgqggg_mmpppp_R(ep,mpc)- Rqgqggg_mmpppp_nf(ep,mpc);
}



template <class T> complex<T> Rqgqggg_ppmmmm_R(const eval_param<T>& ep,const mass_param_coll& mpc)
{

return

- complex<T>(0,1) * (
		(pow(SS(1,2,3),2)*pow(SPB(4,2),2)*SPA(3,4))/
    (complex<T>(2,0)*SPB(3,2)*SPB(4,3)*SPB(5,4)*SPB(6,5)*
      (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))*SPA(2,3)*
      (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))) +
   (SS(1,2,3)*SPA(1,3)*(SPB(4,1)*SPA(1,3) + SPB(4,2)*SPA(2,3))*SPA(3,4))/
    (complex<T>(2,0)*SPB(5,4)*SPB(6,5)*SPA(1,2)*
      (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))*SPA(2,3)*
      (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))) +
   (SPA(1,3)*(SPB(5,1)*SPA(1,3) + SPB(5,2)*SPA(2,3))*SPA(5,6))/
    (complex<T>(2,0)*SPB(5,4)*SPB(6,5)*SPA(1,2)*
      (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))*SPA(2,3)) +
   ((SPB(3,2)*SPA(1,3) + SPB(4,2)*SPA(1,4))*
      (-(SPB(3,2)*SPA(3,6)) - SPB(4,2)*SPA(4,6))*SPA(5,6))/
    (complex<T>(2,0)*SS(2,3,4)*SPB(4,3)*SPB(6,5)*
      (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))*SPA(1,6)) -
   (pow(SPB(5,2),2)*pow(SPA(5,6),2))/
    (complex<T>(2,0)*SS(2,3,4)*SPB(4,3)*SPB(5,4)*SPB(6,5)*SPA(1,6)) -
   (pow(SPB(2,1),2)*pow(SPB(4,2),2)*pow(SPA(4,5),2))/
    (complex<T>(2,0)*SPB(6,1)*SPB(3,2)*SPB(4,3)*SPB(5,4)*
      (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))*
      (SPB(2,1)*SPA(1,5) + SPB(6,2)*SPA(5,6))) +
   (pow(SPB(2,1),2)*SPA(3,5)*SPA(4,5)*
      (SPB(2,1)*SPA(1,3) + SPB(6,2)*SPA(3,6)))/
    (complex<T>(2,0)*SS(3,4,5)*SPB(6,1)*SPB(5,4)*
      (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))*
      (SPB(2,1)*SPA(1,5) + SPB(6,2)*SPA(5,6))) -
   (pow(SPB(2,1),2)*(pow(-(SPB(5,1)*SPB(4,2)*SPA(4,5)) -
           SPB(6,1)*SPB(4,2)*SPA(4,6),2)/
         (S(2,3)*(S(2,3) - SS(1,5,6))*pow(SPB(2,1),2)) +
        (pow(SPB(6,1),2)*pow(SPB(5,2),2)*pow(SPA(5,6),2))/
         (SS(1,5,6)*(-S(1,6) + SS(1,5,6))*pow(SPB(2,1),2))))/
    (complex<T>(2,0)*SPB(6,1)*SPB(4,3)*SPB(5,4)*SPB(6,5))

    )

    ;
}

template <class T> complex<T> Rqgqggg_ppmmmm_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{

return - complex<T>(0,1) * (
		((pow(SPB(2,1),2)*SPA(3,4)*SPA(3,5))/
		      (SS(3,4,5)*SPB(6,1)*SPB(5,4)*
		        (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))) -
		     ((SPB(4,1)*SPA(1,3) + SPB(4,2)*SPA(2,3))*
		        pow(SPB(5,1)*SPA(1,3) + SPB(5,2)*SPA(2,3),2)*SPA(4,5))/
		      (SS(1,2,3)*pow(SPB(5,4),2)*SPB(6,5)*SPA(1,2)*SPA(2,3)*
		        (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))) +
		     (SPA(3,4)*SPA(3,6))/(SPB(5,4)*SPB(6,5)*SPA(1,2)*SPA(2,3)) -
		     ((-(SPB(3,2)*SPA(3,5)) - SPB(4,2)*SPA(4,5))*
		        (-(SPB(3,2)*SPA(3,6)) - SPB(4,2)*SPA(4,6)))/
		      (SS(2,3,4)*SPB(4,3)*SPB(6,5)*
		        (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))) -
		     ((-(SPB(5,2)*SPA(1,2)) - SPB(5,3)*SPA(1,3))*
		        (SPB(5,1)*SPA(1,3) + SPB(5,2)*SPA(2,3))*
		        (SPB(6,1)*SPA(1,3) + SPB(6,2)*SPA(2,3))*SPA(5,6))/
		      (SS(1,2,3)*SPB(5,4)*pow(SPB(6,5),2)*SPA(1,2)*
		        (-(SPB(4,2)*SPA(1,2)) - SPB(4,3)*SPA(1,3))*SPA(2,3)))/complex<T>(3,0)
)
;

}

template <class T> complex<T> Rqgqggg_ppmmmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
return 	-Rqgqggg_ppmmmm_R(ep,mpc)- Rqgqggg_ppmmmm_nf(ep,mpc);
}


template <class T> complex<T> Rqgqggg_mmppmm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
// added by Harald { qbm m qp p m m }_L
return -((complex<T>(0,1)*SPA(2,1)*SPA(5,3)*(-(SPA(5,3)*SPB(1,3))-
SPA(5,4)*SPB(1,4))*((SPA(6,5)*SPB(4,5)*pow(SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6),2))/(complex<T>(3,0)*SPA(2,1)*(-(SPA(5,3)*SPB(1,3))-
SPA(5,4)*SPB(1,4))*SPB(5,6)*pow(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6),2))+
(SPA(6,5)*SPB(4,6)*(SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6)))/(complex<T>(2,0)*SPA(2,1)*SPA(5,3)*SPB(1,6)*SPB(5,6)*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6)))))/(SPA(5,4)*SPB(1,2)*(-(SPA(5,3)*SPB(2,3))-
SPA(5,4)*SPB(2,4)))-
((complex<T>(0,-1)*SPB(1,4)*pow(SPB(3,4),2))/(complex<T>(2,0)*SPB(1,2)*SPB(1,6)*SPB(4,5)*SPB(5,6))-
(complex<T>(0,1)*SPA(2,1)*SPB(3,4)*((pow(SPA(6,5),2)*SPB(2,6)*SPB(3,4)*SPB(4,6))/(complex<T>(2,0)*SPA(2,1)*SPB(1,6)*(-(SPA(5,3)*SPB(2,3))-
SPA(5,4)*SPB(2,4))*SPB(5,6))-
(pow(SPA(6,5),2)*SPB(2,4)*SPB(2,6)*SPB(3,4)*SPB(4,5))/(SPA(2,1)*pow(SPB(1,2),2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*SPB(5,6))-
(pow(SPA(6,5),3)*SPB(2,4)*pow(SPB(2,6),2)*SPB(3,4)*SPB(4,5)*(-(complex<T>(2,0)*SPA(5,4)*SPB(4,5))-
(complex<T>(2,0)*SPA(5,3)*SPB(2,3)*SPB(4,5))/SPB(2,4)-
SPA(6,4)*SPB(4,6)-
(SPA(6,3)*SPB(2,3)*SPB(4,6))/SPB(2,4)-
SPA(6,5)*SPB(5,6)))/(complex<T>(3,0)*SPA(2,1)*pow(SPB(1,2),2)*(-(SPA(5,3)*SPB(2,3))-
SPA(5,4)*SPB(2,4))*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*SPB(5,6)*(SPA(6,4)*SPB(4,6)+
(SPA(6,3)*SPB(2,3)*SPB(4,6))/SPB(2,4)+
SPA(6,5)*SPB(5,6)))))/(SPB(2,4)*(SPA(6,4)*SPB(4,6)+
(SPA(6,3)*SPB(2,3)*SPB(4,6))/SPB(2,4)+
SPA(6,5)*SPB(5,6))))/SPB(2,3)+
(pow(SPB(3,4),2)*((complex<T>(0,-1)*SPA(3,2)*SPB(1,4)*pow(SPB(2,4),2))/(complex<T>(2,0)*SPB(1,2)*SPB(1,6)*SPB(4,5)*SPB(5,6))+
complex<T>(0,1)*((-((SPA(3,2)*pow(SPA(6,5),2)*pow(SPB(2,4),2)*SPB(2,6)*SPB(4,6))/(SPB(1,6)*pow(SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6),2)*SPB(4,5)*(complex<T>(1,0)-
(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6)))/((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(4,5)))*SPB(5,6)))+
(SPA(3,2)*SPA(5,1)*SPA(6,5)*pow(SPB(2,4),2))/((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(5,6)*SS(1,5,6)))/complex<T>(2,0)+
((complex<T>(2,0)*SPA(3,2)*SPA(6,5)*pow(SPB(2,4),3))/(pow(SPB(1,2),2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*SPB(5,6))+
(complex<T>(3,0)*SPA(3,2)*pow(SPA(6,5),2)*pow(SPB(2,4),3)*SPB(2,6)*SPB(4,5))/(pow(SPB(1,2),3)*pow(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6),2)*(complex<T>(1,0)-
((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(4,5))/(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))))*SPB(5,6))+
(SPA(3,2)*pow(SPA(6,5),3)*pow(SPB(2,4),3)*pow(SPB(2,6),2)*pow(SPB(4,5),2)*(-(((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(4,5))/(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))))+
(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6)))/((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(4,5))))/(pow(SPB(1,2),4)*pow(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6),3)*pow(complex<T>(1,0)-
((SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*SPB(4,5))/(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))),3)*SPB(5,6))-
(SPA(3,2)*SPA(5,1)*pow(SPB(2,4),2)*(SPA(6,1)*SPB(1,4)-
SPA(6,5)*SPB(4,5)))/(SPB(1,2)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*SPB(5,6)*SS(1,5,6))-
(SPA(3,2)*SPA(5,1)*SPA(6,5)*pow(SPB(2,4),2)*SPB(2,5)*(SPA(5,1)*SPB(1,4)+
SPA(6,5)*SPB(4,6)))/(SPB(1,2)*(SPA(5,1)*SPB(1,2)+
SPA(6,5)*SPB(2,6))*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*SPB(5,6)*SS(1,5,6)))/complex<T>(3,0))))/(SPA(3,2)*SPB(2,3)*pow(SPB(2,4),2))-
((complex<T>(0,1)*(-((SPA(5,3)*SPA(6,5)*(-(SPA(5,3)*SPB(1,3))-
SPA(5,4)*SPB(1,4))*(SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6)))/(SPB(1,2)*SPB(1,6)*(-(SPA(5,3)*SPB(2,3))-
SPA(5,4)*SPB(2,4))*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6))))-
(SPA(2,1)*SPA(5,3)*pow(-(SPA(5,3)*SPB(1,3))-
SPA(5,4)*SPB(1,4),2))/(SPB(1,2)*SPB(1,6)*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6))*SS(3,4,5))))/complex<T>(2,0)-
(complex<T>(0,1)*pow(SPA(5,3),2)*SPA(6,5)*(SPA(5,3)*SPB(3,6)+
SPA(5,4)*SPB(4,6))*SS(3,4,5))/(complex<T>(3,0)*SPB(1,2)*(-(SPA(5,3)*SPB(2,3))-
SPA(5,4)*SPB(2,4))*pow(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6),2)))/(SPA(4,3)*SPA(5,4))-
(complex<T>(0,1)*pow(SPB(1,3),2)*SPB(1,4)*pow(-(SPA(2,1)*SPB(2,4))-
SPA(3,1)*SPB(3,4),2))/(complex<T>(2,0)*SPB(1,2)*(SPA(2,1)*SPB(1,2)+
SPA(3,1)*SPB(1,3))*SPB(1,6)*SPB(2,3)*SPB(4,5)*SPB(5,6)*SS(4,5,6))+
(complex<T>(0,1)*((pow(SPA(5,2)*SPB(4,5)+
SPA(6,2)*SPB(4,6),2)*(SPA(5,3)*SPB(4,5)+
SPA(6,3)*SPB(4,6)))/(complex<T>(2,0)*SPA(3,2)*SPB(1,2)*SPB(4,5)*SPB(5,6)*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6)))+
(pow(SPA(2,1),2)*SPA(6,5)*pow(SPB(4,6),3)*(-(SPA(4,3)*SPB(4,5))+
SPA(6,3)*SPB(5,6)))/(complex<T>(3,0)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*pow(SPB(5,6),2)*pow(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6),2))))/SS(4,5,6)+
((complex<T>(0,-1)*SPB(1,4)*pow(SPA(5,2)*SPB(4,5)+
SPA(6,2)*SPB(4,6),2))/(complex<T>(2,0)*SPA(3,2)*SPB(1,2)*SPB(1,6)*SPB(4,5)*SPB(5,6))-
(complex<T>(0,1)*(SPA(5,3)*SPB(4,5)+
SPA(6,3)*SPB(4,6))*((pow(SPA(6,5),2)*SPB(3,4)*SPB(4,5)*((SPA(5,1)*SPB(1,6)+
SPA(5,2)*SPB(2,6))*SPB(4,5)+
(SPA(6,1)*SPB(1,6)+
SPA(6,2)*SPB(2,6))*SPB(4,6)))/SPB(5,6)+
(pow(SPA(6,5),2)*SPB(4,5)*pow((SPA(5,1)*SPB(1,6)+
SPA(5,2)*SPB(2,6))*SPB(4,5)+
(SPA(6,1)*SPB(1,6)+
SPA(6,2)*SPB(2,6))*SPB(4,6),2))/(complex<T>(3,0)*SPB(5,6)*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6)))))/(SPA(6,5)*pow(SPB(1,2),2)*SPB(4,5)*(SPA(5,1)*SPB(4,5)+
SPA(6,1)*SPB(4,6))*(-(SPA(4,3)*SPB(4,6))-
SPA(5,3)*SPB(5,6))))/SS(4,5,6)+
(complex<T>(0,1)*SPA(2,1)*SPB(1,3)*(-((pow(SPA(6,5),2)*SPB(3,4)*SPB(3,6)*SPB(4,5))/(SPA(2,1)*SPB(1,3)*SPB(5,6)*(-(SPA(6,4)*SPB(4,6))-
SPA(6,5)*SPB(5,6))))-
(pow(SPA(6,5),3)*pow(SPB(3,6),2)*SPB(4,5)*(-(complex<T>(2,0)*SPA(5,4)*SPB(4,5))-
SPA(6,4)*SPB(4,6)-
SPA(6,5)*SPB(5,6)))/(complex<T>(3,0)*SPA(2,1)*SPA(5,4)*SPB(1,3)*SPB(5,6)*pow(SPA(6,4)*SPB(4,6)+
SPA(6,5)*SPB(5,6),2))-
(pow(SPA(6,5),2)*SPB(3,6)*SPB(4,6)*SS(4,5,6))/(complex<T>(2,0)*SPA(2,1)*SPA(5,4)*SPB(1,6)*SPB(5,6)*(-(SPA(6,4)*SPB(4,6))-
SPA(6,5)*SPB(5,6)))))/(SPB(1,2)*SPB(2,3)*SS(4,5,6)));
}


template <class T> complex<T> Rqgqggg_mppppm_L(const eval_param<T>& ep,const mass_param_coll& mpc)
{
// added by Harald { qbm p qp p p m }_L
return
(complex<T>(0,-1)*SPB(3,2)*SPB(5,1)*(SPA(1,3)*SPB(5,1) +
SPA(3,6)*SPB(6,5))*((SPA(4,6)*SPB(5,4)*(SPA(1,4)*SPB(5,1) +
SPA(4,6)*SPB(6,5)))/(complex<T>(2,0)*SPA(3,4)*SPA(4,5)*SPB(3,2)*SPB(5,1)*(SPA
(4,5)*SPB(5,1) + SPA(4,6)*SPB(6,1))) +
(SPA(5,6)*SPB(5,4)*pow(SPA(1,4)*SPB(5,1) +
SPA(4,6)*SPB(6,5),2))/(complex<T>(3,0)*SPA(4,5)*SPB(3,2)*pow(SPA(4,5)*SPB(5,1
) + SPA(4,6)*SPB(6,1),2)*(SPA(1,3)*SPB(5,1) +
SPA(3,6)*SPB(6,5)))))/(SPA(2,3)*SPB(6,5)*(SPA(1,2)*SPB(5,1) +
SPA(2,6)*SPB(6,5))) +
((complex<T>(0,1)*pow(SPA(1,6),2)*SPA(3,6))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4
)*SPA(4,5)*SPA(5,6)) +
(complex<T>(0,1)*SPA(1,6)*SPB(3,2)*((SPA(1,6)*SPA(2,4)*SPA(2,6)*SPA(5,6)*pow(
SPB(5,4),2))/(pow(SPA(2,3),2)*SPA(4,5)*SPB(3,2)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))) -
(SPA(1,6)*SPA(2,4)*SPA(4,6)*pow(SPB(5,4),2))/(complex<T>(2,0)*SPA(3,4)*SPA(4,
5)*SPB(3,2)*(SPA(1,2)*SPB(5,1) + SPA(2,6)*SPB(6,5))) -
(SPA(1,6)*pow(SPA(2,4),2)*SPA(2,6)*SPA(5,6)*pow(SPB(5,4),3)*(-((SPA(1,2)*SPA(
4,6)*SPB(4,1))/SPA(2,6)) -
(complex<T>(2,0)*SPA(1,2)*SPA(5,6)*SPB(5,1))/SPA(2,6) - SPA(4,5)*SPB(5,4) -
SPA(4,6)*SPB(6,4)
-complex<T>(2,0)*SPA(5,6)*SPB(6,5)))/(complex<T>(3,0)*pow(SPA(2,3),2)*SPA(4,5)
*SPB(3,2)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))*((SPA(1,2)*SPA(4,6)*SPB(4,1))/SPA(2,6) + SPA(4,5)*SPB(5,4)
+ SPA(4,6)*SPB(6,4))*(SPA(1,2)*SPB(5,1) +
SPA(2,6)*SPB(6,5)))))/(SPA(2,6)*((SPA(1,2)*SPA(4,6)*SPB(4,1))/SPA(2,6) +
SPA(4,5)*SPB(5,4) + SPA(4,6)*SPB(6,4))))/SPA(1,2) +
((complex<T>(0,1)*(-((SPB(5,1)*SPB(5,4)*(SPA(1,3)*SPB(5,1) +
SPA(3,6)*SPB(6,5))*(SPA(1,4)*SPB(5,1) +
SPA(4,6)*SPB(6,5)))/(SPA(2,3)*SPA(3,4)*(SPA(4,5)*SPB(5,1) +
SPA(4,6)*SPB(6,1))*(SPA(1,2)*SPB(5,1) + SPA(2,6)*SPB(6,5)))) -
(SPB(3,2)*SPB(5,1)*pow(SPA(1,3)*SPB(5,1) +
SPA(3,6)*SPB(6,5),2))/(SPA(2,3)*SPA(3,4)*(SPA(4,5)*SPB(5,1) +
SPA(4,6)*SPB(6,1))*SS(1,5,6))))/complex<T>(2,0) -
(complex<T>(0,1)*pow(SPB(5,1),2)*SPB(5,4)*(SPA(1,4)*SPB(5,1) +
SPA(4,6)*SPB(6,5))*SS(1,5,6))/(complex<T>(3,0)*SPA(2,3)*pow(SPA(4,5)*SPB(5,1)
+ SPA(4,6)*SPB(6,1),2)*(SPA(1,2)*SPB(5,1) +
SPA(2,6)*SPB(6,5))))/(SPB(6,1)*SPB(6,5)) +
(pow(SPA(1,6),2)*((complex<T>(0,-1)*pow(SPA(2,6),2)*SPA(3,6)*SPB(2,1))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6)) +
complex<T>(0,1)*((-((SPA(2,4)*pow(SPA(2,6),2)*SPA(4,6)*SPB(2,1)*pow(SPB(5,4),2))/(SPA(3,4)*SPA(4,5)*SPA(5,6)*pow(-(SPA(2,3)*SPB(5,3))-SPA(2,4)*SPB(5,4),2)*(complex<T>(1,0) - (SPA(2,3)*(-(SPA(4,6)*SPB(4,3))-SPA(5,6)*SPB(5,3)))/(SPA(5,6)*(-(SPA(2,3)*SPB(5,3)) - SPA(2,4)*SPB(5,4))))))-(pow(SPA(2,6),2)*SPB(2,1)*SPB(5,3)*SPB(5,4))/(SPA(4,5)*(-(SPA(2,3)*SPB(5,3)) - SPA(2,4)*SPB(5,4))*SS(3,4,5)))/complex<T>(2,0) + ((-complex<T>(2,0)*pow(SPA(2,6),3)*SPB(2,1)*SPB(5,4))/(pow(SPA(2,3),2)*SPA(4,5)*(-(SPA(4,6)
*SPB(4,3)) - SPA(5,6)*SPB(5,3))) +
(complex<T>(3,0)*SPA(2,4)*pow(SPA(2,6),3)*SPA(5,6)*SPB(2,1)*pow(SPB(5,4),2))/
(pow(SPA(2,3),3)*SPA(4,5)*pow(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3),2)*(complex<T>(1,0) - (SPA(5,6)*(-(SPA(2,3)*SPB(5,3)) -
SPA(2,4)*SPB(5,4)))/(SPA(2,3)*(-(SPA(4,6)*SPB(4,3)) - SPA(5,6)*SPB(5,3))))) -
(pow(SPA(2,4),2)*pow(SPA(2,6),3)*pow(SPA(5,6),2)*SPB(2,1)*pow(SPB(5,4),3)*((SPA(2,3)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3)))/(SPA(5,6)*(-(SPA(2,3)*SPB(5,3)) - SPA(2,4)*SPB(5,4))) -
(SPA(5,6)*(-(SPA(2,3)*SPB(5,3))
-SPA(2,4)*SPB(5,4)))/(SPA(2,3)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3)))))/(pow(SPA(2,3),4)*SPA(4,5)*pow(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3),3)*pow(complex<T>(1,0) - (SPA(5,6)*(-(SPA(2,3)*SPB(5,3)) -
SPA(2,4)*SPB(5,4)))/(SPA(2,3)*(-(SPA(4,6)*SPB(4,3)) - SPA(5,6)*SPB(5,3))),3))
- (SPA(2,5)*pow(SPA(2,6),2)*SPB(2,1)*SPB(5,3)*SPB(5,4)*(SPA(3,6)*SPB(5,3) +
SPA(4,6)*SPB(5,4)))/(SPA(2,3)*SPA(4,5)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))*(-(SPA(2,3)*SPB(5,3)) - SPA(2,4)*SPB(5,4))*SS(3,4,5)) +
(pow(SPA(2,6),2)*SPB(2,1)*SPB(5,3)*(SPA(3,6)*SPB(4,3) -
SPA(5,6)*SPB(5,4)))/(SPA(2,3)*SPA(4,5)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))*SS(3,4,5)))/complex<T>(3,0))))/(SPA(1,2)*pow(SPA(2,6),2)*SPB(2,1)) - (complex<T>(0,1)*pow(SPA(1,3),2)*SPA(3,6)*pow(SPA(1,6)*SPB(3,1) +
SPA(2,6)*SPB(3,2),2))/(complex<T>(2,0)*SPA(1,2)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6)*(SPA(1,3)*SPB(3,1) + SPA(2,3)*SPB(3,2))*SS(4,5,6)) +
(complex<T>(0,1)*(((-(SPA(4,6)*SPB(4,1)) -
SPA(5,6)*SPB(5,1))*pow(-(SPA(4,6)*SPB(4,2)) -
SPA(5,6)*SPB(5,2),2))/(complex<T>(2,0)*SPA(2,3)*SPA(4,5)*SPA(5,6)*SPB(2,1)*(SPA(4,5)*SPB(5,1) + SPA(4,6)*SPB(6,1))) +
(pow(SPA(4,6),3)*pow(SPB(3,2),2)*SPB(5,4)*(-(SPA(4,5)*SPB(4,1)) +
SPA(5,6)*SPB(6,1)))/(complex<T>(3,0)*pow(SPA(4,5),2)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))*pow(SPA(4,5)*SPB(5,1) + SPA(4,6)*SPB(6,1),2))))/SS(4,5,6)
+ ((complex<T>(0,-1)*SPA(3,6)*pow(-(SPA(4,6)*SPB(4,2)) -
SPA(5,6)*SPB(5,2),2))/(complex<T>(2,0)*SPA(2,3)*SPA(3,4)*SPA(4,5)*SPA(5,6)*SPB(2,1)) - (complex<T>(0,1)*(-(SPA(4,6)*SPB(4,1)) -
SPA(5,6)*SPB(5,1))*((SPA(1,6)*SPA(5,6)*(-(SPA(4,6)*(SPA(2,4)*SPB(4,2) +
SPA(3,4)*SPB(4,3))) - SPA(5,6)*(SPA(2,4)*SPB(5,2) +
SPA(3,4)*SPB(5,3)))*pow(SPB(5,4),2))/SPA(4,5) +
(SPA(5,6)*pow(-(SPA(4,6)*(SPA(2,4)*SPB(4,2) + SPA(3,4)*SPB(4,3))) -
SPA(5,6)*(SPA(2,4)*SPB(5,2) +
SPA(3,4)*SPB(5,3)),2)*pow(SPB(5,4),2))/(complex<T>(3,0)*SPA(4,5)*(SPA(4,5)*SPB(5,1) +
SPA(4,6)*SPB(6,1)))))/(pow(SPA(2,3),2)*SPA(5,6)*(-(SPA(4,6)*SPB(4,3)) -
SPA(5,6)*SPB(5,3))*SPB(5,4)*(SPA(4,5)*SPB(5,1) +
SPA(4,6)*SPB(6,1))))/SS(4,5,6) +
(complex<T>(0,1)*SPA(1,3)*SPB(3,2)*(-((SPA(1,4)*SPA(1,6)*SPA(5,6)*pow(SPB(5,4
),2))/(SPA(1,3)*SPA(4,5)*SPB(3,2)*(-(SPA(4,5)*SPB(5,4)) -
SPA(4,6)*SPB(6,4)))) -
(pow(SPA(1,4),2)*SPA(5,6)*pow(SPB(5,4),3)*(-(SPA(4,5)*SPB(5,4)) -
SPA(4,6)*SPB(6,4) -
complex<T>(2,0)*SPA(5,6)*SPB(6,5)))/(complex<T>(3,0)*SPA(1,3)*SPA(4,5)*SPB(3,
2)*pow(SPA(4,5)*SPB(5,4) + SPA(4,6)*SPB(6,4),2)*SPB(6,5)) -
(SPA(1,4)*SPA(4,6)*pow(SPB(5,4),2)*SS(4,5,6))/(complex<T>(2,0)*SPA(3,4)*SPA(4
,5)*SPB(3,2)*(-(SPA(4,5)*SPB(5,4)) -
SPA(4,6)*SPB(6,4))*SPB(6,5))))/(SPA(1,2)*SPA(2,3)*SS(4,5,6))
;}


/* *************** table of switch values ************* */

/*
2304: m m m m qm qp
1536: m m m m qp qm
2112: m m m qm m qp
576: m m m qm qp m
3648: m m m qm qp p
2880: m m m qm p qp
1152: m m m qp m qm
384: m m m qp qm m
3456: m m m qp qm p
1920: m m m qp p qm
2496: m m m p qm qp
1728: m m m p qp qm
2064: m m qm m m qp
528: m m qm m qp m
3600: m m qm m qp p
2832: m m qm m p qp
144: m m qm qp m m
3216: m m qm qp m p
912: m m qm qp p m
3984: m m qm qp p p
2256: m m qm p m qp
720: m m qm p qp m
3792: m m qm p qp p
3024: m m qm p p qp
1056: m m qp m m qm
288: m m qp m qm m
3360: m m qp m qm p
1824: m m qp m p qm
96: m m qp qm m m
3168: m m qp qm m p
864: m m qp qm p m
3936: m m qp qm p p
1248: m m qp p m qm
480: m m qp p qm m
3552: m m qp p qm p
2016: m m qp p p qm
2352: m m p m qm qp
1584: m m p m qp qm
2160: m m p qm m qp
624: m m p qm qp m
3696: m m p qm qp p
2928: m m p qm p qp
1200: m m p qp m qm
432: m m p qp qm m
3504: m m p qp qm p
1968: m m p qp p qm
2544: m m p p qm qp
1776: m m p p qp qm
2052: m qm m m m qp
516: m qm m m qp m
3588: m qm m m qp p
2820: m qm m m p qp
132: m qm m qp m m
3204: m qm m qp m p
900: m qm m qp p m
3972: m qm m qp p p
2244: m qm m p m qp
708: m qm m p qp m
3780: m qm m p qp p
3012: m qm m p p qp
36: m qm qp m m m
3108: m qm qp m m p
804: m qm qp m p m
3876: m qm qp m p p
228: m qm qp p m m
3300: m qm qp p m p
996: m qm qp p p m
4068: m qm qp p p p
2100: m qm p m m qp
564: m qm p m qp m
3636: m qm p m qp p
2868: m qm p m p qp
180: m qm p qp m m
3252: m qm p qp m p
948: m qm p qp p m
4020: m qm p qp p p
2292: m qm p p m qp
756: m qm p p qp m
3828: m qm p p qp p
3060: m qm p p p qp
1032: m qp m m m qm
264: m qp m m qm m
3336: m qp m m qm p
1800: m qp m m p qm
72: m qp m qm m m
3144: m qp m qm m p
840: m qp m qm p m
3912: m qp m qm p p
1224: m qp m p m qm
456: m qp m p qm m
3528: m qp m p qm p
1992: m qp m p p qm
24: m qp qm m m m
3096: m qp qm m m p
792: m qp qm m p m
3864: m qp qm m p p
216: m qp qm p m m
3288: m qp qm p m p
984: m qp qm p p m
4056: m qp qm p p p
1080: m qp p m m qm
312: m qp p m qm m
3384: m qp p m qm p
1848: m qp p m p qm
120: m qp p qm m m
3192: m qp p qm m p
888: m qp p qm p m
3960: m qp p qm p p
1272: m qp p p m qm
504: m qp p p qm m
3576: m qp p p qm p
2040: m qp p p p qm
2316: m p m m qm qp
1548: m p m m qp qm
2124: m p m qm m qp
588: m p m qm qp m
3660: m p m qm qp p
2892: m p m qm p qp
1164: m p m qp m qm
396: m p m qp qm m
3468: m p m qp qm p
1932: m p m qp p qm
2508: m p m p qm qp
1740: m p m p qp qm
2076: m p qm m m qp
540: m p qm m qp m
3612: m p qm m qp p
2844: m p qm m p qp
156: m p qm qp m m
3228: m p qm qp m p
924: m p qm qp p m
3996: m p qm qp p p
2268: m p qm p m qp
732: m p qm p qp m
3804: m p qm p qp p
3036: m p qm p p qp
1068: m p qp m m qm
300: m p qp m qm m
3372: m p qp m qm p
1836: m p qp m p qm
108: m p qp qm m m
3180: m p qp qm m p
876: m p qp qm p m
3948: m p qp qm p p
1260: m p qp p m qm
492: m p qp p qm m
3564: m p qp p qm p
2028: m p qp p p qm
2364: m p p m qm qp
1596: m p p m qp qm
2172: m p p qm m qp
636: m p p qm qp m
3708: m p p qm qp p
2940: m p p qm p qp
1212: m p p qp m qm
444: m p p qp qm m
3516: m p p qp qm p
1980: m p p qp p qm
2556: m p p p qm qp
1788: m p p p qp qm
2049: qm m m m m qp
513: qm m m m qp m
3585: qm m m m qp p
2817: qm m m m p qp
129: qm m m qp m m
3201: qm m m qp m p
897: qm m m qp p m
3969: qm m m qp p p
2241: qm m m p m qp
705: qm m m p qp m
3777: qm m m p qp p
3009: qm m m p p qp
33: qm m qp m m m
3105: qm m qp m m p
801: qm m qp m p m
3873: qm m qp m p p
225: qm m qp p m m
3297: qm m qp p m p
993: qm m qp p p m
4065: qm m qp p p p
2097: qm m p m m qp
561: qm m p m qp m
3633: qm m p m qp p
2865: qm m p m p qp
177: qm m p qp m m
3249: qm m p qp m p
945: qm m p qp p m
4017: qm m p qp p p
2289: qm m p p m qp
753: qm m p p qp m
3825: qm m p p qp p
3057: qm m p p p qp
9: qm qp m m m m
3081: qm qp m m m p
777: qm qp m m p m
3849: qm qp m m p p
201: qm qp m p m m
3273: qm qp m p m p
969: qm qp m p p m
4041: qm qp m p p p
57: qm qp p m m m
3129: qm qp p m m p
825: qm qp p m p m
3897: qm qp p m p p
249: qm qp p p m m
3321: qm qp p p m p
1017: qm qp p p p m
4089: qm qp p p p p
2061: qm p m m m qp
525: qm p m m qp m
3597: qm p m m qp p
2829: qm p m m p qp
141: qm p m qp m m
3213: qm p m qp m p
909: qm p m qp p m
3981: qm p m qp p p
2253: qm p m p m qp
717: qm p m p qp m
3789: qm p m p qp p
3021: qm p m p p qp
45: qm p qp m m m
3117: qm p qp m m p
813: qm p qp m p m
3885: qm p qp m p p
237: qm p qp p m m
3309: qm p qp p m p
1005: qm p qp p p m
4077: qm p qp p p p
2109: qm p p m m qp
573: qm p p m qp m
3645: qm p p m qp p
2877: qm p p m p qp
189: qm p p qp m m
3261: qm p p qp m p
957: qm p p qp p m
4029: qm p p qp p p
2301: qm p p p m qp
765: qm p p p qp m
3837: qm p p p qp p
3069: qm p p p p qp
1026: qp m m m m qm
258: qp m m m qm m
3330: qp m m m qm p
1794: qp m m m p qm
66: qp m m qm m m
3138: qp m m qm m p
834: qp m m qm p m
3906: qp m m qm p p
1218: qp m m p m qm
450: qp m m p qm m
3522: qp m m p qm p
1986: qp m m p p qm
18: qp m qm m m m
3090: qp m qm m m p
786: qp m qm m p m
3858: qp m qm m p p
210: qp m qm p m m
3282: qp m qm p m p
978: qp m qm p p m
4050: qp m qm p p p
1074: qp m p m m qm
306: qp m p m qm m
3378: qp m p m qm p
1842: qp m p m p qm
114: qp m p qm m m
3186: qp m p qm m p
882: qp m p qm p m
3954: qp m p qm p p
1266: qp m p p m qm
498: qp m p p qm m
3570: qp m p p qm p
2034: qp m p p p qm
6: qp qm m m m m
3078: qp qm m m m p
774: qp qm m m p m
3846: qp qm m m p p
198: qp qm m p m m
3270: qp qm m p m p
966: qp qm m p p m
4038: qp qm m p p p
54: qp qm p m m m
3126: qp qm p m m p
822: qp qm p m p m
3894: qp qm p m p p
246: qp qm p p m m
3318: qp qm p p m p
1014: qp qm p p p m
4086: qp qm p p p p
1038: qp p m m m qm
270: qp p m m qm m
3342: qp p m m qm p
1806: qp p m m p qm
78: qp p m qm m m
3150: qp p m qm m p
846: qp p m qm p m
3918: qp p m qm p p
1230: qp p m p m qm
462: qp p m p qm m
3534: qp p m p qm p
1998: qp p m p p qm
30: qp p qm m m m
3102: qp p qm m m p
798: qp p qm m p m
3870: qp p qm m p p
222: qp p qm p m m
3294: qp p qm p m p
990: qp p qm p p m
4062: qp p qm p p p
1086: qp p p m m qm
318: qp p p m qm m
3390: qp p p m qm p
1854: qp p p m p qm
126: qp p p qm m m
3198: qp p p qm m p
894: qp p p qm p m
3966: qp p p qm p p
1278: qp p p p m qm
510: qp p p p qm m
3582: qp p p p qm p
2046: qp p p p p qm
2307: p m m m qm qp
1539: p m m m qp qm
2115: p m m qm m qp
579: p m m qm qp m
3651: p m m qm qp p
2883: p m m qm p qp
1155: p m m qp m qm
387: p m m qp qm m
3459: p m m qp qm p
1923: p m m qp p qm
2499: p m m p qm qp
1731: p m m p qp qm
2067: p m qm m m qp
531: p m qm m qp m
3603: p m qm m qp p
2835: p m qm m p qp
147: p m qm qp m m
3219: p m qm qp m p
915: p m qm qp p m
3987: p m qm qp p p
2259: p m qm p m qp
723: p m qm p qp m
3795: p m qm p qp p
3027: p m qm p p qp
1059: p m qp m m qm
291: p m qp m qm m
3363: p m qp m qm p
1827: p m qp m p qm
99: p m qp qm m m
3171: p m qp qm m p
867: p m qp qm p m
3939: p m qp qm p p
1251: p m qp p m qm
483: p m qp p qm m
3555: p m qp p qm p
2019: p m qp p p qm
2355: p m p m qm qp
1587: p m p m qp qm
2163: p m p qm m qp
627: p m p qm qp m
3699: p m p qm qp p
2931: p m p qm p qp
1203: p m p qp m qm
435: p m p qp qm m
3507: p m p qp qm p
1971: p m p qp p qm
2547: p m p p qm qp
1779: p m p p qp qm
2055: p qm m m m qp
519: p qm m m qp m
3591: p qm m m qp p
2823: p qm m m p qp
135: p qm m qp m m
3207: p qm m qp m p
903: p qm m qp p m
3975: p qm m qp p p
2247: p qm m p m qp
711: p qm m p qp m
3783: p qm m p qp p
3015: p qm m p p qp
39: p qm qp m m m
3111: p qm qp m m p
807: p qm qp m p m
3879: p qm qp m p p
231: p qm qp p m m
3303: p qm qp p m p
999: p qm qp p p m
4071: p qm qp p p p
2103: p qm p m m qp
567: p qm p m qp m
3639: p qm p m qp p
2871: p qm p m p qp
183: p qm p qp m m
3255: p qm p qp m p
951: p qm p qp p m
4023: p qm p qp p p
2295: p qm p p m qp
759: p qm p p qp m
3831: p qm p p qp p
3063: p qm p p p qp
1035: p qp m m m qm
267: p qp m m qm m
3339: p qp m m qm p
1803: p qp m m p qm
75: p qp m qm m m
3147: p qp m qm m p
843: p qp m qm p m
3915: p qp m qm p p
1227: p qp m p m qm
459: p qp m p qm m
3531: p qp m p qm p
1995: p qp m p p qm
27: p qp qm m m m
3099: p qp qm m m p
795: p qp qm m p m
3867: p qp qm m p p
219: p qp qm p m m
3291: p qp qm p m p
987: p qp qm p p m
4059: p qp qm p p p
1083: p qp p m m qm
315: p qp p m qm m
3387: p qp p m qm p
1851: p qp p m p qm
123: p qp p qm m m
3195: p qp p qm m p
891: p qp p qm p m
3963: p qp p qm p p
1275: p qp p p m qm
507: p qp p p qm m
3579: p qp p p qm p
2043: p qp p p p qm
2319: p p m m qm qp
1551: p p m m qp qm
2127: p p m qm m qp
591: p p m qm qp m
3663: p p m qm qp p
2895: p p m qm p qp
1167: p p m qp m qm
399: p p m qp qm m
3471: p p m qp qm p
1935: p p m qp p qm
2511: p p m p qm qp
1743: p p m p qp qm
2079: p p qm m m qp
543: p p qm m qp m
3615: p p qm m qp p
2847: p p qm m p qp
159: p p qm qp m m
3231: p p qm qp m p
927: p p qm qp p m
3999: p p qm qp p p
2271: p p qm p m qp
735: p p qm p qp m
3807: p p qm p qp p
3039: p p qm p p qp
1071: p p qp m m qm
303: p p qp m qm m
3375: p p qp m qm p
1839: p p qp m p qm
111: p p qp qm m m
3183: p p qp qm m p
879: p p qp qm p m
3951: p p qp qm p p
1263: p p qp p m qm
495: p p qp p qm m
3567: p p qp p qm p
2031: p p qp p p qm
2367: p p p m qm qp
1599: p p p m qp qm
2175: p p p qm m qp
639: p p p qm qp m
3711: p p p qm qp p
2943: p p p qm p qp
1215: p p p qp m qm
447: p p p qp qm m
3519: p p p qp qm p
1983: p p p qp p qm
2559: p p p p qm qp
1791: p p p p qp qm
*/


template <class T> complex<T> R2q4g_L_4089(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mppppp_L(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_R_4089(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mppppp_R(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_nf_4089(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mppppp_nf(ep,mpc);}   // qbm qp p p p p

template <class T> complex<T> R2q4g_L_4086(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmpppp_L(ep,mpc);}   // qbp qm p p p p
template <class T> complex<T> R2q4g_R_4086(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmpppp_R(ep,mpc);}   // qbp qm p p p p
template <class T> complex<T> R2q4g_nf_4086(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmpppp_nf(ep,mpc);}   // qbp qm p p p p


template <class T> complex<T> R2q4g_L_6(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmmmmm_L(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_R_6(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmmmmm_R(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_nf_6(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_pmmmmm_nf(ep,mpc);}   // qbm qp p p p p

template <class T> complex<T> R2q4g_L_9(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mpmmmm_L(ep,mpc);}   // qbp qm p p p p
template <class T> complex<T> R2q4g_R_9(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mpmmmm_R(ep,mpc);}   // qbp qm p p p p
template <class T> complex<T> R2q4g_nf_9(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q4g_mpmmmm_nf(ep,mpc);}   // qbp qm p p p p

template <class T> complex<T> R2q4g_L_4065(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpppp_L(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_R_4065(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpppp_R(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_nf_4065(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpppp_nf(ep,mpc);}   // qbm qp p p p p

template <class T> complex<T> R2q4g_L_30(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmmmm_L(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_R_30(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmmmm_R(ep,mpc);}   // qbm qp p p p p
template <class T> complex<T> R2q4g_nf_30(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmmmm_nf(ep,mpc);}   // qbm qp p p p p

template <class T> complex<T> R2q4g_L_4077(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mppppp_sl(ep,mpc);}   // qbm p qp p p p
template <class T> complex<T> R2q4g_R_4077(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mppppp_R(ep,mpc);}   // qbm p qp p p p
template <class T> complex<T> R2q4g_nf_4077(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mppppp_nf(ep,mpc);}   // qbm p qp p p p

template <class T> complex<T> R2q4g_L_18(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_pmmmmm_L(ep,mpc);}   // qbp m qm m m m
template <class T> complex<T> R2q4g_R_18(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_pmmmmm_R(ep,mpc);}   // qbp m qm m m m
template <class T> complex<T> R2q4g_nf_18(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_pmmmmm_nf(ep,mpc);}   // qbp m qm m m m

template <class T> complex<T> R2q4g_L_33(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpmmm_L(ep,mpc);}   // qbm m qp m m m
template <class T> complex<T> R2q4g_R_33(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpmmm_R(ep,mpc);}   // qbm m qp m m m
template <class T> complex<T> R2q4g_nf_33(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmpmmm_nf(ep,mpc);}   // qbm m qp m m m

template <class T> complex<T> R2q4g_L_4062(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmppp_L(ep,mpc);}   // qbm p qp p p p
template <class T> complex<T> R2q4g_R_4062(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmppp_R(ep,mpc);}   // qbm p qp p p p
template <class T> complex<T> R2q4g_nf_4062(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_ppmppp_nf(ep,mpc);}   // qbm p qp p p p

// added by Harald
template <class T> complex<T> R2q4g_L_1005(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mppppm_L(ep,mpc);}   // qbm p qp p p m
template <class T> complex<T> R2q4g_L_225(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mmppmm_L(ep,mpc);}   // qbm m qp p m m

//template <class T> complex<T> R2q4g_L_4077(const eval_param<T>& ep,const mass_param_coll& mpc){return Rqgqggg_mppppp_sl(ep,mpc);}   // qbm p qp p p p
//template <class T> complex<T> R2q4g_L_4023(const eval_param<T>& ep,const mass_param_coll& mpc){vector<int> new_mpc; rotate_copy(mpc.begin(),mpc.begin()+1,mpc.end(),back_inserter(new_mpc)); return Rqgqggg_mppppp_sl(ep,new_mpc);}
//template <class T> complex<T> R2q4g_L_3807(const eval_param<T>& ep,const mass_param_coll& mpc){vector<int> new_mpc; rotate_copy(mpc.begin(),mpc.begin()+2,mpc.end(),back_inserter(new_mpc)); return Rqgqggg_mppppp_sl(ep,new_mpc);}
//template <class T> complex<T> R2q4g_L_2943(const eval_param<T>& ep,const mass_param_coll& mpc){vector<int> new_mpc; rotate_copy(mpc.begin(),mpc.begin()+3,mpc.end(),back_inserter(new_mpc)); return Rqgqggg_mppppp_sl(ep,new_mpc);}
//template <class T> complex<T> R2q4g_L_3582(const eval_param<T>& ep,const mass_param_coll& mpc){vector<int> new_mpc; rotate_copy(mpc.begin(),mpc.begin()+4,mpc.end(),back_inserter(new_mpc)); return Rqgqggg_mppppp_sl(ep,new_mpc);}
//template <class T> complex<T> R2q4g_L_2043(const eval_param<T>& ep,const mass_param_coll& mpc){vector<int> new_mpc; rotate_copy(mpc.begin(),mpc.begin()+4,mpc.end(),back_inserter(new_mpc)); return Rqgqggg_mppppp_sl(ep,new_mpc);}

#define _CASE_R2q4g_L_Ptr_eval(K) case K : return &R2q4g_L_ ## K

template <class T> complex<T> (*R2q4g_L_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {


	switch (hc) {

	_CASE_R2q4g_L_Ptr_eval(4077);
//	_CASE_R2q4g_L_Ptr_eval(4023);
//	_CASE_R2q4g_L_Ptr_eval(3807);
//	_CASE_R2q4g_L_Ptr_eval(2943);
//	_CASE_R2q4g_L_Ptr_eval(3582);
//	_CASE_R2q4g_L_Ptr_eval(2043);
	_CASE_R2q4g_L_Ptr_eval(18);
	_CASE_R2q4g_L_Ptr_eval(33);
	_CASE_R2q4g_L_Ptr_eval(4062);

	_CASE_R2q4g_L_Ptr_eval(4089);
	_CASE_R2q4g_L_Ptr_eval(6);
	_CASE_R2q4g_L_Ptr_eval(4086);
	_CASE_R2q4g_L_Ptr_eval(9);
	_CASE_R2q4g_L_Ptr_eval(4065);
	_CASE_R2q4g_L_Ptr_eval(30);

	_CASE_R2q4g_L_Ptr_eval(1005);
	_CASE_R2q4g_L_Ptr_eval(225);


	default: return 0;
	}
}

#define _CASE_R2q4g_R_Ptr_eval(K) case K : return &R2q4g_R_ ## K

template <class T> complex<T> (*R2q4g_SLC_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {


	switch (hc) {

	_CASE_R2q4g_R_Ptr_eval(4089);
	_CASE_R2q4g_R_Ptr_eval(6);
	_CASE_R2q4g_R_Ptr_eval(4086);
	_CASE_R2q4g_R_Ptr_eval(9);
	_CASE_R2q4g_R_Ptr_eval(4065);
	_CASE_R2q4g_R_Ptr_eval(30);

	_CASE_R2q4g_R_Ptr_eval(4077);
	_CASE_R2q4g_R_Ptr_eval(18);
	_CASE_R2q4g_R_Ptr_eval(33);
	_CASE_R2q4g_R_Ptr_eval(4062);


	default: return 0;
	}
}
#define _CASE_R2q4g_nf_Ptr_eval(K) case K : return &R2q4g_nf_ ## K

template <class T> complex<T> (*R2q4g_nf_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {


	switch (hc) {

	_CASE_R2q4g_nf_Ptr_eval(4089);
	_CASE_R2q4g_nf_Ptr_eval(6);
	_CASE_R2q4g_nf_Ptr_eval(4086);
	_CASE_R2q4g_nf_Ptr_eval(9);
	_CASE_R2q4g_nf_Ptr_eval(4065);
	_CASE_R2q4g_nf_Ptr_eval(30);

	_CASE_R2q4g_nf_Ptr_eval(4077);
	_CASE_R2q4g_nf_Ptr_eval(18);
	_CASE_R2q4g_nf_Ptr_eval(33);
	_CASE_R2q4g_nf_Ptr_eval(4062);


	default: return 0;
	}
}


template complex<R> (*R2q4g_L_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R2q4g_L_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R2q4g_L_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R2q4g_L_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

 template complex<R> (*R2q4g_SLC_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R2q4g_SLC_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R2q4g_SLC_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R2q4g_SLC_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

  template complex<R> (*R2q4g_nf_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
   template complex<RHP> (*R2q4g_nf_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
   template complex<RVHP> (*R2q4g_nf_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

   template complex<RGMP> (*R2q4g_nf_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif
}

