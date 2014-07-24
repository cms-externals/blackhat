/*
 * A0_2q2l2Q_eval.cpp
 *
 *  Created on: Sep 5, 2008
 *      Author_eval: daniel
 */


#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include <iostream>

using namespace std;

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1,k-1)

namespace BH {

template <class T> complex<T> A2q2l2G_qbllbqGBG_mmppmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return

	(complex<T>(0,-1)*(SPA(1,5)*SPB(4,3)*
	        (SPA(2,3)*SPB(6,3) + SPA(2,4)*SPB(6,4))*SS(1,2,3) -
	       SPA(1,2)*(SPA(1,5)*SPB(3,1) + SPA(2,5)*SPB(3,2))*SPB(6,4)*
	        SS(2,3,4)))/(S(2,3)*S(5,6)*SS(1,2,3)*SS(2,3,4));


}

template <class T> complex<T> A2q2l2G_qbllbqGBG_pmpmmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return
-	(complex<T>(0,1)*(SPA(2,4)*(SPA(2,5)*SPB(3,2) - SPA(4,5)*SPB(4,3))*SPB(6,1)*
	        SS(1,2,3) + SPA(4,5)*SPB(3,1)*
	        (SPA(1,2)*SPB(6,1) - SPA(2,3)*SPB(6,3))*SS(2,3,4)))/
	   (S(2,3)*S(5,6)*SS(1,2,3)*SS(2,3,4));
}

template <class T> complex<T> A2q2l2G_qbllbqGBG_mmpppm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return
	- (complex<T>(0,-1)*(SPA(1,6)*SPB(4,3)*
	        (SPA(2,3)*SPB(5,3) + SPA(2,4)*SPB(5,4))*SS(1,2,3) -
	       SPA(1,2)*(SPA(1,6)*SPB(3,1) + SPA(2,6)*SPB(3,2))*SPB(5,4)*
	        SS(2,3,4)))/(S(2,3)*S(5,6)*SS(1,2,3)*SS(2,3,4));
}

template <class T> complex<T> A2q2l2G_qbllbqGBG_pmpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return
	(complex<T>(0,1)*(SPA(2,4)*(SPA(2,6)*SPB(3,2) - SPA(4,6)*SPB(4,3))*SPB(5,1)*
	        SS(1,2,3) + SPA(4,6)*SPB(3,1)*
	        (SPA(1,2)*SPB(5,1) - SPA(2,3)*SPB(5,3))*SS(2,3,4)))/
	   (S(2,3)*S(5,6)*SS(1,2,3)*SS(2,3,4));

}

template <class T> complex<T> A2q2l2G_qbGBGllbq_mpmmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return
	(complex<T>(0,1)*(SS(2,3,6)*SPA(1,3)*
	        (-(SPA(1,4)*SPB(2,1)) + SPA(3,4)*SPB(3,2))*SPB(6,5) +
	       SS(1,2,3)*SPA(1,4)*SPB(6,2)*(SPA(2,3)*SPB(5,2) + SPA(3,6)*SPB(6,5))
	       ))/(S(2,3)*S(4,5)*SS(1,2,3)*SS(2,3,6));
}

template <class T> complex<T> A2q2l2G_qbGBGllbq_mmpmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return(complex<T>(0,1)*(SPA(1,4)*SPB(6,3)*(SPA(2,3)*SPB(5,3) - SPA(2,6)*SPB(6,5))*
        SS(1,2,3) + SPA(1,2)*(SPA(1,4)*SPB(3,1) + SPA(2,4)*SPB(3,2))*
        SPB(6,5)*SS(2,3,6)))/(S(2,3)*S(4,5)*SS(1,2,3)*SS(2,3,6));
}

template <class T> complex<T> A2q2l2G_qbGBGllbq_pmpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return(complex<T>(0,1)*(SPA(2,6)*SPB(5,1)*
	        (SPA(2,4)*SPB(3,2) + SPA(4,6)*SPB(6,3))*SS(1,2,3) +
	       SPA(4,6)*SPB(3,1)*(-(SPA(1,2)*SPB(5,1)) + SPA(2,3)*SPB(5,3))*
	        SS(2,3,6)))/(S(2,3)*S(4,5)*SS(1,2,3)*SS(2,3,6));
}
template <class T> complex<T> A2q2l2G_qbGBGllbq_ppmmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return (complex<T>(0,1)*(SPA(3,6)*SPB(5,1)*
	        (SPA(3,4)*SPB(3,2) - SPA(4,6)*SPB(6,2))*SS(1,2,3) +
	       SPA(4,6)*SPB(2,1)*(SPA(1,3)*SPB(5,1) + SPA(2,3)*SPB(5,2))*
	        SS(2,3,6)))/(S(2,3)*S(4,5)*SS(1,2,3)*SS(2,3,6))
	;
}


template <class T> complex<T> A2q2l2G_GqqbllbGB_mpmmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return
	(complex<T>(0,1)*(SS(6,1,2)*SPA(3,1)*
	        (-(SPA(3,4)*SPB(6,3)) + SPA(1,4)*SPB(1,6))*SPB(2,5) +
	       SS(3,6,1)*SPA(3,4)*SPB(2,6)*(SPA(6,1)*SPB(5,6) + SPA(1,2)*SPB(2,5))
	       ))/(S(6,1)*S(4,5)*SS(3,6,1)*SS(6,1,2));
}

template <class T> complex<T> A2q2l2G_GqqbllbGB_ppmmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return(complex<T>(0,1)*(SPA(3,4)*SPB(2,1)*(SPA(6,1)*SPB(5,1) - SPA(6,2)*SPB(2,5))*
        SS(3,6,1) + SPA(3,6)*(SPA(3,4)*SPB(1,3) + SPA(6,4)*SPB(1,6))*
        SPB(2,5)*SS(6,1,2)))/(S(6,1)*S(4,5)*SS(3,6,1)*SS(6,1,2));
}

template <class T> complex<T> A2q2l2G_GqqbllbGB_pmpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return(complex<T>(0,1)*(SPA(6,2)*SPB(5,3)*
	        (SPA(6,4)*SPB(1,6) + SPA(4,2)*SPB(2,1))*SS(3,6,1) +
	       SPA(4,2)*SPB(1,3)*(-(SPA(3,6)*SPB(5,3)) + SPA(6,1)*SPB(5,1))*
	        SS(6,1,2)))/(S(6,1)*S(4,5)*SS(3,6,1)*SS(6,1,2));
}
template <class T> complex<T> A2q2l2G_GqqbllbGB_mmpmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return (complex<T>(0,1)*(SPA(1,2)*SPB(5,3)*
	        (SPA(1,4)*SPB(1,6) - SPA(4,2)*SPB(2,6))*SS(3,6,1) +
	       SPA(4,2)*SPB(6,3)*(SPA(3,1)*SPB(5,3) + SPA(6,1)*SPB(5,6))*
	        SS(6,1,2)))/(S(6,1)*S(4,5)*SS(3,6,1)*SS(6,1,2))
	;
}

template <class T> complex<T>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {

   case 107504: return &A2q2l2G_qbllbqGBG_mmppmp_eval;
   case 78832: return &A2q2l2G_qbllbqGBG_mmpppm_eval;
   case 106993: return &A2q2l2G_qbllbqGBG_pmpmmp_eval;
   case  78321: return &A2q2l2G_qbllbqGBG_pmpmpm_eval;

   case  64664: return &A2q2l2G_qbGBGllbq_mpmmpp_eval;
   case  64720: return &A2q2l2G_qbGBGllbq_mmpmpp_eval;
   case  31897: return &A2q2l2G_qbGBGllbq_ppmmpm_eval;
   case  31953: return &A2q2l2G_qbGBGllbq_pmpmpm_eval;

   case  130058: return &A2q2l2G_GqqbllbGB_mpmmpp_eval;
   case  130114: return &A2q2l2G_GqqbllbGB_mmpmpp_eval;
   case  97291: return &A2q2l2G_GqqbllbGB_ppmmpm_eval;
   case  97347: return &A2q2l2G_GqqbllbGB_pmpmpm_eval;

default: return 0;

}
}


template complex<R>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP
template complex<RGMP>  (*A2q2l2G_Tree_Ptr_eval(long hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif

}


/* *************** table of switch values ************* */

/*
1237: qmqpemep
637: qmemepqp
1232: qpqmemep
422: qpemepqm
502: emepqmqp
322: emepqpqm
947: epqmqpem
917: epqpqmem
*/
