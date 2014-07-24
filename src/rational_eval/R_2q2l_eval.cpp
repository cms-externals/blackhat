

#include "R_2q2l_eval.h"
#include "eval_param.h"

#include "BH_error.h"

using namespace std;
namespace BH  {


#define _VERBOSE 0


#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)

template <int i1, int i2,int i3,int i4,class T> complex<T> R2q2l_pmmp(const eval_param<T>& ep,const mass_param_coll& mpc){

#if _VERBOSE
	_MESSAGE("R2q2l : +--+");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*pow(SPA(i2,i3),2))/(SPA(i1,i2)*
			SPA(i3,i4));
}

template <int i1, int i2,int i3,int i4,class T> complex<T> R2q2l_mpmp(const eval_param<T>& ep,const mass_param_coll& mpc){

#if _VERBOSE
	_MESSAGE("R2q2l : -+-+");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*pow(SPA(i1,i3),2))/(SPA(i1,i2)*
			SPA(i3,i4));
}

// obtained by changing <ab> --> [ba]

template <class T> complex<T> R2q2l_mppm(const eval_param<T>& ep,const mass_param_coll& mpc){

#if _VERBOSE
	_MESSAGE("R2q2l : +--+");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*pow(SPB(3,2),2))/(SPB(2,1)*	SPB(4,3));
}

template <class T> complex<T> R2q2l_pmpm(const eval_param<T>& ep,const mass_param_coll& mpc){

#if _VERBOSE
	_MESSAGE("R2q2l : -+-+");
#endif
	return (complex<T>(0,1)/complex<T>(2,0)*pow(SPB(3,1),2))/(SPB(2,1)*	SPB(4,3));
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

template <class T> complex<T> R2q2l_L_1232(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_pmmp<1,2,3,4,T>(ep,mpc);}   // qp qm em ep
template <class T> complex<T> R2q2l_L_917 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_pmmp<2,3,4,1,T>(ep,mpc);}   // ep qp qm em
template <class T> complex<T> R2q2l_L_322 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_pmmp<3,4,1,2,T>(ep,mpc);}   // em ep qp qm
template <class T> complex<T> R2q2l_L_637 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_pmmp<4,1,2,3,T>(ep,mpc);}   // qm em ep qp

template <class T> complex<T> R2q2l_L_1237(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_mpmp<1,2,3,4,T>(ep,mpc);}   // qm qp em ep
template <class T> complex<T> R2q2l_L_947 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_mpmp<2,3,4,1,T>(ep,mpc);}   // ep qm qp em
template <class T> complex<T> R2q2l_L_502 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_mpmp<3,4,1,2,T>(ep,mpc);}   // em ep qm qp
template <class T> complex<T> R2q2l_L_422 (const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_mpmp<4,1,2,3,T>(ep,mpc);}   // qp em ep qm

template <class T> complex<T> R2q2l_L_1057(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_mppm<T>(ep,mpc);}   // qm qp ep em
template <class T> complex<T> R2q2l_L_1052(const eval_param<T>& ep,const mass_param_coll& mpc){return R2q2l_pmpm<T>(ep,mpc);}   // qp qm ep em


#define _CASE_R2q2l_L_Ptr_eval(K) case K : return &R2q2l_L_ ## K

template <class T> complex<T> (*R2q2l_L_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
	_CASE_R2q2l_L_Ptr_eval(1232);
	_CASE_R2q2l_L_Ptr_eval(917);
	_CASE_R2q2l_L_Ptr_eval(322);
	_CASE_R2q2l_L_Ptr_eval(637);
	_CASE_R2q2l_L_Ptr_eval(1237);
	_CASE_R2q2l_L_Ptr_eval(947);
	_CASE_R2q2l_L_Ptr_eval(502);
	_CASE_R2q2l_L_Ptr_eval(422);
	_CASE_R2q2l_L_Ptr_eval(1057);
	_CASE_R2q2l_L_Ptr_eval(1052);


	default: return 0;
	}
}









template complex<R> (*R2q2l_L_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R2q2l_L_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R2q2l_L_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R2q2l_L_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif


 template <class T> complex<T> R2q2l_nf(const eval_param<T>& ep,const mass_param_coll& mpc){

 #if _VERBOSE
 	_MESSAGE("R2q2l : -+-+");
 #endif
 	return complex<T>(0,0);
 }

 template <class T> complex<T> (*R2q2l_nf_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	 return &R2q2l_nf;

 }



 template complex<R> (*R2q2l_nf_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R2q2l_nf_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R2q2l_nf_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R2q2l_nf_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif





}

