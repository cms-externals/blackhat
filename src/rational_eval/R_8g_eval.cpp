

#include "R_8g_eval.h"
#include "eval_param.h"
#include "R_alln_eval.h"
 using namespace std;
namespace BH   {

int helcode_g(const process& pro);
//#define _VERBOSE 1


template <class T> complex<T> R8g254(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);/*Rallmp(ep,mpc);*/}  // 254
template <class T> complex<T> R8g1(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);/*Rallpm(ep,mpc);*/}   //1
template <class T> complex<T> R8g0(const eval_param<T>& ep,const mass_param_coll& mpc){return Rallm(ep,mpc);}   //0
template <class T> complex<T> R8g255(const eval_param<T>& ep,const mass_param_coll& mpc){return Rallp(ep,mpc);}   //255

//#endif

template <class T> complex<T> (*R8g_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
//#if _FAST_COMPILE == 0
	case 0  : return &R8g0;
	case 1  : return 0;//&R8g1;
	case 254: return 0;//&R8g254;
	case 255: return &R8g255;
//#endif
	default : return 0;
	}
}
template <class T> complex<T> R8g(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

switch (helcode_g(p)) {
//#if _FAST_COMPILE == 0
case 0  : return R8g0(ep,mpc);
case 1  : return R8g1(ep,mpc);
case 254: return R8g254(ep,mpc);
case 255: return R8g255(ep,mpc);
//#endif
default: {_WARNING3("using unknown known rational term for ",p," returned 0;"); return complex<T>(0.,0.);};

}
}




template complex<R> R8g(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
template complex<RHP> R8g(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
template complex<RVHP> R8g(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);

template complex<R> (*R8g_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R8g_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R8g_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R8g_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

}

