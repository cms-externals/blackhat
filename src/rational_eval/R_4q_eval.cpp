

#include "R_4g_eval.h"
#include "eval_param.h" 

using namespace std;
namespace BH  {


#define _ONLY_X_PART 1
#define X 2.
#define _VERBOSE 0
template <class T> complex<T> R4q0(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4q : ----");
#endif
	return complex<T>(0,0);
}

#define _CASE_R4q_L_Ptr_eval(K) case K : return &R4q_L_ ## K

template <class T> complex<T> (*R4q_L_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {


default: return &R4q0;
	}
}




template complex<R> (*R4q_L_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R4q_L_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R4q_L_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R4q_L_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

}

