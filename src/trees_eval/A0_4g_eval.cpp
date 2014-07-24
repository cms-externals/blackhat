

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"

using namespace std;

 namespace BH {

template <class T> complex<T> A4g0_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }


template <class T> complex<T> A4g12_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,3)*ep.spa(1,2)*ep.spa(2,3)); }

template <class T> complex<T> A4g10_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,2),4))/(ep.spa(0,1)*ep.spa(0,3)*ep.spa(1,2)*ep.spa(2,3)); }

template <class T> complex<T> A4g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)); }

template <class T> complex<T> A4g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,2),3))/(ep.spa(0,1)*ep.spa(0,3)*ep.spa(2,3)); }

template <class T> complex<T> A4g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,3),4))/(ep.spa(0,1)*ep.spa(0,3)*ep.spa(1,2)*ep.spa(2,3)); }

template <class T> complex<T> A4g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,3),3))/(ep.spa(0,1)*ep.spa(0,3)*ep.spa(1,2)); }


template <class T> complex<T> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 12 : return &A4g12_eval;
	case 10 : return &A4g10_eval;
	case 6 : return &A4g6_eval;
	case 9 : return &A4g9_eval;
	case 5 : return &A4g5_eval;
	case 3 : return &A4g3_eval;
	default: return &A4g0_eval;}
}


template complex<R> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&) ;
 template complex<RHP> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&) ;
 template complex<RVHP> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*A4g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&) ;
#endif



}

