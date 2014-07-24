

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"

using namespace std;

namespace BH {


template <class T> complex<T> A5g0_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }


template <class T> complex<T> A5g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(1,0),3))/(ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(2,0),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(2,1),3))/(ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g7_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(3,4),3))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)); }

template <class T> complex<T> A5g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,0),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g10_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,1),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g11_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,4),4))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g12_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,2),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g13_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,4),4))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g14_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,4),3))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g17_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,0),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)); }

template <class T> complex<T> A5g18_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,1),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g19_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,3),3))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4)); }

template <class T> complex<T> A5g20_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,2),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3)); }

template <class T> complex<T> A5g21_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,3),4))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g22_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,3),4))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g24_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,3),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)); }

template <class T> complex<T> A5g25_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,2),3))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g26_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,2),4))/(ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A5g28_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 3 : return &A5g3_eval;
	case 5 : return &A5g5_eval;
	case 6 : return &A5g6_eval;
	case 7 : return &A5g7_eval;
	case 9 : return &A5g9_eval;
	case 10 : return &A5g10_eval;
	case 11 : return &A5g11_eval;
	case 12 : return &A5g12_eval;
	case 13 : return &A5g13_eval;
	case 14 : return &A5g14_eval;
	case 17 : return &A5g17_eval;
	case 18 : return &A5g18_eval;
	case 19 : return &A5g19_eval;
	case 20 : return &A5g20_eval;
	case 21 : return &A5g21_eval;
	case 22 : return &A5g22_eval;
	case 24 : return &A5g24_eval;
	case 25 : return &A5g25_eval;
	case 26 : return &A5g26_eval;
	case 28 : return &A5g28_eval;
	default: return &A5g0_eval;
	}
}



template complex<R> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&) ;
 template complex<RHP> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&) ;
 template complex<RVHP> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*A5g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&) ;
#endif


}

