#define SPA(i,j) ep.spa(ind.at(i-1),ep.spa(ind.at(j-1)))

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"


using namespace std;

namespace BH {

template <class T> complex<T> A3g0_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A3g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spa(0,1),3))/(ep.spa(0,2)*ep.spa(1,2)); }

template <class T> complex<T> A3g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*ep.spa(1,2)); }

template <class T> complex<T> A3g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spb(2,1),3))/(ep.spb(1,0)*ep.spb(2,0)); }

template <class T> complex<T> A3g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spa(1,2),3))/(ep.spa(0,1)*ep.spa(0,2)); }

template <class T> complex<T> A3g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*ep.spb(2,1)); }

template <class T> complex<T> A3g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0.,-1.)*pow(ep.spb(1,0),3))/(ep.spb(2,0)*ep.spb(2,1)); }






template <class T> complex<T> A3g_eval(int hc,const eval_param<T>& ep, const mass_param_coll& masses){
switch (hc) {
case 4 : return A3g4_eval(ep,masses);
case 2 : return A3g2_eval(ep,masses);
case 6 : return A3g6_eval(ep,masses);
case 1 : return A3g1_eval(ep,masses);
case 5 : return A3g5_eval(ep,masses);
case 3 : return A3g3_eval(ep,masses);

default: return complex<T>(0.,0.);

}
}

template <class T> complex<T> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	 switch (hc) {
	 case 4 : return &A3g4_eval;
	 case 2 : return &A3g2_eval;
	 case 6 : return &A3g6_eval;
	 case 1 : return &A3g1_eval;
	 case 5 : return &A3g5_eval;
	 case 3 : return &A3g3_eval;
	 default: return &A3g0_eval;
	 }
}

template complex<R> A3g_eval(int hc,const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP> A3g_eval(int hc,const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP> A3g_eval(int hc,const eval_param<RVHP>& ep, const mass_param_coll& masses);

template complex<R> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&) ;
 template complex<RHP> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&) ;
 template complex<RVHP> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*A3g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&) ;
#endif



}


