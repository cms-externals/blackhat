/*
 * A0_2q2Q2l_eval.cpp
 *
 * created with Harald's scripts from trees computed by Lance
 * notes:
 * *) trees are manifest dual super conformal expressions
 * *) speed-up from idendifying commone sumbexpr. still possible
 *  
 *
 */

/* Implementation of super dual conformal trees  */

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"

using namespace std;


#define SPA(i,j) ep.spa(i,j)
#define SPB(i,j) ep.spb(i,j)

/* 
 * more spab & spaa  macros 
*/
namespace BH {

template<class T> static inline complex<T> square(complex<T> x)
{return(x*x);}
template<class T> static inline complex<T> cube(complex<T> x)
{return(x*x*x);}



//template <class T> complex<T>  ZeroF(const eval_param<T>& ep, const mass_param_coll& masses){return complex<T>(0,0);}

/*
 *
 *
 * The 2q2Q quarks  2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2Q2l_qmq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb24)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa13*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2Q2l_qmqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb34)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa12*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2Q2l_qpq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb34)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa12*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2Q2l_qpqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb24)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa13*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2Q2l_qmq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb24)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa13*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2Q2l_qmqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb34)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa12*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2Q2l_qpq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb34)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa12*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2Q2l_qpqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);

return(
((pow2_spa15*pow2_spb24)/(s234*spa56*spab_1_23_4*spb23)
-
(pow2_spa13*pow2_spb46)/(s123*spa23*spab_1_23_4*spb56))*complex<T>(0,1)
);
}




template <class T> complex<T> (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 225433:	 return &A2q2Q2l_qpq2pqb2mqbmlplbm_eval<0,1,2,3,4,5>;//qp q2p qb2m qbm lp lbm
case 225489:	 return &A2q2Q2l_qpqb2mq2pqbmlplbm_eval<0,1,2,3,4,5>;//qp qb2m q2p qbm lp lbm
case 225944:	 return &A2q2Q2l_qmq2pqb2mqbplplbm_eval<0,1,2,3,4,5>;//qm q2p qb2m qbp lp lbm
case 226000:	 return &A2q2Q2l_qmqb2mq2pqbplplbm_eval<0,1,2,3,4,5>;//qm qb2m q2p qbp lp lbm
case 254105:	 return &A2q2Q2l_qpq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5>;//qp q2p qb2m qbm lm lbp
case 254161:	 return &A2q2Q2l_qpqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5>;//qp qb2m q2p qbm lm lbp
case 254616:	 return &A2q2Q2l_qmq2pqb2mqbplmlbp_eval<0,1,2,3,4,5>;//qm q2p qb2m qbp lm lbp
case 254672:	 return &A2q2Q2l_qmqb2mq2pqbplmlbp_eval<0,1,2,3,4,5>;//qm qb2m q2p qbp lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q2Q2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A2q2Q2l_Tree_Ptr_eval(long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif
}
