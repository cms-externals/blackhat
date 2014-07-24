

#include "R_4g_eval.h"
#include "eval_param.h"
#include "tree_amp.h"

using namespace std;
namespace BH  {

int helcode_g(const process& pro);

#define _ONLY_X_PART 1
#define X 2
#define _VERBOSE 0


template <class T> complex<T> R4g0(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : ----");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	-complex<T>(0,1)/complex<T>(3,0)*ep.spa(0,1)*ep.spa(2,3)/(ep.spb(0,1)*ep.spb(2,3))
#endif
;
}

template <class T> complex<T> R4g15(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : ++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	-complex<T>(0,1)/complex<T>(3,0)*ep.spb(0,1)*ep.spb(2,3)/(ep.spa(0,1)*ep.spa(2,3))
#endif
;
}

template <int i1, int i2, class T> complex<T> R4g2m2p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : --++");
#endif
	return complex<T>(X,0)/complex<T>(9,0)*A04g_eval<i1,i2>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g2m2p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : -+-+");
#endif
	return complex<T>(X,0)/complex<T>(9,0)*A04g_eval<i1,i3>(ep,mpc)
#if _ONLY_X_PART
	-complex<T>(0,1)*pow(ep.spa(i1,i3),2)*ep.spb(i1,i2)*ep.spb(i2,i3)/(ep.spa(i3,i4)*ep.spa(i4,i1)*pow(ep.spb(i1,i3),2))
#endif
;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : +---");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	+complex<T>(0,1)/complex<T>(3,0)*ep.spb(i2,i4)*pow(ep.spa(i2,i4),3)/(ep.spa(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i1))
#endif
	;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : -+++");
#endif

	return complex<T>(0,0)
	#if _ONLY_X_PART
	+ complex<T>(0,1)/complex<T>(3,0)*ep.spa(i2,i4)*pow(ep.spb(i2,i4),3)/(ep.spb(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i1))
#endif
;
}

//------------------  NF   -----------------------------
template <class T> complex<T> R4g0_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: ----");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	+complex<T>(0,1)/complex<T>(3,0)*ep.spa(0,1)*ep.spa(2,3)/(ep.spb(0,1)*ep.spb(2,3))
#endif
;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g2m2p_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: -+-+");
#endif

	return complex<T>(-T(2)/T(9),0)*A04g_eval<i1,i3>(ep,mpc)
	#if _ONLY_X_PART

+(ep.s(i1,i2)*ep.s(i2,i3)/ep.s(i1,i3)/ep.s(i1,i3))*A04g_eval<i1,i3>(ep,mpc);
#endif
;
}

template <int i1,int i2,class T> complex<T> R4g2m2p_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: --++");
#endif

	return complex<T>(-T(2)/T(9),0)*A04g_eval<i1,i2>(ep,mpc);
;
}

template <class T> complex<T> R4g15_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: ++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	+complex<T>(0,1)/complex<T>(3,0)*ep.spb(0,1)*ep.spb(2,3)/(ep.spa(0,1)*ep.spa(2,3))
#endif
;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1p_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: +---");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	-complex<T>(0,1)/complex<T>(3,0)*ep.spb(i2,i4)*pow(ep.spa(i2,i4),3)/(ep.spa(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i1))
#endif
	;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1m_nf(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 nf: -+++");
#endif

	return complex<T>(0,0)
	#if _ONLY_X_PART
	- complex<T>(0,1)/complex<T>(3,0)*ep.spa(i2,i4)*pow(ep.spb(i2,i4),3)/(ep.spb(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i1))
#endif
;
}

//----------  LT  --------------------------------------------

template <class T> complex<T> R4g0_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT: ----");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	-complex<T>(0,1)/complex<T>(3,0)*ep.spa(0,1)*ep.spa(2,3)/(ep.spb(0,1)*ep.spb(2,3))
#endif
;
}

template <class T> complex<T> R4g15_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT : ++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
	-complex<T>(0,1)/complex<T>(3,0)*ep.spb(0,1)*ep.spb(2,3)/(ep.spa(0,1)*ep.spa(2,3))
#endif
;
}

template <int i1, int i2, class T> complex<T> R4g2m2p_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT: --++");
#endif
	return complex<T>(1,0)/complex<T>(2,0)*A04g_eval<i1,i2>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g2m2p_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT: -+-+");
#endif
	return complex<T>(0,0);//(complex<T>(-1,0)/complex<T>(2,0) - ep.s(i1,i2)*ep.s(i2,i3)/ep.s(i1,i3)/ep.s(i1,i3))*A04g_eval<i1,i3>(ep,mpc);//R4g2m2p<i1,i2,i3,i4,T>(ep,mpc);
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1p_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT: +---");
#endif
	return complex<T>(0,0)
#if _ONLY_X_PART
	+complex<T>(0,1)/complex<T>(3,0)*ep.spb(i2,i4)*pow(ep.spa(i2,i4),3)/(ep.spa(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i1))
#endif
	;
}

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1m_LT(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 : -+++");
#endif

	return complex<T>(0,0)
	#if _ONLY_X_PART
	+ complex<T>(0,1)/complex<T>(3,0)*ep.spa(i2,i4)*pow(ep.spb(i2,i4),3)/(ep.spb(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i1))
#endif
;
}
//
////----------  RT  --------------------------------------------
//
//template <class T> complex<T> R4g0_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 RT: ----");
//#endif
//	return complex<T>(0,0)
//#if _ONLY_X_PART
//	-complex<T>(0,1)/complex<T>(3,0)*ep.spa(0,1)*ep.spa(2,3)/(ep.spb(0,1)*ep.spb(2,3))
//#endif
//;
//}
//
//template <class T> complex<T> R4g15_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 RT : ++++");
//#endif
//	return complex<T>(0,0)
//	#if _ONLY_X_PART
//	-complex<T>(0,1)/complex<T>(3,0)*ep.spb(0,1)*ep.spb(2,3)/(ep.spa(0,1)*ep.spa(2,3))
//#endif
//;
//}
//
//template <int i1, int i2, class T> complex<T> R4g2m2p_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 RT: --++");
//#endif
//	return complex<T>(1,0)/complex<T>(2,0)*A04g_eval<i1,i2>(ep,mpc);
//}
//
//template <int i1, int i2, int i3, int i4, class T> complex<T> R4g2m2p_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 RT: -+-+");
//#endif
//	return (complex<T>(-1,0)/complex<T>(2,0) - ep.s(i1,i2)*ep.s(i2,i3)/ep.s(i1,i3)/ep.s(i1,i3))*A04g_eval<i1,i3>(ep,mpc);//R4g2m2p<i1,i2,i3,i4,T>(ep,mpc);
//}
//
//template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1p_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 RT: +---");
//#endif
//	return complex<T>(0,0)
//#if _ONLY_X_PART
//	+complex<T>(0,1)/complex<T>(3,0)*ep.spb(i2,i4)*pow(ep.spa(i2,i4),3)/(ep.spa(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spa(i4,i1))
//#endif
//	;
//}
//
//template <int i1, int i2, int i3, int i4, class T> complex<T> R4g1m_RT(const eval_param<T>& ep,const mass_param_coll& mpc)
//{
//#if _VERBOSE
//	_MESSAGE("R4 : -+++");
//#endif
//
//	return complex<T>(0,0)
//	#if _ONLY_X_PART
//	+ complex<T>(0,1)/complex<T>(3,0)*ep.spa(i2,i4)*pow(ep.spb(i2,i4),3)/(ep.spb(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spb(i4,i1))
//#endif
//;
//}
//

template <int i1, int i2, int i3, int i4, class T> complex<T> R4g2m2p_SLC(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R4 LT: -+-+");
#endif
	return (complex<T>(-1,0)/complex<T>(2,0) )*A04g_eval<i1,i3>(ep,mpc);//R4g2m2p<i1,i2,i3,i4,T>(ep,mpc);
}
//R4g0 and R4g15 already in rat_ampl;

template <class T> complex<T> R4g1(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p<0,1,2,3>(ep,mpc);}
template <class T> complex<T> R4g2(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p<1,2,3,0>(ep,mpc);}
template <class T> complex<T> R4g3(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<2,3>(ep,mpc);}
template <class T> complex<T> R4g4(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p<2,3,0,1>(ep,mpc);}
template <class T> complex<T> R4g5(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<3,0,1,2>(ep,mpc);}
template <class T> complex<T> R4g6(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<0,3>(ep,mpc);}
template <class T> complex<T> R4g7(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m<3,0,1,2>(ep,mpc);}
template <class T> complex<T> R4g8(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p<3,0,1,2>(ep,mpc);}
template <class T> complex<T> R4g9(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<1,2>(ep,mpc);}
template <class T> complex<T> R4g10(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<0,1,2,3>(ep,mpc);}
template <class T> complex<T> R4g11(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m<2,3,0,1>(ep,mpc);}
template <class T> complex<T> R4g12(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p<0,1>(ep,mpc);}
template <class T> complex<T> R4g13(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m<1,2,3,0>(ep,mpc);}
template <class T> complex<T> R4g14(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m<0,1,2,3>(ep,mpc);}


template <class T> complex<T> (*R4g_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
	case 0 : return &R4g0;
	case 1 : return &R4g1;
	case 2 : return &R4g2;
	case 3 : return &R4g3;
	case 4 : return &R4g4;
	case 5 : return &R4g5;
	case 6 : return &R4g6;
	case 7 : return &R4g7;
	case 8 : return &R4g8;
	case 9 : return &R4g9;
	case 10 : return &R4g10;
	case 11 : return &R4g11;
	case 12 : return &R4g12;
	case 13 : return &R4g13;
	case 14 : return &R4g14;
	case 15 : return &R4g15;
	}
}

template <class T> complex<T> R4g(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

switch (helcode_g(p)) {
case 0 : return R4g0(ep,mpc);
case 1 : return R4g1(ep,mpc);
case 2 : return R4g2(ep,mpc);
case 3 : return R4g3(ep,mpc);
case 4 : return R4g4(ep,mpc);
case 5 : return R4g5(ep,mpc);
case 6 : return R4g6(ep,mpc);
case 7 : return R4g7(ep,mpc);
case 8 : return R4g8(ep,mpc);
case 9 : return R4g9(ep,mpc);
case 10 : return R4g10(ep,mpc);
case 11 : return R4g11(ep,mpc);
case 12 : return R4g12(ep,mpc);
case 13 : return R4g13(ep,mpc);
case 14 : return R4g14(ep,mpc);
case 15 : return R4g15(ep,mpc);
default: return complex<T>(0,0);

}
}






// --- NF ------


// template <class T> complex<T> R4g1_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g2_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g3_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m_nf<2,3>(ep,mpc);}
// template <class T> complex<T> R4g4_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g5_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<3,0,1,2>(ep,mpc);}
// template <class T> complex<T> R4g6_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m_nf<0,3>(ep,mpc);}
// template <class T> complex<T> R4g7_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g8_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g9_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m_nf<1,2>(ep,mpc);}
// template <class T> complex<T> R4g10_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<0,1,2,3>(ep,mpc);}
// template <class T> complex<T> R4g11_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g12_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m_nf<0,1>(ep,mpc);}
// template <class T> complex<T> R4g13_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g14_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}
// template <class T> complex<T> R4g15_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g0_nf(ep,mpc);}

 template <class T> complex<T> R4g1_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_nf<0,1,2,3>(ep,mpc);}
 template <class T> complex<T> R4g2_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_nf<1,2,3,0>(ep,mpc);}
 template <class T> complex<T> R4g3_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<2,3>(ep,mpc);}
 template <class T> complex<T> R4g4_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_nf<2,3,0,1>(ep,mpc);}
 template <class T> complex<T> R4g5_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g6_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<0,3>(ep,mpc);}
 template <class T> complex<T> R4g7_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_nf<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g8_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_nf<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g9_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<1,2>(ep,mpc);}
 template <class T> complex<T> R4g10_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<0,1,2,3>(ep,mpc);}
 template <class T> complex<T> R4g11_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_nf<2,3,0,1>(ep,mpc);}
 template <class T> complex<T> R4g12_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_nf<0,1>(ep,mpc);}
 template <class T> complex<T> R4g13_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_nf<1,2,3,0>(ep,mpc);}
 template <class T> complex<T> R4g14_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_nf<0,1,2,3>(ep,mpc);}



 template <class T> complex<T> (*R4g_nf_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
 	switch (hc) {
 	case 0 : return &R4g0_nf;
 	case 1 : return &R4g1_nf;
 	case 2 : return &R4g2_nf;
 	case 3 : return &R4g3_nf;
 	case 4 : return &R4g4_nf;
 	case 5 : return &R4g5_nf;
 	case 6 : return &R4g6_nf;
 	case 7 : return &R4g7_nf;
 	case 8 : return &R4g8_nf;
 	case 9 : return &R4g9_nf;
 	case 10 : return &R4g10_nf;
 	case 11 : return &R4g11_nf;
 	case 12 : return &R4g12_nf;
 	case 13 : return &R4g13_nf;
 	case 14 : return &R4g14_nf;
 	case 15 : return &R4g15_nf;
 }
}

 template <class T> complex<T> R4g_nf(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

 switch (helcode_g(p)) {
 case 0 : return R4g0(ep,mpc);
 case 1 : return R4g1(ep,mpc);
 case 2 : return R4g2(ep,mpc);
 case 3 : return R4g3(ep,mpc);
 case 4 : return R4g4(ep,mpc);
 case 5 : return R4g5(ep,mpc);
 case 6 : return R4g6(ep,mpc);
 case 7 : return R4g7(ep,mpc);
 case 8 : return R4g8(ep,mpc);
 case 9 : return R4g9(ep,mpc);
 case 10 : return R4g10(ep,mpc);
 case 11 : return R4g11(ep,mpc);
 case 12 : return R4g12(ep,mpc);
 case 13 : return R4g13(ep,mpc);
 case 14 : return R4g14(ep,mpc);
 case 15 : return R4g15(ep,mpc);
 default: return complex<T>(0,0);

 }
 }


 template <class T> complex<T> R4g1_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_LT<0,1,2,3>(ep,mpc);}
 template <class T> complex<T> R4g2_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_LT<1,2,3,0>(ep,mpc);}
 template <class T> complex<T> R4g3_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<2,3>(ep,mpc);}
 template <class T> complex<T> R4g4_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_LT<2,3,0,1>(ep,mpc);}
 template <class T> complex<T> R4g5_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g6_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<0,3>(ep,mpc);}
 template <class T> complex<T> R4g7_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_LT<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g8_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1p_LT<3,0,1,2>(ep,mpc);}
 template <class T> complex<T> R4g9_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<1,2>(ep,mpc);}
 template <class T> complex<T> R4g10_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<0,1,2,3>(ep,mpc);}
 template <class T> complex<T> R4g11_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_LT<2,3,0,1>(ep,mpc);}
 template <class T> complex<T> R4g12_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g2m2p_LT<0,1>(ep,mpc);}
 template <class T> complex<T> R4g13_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_LT<1,2,3,0>(ep,mpc);}
 template <class T> complex<T> R4g14_LT(const eval_param<T>& ep,const mass_param_coll& mpc){return R4g1m_LT<0,1,2,3>(ep,mpc);}


 template <class T> complex<T> (*R4g_LT_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
 	switch (hc) {
 	case 0 : return &R4g0_LT;
 	case 1 : return &R4g1_LT;
 	case 2 : return &R4g2_LT;
 	case 3 : return &R4g3_LT;
 	case 4 : return &R4g4_LT;
 	case 5 : return &R4g5_LT;
 	case 6 : return &R4g6_LT;
 	case 7 : return &R4g7_LT;
 	case 8 : return &R4g8_LT;
 	case 9 : return &R4g9_LT;
 	case 10 : return &R4g10_LT;
 	case 11 : return &R4g11_LT;
 	case 12 : return &R4g12_LT;
 	case 13 : return &R4g13_LT;
 	case 14 : return &R4g14_LT;
 	case 15 : return &R4g15_LT;
 	}
 }

 template <class T> complex<T> R4g_LT(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

 switch (helcode_g(p)) {
 case 0 : return R4g0_LT(ep,mpc);
 case 1 : return R4g1_LT(ep,mpc);
 case 2 : return R4g2_LT(ep,mpc);
 case 3 : return R4g3_LT(ep,mpc);
 case 4 : return R4g4_LT(ep,mpc);
 case 5 : return R4g5_LT(ep,mpc);
 case 6 : return R4g6_LT(ep,mpc);
 case 7 : return R4g7_LT(ep,mpc);
 case 8 : return R4g8_LT(ep,mpc);
 case 9 : return R4g9_LT(ep,mpc);
 case 10 : return R4g10_LT(ep,mpc);
 case 11 : return R4g11_LT(ep,mpc);
 case 12 : return R4g12_LT(ep,mpc);
 case 13 : return R4g13_LT(ep,mpc);
 case 14 : return R4g14_LT(ep,mpc);
 case 15 : return R4g15_LT(ep,mpc);
 default: return complex<T>(0,0);

 }
 }

 template <class T> complex<T> R4g1_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g1p_LT<0,1,2,3>(ep,mpc);}
  template <class T> complex<T> R4g2_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //complex<T>(0,0);} //R4g1p_LT<1,2,3,0>(ep,mpc);}
  template <class T> complex<T> R4g3_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //-R4g2m2p_LT<2,3>(ep,mpc);}
  template <class T> complex<T> R4g4_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g1p_LT<2,3,0,1>(ep,mpc);}
  template <class T> complex<T> R4g5_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g2m2p_SLC<3,0,1,2>(ep,mpc);}
  template <class T> complex<T> R4g6_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //-R4g2m2p_LT<0,3>(ep,mpc);}
  template <class T> complex<T> R4g7_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g1m_LT<3,0,1,2>(ep,mpc);}
  template <class T> complex<T> R4g8_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} // R4g1p_LT<3,0,1,2>(ep,mpc);}
  template <class T> complex<T> R4g9_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //-R4g2m2p_LT<1,2>(ep,mpc);}
  template <class T> complex<T> R4g10_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} // R4g2m2p_SLC<0,1,2,3>(ep,mpc);}
  template <class T> complex<T> R4g11_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g1m_LT<2,3,0,1>(ep,mpc);}
  template <class T> complex<T> R4g12_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //-R4g2m2p_LT<0,1>(ep,mpc);}
  template <class T> complex<T> R4g13_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} // R4g1m_LT<1,2,3,0>(ep,mpc);}
  template <class T> complex<T> R4g14_SLC(const eval_param<T>& ep,const mass_param_coll& mpc){return complex<T>(0,0);} //R4g1m_LT<0,1,2,3>(ep,mpc);}


  template <class T> complex<T> (*R4g_SLC_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
  	switch (hc) {
  	case 0 : return &R4g0_LT;
  	case 1 : return &R4g1_SLC;
  	case 2 : return &R4g2_SLC;
  	case 3 : return &R4g3_SLC;
  	case 4 : return &R4g4_SLC;
  	case 5 : return &R4g5_SLC;
  	case 6 : return &R4g6_SLC;
  	case 7 : return &R4g7_SLC;
  	case 8 : return &R4g8_SLC;
  	case 9 : return &R4g9_SLC;
  	case 10 : return &R4g10_SLC;
  	case 11 : return &R4g11_SLC;
  	case 12 : return &R4g12_SLC;
  	case 13 : return &R4g13_SLC;
  	case 14 : return &R4g14_SLC;
  	case 15 : return &R4g15_LT;
  	}
  }

  template <class T> complex<T> R4g_SLC(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

  switch (helcode_g(p)) {
  case 0 : return R4g0_LT(ep,mpc);
  case 1 : return R4g1_SLC(ep,mpc);
  case 2 : return R4g2_SLC(ep,mpc);
  case 3 : return R4g3_SLC(ep,mpc);
  case 4 : return R4g4_SLC(ep,mpc);
  case 5 : return R4g5_SLC(ep,mpc);
  case 6 : return R4g6_SLC(ep,mpc);
  case 7 : return R4g7_SLC(ep,mpc);
  case 8 : return R4g8_SLC(ep,mpc);
  case 9 : return R4g9_SLC(ep,mpc);
  case 10 : return R4g10_SLC(ep,mpc);
  case 11 : return R4g11_SLC(ep,mpc);
  case 12 : return R4g12_SLC(ep,mpc);
  case 13 : return R4g13_SLC(ep,mpc);
  case 14 : return R4g14_SLC(ep,mpc);
  case 15 : return R4g15_LT(ep,mpc);
  default: return complex<T>(0,0);

  }
  }


 template complex<R> R4g(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
 template complex<RHP> R4g(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
 template complex<RVHP> R4g(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);

 template complex<R> (*R4g_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R4g_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R4g_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R4g_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif


 template complex<R> R4g_nf(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
 template complex<RHP> R4g_nf(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
 template complex<RVHP> R4g_nf(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);

 template complex<R> (*R4g_nf_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
  template complex<RHP> (*R4g_nf_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
  template complex<RVHP> (*R4g_nf_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

  template complex<RGMP> (*R4g_nf_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif




  template complex<R> R4g_LT(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
  template complex<RHP> R4g_LT(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
  template complex<RVHP> R4g_LT(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);

  template complex<R> (*R4g_LT_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
   template complex<RHP> (*R4g_LT_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
   template complex<RVHP> (*R4g_LT_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

   template complex<RGMP> (*R4g_LT_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

   template complex<R> R4g_SLC(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
    template complex<RHP> R4g_SLC(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
    template complex<RVHP> R4g_SLC(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);

    template complex<R> (*R4g_SLC_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
     template complex<RHP> (*R4g_SLC_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
     template complex<RVHP> (*R4g_SLC_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

     template complex<RGMP> (*R4g_SLC_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif


}

