/*
 * A0_2q2g2l_eval.cpp
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
 * The 2 quarks 2g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2g2l_qmmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spab_5_12_3=SPA(i1,i5)*spb31+SPA(i2,i5)*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_34_6=-(spa23*SPB(i6,i3))-spa24*SPB(i6,i4);

return(
complex<T>(0,1)*(-((pow(spab_5_12_3,2)*SPB(i1,i3))/(spab_4_23_1*(spa23*spb32+spb21*SPA(i1,i2
)+
spb31*SPA(i1,i3))*SPA(i5,i6)*SPB(i1,i2)*SPB(i2,i3))
)-
(spa24*pow(spab_2_34_6,2))/(spa23*spa34*spab_4_23_1*(spa23*spb32+spa24*SPB(i4,i2
)+
spa34*SPB(i4,i3))*SPB(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2g2l_qmmpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spab_5_12_3=SPA(i1,i5)*spb31+SPA(i2,i5)*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_34_6=-(spa23*SPB(i6,i3))-spa24*SPB(i6,i4);

return(
complex<T>(0,-1)*(-((pow(spab_5_12_3,2)*SPB(i1,i3))/(spab_4_23_1*(spa23*spb32+spb21*SPA(i1,i2
)+
spb31*SPA(i1,i3))*SPA(i5,i6)*SPB(i1,i2)*SPB(i2,i3))
)-
(spa24*pow(spab_2_34_6,2))/(spa23*spa34*spab_4_23_1*(spa23*spb32+spa24*SPB(i4,i2
)+
spa34*SPB(i4,i3))*SPB(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2g2l_qmpmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;

return(
complex<T>(0,1)*(-((pow(SPA(i1,i5),2)*pow(SPB(i2,i4),3))/(spab_1_23_4*(spa23*spb32+spb42*SPA(i2,i4
)+
spb43*SPA(i3,i4))*SPA(i5,i6)*SPB(i2,i3)*SPB(i3,i4))
)-
(pow(spa13,3)*pow(SPB(i4,i6),2))/(spa12*spa23*spab_1_23_4*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1))*SPB(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2g2l_qmpmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;

return(
complex<T>(0,-1)*(-((pow(SPA(i1,i5),2)*pow(SPB(i2,i4),3))/(spab_1_23_4*(spa23*spb32+spb42*SPA(i2,i4
)+
spb43*SPA(i3,i4))*SPA(i5,i6)*SPB(i2,i3)*SPB(i3,i4))
)-
(pow(spa13,3)*pow(SPB(i4,i6),2))/(spa12*spa23*spab_1_23_4*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1))*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2g2l_qpmpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;

return(
complex<T>(0,-1)*(-((pow(SPA(i1,i5),2)*pow(SPB(i2,i4),3))/(spab_1_23_4*(spa23*spb32+spb42*SPA(i2,i4
)+
spb43*SPA(i3,i4))*SPA(i5,i6)*SPB(i2,i3)*SPB(i3,i4))
)-
(pow(spa13,3)*pow(SPB(i4,i6),2))/(spa12*spa23*spab_1_23_4*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1))*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2g2l_qpmpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;

return(
complex<T>(0,1)*(-((pow(SPA(i1,i5),2)*pow(SPB(i2,i4),3))/(spab_1_23_4*(spa23*spb32+spb42*SPA(i2,i4
)+
spb43*SPA(i3,i4))*SPA(i5,i6)*SPB(i2,i3)*SPB(i3,i4))
)-
(pow(spa13,3)*pow(SPB(i4,i6),2))/(spa12*spa23*spab_1_23_4*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1))*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2g2l_qppmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spab_5_12_3=SPA(i1,i5)*spb31+SPA(i2,i5)*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_34_6=-(spa23*SPB(i6,i3))-spa24*SPB(i6,i4);

return(
complex<T>(0,-1)*(-((pow(spab_5_12_3,2)*SPB(i1,i3))/(spab_4_23_1*(spa23*spb32+spb21*SPA(i1,i2
)+
spb31*SPA(i1,i3))*SPA(i5,i6)*SPB(i1,i2)*SPB(i2,i3))
)-
(spa24*pow(spab_2_34_6,2))/(spa23*spa34*spab_4_23_1*(spa23*spb32+spa24*SPB(i4,i2
)+
spa34*SPB(i4,i3))*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2g2l_qppmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spab_5_12_3=SPA(i1,i5)*spb31+SPA(i2,i5)*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_34_6=-(spa23*SPB(i6,i3))-spa24*SPB(i6,i4);

return(
complex<T>(0,1)*(-((pow(spab_5_12_3,2)*SPB(i1,i3))/(spab_4_23_1*(spa23*spb32+spb21*SPA(i1,i2
)+
spb31*SPA(i1,i3))*SPA(i5,i6)*SPB(i1,i2)*SPB(i2,i3))
)-
(spa24*pow(spab_2_34_6,2))/(spa23*spa34*spab_4_23_1*(spa23*spb32+spa24*SPB(i4,i2
)+
spa34*SPB(i4,i3))*SPB(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2g2l_qmppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i5),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T>  A2q2g2l_qmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPB(i4,i6),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2g2l_qpppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i5),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2g2l_qmppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i5),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, class T> complex<T>  A2q2g2l_qpmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPB(i4,i6),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i6)))
);
}


template <int i1, int i2, int i3, int i4, int i6, int i5, class T> complex<T>  A2q2g2l_qmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPB(i4,i6),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2g2l_qpppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i5),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i5,i6)))
);
}


template <int i4, int i3, int i2, int i1, int i6, int i5, class T> complex<T>  A2q2g2l_qpmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPB(i4,i6),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i5,i6)))
);
}




template <class T> complex<T> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

//case 11765:	 return &A2q2g2l_qpmpqbmlplbm_eval<2,3,4,5,0,1>;//lp lbm qp m p qbm
case 11770:	 return &A2q2g2l_qpmpqbmlmlbp_eval<2,3,4,5,0,1>;//lm lbp qp m p qbm
//case 12413:	 return &A2q2g2l_qpppqbmlplbm_eval<2,3,4,5,0,1>;//lp lbm qp p p qbm
case 12418:	 return &A2q2g2l_qpppqbmlmlbp_eval<2,3,4,5,0,1>;//lm lbp qp p p qbm
//case 14112:	 return &A2q2g2l_qmmmqbplplbm_eval<5,0,1,2,3,4>;//m m qbp lp lbm qm
//case 14115:	 return &A2q2g2l_qmpmqbplplbm_eval<5,0,1,2,3,4>;//p m qbp lp lbm qm
//case 14130:	 return &A2q2g2l_qmmpqbplplbm_eval<5,0,1,2,3,4>;//m p qbp lp lbm qm
//case 14133:	 return &A2q2g2l_qmppqbplplbm_eval<5,0,1,2,3,4>;//p p qbp lp lbm qm
case 15192:	 return &A2q2g2l_qmmmqbplmlbp_eval<5,0,1,2,3,4>;//m m qbp lm lbp qm
case 15195:	 return &A2q2g2l_qmpmqbplmlbp_eval<5,0,1,2,3,4>;//p m qbp lm lbp qm
case 15210:	 return &A2q2g2l_qmmpqbplmlbp_eval<5,0,1,2,3,4>;//m p qbp lm lbp qm
case 15213:	 return &A2q2g2l_qmppqbplmlbp_eval<5,0,1,2,3,4>;//p p qbp lm lbp qm
//case 15617:	 return &A2q2g2l_qmmmqbplplbm_eval<2,3,4,5,0,1>;//lp lbm qm m m qbp
case 15622:	 return &A2q2g2l_qmmmqbplmlbp_eval<2,3,4,5,0,1>;//lm lbp qm m m qbp
//case 16265:	 return &A2q2g2l_qmpmqbplplbm_eval<2,3,4,5,0,1>;//lp lbm qm p m qbp
case 16270:	 return &A2q2g2l_qmpmqbplmlbp_eval<2,3,4,5,0,1>;//lm lbp qm p m qbp
//case 19505:	 return &A2q2g2l_qmmpqbplplbm_eval<2,3,4,5,0,1>;//lp lbm qm m p qbp
case 19510:	 return &A2q2g2l_qmmpqbplmlbp_eval<2,3,4,5,0,1>;//lm lbp qm m p qbp
//case 20153:	 return &A2q2g2l_qmppqbplplbm_eval<2,3,4,5,0,1>;//lp lbm qm p p qbp
case 20158:	 return &A2q2g2l_qmppqbplmlbp_eval<2,3,4,5,0,1>;//lm lbp qm p p qbp
//case 21852:	 return &A2q2g2l_qpmmqbmlplbm_eval<5,0,1,2,3,4>;//m m qbm lp lbm qp
//case 21855:	 return &A2q2g2l_qppmqbmlplbm_eval<5,0,1,2,3,4>;//p m qbm lp lbm qp
//case 21870:	 return &A2q2g2l_qpmpqbmlplbm_eval<5,0,1,2,3,4>;//m p qbm lp lbm qp
//case 21873:	 return &A2q2g2l_qpppqbmlplbm_eval<5,0,1,2,3,4>;//p p qbm lp lbm qp
case 22932:	 return &A2q2g2l_qpmmqbmlmlbp_eval<5,0,1,2,3,4>;//m m qbm lm lbp qp
case 22935:	 return &A2q2g2l_qppmqbmlmlbp_eval<5,0,1,2,3,4>;//p m qbm lm lbp qp
case 22950:	 return &A2q2g2l_qpmpqbmlmlbp_eval<5,0,1,2,3,4>;//m p qbm lm lbp qp
case 22953:	 return &A2q2g2l_qpppqbmlmlbp_eval<5,0,1,2,3,4>;//p p qbm lm lbp qp
//case 2352:	 return &A2q2g2l_qmmmqbplplbm_eval<4,5,0,1,2,3>;//m qbp lp lbm qm m
//case 2355:	 return &A2q2g2l_qmmpqbplplbm_eval<4,5,0,1,2,3>;//p qbp lp lbm qm m
//case 23720:	 return &A2q2g2l_qmmpqbplplbm_eval<3,4,5,0,1,2>;//qbp lp lbm qm m p
case 23750:	 return &A2q2g2l_qmmpqbplmlbp_eval<3,4,5,0,1,2>;//qbp lm lbp qm m p
//case 23935:	 return &A2q2g2l_qpmpqbmlplbm_eval<3,4,5,0,1,2>;//qbm lp lbm qp m p
case 23965:	 return &A2q2g2l_qpmpqbmlmlbp_eval<3,4,5,0,1,2>;//qbm lm lbp qp m p
case 2532:	 return &A2q2g2l_qmmmqbplmlbp_eval<4,5,0,1,2,3>;//m qbp lm lbp qm m
case 2535:	 return &A2q2g2l_qmmpqbplmlbp_eval<4,5,0,1,2,3>;//p qbp lm lbp qm m
//case 25680:	 return &A2q2g2l_qmpmqbplplbm_eval<4,5,0,1,2,3>;//m qbp lp lbm qm p
//case 25683:	 return &A2q2g2l_qmppqbplplbm_eval<4,5,0,1,2,3>;//p qbp lp lbm qm p
case 25860:	 return &A2q2g2l_qmpmqbplmlbp_eval<4,5,0,1,2,3>;//m qbp lm lbp qm p
case 25863:	 return &A2q2g2l_qmppqbplmlbp_eval<4,5,0,1,2,3>;//p qbp lm lbp qm p
//case 26970:	 return &A2q2g2l_qppmqbmlplbm_eval<4,5,0,1,2,3>;//m qbm lp lbm qp p
//case 26973:	 return &A2q2g2l_qpppqbmlplbm_eval<4,5,0,1,2,3>;//p qbm lp lbm qp p
case 27150:	 return &A2q2g2l_qppmqbmlmlbp_eval<4,5,0,1,2,3>;//m qbm lm lbp qp p
case 27153:	 return &A2q2g2l_qpppqbmlmlbp_eval<4,5,0,1,2,3>;//p qbm lm lbp qp p
//case 27608:	 return &A2q2g2l_qmppqbplplbm_eval<3,4,5,0,1,2>;//qbp lp lbm qm p p
case 27638:	 return &A2q2g2l_qmppqbplmlbp_eval<3,4,5,0,1,2>;//qbp lm lbp qm p p
//case 27823:	 return &A2q2g2l_qpppqbmlplbm_eval<3,4,5,0,1,2>;//qbm lp lbm qp p p
case 27853:	 return &A2q2g2l_qpppqbmlmlbp_eval<3,4,5,0,1,2>;//qbm lm lbp qp p p
case 32417:	 return &A2q2g2l_qpmmqbmlmlbp_eval<1,2,3,4,5,0>;//lbp qp m m qbm lm
case 32525:	 return &A2q2g2l_qppmqbmlmlbp_eval<1,2,3,4,5,0>;//lbp qp p m qbm lm
case 33065:	 return &A2q2g2l_qpmpqbmlmlbp_eval<1,2,3,4,5,0>;//lbp qp m p qbm lm
case 33173:	 return &A2q2g2l_qpppqbmlmlbp_eval<1,2,3,4,5,0>;//lbp qp p p qbm lm
case 33707:	 return &A2q2g2l_qmmmqbplmlbp_eval<1,2,3,4,5,0>;//lbp qm m m qbp lm
case 33815:	 return &A2q2g2l_qmpmqbplmlbp_eval<1,2,3,4,5,0>;//lbp qm p m qbp lm
case 34355:	 return &A2q2g2l_qmmpqbplmlbp_eval<1,2,3,4,5,0>;//lbp qm m p qbp lm
case 34463:	 return &A2q2g2l_qmppqbplmlbp_eval<1,2,3,4,5,0>;//lbp qm p p qbp lm
//case 3642:	 return &A2q2g2l_qpmmqbmlplbm_eval<4,5,0,1,2,3>;//m qbm lp lbm qp m
//case 3645:	 return &A2q2g2l_qpmpqbmlplbm_eval<4,5,0,1,2,3>;//p qbm lp lbm qp m
//case 37802:	 return &A2q2g2l_qpmmqbmlplbm_eval<0,1,2,3,4,5>;//qp m m qbm lp lbm
//case 37820:	 return &A2q2g2l_qppmqbmlplbm_eval<0,1,2,3,4,5>;//qp p m qbm lp lbm
//case 37910:	 return &A2q2g2l_qpmpqbmlplbm_eval<0,1,2,3,4,5>;//qp m p qbm lp lbm
//case 37928:	 return &A2q2g2l_qpppqbmlplbm_eval<0,1,2,3,4,5>;//qp p p qbm lp lbm
//case 38017:	 return &A2q2g2l_qmmmqbplplbm_eval<0,1,2,3,4,5>;//qm m m qbp lp lbm
//case 38035:	 return &A2q2g2l_qmpmqbplplbm_eval<0,1,2,3,4,5>;//qm p m qbp lp lbm
//case 38125:	 return &A2q2g2l_qmmpqbplplbm_eval<0,1,2,3,4,5>;//qm m p qbp lp lbm
//case 38143:	 return &A2q2g2l_qmppqbplplbm_eval<0,1,2,3,4,5>;//qm p p qbp lp lbm
case 3822:	 return &A2q2g2l_qpmmqbmlmlbp_eval<4,5,0,1,2,3>;//m qbm lm lbp qp m
case 3825:	 return &A2q2g2l_qpmpqbmlmlbp_eval<4,5,0,1,2,3>;//p qbm lm lbp qp m
//case 392:	 return &A2q2g2l_qmmmqbplplbm_eval<3,4,5,0,1,2>;//qbp lp lbm qm m m
//case 40192:	 return &A2q2g2l_qpmmqbmlplbm_eval<1,2,3,4,5,0>;//lbm qp m m qbm lp
//case 40300:	 return &A2q2g2l_qppmqbmlplbm_eval<1,2,3,4,5,0>;//lbm qp p m qbm lp
//case 40840:	 return &A2q2g2l_qpmpqbmlplbm_eval<1,2,3,4,5,0>;//lbm qp m p qbm lp
//case 40948:	 return &A2q2g2l_qpppqbmlplbm_eval<1,2,3,4,5,0>;//lbm qp p p qbm lp
//case 41482:	 return &A2q2g2l_qmmmqbplplbm_eval<1,2,3,4,5,0>;//lbm qm m m qbp lp
//case 41590:	 return &A2q2g2l_qmpmqbplplbm_eval<1,2,3,4,5,0>;//lbm qm p m qbp lp
//case 42130:	 return &A2q2g2l_qmmpqbplplbm_eval<1,2,3,4,5,0>;//lbm qm m p qbp lp
//case 42238:	 return &A2q2g2l_qmppqbplplbm_eval<1,2,3,4,5,0>;//lbm qm p p qbp lp
case 422:	 return &A2q2g2l_qmmmqbplmlbp_eval<3,4,5,0,1,2>;//qbp lm lbp qm m m
//case 4280:	 return &A2q2g2l_qmpmqbplplbm_eval<3,4,5,0,1,2>;//qbp lp lbm qm p m
case 4310:	 return &A2q2g2l_qmpmqbplmlbp_eval<3,4,5,0,1,2>;//qbp lm lbp qm p m
case 44282:	 return &A2q2g2l_qpmmqbmlmlbp_eval<0,1,2,3,4,5>;//qp m m qbm lm lbp
case 44300:	 return &A2q2g2l_qppmqbmlmlbp_eval<0,1,2,3,4,5>;//qp p m qbm lm lbp
case 44390:	 return &A2q2g2l_qpmpqbmlmlbp_eval<0,1,2,3,4,5>;//qp m p qbm lm lbp
case 44408:	 return &A2q2g2l_qpppqbmlmlbp_eval<0,1,2,3,4,5>;//qp p p qbm lm lbp
case 44497:	 return &A2q2g2l_qmmmqbplmlbp_eval<0,1,2,3,4,5>;//qm m m qbp lm lbp
case 44515:	 return &A2q2g2l_qmpmqbplmlbp_eval<0,1,2,3,4,5>;//qm p m qbp lm lbp
case 44605:	 return &A2q2g2l_qmmpqbplmlbp_eval<0,1,2,3,4,5>;//qm m p qbp lm lbp
case 44623:	 return &A2q2g2l_qmppqbplmlbp_eval<0,1,2,3,4,5>;//qm p p qbp lm lbp
//case 4495:	 return &A2q2g2l_qppmqbmlplbm_eval<3,4,5,0,1,2>;//qbm lp lbm qp p m
case 4525:	 return &A2q2g2l_qppmqbmlmlbp_eval<3,4,5,0,1,2>;//qbm lm lbp qp p m
//case 607:	 return &A2q2g2l_qpmmqbmlplbm_eval<3,4,5,0,1,2>;//qbm lp lbm qp m m
case 637:	 return &A2q2g2l_qpmmqbmlmlbp_eval<3,4,5,0,1,2>;//qbm lm lbp qp m m
//case 7877:	 return &A2q2g2l_qpmmqbmlplbm_eval<2,3,4,5,0,1>;//lp lbm qp m m qbm
case 7882:	 return &A2q2g2l_qpmmqbmlmlbp_eval<2,3,4,5,0,1>;//lm lbp qp m m qbm
//case 8525:	 return &A2q2g2l_qppmqbmlplbm_eval<2,3,4,5,0,1>;//lp lbm qp p m qbm
case 8530:	 return &A2q2g2l_qppmqbmlmlbp_eval<2,3,4,5,0,1>;//lm lbp qp p m qbm
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q2g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}


template complex<R> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2q2g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}

