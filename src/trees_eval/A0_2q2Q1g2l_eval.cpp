/*
 * A0_2q2Q1g2l_eval.cpp
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
 * The 2q2Q quarks 1g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_67_1=spa34*(spa16*spb64+spa17*spb74)+spa35*(spa16*spb65+spa17*spb75);
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));
const complex<T> pow2_spab_3_45_2=(pow(spab_3_45_2,2));

return(
(-((pow2_spa13*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))
-
(pow2_spa16*pow2_spab_3_45_2*spa35)/(s167*spa34*spa45*spa67*spaa_3_45_67_1*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_3_12_67_5*spaa_3_45_67_1*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2ppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa14*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmpq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5*spab_3_24_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5*spab_2_34_5)+(pow2_spa16*pow2_spb35)/(s345*spa12*spa67*spab_2_34_5*spb34)+(pow2_spa14*pow2_spb57*spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mqbpplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_2_45_7=-(spa24*spb74)-spa25*spb75;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_2_13_67_5=-(spa12*(spa56*spb61+spa57*spb71))+spa23*(spa56*spb63+spa57*spb73);
const complex<T> spaa_1_67_45_2=spa16*(-(spa24*spb64)-spa25*spb65)+spa17*(-(spa24*spb74)-spa25*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s245=(spa24*spb42+spa25*spb52+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_2_45_7=(pow(spab_2_45_7,2));

return(
(-((pow2_spa16*pow(s245,2))/(s167*spa45*spa67*spaa_1_67_45_2*spab_5_24_3))+(pow2_spa13*pow2_spab_6_57_4*spab_2_13_4)/(s123*spa23*spa67*spaa_2_13_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb24*spab_1_24_3)/(s234*spa15*spa67*spab_1_23_4*spab_5_24_3*spb23)+(pow2_spa13*pow2_spab_2_45_7)/(spa23*spa45*spaa_1_67_45_2*spaa_2_13_67_5*spb67)
-
(pow2_spa13*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25*spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spab_1_34_2*pow2_spb57*spab_1_24_3)/(s234*s567*spab_1_23_4*spab_1_67_5*spb23*spb34*spb67)
-
(pow2_spa13*pow2_spb57)/(s123*spa23*spab_1_23_4*spb45*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2pmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_34_2*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmmq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spbb_3_12_67_5=spb31*(spa16*spb65+spa17*spb75)+spb32*(spa26*spb65+spa27*spb75);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spb57*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_4_12_3*pow2_spb57*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_3_12_67_5)
-
(pow2_spab_2_16_7*pow2_spb35*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spab_6_12_3*pow2_spb35*spb13)/(spa67*spb12*spb23*spb34*spbb_3_12_67_5*spbb_3_45_67_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mqbpmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spbb_4_67_15_3=-((spa16*spb31-spa56*spb53)*spb64)-(spa17*spb31-spa57*spb53)*spb74;
const complex<T> spbb_3_24_67_5=spb32*(spa26*spb65+spa27*spb75)-spb43*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_3=spa16*spb31-spa56*spb53;
const complex<T> spab_2_13_5=spa12*spb51-spa23*spb53;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s135=(spa13*spb31+spa15*spb51+spa35*spb53);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_3=(pow(spab_6_15_3,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb24)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa13*pow2_spb47*spab_2_13_4)/(s123*spa23*spab_1_23_4*spab_2_13_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb24*spab_1_24_3)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_3_24_67_5)+(pow2_spb47*pow(s135,2))/(s467*spab_2_13_5*spb15*spb67*spbb_4_67_15_3)+(pow2_spab_6_15_3*pow2_spb24)/(spa67*spb15*spb23*spbb_3_24_67_5*spbb_4_67_15_3))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
(-((pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))+(pow2_spa16*pow(s345,2)*spa35)/(s167*spa34*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow2_spb34*spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mpq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_2_34_5=(pow(spab_2_34_5,2));

return(
((pow2_spa16*pow2_spab_2_34_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa12*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmpqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));

return(
((pow2_spa16*pow2_spb45*pow3_spa13)/(spa12*spa23*spa67*spaa_1_67_45_3*spab_1_23_4*spab_1_67_5)
-
(pow2_spa16*pow3_spab_3_45_2)/(s167*s345*spa34*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb34)+(pow2_spb57*pow3_spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pqbpplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
((pow2_spa16*pow(s345,2))/(s167*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*spa23*spa67*spaa_3_12_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb34*spab_1_34_2)/(s234*spa15*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7)/(spa23*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67)
-
(pow2_spa12*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));

return(
(-((pow2_spa16*pow3_spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spb57*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb67)
-
(pow2_spb57*pow3_spab_4_12_3)/(s123*s567*spab_4_23_1*spb23*spb67*spbb_5_67_12_3)+(pow2_spa12*pow2_spb57*pow3_spb35)/(spab_1_67_5*spab_2_34_5*spb34*spb45*spb67*spbb_5_67_12_3))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_23_4=(pow(spab_1_23_4,2));

return(
(-((pow2_spa16*pow2_spb45)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_23_4*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmmqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa23=(pow(spa23,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spa23*pow2_spb57*spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_2_16_7*pow2_spb45*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spb57*pow(s123,2)*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_5_67_12_3)
-
(pow2_spab_6_12_3*pow2_spb45*spb13)/(spa67*spb12*spb23*spb34*spbb_3_45_67_1*spbb_5_67_12_3))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pqbpmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_4_67_15_2=-((spa16*spb21-spa56*spb52)*spb64)-(spa17*spb21-spa57*spb52)*spb74;
const complex<T> spbb_2_34_67_5=-(spb32*(spa36*spb65+spa37*spb75))-spb42*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_2=spa16*spb21-spa56*spb52;
const complex<T> spab_3_12_5=spa13*spb51+spa23*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s125=(spa12*spb21+spa15*spb51+spa25*spb52);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_2=(pow(spab_6_15_2,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb34)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa12*pow2_spb47*spab_3_12_4)/(s123*spa23*spab_1_23_4*spab_3_12_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb34*spab_1_34_2)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_2_34_67_5)
-
(pow2_spb47*pow(s125,2))/(s467*spab_3_12_5*spb15*spb67*spbb_4_67_15_2)+(pow2_spab_6_15_2*pow2_spb34)/(spa67*spb15*spb23*spbb_2_34_67_5*spbb_4_67_15_2))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));

return(
((pow2_spa16*pow2_spb45*pow3_spa13)/(spa12*spa23*spa67*spaa_1_67_45_3*spab_1_23_4*spab_1_67_5)
-
(pow2_spa16*pow3_spab_3_45_2)/(s167*s345*spa34*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb34)+(pow2_spb57*pow3_spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2ppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_2_34_5=(pow(spab_2_34_5,2));

return(
((pow2_spa16*pow2_spab_2_34_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa12*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qppq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
(-((pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))+(pow2_spa16*pow(s345,2)*spa35)/(s167*spa34*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow2_spb34*spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mqbmplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
((pow2_spa16*pow(s345,2))/(s167*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*spa23*spa67*spaa_3_12_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb34*spab_1_34_2)/(s234*spa15*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7)/(spa23*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67)
-
(pow2_spa12*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa23=(pow(spa23,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spa23*pow2_spb57*spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_2_16_7*pow2_spb45*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spb57*pow(s123,2)*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_5_67_12_3)
-
(pow2_spab_6_12_3*pow2_spb45*spb13)/(spa67*spb12*spb23*spb34*spbb_3_45_67_1*spbb_5_67_12_3))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2pmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_23_4=(pow(spab_1_23_4,2));

return(
(-((pow2_spa16*pow2_spb45)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_23_4*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpmq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));

return(
(-((pow2_spa16*pow3_spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spb57*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb67)
-
(pow2_spb57*pow3_spab_4_12_3)/(s123*s567*spab_4_23_1*spb23*spb67*spbb_5_67_12_3)+(pow2_spa12*pow2_spb57*pow3_spb35)/(spab_1_67_5*spab_2_34_5*spb34*spb45*spb67*spbb_5_67_12_3))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mqbmmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_4_67_15_2=-((spa16*spb21-spa56*spb52)*spb64)-(spa17*spb21-spa57*spb52)*spb74;
const complex<T> spbb_2_34_67_5=-(spb32*(spa36*spb65+spa37*spb75))-spb42*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_2=spa16*spb21-spa56*spb52;
const complex<T> spab_3_12_5=spa13*spb51+spa23*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s125=(spa12*spb21+spa15*spb51+spa25*spb52);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_2=(pow(spab_6_15_2,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb34)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa12*pow2_spb47*spab_3_12_4)/(s123*spa23*spab_1_23_4*spab_3_12_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb34*spab_1_34_2)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_2_34_67_5)
-
(pow2_spb47*pow(s125,2))/(s467*spab_3_12_5*spb15*spb67*spbb_4_67_15_2)+(pow2_spab_6_15_2*pow2_spb34)/(spa67*spb15*spb23*spbb_2_34_67_5*spbb_4_67_15_2))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5*spab_3_24_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5*spab_2_34_5)+(pow2_spa16*pow2_spb35)/(s345*spa12*spa67*spab_2_34_5*spb34)+(pow2_spa14*pow2_spb57*spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mpq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa14*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qppqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_67_1=spa34*(spa16*spb64+spa17*spb74)+spa35*(spa16*spb65+spa17*spb75);
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));
const complex<T> pow2_spab_3_45_2=(pow(spab_3_45_2,2));

return(
(-((pow2_spa13*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))
-
(pow2_spa16*pow2_spab_3_45_2*spa35)/(s167*spa34*spa45*spa67*spaa_3_45_67_1*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_3_12_67_5*spaa_3_45_67_1*spb67))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pqbmplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_2_45_7=-(spa24*spb74)-spa25*spb75;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_2_13_67_5=-(spa12*(spa56*spb61+spa57*spb71))+spa23*(spa56*spb63+spa57*spb73);
const complex<T> spaa_1_67_45_2=spa16*(-(spa24*spb64)-spa25*spb65)+spa17*(-(spa24*spb74)-spa25*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s245=(spa24*spb42+spa25*spb52+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_2_45_7=(pow(spab_2_45_7,2));

return(
(-((pow2_spa16*pow(s245,2))/(s167*spa45*spa67*spaa_1_67_45_2*spab_5_24_3))+(pow2_spa13*pow2_spab_6_57_4*spab_2_13_4)/(s123*spa23*spa67*spaa_2_13_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb24*spab_1_24_3)/(s234*spa15*spa67*spab_1_23_4*spab_5_24_3*spb23)+(pow2_spa13*pow2_spab_2_45_7)/(spa23*spa45*spaa_1_67_45_2*spaa_2_13_67_5*spb67)
-
(pow2_spa13*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spbb_3_12_67_5=spb31*(spa16*spb65+spa17*spb75)+spb32*(spa26*spb65+spa27*spb75);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spb57*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_4_12_3*pow2_spb57*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_3_12_67_5)
-
(pow2_spab_2_16_7*pow2_spb35*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spab_6_12_3*pow2_spb35*spb13)/(spa67*spb12*spb23*spb34*spbb_3_12_67_5*spbb_3_45_67_1))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_34_2*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpmqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25*spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spab_1_34_2*pow2_spb57*spab_1_24_3)/(s234*s567*spab_1_23_4*spab_1_67_5*spb23*spb34*spb67)
-
(pow2_spa13*pow2_spb57)/(s123*spa23*spab_1_23_4*spb45*spb67))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i6, int i7, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pqbmmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spbb_4_67_15_3=-((spa16*spb31-spa56*spb53)*spb64)-(spa17*spb31-spa57*spb53)*spb74;
const complex<T> spbb_3_24_67_5=spb32*(spa26*spb65+spa27*spb75)-spb43*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_3=spa16*spb31-spa56*spb53;
const complex<T> spab_2_13_5=spa12*spb51-spa23*spb53;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s135=(spa13*spb31+spa15*spb51+spa35*spb53);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_3=(pow(spab_6_15_3,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb24)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa13*pow2_spb47*spab_2_13_4)/(s123*spa23*spab_1_23_4*spab_2_13_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb24*spab_1_24_3)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_3_24_67_5)+(pow2_spb47*pow(s135,2))/(s467*spab_2_13_5*spb15*spb67*spbb_4_67_15_3)+(pow2_spab_6_15_3*pow2_spb24)/(spa67*spb15*spb23*spbb_3_24_67_5*spbb_4_67_15_3))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_67_1=spa34*(spa16*spb64+spa17*spb74)+spa35*(spa16*spb65+spa17*spb75);
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));
const complex<T> pow2_spab_3_45_2=(pow(spab_3_45_2,2));

return(
(-((pow2_spa13*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))
-
(pow2_spa16*pow2_spab_3_45_2*spa35)/(s167*spa34*spa45*spa67*spaa_3_45_67_1*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_3_12_67_5*spaa_3_45_67_1*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2ppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa14*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmpq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5*spab_3_24_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5*spab_2_34_5)+(pow2_spa16*pow2_spb35)/(s345*spa12*spa67*spab_2_34_5*spb34)+(pow2_spa14*pow2_spb57*spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mqbpplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_2_45_7=-(spa24*spb74)-spa25*spb75;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_2_13_67_5=-(spa12*(spa56*spb61+spa57*spb71))+spa23*(spa56*spb63+spa57*spb73);
const complex<T> spaa_1_67_45_2=spa16*(-(spa24*spb64)-spa25*spb65)+spa17*(-(spa24*spb74)-spa25*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s245=(spa24*spb42+spa25*spb52+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_2_45_7=(pow(spab_2_45_7,2));

return(
(-((pow2_spa16*pow(s245,2))/(s167*spa45*spa67*spaa_1_67_45_2*spab_5_24_3))+(pow2_spa13*pow2_spab_6_57_4*spab_2_13_4)/(s123*spa23*spa67*spaa_2_13_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb24*spab_1_24_3)/(s234*spa15*spa67*spab_1_23_4*spab_5_24_3*spb23)+(pow2_spa13*pow2_spab_2_45_7)/(spa23*spa45*spaa_1_67_45_2*spaa_2_13_67_5*spb67)
-
(pow2_spa13*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25*spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spab_1_34_2*pow2_spb57*spab_1_24_3)/(s234*s567*spab_1_23_4*spab_1_67_5*spb23*spb34*spb67)
-
(pow2_spa13*pow2_spb57)/(s123*spa23*spab_1_23_4*spb45*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2pmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_34_2*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmmq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spbb_3_12_67_5=spb31*(spa16*spb65+spa17*spb75)+spb32*(spa26*spb65+spa27*spb75);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spb57*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_4_12_3*pow2_spb57*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_3_12_67_5)
-
(pow2_spab_2_16_7*pow2_spb35*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spab_6_12_3*pow2_spb35*spb13)/(spa67*spb12*spb23*spb34*spbb_3_12_67_5*spbb_3_45_67_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmq2pqb2mqbpmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spbb_4_67_15_3=-((spa16*spb31-spa56*spb53)*spb64)-(spa17*spb31-spa57*spb53)*spb74;
const complex<T> spbb_3_24_67_5=spb32*(spa26*spb65+spa27*spb75)-spb43*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_3=spa16*spb31-spa56*spb53;
const complex<T> spab_2_13_5=spa12*spb51-spa23*spb53;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s135=(spa13*spb31+spa15*spb51+spa35*spb53);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_3=(pow(spab_6_15_3,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb24)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa13*pow2_spb47*spab_2_13_4)/(s123*spa23*spab_1_23_4*spab_2_13_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb24*spab_1_24_3)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_3_24_67_5)+(pow2_spb47*pow(s135,2))/(s467*spab_2_13_5*spb15*spb67*spbb_4_67_15_3)+(pow2_spab_6_15_3*pow2_spb24)/(spa67*spb15*spb23*spbb_3_24_67_5*spbb_4_67_15_3))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
(-((pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))+(pow2_spa16*pow(s345,2)*spa35)/(s167*spa34*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow2_spb34*spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mpq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_2_34_5=(pow(spab_2_34_5,2));

return(
((pow2_spa16*pow2_spab_2_34_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa12*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmpqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));

return(
((pow2_spa16*pow2_spb45*pow3_spa13)/(spa12*spa23*spa67*spaa_1_67_45_3*spab_1_23_4*spab_1_67_5)
-
(pow2_spa16*pow3_spab_3_45_2)/(s167*s345*spa34*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb34)+(pow2_spb57*pow3_spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pqbpplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
((pow2_spa16*pow(s345,2))/(s167*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*spa23*spa67*spaa_3_12_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb34*spab_1_34_2)/(s234*spa15*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7)/(spa23*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67)
-
(pow2_spa12*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));

return(
(-((pow2_spa16*pow3_spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spb57*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb67)
-
(pow2_spb57*pow3_spab_4_12_3)/(s123*s567*spab_4_23_1*spb23*spb67*spbb_5_67_12_3)+(pow2_spa12*pow2_spb57*pow3_spb35)/(spab_1_67_5*spab_2_34_5*spb34*spb45*spb67*spbb_5_67_12_3))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_23_4=(pow(spab_1_23_4,2));

return(
(-((pow2_spa16*pow2_spb45)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_23_4*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmmqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa23=(pow(spa23,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spa23*pow2_spb57*spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_2_16_7*pow2_spb45*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spb57*pow(s123,2)*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_5_67_12_3)
-
(pow2_spab_6_12_3*pow2_spb45*spb13)/(spa67*spb12*spb23*spb34*spbb_3_45_67_1*spbb_5_67_12_3))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qmqb2mq2pqbpmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_4_67_15_2=-((spa16*spb21-spa56*spb52)*spb64)-(spa17*spb21-spa57*spb52)*spb74;
const complex<T> spbb_2_34_67_5=-(spb32*(spa36*spb65+spa37*spb75))-spb42*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_2=spa16*spb21-spa56*spb52;
const complex<T> spab_3_12_5=spa13*spb51+spa23*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s125=(spa12*spb21+spa15*spb51+spa25*spb52);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_2=(pow(spab_6_15_2,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb34)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa12*pow2_spb47*spab_3_12_4)/(s123*spa23*spab_1_23_4*spab_3_12_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb34*spab_1_34_2)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_2_34_67_5)
-
(pow2_spb47*pow(s125,2))/(s467*spab_3_12_5*spb15*spb67*spbb_4_67_15_2)+(pow2_spab_6_15_2*pow2_spb34)/(spa67*spb15*spb23*spbb_2_34_67_5*spbb_4_67_15_2))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));

return(
((pow2_spa16*pow2_spb45*pow3_spa13)/(spa12*spa23*spa67*spaa_1_67_45_3*spab_1_23_4*spab_1_67_5)
-
(pow2_spa16*pow3_spab_3_45_2)/(s167*s345*spa34*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb34)+(pow2_spb57*pow3_spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2ppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_2_34_5=(pow(spab_2_34_5,2));

return(
((pow2_spa16*pow2_spab_2_34_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa12*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qppq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
(-((pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))+(pow2_spa16*pow(s345,2)*spa35)/(s167*spa34*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa16*pow2_spb34*spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mqbmplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> spaa_1_67_45_3=spa16*(-(spa34*spb64)-spa35*spb65)+spa17*(-(spa34*spb74)-spa35*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));

return(
((pow2_spa16*pow(s345,2))/(s167*spa45*spa67*spaa_1_67_45_3*spab_5_34_2)+(pow2_spa12*pow2_spab_6_57_4*spab_3_12_4)/(s123*spa23*spa67*spaa_3_12_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb34*spab_1_34_2)/(s234*spa15*spa67*spab_1_23_4*spab_5_34_2*spb23)+(pow2_spa12*pow2_spab_3_45_7)/(spa23*spa45*spaa_1_67_45_3*spaa_3_12_67_5*spb67)
-
(pow2_spa12*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa23=(pow(spa23,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spa23*pow2_spb57*spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_2_16_7*pow2_spb45*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spb57*pow(s123,2)*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_5_67_12_3)
-
(pow2_spab_6_12_3*pow2_spb45*spb13)/(spa67*spb12*spb23*spb34*spbb_3_45_67_1*spbb_5_67_12_3))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2pmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb45=(pow(spb45,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_23_4=(pow(spab_1_23_4,2));

return(
(-((pow2_spa16*pow2_spb45)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_23_4*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpmq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_5_67_12_3=-((spa16*spb31+spa26*spb32)*spb65)-(spa17*spb31+spa27*spb32)*spb75;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));

return(
(-((pow2_spa16*pow3_spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spb57*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb67)
-
(pow2_spb57*pow3_spab_4_12_3)/(s123*s567*spab_4_23_1*spb23*spb67*spbb_5_67_12_3)+(pow2_spa12*pow2_spb57*pow3_spb35)/(spab_1_67_5*spab_2_34_5*spb34*spb45*spb67*spbb_5_67_12_3))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpq2pqb2mqbmmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb34=(pow(spb34,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spbb_4_67_15_2=-((spa16*spb21-spa56*spb52)*spb64)-(spa17*spb21-spa57*spb52)*spb74;
const complex<T> spbb_2_34_67_5=-(spb32*(spa36*spb65+spa37*spb75))-spb42*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_2=spa16*spb21-spa56*spb52;
const complex<T> spab_3_12_5=spa13*spb51+spa23*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s125=(spa12*spb21+spa15*spb51+spa25*spb52);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_2=(pow(spab_6_15_2,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb34)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa12*pow2_spb47*spab_3_12_4)/(s123*spa23*spab_1_23_4*spab_3_12_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb34*spab_1_34_2)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_2_34_67_5)
-
(pow2_spb47*pow(s125,2))/(s467*spab_3_12_5*spb15*spb67*spbb_4_67_15_2)+(pow2_spab_6_15_2*pow2_spb34)/(spa67*spb15*spb23*spbb_2_34_67_5*spbb_4_67_15_2))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5*spab_3_24_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5*spab_2_34_5)+(pow2_spa16*pow2_spb35)/(s345*spa12*spa67*spab_2_34_5*spb34)+(pow2_spa14*pow2_spb57*spa13)/(s567*spa12*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mpq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_4_23_5=(pow(spab_4_23_5,2));

return(
((pow2_spa16*pow2_spab_4_23_5)/(s167*s234*spa23*spa34*spa67*spab_1_67_5)+(pow2_spa14*pow2_spb57)/(s567*spa23*spa34*spab_1_67_5*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qppqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_67_1=spa34*(spa16*spb64+spa17*spb74)+spa35*(spa16*spb65+spa17*spb75);
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_3_45_7=(pow(spab_3_45_7,2));
const complex<T> pow2_spab_3_45_2=(pow(spab_3_45_2,2));

return(
(-((pow2_spa13*pow2_spab_6_57_4*spab_3_12_4)/(s123*s567*spa23*spa67*spaa_3_12_67_5*spab_1_23_4))
-
(pow2_spa16*pow2_spab_3_45_2*spa35)/(s167*spa34*spa45*spa67*spaa_3_45_67_1*spab_5_34_2)+(pow2_spa16*pow3_spb24)/(s234*spa67*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_45_7*spa35)/(spa23*spa34*spa45*spaa_3_12_67_5*spaa_3_45_67_1*spb67))*complex<T>(0,-1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pqbmplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_2_45_7=-(spa24*spb74)-spa25*spb75;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_2_13_67_5=-(spa12*(spa56*spb61+spa57*spb71))+spa23*(spa56*spb63+spa57*spb73);
const complex<T> spaa_1_67_45_2=spa16*(-(spa24*spb64)-spa25*spb65)+spa17*(-(spa24*spb74)-spa25*spb75);
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s245=(spa24*spb42+spa25*spb52+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_57_4=(pow(spab_6_57_4,2));
const complex<T> pow2_spab_2_45_7=(pow(spab_2_45_7,2));

return(
(-((pow2_spa16*pow(s245,2))/(s167*spa45*spa67*spaa_1_67_45_2*spab_5_24_3))+(pow2_spa13*pow2_spab_6_57_4*spab_2_13_4)/(s123*spa23*spa67*spaa_2_13_67_5*spab_1_23_4*spab_5_67_4)
-
(pow2_spa16*pow2_spb24*spab_1_24_3)/(s234*spa15*spa67*spab_1_23_4*spab_5_24_3*spb23)+(pow2_spa13*pow2_spab_2_45_7)/(spa23*spa45*spaa_1_67_45_2*spaa_2_13_67_5*spb67)
-
(pow2_spa13*pow2_spb47)/(s467*spa15*spa23*spab_5_67_4*spb67))*complex<T>(0,-1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb35=(pow(spb35,2));
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spbb_3_12_67_5=spb31*(spa16*spb65+spa17*spb75)+spb32*(spa26*spb65+spa27*spb75);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_6_12_3=(pow(spab_6_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_16_7=(pow(spab_2_16_7,2));

return(
(-((pow2_spb57*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb67))
-
(pow2_spab_4_12_3*pow2_spb57*spb13)/(s567*spab_4_23_1*spb12*spb23*spb67*spbb_3_12_67_5)
-
(pow2_spab_2_16_7*pow2_spb35*spab_2_45_3)/(s167*s345*spab_2_34_5*spb34*spb67*spbb_3_45_67_1)+(pow2_spab_6_12_3*pow2_spb35*spb13)/(spa67*spb12*spb23*spb34*spbb_3_12_67_5*spbb_3_45_67_1))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25)/(s167*spa67*spab_1_67_5*spb23*spb34))
-
(pow2_spab_1_34_2*pow2_spb57)/(s234*s567*spab_1_67_5*spb23*spb34*spb67))*complex<T>(0,1)
);
}


template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpmqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spb25=(pow(spb25,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa16*pow2_spb25*spb35)/(s167*spa67*spab_1_67_5*spb23*spb34*spb45))
-
(pow2_spab_1_34_2*pow2_spb57*spab_1_24_3)/(s234*s567*spab_1_23_4*spab_1_67_5*spb23*spb34*spb67)
-
(pow2_spa13*pow2_spb57)/(s123*spa23*spab_1_23_4*spb45*spb67))*complex<T>(0,1)
);
}


template <int i4, int i3, int i2, int i1, int i5, int i7, int i6, class T> complex<T>  A2q2Q1g2l_qpqb2mq2pqbmmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> pow2_spb24=(pow(spb24,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spbb_4_67_15_3=-((spa16*spb31-spa56*spb53)*spb64)-(spa17*spb31-spa57*spb53)*spb74;
const complex<T> spbb_3_24_67_5=spb32*(spa26*spb65+spa27*spb75)-spb43*(spa46*spb65+spa47*spb75);
const complex<T> spab_6_15_3=spa16*spb31-spa56*spb53;
const complex<T> spab_2_13_5=spa12*spb51-spa23*spb53;
const complex<T> spab_2_13_4=spa12*spb41-spa23*spb43;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_56_7=-(spa15*spb75)-spa16*spb76;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s467=(spa46*spb64+spa47*spb74+spa67*spb76);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s167=(spa16*spb61+spa17*spb71+spa67*spb76);
const complex<T> s135=(spa13*spb31+spa15*spb51+spa35*spb53);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> pow2_spab_6_15_3=(pow(spab_6_15_3,2));
const complex<T> pow2_spab_1_56_7=(pow(spab_1_56_7,2));

return(
(-((pow2_spa16*pow2_spb24)/(s167*spa67*spab_1_67_5*spb23*spb45))
-
(pow2_spa13*pow2_spb47*spab_2_13_4)/(s123*spa23*spab_1_23_4*spab_2_13_5*spb45*spb67)+(pow2_spab_1_56_7*pow2_spb24*spab_1_24_3)/(s234*spab_1_23_4*spab_1_67_5*spb23*spb67*spbb_3_24_67_5)+(pow2_spb47*pow(s135,2))/(s467*spab_2_13_5*spb15*spb67*spbb_4_67_15_3)+(pow2_spab_6_15_3*pow2_spb24)/(spa67*spb15*spb23*spbb_3_24_67_5*spbb_4_67_15_3))*complex<T>(0,1)
);
}




template <class T> complex<T> (*A2q2Q1g2l_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 1803489:	 return &A2q2Q1g2l_qpmq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m q2p qb2m qbm lp lbm
case 1803497:	 return &A2q2Q1g2l_qppq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p q2p qb2m qbm lp lbm
case 1803545:	 return &A2q2Q1g2l_qpq2pmqb2mqbmlplbm_eval<0,1,2,3,4,5,6>;//qp q2p m qb2m qbm lp lbm
case 1803609:	 return &A2q2Q1g2l_qpq2ppqb2mqbmlplbm_eval<0,1,2,3,4,5,6>;//qp q2p p qb2m qbm lp lbm
case 1803937:	 return &A2q2Q1g2l_qpmqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6>;//qp m qb2m q2p qbm lp lbm
case 1803945:	 return &A2q2Q1g2l_qppqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6>;//qp p qb2m q2p qbm lp lbm
case 1804049:	 return &A2q2Q1g2l_qpqb2mmq2pqbmlplbm_eval<0,1,2,3,4,5,6>;//qp qb2m m q2p qbm lp lbm
case 1804113:	 return &A2q2Q1g2l_qpqb2mpq2pqbmlplbm_eval<0,1,2,3,4,5,6>;//qp qb2m p q2p qbm lp lbm
case 1804441:	 return &A2q2Q1g2l_qpq2pqb2mmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp q2p qb2m m qbm lp lbm
case 1804497:	 return &A2q2Q1g2l_qpqb2mq2pmqbmlplbm_eval<0,1,2,3,4,5,6>;//qp qb2m q2p m qbm lp lbm
case 1804953:	 return &A2q2Q1g2l_qpq2pqb2mpqbmlplbm_eval<0,1,2,3,4,5,6>;//qp q2p qb2m p qbm lp lbm
case 1805009:	 return &A2q2Q1g2l_qpqb2mq2ppqbmlplbm_eval<0,1,2,3,4,5,6>;//qp qb2m q2p p qbm lp lbm
case 1807584:	 return &A2q2Q1g2l_qmmq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6>;//qm m q2p qb2m qbp lp lbm
case 1807592:	 return &A2q2Q1g2l_qmpq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6>;//qm p q2p qb2m qbp lp lbm
case 1807640:	 return &A2q2Q1g2l_qmq2pmqb2mqbplplbm_eval<0,1,2,3,4,5,6>;//qm q2p m qb2m qbp lp lbm
case 1807704:	 return &A2q2Q1g2l_qmq2ppqb2mqbplplbm_eval<0,1,2,3,4,5,6>;//qm q2p p qb2m qbp lp lbm
case 1808032:	 return &A2q2Q1g2l_qmmqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6>;//qm m qb2m q2p qbp lp lbm
case 1808040:	 return &A2q2Q1g2l_qmpqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6>;//qm p qb2m q2p qbp lp lbm
case 1808144:	 return &A2q2Q1g2l_qmqb2mmq2pqbplplbm_eval<0,1,2,3,4,5,6>;//qm qb2m m q2p qbp lp lbm
case 1808208:	 return &A2q2Q1g2l_qmqb2mpq2pqbplplbm_eval<0,1,2,3,4,5,6>;//qm qb2m p q2p qbp lp lbm
case 1808536:	 return &A2q2Q1g2l_qmq2pqb2mmqbplplbm_eval<0,1,2,3,4,5,6>;//qm q2p qb2m m qbp lp lbm
case 1808592:	 return &A2q2Q1g2l_qmqb2mq2pmqbplplbm_eval<0,1,2,3,4,5,6>;//qm qb2m q2p m qbp lp lbm
case 1809048:	 return &A2q2Q1g2l_qmq2pqb2mpqbplplbm_eval<0,1,2,3,4,5,6>;//qm q2p qb2m p qbp lp lbm
case 1809104:	 return &A2q2Q1g2l_qmqb2mq2ppqbplplbm_eval<0,1,2,3,4,5,6>;//qm qb2m q2p p qbp lp lbm
case 1818777:	 return &A2q2Q1g2l_qpq2pqb2mqbmmlplbm_eval<0,1,2,3,4,5,6>;//qp q2p qb2m qbm m lp lbm
case 1818833:	 return &A2q2Q1g2l_qpqb2mq2pqbmmlplbm_eval<0,1,2,3,4,5,6>;//qp qb2m q2p qbm m lp lbm
case 1819288:	 return &A2q2Q1g2l_qmq2pqb2mqbpmlplbm_eval<0,1,2,3,4,5,6>;//qm q2p qb2m qbp m lp lbm
case 1819344:	 return &A2q2Q1g2l_qmqb2mq2pqbpmlplbm_eval<0,1,2,3,4,5,6>;//qm qb2m q2p qbp m lp lbm
case 1822873:	 return &A2q2Q1g2l_qpq2pqb2mqbmplplbm_eval<0,1,2,3,4,5,6>;//qp q2p qb2m qbm p lp lbm
case 1822929:	 return &A2q2Q1g2l_qpqb2mq2pqbmplplbm_eval<0,1,2,3,4,5,6>;//qp qb2m q2p qbm p lp lbm
case 1823384:	 return &A2q2Q1g2l_qmq2pqb2mqbpplplbm_eval<0,1,2,3,4,5,6>;//qm q2p qb2m qbp p lp lbm
case 1823440:	 return &A2q2Q1g2l_qmqb2mq2pqbpplplbm_eval<0,1,2,3,4,5,6>;//qm qb2m q2p qbp p lp lbm
case 2032865:	 return &A2q2Q1g2l_qpmq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m q2p qb2m qbm lm lbp
case 2032873:	 return &A2q2Q1g2l_qppq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p q2p qb2m qbm lm lbp
case 2032921:	 return &A2q2Q1g2l_qpq2pmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp q2p m qb2m qbm lm lbp
case 2032985:	 return &A2q2Q1g2l_qpq2ppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp q2p p qb2m qbm lm lbp
case 2033313:	 return &A2q2Q1g2l_qpmqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp m qb2m q2p qbm lm lbp
case 2033321:	 return &A2q2Q1g2l_qppqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp p qb2m q2p qbm lm lbp
case 2033425:	 return &A2q2Q1g2l_qpqb2mmq2pqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m m q2p qbm lm lbp
case 2033489:	 return &A2q2Q1g2l_qpqb2mpq2pqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m p q2p qbm lm lbp
case 2033817:	 return &A2q2Q1g2l_qpq2pqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp q2p qb2m m qbm lm lbp
case 2033873:	 return &A2q2Q1g2l_qpqb2mq2pmqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m q2p m qbm lm lbp
case 2034329:	 return &A2q2Q1g2l_qpq2pqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp q2p qb2m p qbm lm lbp
case 2034385:	 return &A2q2Q1g2l_qpqb2mq2ppqbmlmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m q2p p qbm lm lbp
case 2036960:	 return &A2q2Q1g2l_qmmq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m q2p qb2m qbp lm lbp
case 2036968:	 return &A2q2Q1g2l_qmpq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p q2p qb2m qbp lm lbp
case 2037016:	 return &A2q2Q1g2l_qmq2pmqb2mqbplmlbp_eval<0,1,2,3,4,5,6>;//qm q2p m qb2m qbp lm lbp
case 2037080:	 return &A2q2Q1g2l_qmq2ppqb2mqbplmlbp_eval<0,1,2,3,4,5,6>;//qm q2p p qb2m qbp lm lbp
case 2037408:	 return &A2q2Q1g2l_qmmqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6>;//qm m qb2m q2p qbp lm lbp
case 2037416:	 return &A2q2Q1g2l_qmpqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6>;//qm p qb2m q2p qbp lm lbp
case 2037520:	 return &A2q2Q1g2l_qmqb2mmq2pqbplmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m m q2p qbp lm lbp
case 2037584:	 return &A2q2Q1g2l_qmqb2mpq2pqbplmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m p q2p qbp lm lbp
case 2037912:	 return &A2q2Q1g2l_qmq2pqb2mmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm q2p qb2m m qbp lm lbp
case 2037968:	 return &A2q2Q1g2l_qmqb2mq2pmqbplmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m q2p m qbp lm lbp
case 2038424:	 return &A2q2Q1g2l_qmq2pqb2mpqbplmlbp_eval<0,1,2,3,4,5,6>;//qm q2p qb2m p qbp lm lbp
case 2038480:	 return &A2q2Q1g2l_qmqb2mq2ppqbplmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m q2p p qbp lm lbp
case 2048153:	 return &A2q2Q1g2l_qpq2pqb2mqbmmlmlbp_eval<0,1,2,3,4,5,6>;//qp q2p qb2m qbm m lm lbp
case 2048209:	 return &A2q2Q1g2l_qpqb2mq2pqbmmlmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m q2p qbm m lm lbp
case 2048664:	 return &A2q2Q1g2l_qmq2pqb2mqbpmlmlbp_eval<0,1,2,3,4,5,6>;//qm q2p qb2m qbp m lm lbp
case 2048720:	 return &A2q2Q1g2l_qmqb2mq2pqbpmlmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m q2p qbp m lm lbp
case 2052249:	 return &A2q2Q1g2l_qpq2pqb2mqbmplmlbp_eval<0,1,2,3,4,5,6>;//qp q2p qb2m qbm p lm lbp
case 2052305:	 return &A2q2Q1g2l_qpqb2mq2pqbmplmlbp_eval<0,1,2,3,4,5,6>;//qp qb2m q2p qbm p lm lbp
case 2052760:	 return &A2q2Q1g2l_qmq2pqb2mqbpplmlbp_eval<0,1,2,3,4,5,6>;//qm q2p qb2m qbp p lm lbp
case 2052816:	 return &A2q2Q1g2l_qmqb2mq2pqbpplmlbp_eval<0,1,2,3,4,5,6>;//qm qb2m q2p qbp p lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q2Q1g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q2Q1g2l_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2Q1g2l_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2Q1g2l_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);


#if BH_USE_GMP

template complex<RGMP> (*A2q2Q1g2l_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);

#endif
}



