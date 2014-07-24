/*
 * A0_2q2Q2g2l_eval.cpp
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
 * The 2q2Q quarks 2g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2pqb2mppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spab_3_456_2=(pow(spab_3_456_2,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
((pow2_spa13*pow2_spaa_3_45_68_7*spaa_3_12_45_3)/(s678*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1)
-
(pow2_spa13*pow2_spab_7_568_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)+(pow2_spa17*pow2_spab_3_456_2*spa36)/(s178*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)
-
(pow2_spa17*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_456_8*spa36)/(spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2ppqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5)+(pow2_spa14*pow2_spab_7_68_5*spab_4_123_5)/(s1234*s678*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa14*pow2_spab_4_56_8*spa46)/(spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2pppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa15*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmpq2pqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spab_4_56_3=(pow(spab_4_56_3,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46*spaa_4_56_178_3)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa14*pow2_spab_7_68_5*spa13*spab_4_123_5)/(s1234*s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa17*pow2_spab_4_56_3*spa46)/(s3456*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_56_8*spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmpq2ppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_178_6=spa13*spb61+spa37*spb76+spa38*spb86;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6)+(pow2_spa17*pow2_spab_5_234_6*spab_3_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa15*pow2_spb68*spa13)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmppq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_178_6=spa14*spb61+spa47*spb76+spa48*spb86;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6*spab_4_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa17*pow2_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6*spab_3_45_6)+(pow2_spa17*pow2_spb46)/(s456*spa12*spa23*spa78*spab_3_45_6*spb45)+(pow2_spa15*pow2_spb68*spa14)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2pqb2mmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_3=spa16*spb63+spa17*spb73+spa18*spb83;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb36)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_3)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68*spab_1_24_3)/(s1234*s234*spab_1_23_4*spab_1_678_5*spb23*spb34*spb56*spb78)+(pow2_spa13*pow2_spb68)/(s123*spa23*spab_1_23_4*spb45*spb56*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2pmqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_4=spa16*spb64+spa17*spb74+spa18*spb84;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb46)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_4)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68)/(s1234*s234*spab_1_678_5*spb23*spb34*spb56*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmq2pmmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
((pow2_spa17*pow2_spb26)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_345_2*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmq2pqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spbb_3_12_678_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)+(pow2_spab_4_12_3*pow2_spb68*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5)
-
(pow2_spab_2_17_8*pow2_spb36*spab_2_456_3*spb46)/(s178*s3456*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow2_spb36*spb13*spb46)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13*spbb_3_12_678_4)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmq2pmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
(-((pow2_spab_2_17_8*pow2_spb36*spab_2_456_3)/(s178*s3456*spab_2_178_6*spb34*spb45*spb78*spbb_3_456_78_1))
-
(pow2_spab_7_12_3*pow2_spb36*spb13)/(spa78*spb12*spb23*spb34*spb45*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3)/(s2345*s345*spab_2_178_6*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmmq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
-
(pow2_spab_5_123_4*pow2_spb68*spb14)/(s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)+(pow2_spb68*pow3_spab_5_23_4)/(s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)+(pow2_spab_3_127_8*pow2_spb46*spab_3_56_4)/(s3456*s456*spab_3_45_6*spb12*spb45*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow2_spb46*spb14)/(spa78*spb12*spb23*spb34*spb45*spbb_4_123_78_6*spbb_4_56_78_1)+(pow2_spb46*pow2_spbb_4_23_17_8*spbb_4_56_23_4)/(s178*spb23*spb34*spb45*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mq2pppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_678_1=spa36*spb61+spa37*spb71+spa38*spb81;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-((pow2_spaa_2_34_568_7*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_3_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+pow2_spab_7_12_3/(s123*spa45*spa56*spa78*spab_4_23_1*spb23)+(pow2_spab_2_17_8*spa36)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mpq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_678_1=spa46*spb61+spa47*spb71+spa48*spb81;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-(pow2_spaa_2_34_568_7/(s234*s5678*spa23*spa34*spa56*spa78*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_4_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+(pow2_spab_2_17_8*spa46)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mppq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_345_6=(pow(spab_2_345_6,2));

return(
((pow2_spa17*pow2_spab_2_345_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa12*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmpqb2mq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_678_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13*spaa_3_12_678_4)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1))+(pow2_spab_7_568_4*pow3_spa13)/(s5678*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)
-
(pow2_spa17*pow3_spab_3_456_2*spa46)/(s178*s3456*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_456_8*pow3_spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmpqb2mpq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_78_6*spaa_3_45_678_1))
-
(pow2_spa17*pow3_spab_3_456_2)/(s178*s3456*spa34*spa45*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_6_178_2)
-
(pow2_spab_3_456_8*pow3_spa13)/(spa12*spa23*spa34*spa45*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmppqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
(-((pow2_spa17*pow3_spaa_4_23_56_4)/(s178*spa23*spa34*spa45*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spab_7_68_5*pow3_spa14)/(s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)+(pow2_spa17*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa17*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_56_8*pow3_spa14)/(spa12*spa23*spa34*spa45*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mq2pmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=(pow(spb36,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb78)
-
(pow2_spb68*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb78*spbb_3_12_678_5)+(pow2_spab_2_17_8*pow3_spb36)/(s178*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow3_spb36)/(spa78*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)
-
(pow2_spb68*pow3_spab_2_45_3)/(s2345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow3_spbb_3_45_12_3)/(s678*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mmq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));

return(
((pow2_spb68*pow3_spab_5_123_4)/(s1234*s678*spab_5_678_1*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4)/(s234*s2345*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_7_123_4*pow3_spb46)/(spa78*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmqb2mmmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb56=(pow(spb56,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_234_5=(pow(spab_1_234_5,2));

return(
((pow2_spa17*pow2_spb56)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_234_5*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmqb2mq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb78)+(pow2_spb68*pow3_spab_5_123_4*spb13)/(s1234*s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_3_127_8*pow3_spb46)/(s3456*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow3_spb46*spb13)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46*spbb_4_56_178_3)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmqb2mmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_3=spa16*spb31+spa67*spb73+spa68*spb83;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));

return(
((pow2_spab_7_68_5*spb13)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)
-
pow2_spbb_5_34_127_8/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_3)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qmmmqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_4=spa16*spb41+spa67*spb74+spa68*spb84;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
((pow2_spab_7_68_5*spb14)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)+pow2_spab_4_56_8/(s456*spa45*spab_6_45_3*spb12*spb23*spb78)
-
(pow2_spbb_5_34_127_8*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_4)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2pqb2mppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
(-((pow2_spa17*pow3_spaa_4_23_56_4)/(s178*spa23*spa34*spa45*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spab_7_68_5*pow3_spa14)/(s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)+(pow2_spa17*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa17*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_56_8*pow3_spa14)/(spa12*spa23*spa34*spa45*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2ppqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_78_6*spaa_3_45_678_1))
-
(pow2_spa17*pow3_spab_3_456_2)/(s178*s3456*spa34*spa45*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_6_178_2)
-
(pow2_spab_3_456_8*pow3_spa13)/(spa12*spa23*spa34*spa45*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2pppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_345_6=(pow(spab_2_345_6,2));

return(
((pow2_spa17*pow2_spab_2_345_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa12*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qppq2pqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_678_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13*spaa_3_12_678_4)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1))+(pow2_spab_7_568_4*pow3_spa13)/(s5678*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)
-
(pow2_spa17*pow3_spab_3_456_2*spa46)/(s178*s3456*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_456_8*pow3_spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qppq2ppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_678_1=spa46*spb61+spa47*spb71+spa48*spb81;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-(pow2_spaa_2_34_568_7/(s234*s5678*spa23*spa34*spa56*spa78*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_4_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+(pow2_spab_2_17_8*spa46)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpppq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_678_1=spa36*spb61+spa37*spb71+spa38*spb81;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-((pow2_spaa_2_34_568_7*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_3_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+pow2_spab_7_12_3/(s123*spa45*spa56*spa78*spab_4_23_1*spb23)+(pow2_spab_2_17_8*spa36)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2pqb2mmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_4=spa16*spb41+spa67*spb74+spa68*spb84;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
((pow2_spab_7_68_5*spb14)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)+pow2_spab_4_56_8/(s456*spa45*spab_6_45_3*spb12*spb23*spb78)
-
(pow2_spbb_5_34_127_8*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_4)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2pmqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_3=spa16*spb31+spa67*spb73+spa68*spb83;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));

return(
((pow2_spab_7_68_5*spb13)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)
-
pow2_spbb_5_34_127_8/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_3)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpq2pmmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb56=(pow(spb56,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_234_5=(pow(spab_1_234_5,2));

return(
((pow2_spa17*pow2_spb56)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_234_5*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmq2pqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb78)+(pow2_spb68*pow3_spab_5_123_4*spb13)/(s1234*s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_3_127_8*pow3_spb46)/(s3456*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow3_spb46*spb13)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46*spbb_4_56_178_3)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmq2pmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));

return(
((pow2_spb68*pow3_spab_5_123_4)/(s1234*s678*spab_5_678_1*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4)/(s234*s2345*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_7_123_4*pow3_spb46)/(spa78*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmmq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=(pow(spb36,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb78)
-
(pow2_spb68*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb78*spbb_3_12_678_5)+(pow2_spab_2_17_8*pow3_spb36)/(s178*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow3_spb36)/(spa78*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)
-
(pow2_spb68*pow3_spab_2_45_3)/(s2345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow3_spbb_3_45_12_3)/(s678*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mq2pppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_178_6=spa14*spb61+spa47*spb76+spa48*spb86;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6*spab_4_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa17*pow2_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6*spab_3_45_6)+(pow2_spa17*pow2_spb46)/(s456*spa12*spa23*spa78*spab_3_45_6*spb45)+(pow2_spa15*pow2_spb68*spa14)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mpq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_178_6=spa13*spb61+spa37*spb76+spa38*spb86;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6)+(pow2_spa17*pow2_spab_5_234_6*spab_3_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa15*pow2_spb68*spa13)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mppq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa15*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qppqb2mq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spab_4_56_3=(pow(spab_4_56_3,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46*spaa_4_56_178_3)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa14*pow2_spab_7_68_5*spa13*spab_4_123_5)/(s1234*s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa17*pow2_spab_4_56_3*spa46)/(s3456*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_56_8*spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qppqb2mpq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5)+(pow2_spa14*pow2_spab_7_68_5*spab_4_123_5)/(s1234*s678*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa14*pow2_spab_4_56_8*spa46)/(spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpppqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spab_3_456_2=(pow(spab_3_456_2,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
((pow2_spa13*pow2_spaa_3_45_68_7*spaa_3_12_45_3)/(s678*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1)
-
(pow2_spa13*pow2_spab_7_568_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)+(pow2_spa17*pow2_spab_3_456_2*spa36)/(s178*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)
-
(pow2_spa17*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_456_8*spa36)/(spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mq2pmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
-
(pow2_spab_5_123_4*pow2_spb68*spb14)/(s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)+(pow2_spb68*pow3_spab_5_23_4)/(s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)+(pow2_spab_3_127_8*pow2_spb46*spab_3_56_4)/(s3456*s456*spab_3_45_6*spb12*spb45*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow2_spb46*spb14)/(spa78*spb12*spb23*spb34*spb45*spbb_4_123_78_6*spbb_4_56_78_1)+(pow2_spb46*pow2_spbb_4_23_17_8*spbb_4_56_23_4)/(s178*spb23*spb34*spb45*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mmq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
(-((pow2_spab_2_17_8*pow2_spb36*spab_2_456_3)/(s178*s3456*spab_2_178_6*spb34*spb45*spb78*spbb_3_456_78_1))
-
(pow2_spab_7_12_3*pow2_spb36*spb13)/(spa78*spb12*spb23*spb34*spb45*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3)/(s2345*s345*spab_2_178_6*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpqb2mmmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
((pow2_spa17*pow2_spb26)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_345_2*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmqb2mq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spbb_3_12_678_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)+(pow2_spab_4_12_3*pow2_spb68*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5)
-
(pow2_spab_2_17_8*pow2_spb36*spab_2_456_3*spb46)/(s178*s3456*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow2_spb36*spb13*spb46)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13*spbb_3_12_678_4)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmqb2mmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_4=spa16*spb64+spa17*spb74+spa18*spb84;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb46)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_4)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68)/(s1234*s234*spab_1_678_5*spb23*spb34*spb56*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q2Q2g2l_qpmmqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_3=spa16*spb63+spa17*spb73+spa18*spb83;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb36)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_3)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68*spab_1_24_3)/(s1234*s234*spab_1_23_4*spab_1_678_5*spb23*spb34*spb56*spb78)+(pow2_spa13*pow2_spb68)/(s123*spa23*spab_1_23_4*spb45*spb56*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2pqb2mppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spab_3_456_2=(pow(spab_3_456_2,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
((pow2_spa13*pow2_spaa_3_45_68_7*spaa_3_12_45_3)/(s678*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1)
-
(pow2_spa13*pow2_spab_7_568_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)+(pow2_spa17*pow2_spab_3_456_2*spa36)/(s178*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)
-
(pow2_spa17*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_456_8*spa36)/(spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2ppqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5)+(pow2_spa14*pow2_spab_7_68_5*spab_4_123_5)/(s1234*s678*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa14*pow2_spab_4_56_8*spa46)/(spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2pppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa15*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmpq2pqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spab_4_56_3=(pow(spab_4_56_3,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46*spaa_4_56_178_3)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa14*pow2_spab_7_68_5*spa13*spab_4_123_5)/(s1234*s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa17*pow2_spab_4_56_3*spa46)/(s3456*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_56_8*spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmpq2ppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_178_6=spa13*spb61+spa37*spb76+spa38*spb86;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6)+(pow2_spa17*pow2_spab_5_234_6*spab_3_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa15*pow2_spb68*spa13)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmppq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_178_6=spa14*spb61+spa47*spb76+spa48*spb86;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6*spab_4_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa17*pow2_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6*spab_3_45_6)+(pow2_spa17*pow2_spb46)/(s456*spa12*spa23*spa78*spab_3_45_6*spb45)+(pow2_spa15*pow2_spb68*spa14)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2pqb2mmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_3=spa16*spb63+spa17*spb73+spa18*spb83;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb36)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_3)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68*spab_1_24_3)/(s1234*s234*spab_1_23_4*spab_1_678_5*spb23*spb34*spb56*spb78)+(pow2_spa13*pow2_spb68)/(s123*spa23*spab_1_23_4*spb45*spb56*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2pmqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_4=spa16*spb64+spa17*spb74+spa18*spb84;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb46)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_4)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68)/(s1234*s234*spab_1_678_5*spb23*spb34*spb56*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmq2pmmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
((pow2_spa17*pow2_spb26)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_345_2*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmq2pqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spbb_3_12_678_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)+(pow2_spab_4_12_3*pow2_spb68*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5)
-
(pow2_spab_2_17_8*pow2_spb36*spab_2_456_3*spb46)/(s178*s3456*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow2_spb36*spb13*spb46)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13*spbb_3_12_678_4)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmq2pmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
(-((pow2_spab_2_17_8*pow2_spb36*spab_2_456_3)/(s178*s3456*spab_2_178_6*spb34*spb45*spb78*spbb_3_456_78_1))
-
(pow2_spab_7_12_3*pow2_spb36*spb13)/(spa78*spb12*spb23*spb34*spb45*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3)/(s2345*s345*spab_2_178_6*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmmq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
-
(pow2_spab_5_123_4*pow2_spb68*spb14)/(s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)+(pow2_spb68*pow3_spab_5_23_4)/(s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)+(pow2_spab_3_127_8*pow2_spb46*spab_3_56_4)/(s3456*s456*spab_3_45_6*spb12*spb45*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow2_spb46*spb14)/(spa78*spb12*spb23*spb34*spb45*spbb_4_123_78_6*spbb_4_56_78_1)+(pow2_spb46*pow2_spbb_4_23_17_8*spbb_4_56_23_4)/(s178*spb23*spb34*spb45*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mq2pppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_678_1=spa36*spb61+spa37*spb71+spa38*spb81;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-((pow2_spaa_2_34_568_7*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_3_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+pow2_spab_7_12_3/(s123*spa45*spa56*spa78*spab_4_23_1*spb23)+(pow2_spab_2_17_8*spa36)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mpq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_678_1=spa46*spb61+spa47*spb71+spa48*spb81;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-(pow2_spaa_2_34_568_7/(s234*s5678*spa23*spa34*spa56*spa78*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_4_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+(pow2_spab_2_17_8*spa46)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mppq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_345_6=(pow(spab_2_345_6,2));

return(
((pow2_spa17*pow2_spab_2_345_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa12*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmpqb2mq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_678_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13*spaa_3_12_678_4)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1))+(pow2_spab_7_568_4*pow3_spa13)/(s5678*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)
-
(pow2_spa17*pow3_spab_3_456_2*spa46)/(s178*s3456*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_456_8*pow3_spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmpqb2mpq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_78_6*spaa_3_45_678_1))
-
(pow2_spa17*pow3_spab_3_456_2)/(s178*s3456*spa34*spa45*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_6_178_2)
-
(pow2_spab_3_456_8*pow3_spa13)/(spa12*spa23*spa34*spa45*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmppqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
(-((pow2_spa17*pow3_spaa_4_23_56_4)/(s178*spa23*spa34*spa45*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spab_7_68_5*pow3_spa14)/(s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)+(pow2_spa17*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa17*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_56_8*pow3_spa14)/(spa12*spa23*spa34*spa45*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mq2pmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=(pow(spb36,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb78)
-
(pow2_spb68*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb78*spbb_3_12_678_5)+(pow2_spab_2_17_8*pow3_spb36)/(s178*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow3_spb36)/(spa78*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)
-
(pow2_spb68*pow3_spab_2_45_3)/(s2345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow3_spbb_3_45_12_3)/(s678*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mmq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));

return(
((pow2_spb68*pow3_spab_5_123_4)/(s1234*s678*spab_5_678_1*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4)/(s234*s2345*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_7_123_4*pow3_spb46)/(spa78*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmqb2mmmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb56=(pow(spb56,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_234_5=(pow(spab_1_234_5,2));

return(
((pow2_spa17*pow2_spb56)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_234_5*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmqb2mq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb78)+(pow2_spb68*pow3_spab_5_123_4*spb13)/(s1234*s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_3_127_8*pow3_spb46)/(s3456*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow3_spb46*spb13)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46*spbb_4_56_178_3)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmqb2mmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_3=spa16*spb31+spa67*spb73+spa68*spb83;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));

return(
((pow2_spab_7_68_5*spb13)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)
-
pow2_spbb_5_34_127_8/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_3)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qmmmqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_4=spa16*spb41+spa67*spb74+spa68*spb84;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
((pow2_spab_7_68_5*spb14)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)+pow2_spab_4_56_8/(s456*spa45*spab_6_45_3*spb12*spb23*spb78)
-
(pow2_spbb_5_34_127_8*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_4)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2pqb2mppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
(-((pow2_spa17*pow3_spaa_4_23_56_4)/(s178*spa23*spa34*spa45*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spab_7_68_5*pow3_spa14)/(s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)+(pow2_spa17*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa17*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_56_8*pow3_spa14)/(spa12*spa23*spa34*spa45*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2ppqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_78_6*spaa_3_45_678_1))
-
(pow2_spa17*pow3_spab_3_456_2)/(s178*s3456*spa34*spa45*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_6_178_2)
-
(pow2_spab_3_456_8*pow3_spa13)/(spa12*spa23*spa34*spa45*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2pppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa12=(pow(spa12,2));
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_345_6=(pow(spab_2_345_6,2));

return(
((pow2_spa17*pow2_spab_2_345_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa12*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qppq2pqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_678_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
(-((pow2_spaa_3_45_68_7*pow3_spa13*spaa_3_12_678_4)/(s678*spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1))+(pow2_spab_7_568_4*pow3_spa13)/(s5678*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)
-
(pow2_spa17*pow3_spab_3_456_2*spa46)/(s178*s3456*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)+(pow2_spa17*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_456_8*pow3_spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qppq2ppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_678_1=spa46*spb61+spa47*spb71+spa48*spb81;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s5678=(spa56*spb65+spa57*spb75+spa67*spb76+spa58*spb85+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-(pow2_spaa_2_34_568_7/(s234*s5678*spa23*spa34*spa56*spa78*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_4_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+(pow2_spab_2_17_8*spa46)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpppq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_678_1=spa36*spb61+spa37*spb71+spa38*spb81;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*spb65)-spa78*(-(spa23*spb83)-spa24*spb84-spa25*spb85);
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));
const complex<T> pow2_spaa_2_345_68_7=(pow(spaa_2_345_68_7,2));
const complex<T> pow2_spaa_2_34_568_7=(pow(spaa_2_34_568_7,2));

return(
(-((pow2_spaa_2_34_568_7*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1))
-
(pow2_spaa_2_345_68_7*spab_3_678_1)/(s2345*s678*spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1)+pow2_spab_7_12_3/(s123*spa45*spa56*spa78*spab_4_23_1*spb23)+(pow2_spab_2_17_8*spa36)/(s178*spa23*spa34*spa45*spa56*spab_6_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2pqb2mmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_4=spa16*spb41+spa67*spb74+spa68*spb84;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));

return(
((pow2_spab_7_68_5*spb14)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)+pow2_spab_4_56_8/(s456*spa45*spab_6_45_3*spb12*spb23*spb78)
-
(pow2_spbb_5_34_127_8*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_4)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2pmqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb81)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_178_3=spa16*spb31+spa67*spb73+spa68*spb83;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spbb_5_34_127_8=(pow(spbb_5_34_127_8,2));
const complex<T> pow2_spbb_5_234_17_8=(pow(spbb_5_234_17_8,2));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));

return(
((pow2_spab_7_68_5*spb13)/(s678*spa78*spab_6_78_1*spb12*spb23*spb34*spb45)
-
pow2_spbb_5_34_127_8/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb78)
-
(pow2_spbb_5_234_17_8*spab_6_178_3)/(s178*s2345*spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpq2pmmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb56=(pow(spb56,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_234_5=(pow(spab_1_234_5,2));

return(
((pow2_spa17*pow2_spb56)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_234_5*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmq2pqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb78)+(pow2_spb68*pow3_spab_5_123_4*spb13)/(s1234*s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_3_127_8*pow3_spb46)/(s3456*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow3_spb46*spb13)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46*spbb_4_56_178_3)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmq2pmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));

return(
((pow2_spb68*pow3_spab_5_123_4)/(s1234*s678*spab_5_678_1*spb23*spb34*spb78*spbb_4_123_78_6)
-
(pow2_spb68*pow3_spab_5_23_4)/(s234*s2345*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)
-
(pow2_spab_7_123_4*pow3_spb46)/(spa78*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1)
-
(pow2_spbb_4_23_17_8*pow3_spb46)/(s178*spb23*spb34*spb45*spb56*spb78*spbb_4_23_178_6*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmmq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=(pow(spb36,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb78)
-
(pow2_spb68*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb78*spbb_3_12_678_5)+(pow2_spab_2_17_8*pow3_spb36)/(s178*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow3_spb36)/(spa78*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)
-
(pow2_spb68*pow3_spab_2_45_3)/(s2345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow3_spbb_3_45_12_3)/(s678*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mq2pppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_178_6=spa14*spb61+spa47*spb76+spa48*spb86;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6*spab_4_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa17*pow2_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6*spab_3_45_6)+(pow2_spa17*pow2_spb46)/(s456*spa12*spa23*spa78*spab_3_45_6*spb45)+(pow2_spa15*pow2_spb68*spa14)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mpq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_178_6=spa13*spb61+spa37*spb76+spa38*spb86;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_34_6=(pow(spab_5_34_6,2));
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa78*spab_2_178_6)+(pow2_spa17*pow2_spab_5_234_6*spab_3_178_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_178_6)+(pow2_spa15*pow2_spb68*spa13)/(s678*spa12*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mppq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_5_234_6=(pow(spab_5_234_6,2));

return(
((pow2_spa17*pow2_spab_5_234_6)/(s178*s2345*spa23*spa34*spa45*spa78*spab_1_78_6)+(pow2_spa15*pow2_spb68)/(s678*spa23*spa34*spa45*spab_1_78_6*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qppqb2mq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spab_4_56_3=(pow(spab_4_56_3,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46*spaa_4_56_178_3)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5*spab_2_34_5)+(pow2_spa14*pow2_spab_7_68_5*spa13*spab_4_123_5)/(s1234*s678*spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa17*pow2_spab_4_56_3*spa46)/(s3456*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3)+(pow2_spa17*pow3_spb35)/(s345*spa12*spa78*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_56_8*spa13*spa46)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qppqb2mpq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow2_spab_7_68_5=(pow(spab_7_68_5,2));
const complex<T> pow2_spab_4_56_8=(pow(spab_4_56_8,2));
const complex<T> pow2_spaa_4_23_56_4=(pow(spaa_4_23_56_4,2));

return(
(-((pow2_spa17*pow2_spaa_4_23_56_4*spa46)/(s178*spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_78_1))
-
(pow2_spa17*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_678_5)+(pow2_spa14*pow2_spab_7_68_5*spab_4_123_5)/(s1234*s678*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_678_5)
-
(pow2_spa14*pow2_spab_4_56_8*spa46)/(spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpppqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_7_568_4=-(spa57*spb54)-spa67*spb64+spa78*spb84;
const complex<T> spab_6_178_2=spa16*spb21+spa67*spb72+spa68*spb82;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_7_568_4=(pow(spab_7_568_4,2));
const complex<T> pow2_spab_3_456_8=(pow(spab_3_456_8,2));
const complex<T> pow2_spab_3_456_2=(pow(spab_3_456_2,2));
const complex<T> pow2_spaa_3_45_68_7=(pow(spaa_3_45_68_7,2));

return(
((pow2_spa13*pow2_spaa_3_45_68_7*spaa_3_12_45_3)/(s678*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1)
-
(pow2_spa13*pow2_spab_7_568_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4)+(pow2_spa17*pow2_spab_3_456_2*spa36)/(s178*spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_178_2)
-
(pow2_spa17*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_178_2)+(pow2_spa17*pow3_spb24)/(s234*spa56*spa78*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_456_8*spa36)/(spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*spb78))*complex<T>(0,1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mq2pmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_678_1=spa56*spb61+spa57*spb71+spa58*spb81;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_127_8=spa13*spb81+spa23*spb82-spa37*spb87;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow2_spbb_4_23_17_8=(pow(spbb_4_23_17_8,2));
const complex<T> pow2_spab_7_123_4=(pow(spab_7_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_127_8=(pow(spab_3_127_8,2));

return(
((pow2_spb68*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
-
(pow2_spab_5_123_4*pow2_spb68*spb14)/(s678*spab_5_678_1*spb12*spb23*spb34*spb78*spbb_4_123_78_6)+(pow2_spb68*pow3_spab_5_23_4)/(s2345*spab_5_34_2*spab_5_678_1*spb23*spb34*spb78*spbb_4_23_178_6)+(pow2_spab_3_127_8*pow2_spb46*spab_3_56_4)/(s3456*s456*spab_3_45_6*spb12*spb45*spb78*spbb_4_56_178_2)
-
(pow2_spab_7_123_4*pow2_spb46*spb14)/(spa78*spb12*spb23*spb34*spb45*spbb_4_123_78_6*spbb_4_56_78_1)+(pow2_spb46*pow2_spbb_4_23_17_8*spbb_4_56_23_4)/(s178*spb23*spb34*spb45*spb78*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mmq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
(-((pow2_spab_2_17_8*pow2_spb36*spab_2_456_3)/(s178*s3456*spab_2_178_6*spb34*spb45*spb78*spbb_3_456_78_1))
-
(pow2_spab_7_12_3*pow2_spb36*spb13)/(spa78*spb12*spb23*spb34*spb45*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3)/(s2345*s345*spab_2_178_6*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpqb2mmmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
((pow2_spa17*pow2_spb26)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45)+(pow2_spab_1_345_2*pow2_spb68)/(s2345*s678*spab_1_78_6*spb23*spb34*spb45*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmqb2mq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb36=(pow(spb36,2));
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spbb_3_12_678_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spab_2_178_6=spa12*spb61+spa27*spb76+spa28*spb86;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_7_12_3=(pow(spab_7_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_17_8=(pow(spab_2_17_8,2));

return(
((pow2_spb68*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)+(pow2_spab_4_12_3*pow2_spb68*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5)
-
(pow2_spab_2_17_8*pow2_spb36*spab_2_456_3*spb46)/(s178*s3456*spab_2_178_6*spb34*spb45*spb56*spb78*spbb_3_456_78_1)
-
(pow2_spab_7_12_3*pow2_spb36*spb13*spb46)/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1)+(pow2_spb68*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_178_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1)
-
(pow2_spb68*pow2_spbb_3_45_12_3*spb13*spbb_3_12_678_4)/(s678*spb12*spb23*spb34*spb45*spb78*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmqb2mmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_4=spa16*spb64+spa17*spb74+spa18*spb84;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb46)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_4)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68)/(s1234*s234*spab_1_678_5*spb23*spb34*spb56*spb78))*complex<T>(0,-1)
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q2Q2g2l_qpmmqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb68=SPB(i6,i8);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=(pow(spb68,2));
const complex<T> pow2_spb26=(pow(spb26,2));
const complex<T> pow2_spa17=(pow(spa17,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_678_5=spa16*spb65+spa17*spb75+spa18*spb85;
const complex<T> spab_1_678_3=spa16*spb63+spa17*spb73+spa18*spb83;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s678=(spa67*spb76+spa68*spb86+spa78*spb87);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s178=(spa17*spb71+spa18*spb81+spa78*spb87);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
((pow2_spa17*pow2_spb26*spb36)/(s178*spa78*spab_1_78_6*spb23*spb34*spb45*spb56)+(pow2_spab_1_345_2*pow2_spb68*spab_1_678_3)/(s2345*s678*spab_1_678_5*spab_1_78_6*spb23*spb34*spb45*spb78)+(pow2_spab_1_34_2*pow2_spb68*spab_1_24_3)/(s1234*s234*spab_1_23_4*spab_1_678_5*spb23*spb34*spb56*spb78)+(pow2_spa13*pow2_spb68)/(s123*spa23*spab_1_23_4*spb45*spb56*spb78))*complex<T>(0,-1)
);
}




template <class T> complex<T> (*A2q2Q2g2l_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 14427937:	 return &A2q2Q2g2l_qpmmq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m m q2p qb2m qbm lp lbm
case 14428009:	 return &A2q2Q2g2l_qpppq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p p q2p qb2m qbm lp lbm
case 14428385:	 return &A2q2Q2g2l_qpmq2pmqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m q2p m qb2m qbm lp lbm
case 14428441:	 return &A2q2Q2g2l_qpq2pmmqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p m m qb2m qbm lp lbm
case 14428905:	 return &A2q2Q2g2l_qppq2ppqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p q2p p qb2m qbm lp lbm
case 14429017:	 return &A2q2Q2g2l_qpq2pppqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p p p qb2m qbm lp lbm
case 14431521:	 return &A2q2Q2g2l_qpmmqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m m qb2m q2p qbm lp lbm
case 14431593:	 return &A2q2Q2g2l_qpppqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p p qb2m q2p qbm lp lbm
case 14432417:	 return &A2q2Q2g2l_qpmqb2mmq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m qb2m m q2p qbm lp lbm
case 14432529:	 return &A2q2Q2g2l_qpqb2mmmq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m m m q2p qbm lp lbm
case 14432937:	 return &A2q2Q2g2l_qppqb2mpq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p qb2m p q2p qbm lp lbm
case 14433105:	 return &A2q2Q2g2l_qpqb2mppq2pqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m p p q2p qbm lp lbm
case 14435553:	 return &A2q2Q2g2l_qpmq2pqb2mmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m q2p qb2m m qbm lp lbm
case 14435609:	 return &A2q2Q2g2l_qpq2pmqb2mmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p m qb2m m qbm lp lbm
case 14436001:	 return &A2q2Q2g2l_qpmqb2mq2pmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m qb2m q2p m qbm lp lbm
case 14436113:	 return &A2q2Q2g2l_qpqb2mmq2pmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m m q2p m qbm lp lbm
case 14436505:	 return &A2q2Q2g2l_qpq2pqb2mmmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p qb2m m m qbm lp lbm
case 14436561:	 return &A2q2Q2g2l_qpqb2mq2pmmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m q2p m m qbm lp lbm
case 14439657:	 return &A2q2Q2g2l_qppq2pqb2mpqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p q2p qb2m p qbm lp lbm
case 14439769:	 return &A2q2Q2g2l_qpq2ppqb2mpqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p p qb2m p qbm lp lbm
case 14440105:	 return &A2q2Q2g2l_qppqb2mq2ppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p qb2m q2p p qbm lp lbm
case 14440273:	 return &A2q2Q2g2l_qpqb2mpq2ppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m p q2p p qbm lp lbm
case 14441113:	 return &A2q2Q2g2l_qpq2pqb2mppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp q2p qb2m p p qbm lp lbm
case 14441169:	 return &A2q2Q2g2l_qpqb2mq2pppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp qb2m q2p p p qbm lp lbm
case 14460704:	 return &A2q2Q2g2l_qmmmq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m m q2p qb2m qbp lp lbm
case 14460776:	 return &A2q2Q2g2l_qmppq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p p q2p qb2m qbp lp lbm
case 14461152:	 return &A2q2Q2g2l_qmmq2pmqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m q2p m qb2m qbp lp lbm
case 14461208:	 return &A2q2Q2g2l_qmq2pmmqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p m m qb2m qbp lp lbm
case 14461672:	 return &A2q2Q2g2l_qmpq2ppqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p q2p p qb2m qbp lp lbm
case 14461784:	 return &A2q2Q2g2l_qmq2pppqb2mqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p p p qb2m qbp lp lbm
case 14464288:	 return &A2q2Q2g2l_qmmmqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m m qb2m q2p qbp lp lbm
case 14464360:	 return &A2q2Q2g2l_qmppqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p p qb2m q2p qbp lp lbm
case 14465184:	 return &A2q2Q2g2l_qmmqb2mmq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m qb2m m q2p qbp lp lbm
case 14465296:	 return &A2q2Q2g2l_qmqb2mmmq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m m m q2p qbp lp lbm
case 14465704:	 return &A2q2Q2g2l_qmpqb2mpq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p qb2m p q2p qbp lp lbm
case 14465872:	 return &A2q2Q2g2l_qmqb2mppq2pqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m p p q2p qbp lp lbm
case 14468320:	 return &A2q2Q2g2l_qmmq2pqb2mmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m q2p qb2m m qbp lp lbm
case 14468376:	 return &A2q2Q2g2l_qmq2pmqb2mmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p m qb2m m qbp lp lbm
case 14468768:	 return &A2q2Q2g2l_qmmqb2mq2pmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m qb2m q2p m qbp lp lbm
case 14468880:	 return &A2q2Q2g2l_qmqb2mmq2pmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m m q2p m qbp lp lbm
case 14469272:	 return &A2q2Q2g2l_qmq2pqb2mmmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p qb2m m m qbp lp lbm
case 14469328:	 return &A2q2Q2g2l_qmqb2mq2pmmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m q2p m m qbp lp lbm
case 14472424:	 return &A2q2Q2g2l_qmpq2pqb2mpqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p q2p qb2m p qbp lp lbm
case 14472536:	 return &A2q2Q2g2l_qmq2ppqb2mpqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p p qb2m p qbp lp lbm
case 14472872:	 return &A2q2Q2g2l_qmpqb2mq2ppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p qb2m q2p p qbp lp lbm
case 14473040:	 return &A2q2Q2g2l_qmqb2mpq2ppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m p q2p p qbp lp lbm
case 14473880:	 return &A2q2Q2g2l_qmq2pqb2mppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm q2p qb2m p p qbp lp lbm
case 14473936:	 return &A2q2Q2g2l_qmqb2mq2pppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm qb2m q2p p p qbp lp lbm
case 16262945:	 return &A2q2Q2g2l_qpmmq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m m q2p qb2m qbm lm lbp
case 16263017:	 return &A2q2Q2g2l_qpppq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p p q2p qb2m qbm lm lbp
case 16263393:	 return &A2q2Q2g2l_qpmq2pmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m q2p m qb2m qbm lm lbp
case 16263449:	 return &A2q2Q2g2l_qpq2pmmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p m m qb2m qbm lm lbp
case 16263913:	 return &A2q2Q2g2l_qppq2ppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p q2p p qb2m qbm lm lbp
case 16264025:	 return &A2q2Q2g2l_qpq2pppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p p p qb2m qbm lm lbp
case 16266529:	 return &A2q2Q2g2l_qpmmqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m m qb2m q2p qbm lm lbp
case 16266601:	 return &A2q2Q2g2l_qpppqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p p qb2m q2p qbm lm lbp
case 16267425:	 return &A2q2Q2g2l_qpmqb2mmq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m qb2m m q2p qbm lm lbp
case 16267537:	 return &A2q2Q2g2l_qpqb2mmmq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m m m q2p qbm lm lbp
case 16267945:	 return &A2q2Q2g2l_qppqb2mpq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p qb2m p q2p qbm lm lbp
case 16268113:	 return &A2q2Q2g2l_qpqb2mppq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m p p q2p qbm lm lbp
case 16270561:	 return &A2q2Q2g2l_qpmq2pqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m q2p qb2m m qbm lm lbp
case 16270617:	 return &A2q2Q2g2l_qpq2pmqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p m qb2m m qbm lm lbp
case 16271009:	 return &A2q2Q2g2l_qpmqb2mq2pmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m qb2m q2p m qbm lm lbp
case 16271121:	 return &A2q2Q2g2l_qpqb2mmq2pmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m m q2p m qbm lm lbp
case 16271513:	 return &A2q2Q2g2l_qpq2pqb2mmmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p qb2m m m qbm lm lbp
case 16271569:	 return &A2q2Q2g2l_qpqb2mq2pmmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m q2p m m qbm lm lbp
case 16274665:	 return &A2q2Q2g2l_qppq2pqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p q2p qb2m p qbm lm lbp
case 16274777:	 return &A2q2Q2g2l_qpq2ppqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p p qb2m p qbm lm lbp
case 16275113:	 return &A2q2Q2g2l_qppqb2mq2ppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p qb2m q2p p qbm lm lbp
case 16275281:	 return &A2q2Q2g2l_qpqb2mpq2ppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m p q2p p qbm lm lbp
case 16276121:	 return &A2q2Q2g2l_qpq2pqb2mppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp q2p qb2m p p qbm lm lbp
case 16276177:	 return &A2q2Q2g2l_qpqb2mq2pppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp qb2m q2p p p qbm lm lbp
case 16295712:	 return &A2q2Q2g2l_qmmmq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m m q2p qb2m qbp lm lbp
case 16295784:	 return &A2q2Q2g2l_qmppq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p p q2p qb2m qbp lm lbp
case 16296160:	 return &A2q2Q2g2l_qmmq2pmqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m q2p m qb2m qbp lm lbp
case 16296216:	 return &A2q2Q2g2l_qmq2pmmqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p m m qb2m qbp lm lbp
case 16296680:	 return &A2q2Q2g2l_qmpq2ppqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p q2p p qb2m qbp lm lbp
case 16296792:	 return &A2q2Q2g2l_qmq2pppqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p p p qb2m qbp lm lbp
case 16299296:	 return &A2q2Q2g2l_qmmmqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m m qb2m q2p qbp lm lbp
case 16299368:	 return &A2q2Q2g2l_qmppqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p p qb2m q2p qbp lm lbp
case 16300192:	 return &A2q2Q2g2l_qmmqb2mmq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m qb2m m q2p qbp lm lbp
case 16300304:	 return &A2q2Q2g2l_qmqb2mmmq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m m m q2p qbp lm lbp
case 16300712:	 return &A2q2Q2g2l_qmpqb2mpq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p qb2m p q2p qbp lm lbp
case 16300880:	 return &A2q2Q2g2l_qmqb2mppq2pqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m p p q2p qbp lm lbp
case 16303328:	 return &A2q2Q2g2l_qmmq2pqb2mmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m q2p qb2m m qbp lm lbp
case 16303384:	 return &A2q2Q2g2l_qmq2pmqb2mmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p m qb2m m qbp lm lbp
case 16303776:	 return &A2q2Q2g2l_qmmqb2mq2pmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m qb2m q2p m qbp lm lbp
case 16303888:	 return &A2q2Q2g2l_qmqb2mmq2pmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m m q2p m qbp lm lbp
case 16304280:	 return &A2q2Q2g2l_qmq2pqb2mmmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p qb2m m m qbp lm lbp
case 16304336:	 return &A2q2Q2g2l_qmqb2mq2pmmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m q2p m m qbp lm lbp
case 16307432:	 return &A2q2Q2g2l_qmpq2pqb2mpqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p q2p qb2m p qbp lm lbp
case 16307544:	 return &A2q2Q2g2l_qmq2ppqb2mpqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p p qb2m p qbp lm lbp
case 16307880:	 return &A2q2Q2g2l_qmpqb2mq2ppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p qb2m q2p p qbp lm lbp
case 16308048:	 return &A2q2Q2g2l_qmqb2mpq2ppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m p q2p p qbp lm lbp
case 16308888:	 return &A2q2Q2g2l_qmq2pqb2mppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm q2p qb2m p p qbp lm lbp
case 16308944:	 return &A2q2Q2g2l_qmqb2mq2pppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm qb2m q2p p p qbp lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q2Q2g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q2Q2g2l_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2Q2g2l_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2Q2g2l_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);
#if BH_USE_GMP
template complex<RGMP> (*A2q2Q2g2l_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif
}
