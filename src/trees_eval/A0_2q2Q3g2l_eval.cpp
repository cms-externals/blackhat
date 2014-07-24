/*
 * A0_2q2Q3g2l_eval.cpp
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
 * The 2q2Q quarks 3g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pqb2mpppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_456_3=-(spa13*(spa34*spb41+spa35*spb51+spa36*spb61))-spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_189_2=(pow(spab_3_189_2,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
((pow2_spa13*pow2_spaa_3_456_79_8*spaa_3_12_456_3)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)+(pow2_spa13*pow2_spaa_3_45_679_8*spaa_3_12_45_3)/(s6789*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa13*pow2_spab_8_123_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)+(pow2_spa18*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow2_spab_3_189_2*spa37)/(s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)
-
(pow2_spa18*pow3_spab_3_456_2)/(s1789*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_128_9*spa37)/(spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2ppqb2mppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
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
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spaa_4_123_56_4)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_8_679_5*spab_4_123_5)/(s1234*s6789*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_4_567_9*spa47)/(spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pppqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_234_67_5*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_89_1))
-
(pow2_spa18*pow3_spab_5_234_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_8_79_6*spab_5_789_6)/(s6789*s789*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa57)/(spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2ppppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7)+(pow2_spa16*pow2_spb79)/(s789*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpq2pqb2mppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa34=(pow(spa34,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_3=spa45*(spa13*spb51+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa38*spb86+spa39*spb96)+spa47*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_56_1789_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spab_4_567_3=(pow(spab_4_567_3,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47*spaa_4_567_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spa13*spaa_4_123_56_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4*spaa_4_56_1789_3)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)
-
(pow2_spa14*pow2_spab_8_679_5*spa13*spab_4_123_5)/(s1234*s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow2_spab_4_567_3*spa47)/(s1289*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)+(pow2_spa18*pow3_spab_4_56_3)/(s3456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_567_9*spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpq2ppqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_245_6=spa23*spb62-spa34*spb64-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_3=spa56*(spa13*spb61+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)+(pow2_spa18*pow3_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_234_6*spab_3_245_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa13*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpq2pppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_3_189_7=spa13*spb71+spa38*spb87+spa39*spb97;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_345_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7)+(pow2_spa18*pow2_spab_6_189_7*spab_3_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa16*pow2_spb79*spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmppq2pqb2mpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_235_6=spa24*spb62+spa34*spb63-spa45*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_4=spa56*(spa14*spb61+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_4=spa56*(spa14*spb61+spa24*spb62+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spab_5_67_4=(pow(spab_5_67_4,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57*spaa_5_67_1289_4)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_4)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spa18*pow3_spab_5_234_6*spab_4_235_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)+(pow2_spa18*pow3_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa14*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa18*pow2_spab_5_67_4*spa57)/(s4567*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb45)
-
(pow2_spa15*pow2_spab_5_67_9*spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmppq2ppqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_4_356_7=spa34*spb73-spa45*spb75-spa46*spb76;
const complex<T> spab_4_189_7=spa14*spb71+spa48*spb87+spa49*spb97;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
(-((pow2_spa18*pow2_spab_6_45_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7))+(pow2_spa18*pow2_spab_6_189_7*spab_4_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_4_356_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)+(pow2_spa16*pow2_spb79*spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpppq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_5_46_7=spa45*spb74-spa56*spb76;
const complex<T> spab_5_346_7=spa35*spb73+spa45*spb74-spa56*spb76;
const complex<T> spab_5_189_7=spa15*spb71+spa58*spb87+spa59*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7*spab_5_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_5_346_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)
-
(pow2_spa18*pow2_spab_6_45_7*spab_5_46_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7)+(pow2_spa18*pow2_spb57)/(s567*spa12*spa23*spa34*spa89*spab_4_56_7*spb56)+(pow2_spa16*pow2_spb79*spa15)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pqb2mmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_3=spa17*spb73+spa18*spb83+spa19*spb93;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_245_3=-(spa12*spb32)+spa14*spb43+spa15*spb53;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb37)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_3)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_245_3)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79*spab_1_24_3)/(s1234*s234*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb89)
-
(pow2_spa13*pow2_spb79)/(s123*spa23*spab_1_23_4*spb45*spb56*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pmqb2mmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_4=spa17*spb74+spa18*spb84+spa19*spb94;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_235_4=-(spa12*spb42)-spa13*spb43+spa15*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb47)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_4)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_235_4)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79)/(s1234*s234*spab_1_234_5*spb23*spb34*spb56*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pmmqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_5=spa17*spb75+spa18*spb85+spa19*spb95;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb57)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_5)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79)/(s2345*s6789*spab_1_789_6*spb23*spb34*spb45*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmq2pmmmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));

return(
(-((pow2_spa18*pow2_spb27)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56))
-
(pow2_spab_1_789_2*pow2_spb79)/(s1789*s789*spab_1_89_7*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmq2pqb2mmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
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
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_4=spb31*(spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spbb_3_12_6789_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_356_4=-(spa23*spb43)+spa25*spb54+spa26*spb64;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))
-
(pow2_spab_4_12_3*pow2_spb79*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb47)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow2_spb37*spb13*spb47)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_356_4)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_4)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13*spbb_3_12_6789_4)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmq2pmqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_5=spb31*(spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_346_5=-(spa23*spb53)-spa24*spb54+spa26*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb57)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_346_5)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3)/(s2345*s345*spab_2_345_6*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmq2pmmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3)/(s1789*s3456*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_89_7*spbb_3_456_789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmq2pqb2mmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_23_1789_5=spb42*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_789_5=spb41*(spa17*spb75+spa18*spb85+spa19*spb95)+spb42*(spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_46_5=-(spa34*spb54)+spa36*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))
-
(pow2_spab_5_123_4*pow2_spb79*spb14)/(s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)+(pow2_spb79*pow3_spab_5_23_4)/(s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)
-
(pow2_spab_3_128_9*pow2_spb47*spab_3_567_4*spb57)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow2_spb47*spb14*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spb57*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4*spab_3_46_5)/(s3456*s456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14*spbb_4_123_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4*spbb_4_23_1789_5)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmq2pmqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spab_3_128_9*pow2_spb47*spab_3_567_4)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_567_189_2))+(pow2_spab_8_123_4*pow2_spb47*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4)/(s3456*s456*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmmq2pqb2mqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_34_5=-((spa36*spb53+spa46*spb54)*spb65)-(spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_234_5=-((spa26*spb52+spa36*spb53+spa46*spb54)*spb65)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_6_789_5=(pow(spab_6_789_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spab_6_789_5*pow2_spb79*spb15)/(s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)
-
(pow2_spb79*pow3_spab_6_234_5)/(s1789*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spb79*pow3_spab_6_34_5)/(s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)
-
(pow2_spab_4_567_9*pow2_spb57*spab_4_67_5)/(s4567*s567*spab_4_56_7*spb12*spb23*spb56*spb89*spbb_5_67_1289_3)
-
(pow2_spb57*pow2_spbb_5_34_128_9*spbb_5_67_34_5)/(s1289*spb12*spb34*spb45*spb56*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow2_spb57*spb15)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_5_1234_89_7*spbb_5_67_89_1)
-
(pow2_spb57*pow2_spbb_5_234_18_9*spbb_5_67_234_5)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mq2ppppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_789_1=spa37*spb71+spa38*spb81+spa39*spb91;
const complex<T> spab_3_245_1=-(spa23*spb21)+spa34*spb41+spa35*spb51;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-((pow2_spaa_2_34_5679_8*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_3_245_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_3_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+pow2_spab_8_12_3/(s123*spa45*spa56*spa67*spa89*spab_4_23_1*spb23)+(pow2_spab_2_18_9*spa37)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mpq2pppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_789_1=spa47*spb71+spa48*spb81+spa49*spb91;
const complex<T> spab_4_235_1=-(spa24*spb21)-spa34*spb31+spa45*spb51;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-(pow2_spaa_2_34_5679_8/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_4_235_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_4_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mppq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_789_1=spa57*spb71+spa58*spb81+spa59*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));

return(
(pow2_spaa_2_345_679_8/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_5_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mpppq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));

return(
(pow2_spaa_2_3456_79_8/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_7_89_1)+pow2_spab_2_18_9/(s189*spa23*spa34*spa45*spa56*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpqb2mq2pppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_356_2=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_4=-(spa13*(spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_6789_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13*spaa_3_12_6789_4)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)+(pow2_spab_8_123_4*pow3_spa13)/(s1234*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)
-
(pow2_spa18*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa47)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_4_356_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_128_9*pow3_spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpqb2mpq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_346_2=-(spa35*spb32)-spa45*spb42+spa56*spb62;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_5=-(spa13*(spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa18*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa57)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_5_346_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpqb2mppq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_89_7*spaa_3_456_789_1))+(pow2_spa18*pow3_spab_3_189_2)/(s1289*s189*spa34*spa45*spa56*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmppqb2mq2ppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_5_46_3=-(spa45*spb43)+spa56*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_23_1789_5=-(spa24*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spa34*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_789_5=-(spa14*(spa57*spb71+spa58*spb81+spa59*spb91))-spa24*(spa57*spb72+spa58*spb82+spa59*spb92)-spa34*(spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14*spaa_4_123_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4*spaa_4_23_1789_5)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spab_8_679_5*pow3_spa14)/(s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)+(pow2_spa18*pow3_spab_4_567_3*spa57)/(s1289*s4567*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3*spab_5_46_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_567_9*pow3_spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmppqb2mpq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_567_3)/(s1289*s4567*spa12*spa45*spa56*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_7_456_3)
-
(pow2_spab_4_567_9*pow3_spa14)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmpppqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa15=(pow(spa15,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_67_4=(pow(spab_5_67_4,3));
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow3_spaa_5_34_67_5=(pow(spaa_5_34_67_5,3));
const complex<T> pow3_spaa_5_234_67_5=(pow(spaa_5_234_67_5,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spa18*pow3_spaa_5_34_67_5)/(s1289*spa12*spa34*spa45*spa56*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow3_spaa_5_234_67_5)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spab_8_79_6*pow3_spa15)/(s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)+(pow2_spa18*pow3_spab_5_234_6)/(s1789*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_34_6)/(s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)+(pow2_spa18*pow3_spab_5_67_4)/(s4567*s567*spa12*spa23*spa56*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb56)
-
(pow2_spab_5_67_9*pow3_spa15)/(spa12*spa23*spa34*spa45*spa56*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mq2pmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
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
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb37=(pow(spb37,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_456_12_3=(pow(spbb_3_456_12_3,3));
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))+(pow2_spb79*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow3_spb37)/(s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow3_spb37)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)+(pow2_spb79*pow3_spab_2_456_3)/(s1789*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow3_spbb_3_456_12_3)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)
-
(pow2_spb79*pow3_spab_2_45_3)/(s2345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow3_spbb_3_45_12_3)/(s6789*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mmq2pmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));

return(
((pow2_spb79*pow3_spab_5_123_4)/(s1234*s6789*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4)/(s234*s2345*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_8_123_4*pow3_spb47)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_89_1)+(pow2_spb79*pow3_spbb_4_56_123_4)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mmmq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5)/(s6789*s789*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5)/(s1789*s2345*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spab_8_679_5*pow3_spb57)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmqb2mmmmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-(pow2_spab_8_79_6/(s789*spa89*spab_7_89_1*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_2345_18_9/(s1789*s189*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmqb2mq2pmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_3=-(spb54*(spa15*spb31+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa68*spb83+spa69*spb93)-spb74*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))+(pow2_spb79*pow3_spab_5_123_4*spb13)/(s1234*s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_3_128_9*pow3_spb47)/(s1289*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow3_spb47*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47*spbb_4_567_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)
-
(pow2_spb79*pow3_spab_3_56_4)/(s3456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow3_spbb_4_56_123_4*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4*spbb_4_56_1789_3)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmqb2mmq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_3=-(spb65*(spa16*spb31+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_245_3=spa26*spb32-spa46*spb43-spa56*spb53;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5*spb13)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_245_3)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5)/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spbb_5_34_128_9*pow3_spb57)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmqb2mmmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_189_3=spa17*spb31+spa78*spb83+spa79*spb93;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb13)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_345_128_9/(s1289*s3456*spab_7_189_2*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_3)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmqb2mq2pmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_4=-(spb65*(spa16*spb41+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_4=-(spb65*(spa16*spb41+spa26*spb42+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_235_4=spa26*spb42+spa36*spb43-spa56*spb54;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa45*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spb79*pow3_spab_6_789_5*spb14)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_235_4)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spab_4_567_9*pow3_spb57)/(s4567*spab_4_56_7*spb12*spb23*spb56*spb67*spb89*spbb_5_67_1289_3)+(pow2_spbb_5_34_128_9*pow3_spb57*spbb_5_67_1289_4)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmqb2mmq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_356_4=spa37*spb43-spa57*spb54-spa67*spb64;
const complex<T> spab_7_189_4=spa17*spb41+spa78*spb84+spa79*spb94;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb14)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))+pow2_spbb_6_45_1238_9/(s456*s4567*spab_7_456_3*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_356_4)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_4)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qmmmmqb2mq2pqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_46_5=spa47*spb54-spa67*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_346_5=spa37*spb53+spa47*spb54-spa67*spb65;
const complex<T> spab_7_189_5=spa17*spb51+spa78*spb85+spa79*spb95;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spab_8_79_6*spb15)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spab_5_67_9/(s567*spa56*spab_7_56_4*spb12*spb23*spb34*spb89)+(pow2_spbb_6_45_1238_9*spab_7_46_5)/(s456*s4567*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_346_5)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_5)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pqb2mpppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa15=(pow(spa15,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_67_4=(pow(spab_5_67_4,3));
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow3_spaa_5_34_67_5=(pow(spaa_5_34_67_5,3));
const complex<T> pow3_spaa_5_234_67_5=(pow(spaa_5_234_67_5,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spa18*pow3_spaa_5_34_67_5)/(s1289*spa12*spa34*spa45*spa56*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow3_spaa_5_234_67_5)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spab_8_79_6*pow3_spa15)/(s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)+(pow2_spa18*pow3_spab_5_234_6)/(s1789*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_34_6)/(s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)+(pow2_spa18*pow3_spab_5_67_4)/(s4567*s567*spa12*spa23*spa56*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb56)
-
(pow2_spab_5_67_9*pow3_spa15)/(spa12*spa23*spa34*spa45*spa56*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2ppqb2mppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_567_3)/(s1289*s4567*spa12*spa45*spa56*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_7_456_3)
-
(pow2_spab_4_567_9*pow3_spa14)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pppqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_89_7*spaa_3_456_789_1))+(pow2_spa18*pow3_spab_3_189_2)/(s1289*s189*spa34*spa45*spa56*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2ppppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));

return(
(pow2_spaa_2_3456_79_8/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_7_89_1)+pow2_spab_2_18_9/(s189*spa23*spa34*spa45*spa56*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppq2pqb2mppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_5_46_3=-(spa45*spb43)+spa56*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_23_1789_5=-(spa24*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spa34*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_789_5=-(spa14*(spa57*spb71+spa58*spb81+spa59*spb91))-spa24*(spa57*spb72+spa58*spb82+spa59*spb92)-spa34*(spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14*spaa_4_123_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4*spaa_4_23_1789_5)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spab_8_679_5*pow3_spa14)/(s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)+(pow2_spa18*pow3_spab_4_567_3*spa57)/(s1289*s4567*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3*spab_5_46_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_567_9*pow3_spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppq2ppqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_346_2=-(spa35*spb32)-spa45*spb42+spa56*spb62;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_5=-(spa13*(spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa18*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa57)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_5_346_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppq2pppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_789_1=spa57*spb71+spa58*spb81+spa59*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));

return(
(pow2_spaa_2_345_679_8/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_5_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpppq2pqb2mpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_356_2=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_4=-(spa13*(spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_6789_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13*spaa_3_12_6789_4)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)+(pow2_spab_8_123_4*pow3_spa13)/(s1234*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)
-
(pow2_spa18*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa47)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_4_356_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_128_9*pow3_spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpppq2ppqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_789_1=spa47*spb71+spa48*spb81+spa49*spb91;
const complex<T> spab_4_235_1=-(spa24*spb21)-spa34*spb31+spa45*spb51;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-(pow2_spaa_2_34_5679_8/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_4_235_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_4_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppppq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_789_1=spa37*spb71+spa38*spb81+spa39*spb91;
const complex<T> spab_3_245_1=-(spa23*spb21)+spa34*spb41+spa35*spb51;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-((pow2_spaa_2_34_5679_8*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_3_245_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_3_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+pow2_spab_8_12_3/(s123*spa45*spa56*spa67*spa89*spab_4_23_1*spb23)+(pow2_spab_2_18_9*spa37)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pqb2mmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_46_5=spa47*spb54-spa67*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_346_5=spa37*spb53+spa47*spb54-spa67*spb65;
const complex<T> spab_7_189_5=spa17*spb51+spa78*spb85+spa79*spb95;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spab_8_79_6*spb15)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spab_5_67_9/(s567*spa56*spab_7_56_4*spb12*spb23*spb34*spb89)+(pow2_spbb_6_45_1238_9*spab_7_46_5)/(s456*s4567*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_346_5)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_5)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pmqb2mmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_356_4=spa37*spb43-spa57*spb54-spa67*spb64;
const complex<T> spab_7_189_4=spa17*spb41+spa78*spb84+spa79*spb94;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb14)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))+pow2_spbb_6_45_1238_9/(s456*s4567*spab_7_456_3*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_356_4)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_4)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pmmqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_189_3=spa17*spb31+spa78*spb83+spa79*spb93;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb13)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_345_128_9/(s1289*s3456*spab_7_189_2*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_3)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpq2pmmmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-(pow2_spab_8_79_6/(s789*spa89*spab_7_89_1*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_2345_18_9/(s1789*s189*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmq2pqb2mmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_4=-(spb65*(spa16*spb41+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_4=-(spb65*(spa16*spb41+spa26*spb42+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_235_4=spa26*spb42+spa36*spb43-spa56*spb54;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa45*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spb79*pow3_spab_6_789_5*spb14)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_235_4)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spab_4_567_9*pow3_spb57)/(s4567*spab_4_56_7*spb12*spb23*spb56*spb67*spb89*spbb_5_67_1289_3)+(pow2_spbb_5_34_128_9*pow3_spb57*spbb_5_67_1289_4)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmq2pmqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_3=-(spb65*(spa16*spb31+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_245_3=spa26*spb32-spa46*spb43-spa56*spb53;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5*spb13)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_245_3)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5)/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spbb_5_34_128_9*pow3_spb57)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmq2pmmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5)/(s6789*s789*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5)/(s1789*s2345*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spab_8_679_5*pow3_spb57)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmq2pqb2mmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_3=-(spb54*(spa15*spb31+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa68*spb83+spa69*spb93)-spb74*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))+(pow2_spb79*pow3_spab_5_123_4*spb13)/(s1234*s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_3_128_9*pow3_spb47)/(s1289*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow3_spb47*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47*spbb_4_567_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)
-
(pow2_spb79*pow3_spab_3_56_4)/(s3456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow3_spbb_4_56_123_4*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4*spbb_4_56_1789_3)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmq2pmqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));

return(
((pow2_spb79*pow3_spab_5_123_4)/(s1234*s6789*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4)/(s234*s2345*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_8_123_4*pow3_spb47)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_89_1)+(pow2_spb79*pow3_spbb_4_56_123_4)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmmq2pqb2mqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
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
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb37=(pow(spb37,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_456_12_3=(pow(spbb_3_456_12_3,3));
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))+(pow2_spb79*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow3_spb37)/(s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow3_spb37)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)+(pow2_spb79*pow3_spab_2_456_3)/(s1789*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow3_spbb_3_456_12_3)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)
-
(pow2_spb79*pow3_spab_2_45_3)/(s2345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow3_spbb_3_45_12_3)/(s6789*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mq2ppppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_5_46_7=spa45*spb74-spa56*spb76;
const complex<T> spab_5_346_7=spa35*spb73+spa45*spb74-spa56*spb76;
const complex<T> spab_5_189_7=spa15*spb71+spa58*spb87+spa59*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7*spab_5_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_5_346_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)
-
(pow2_spa18*pow2_spab_6_45_7*spab_5_46_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7)+(pow2_spa18*pow2_spb57)/(s567*spa12*spa23*spa34*spa89*spab_4_56_7*spb56)+(pow2_spa16*pow2_spb79*spa15)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mpq2pppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_4_356_7=spa34*spb73-spa45*spb75-spa46*spb76;
const complex<T> spab_4_189_7=spa14*spb71+spa48*spb87+spa49*spb97;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
(-((pow2_spa18*pow2_spab_6_45_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7))+(pow2_spa18*pow2_spab_6_189_7*spab_4_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_4_356_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)+(pow2_spa16*pow2_spb79*spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mppq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_3_189_7=spa13*spb71+spa38*spb87+spa39*spb97;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_345_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7)+(pow2_spa18*pow2_spab_6_189_7*spab_3_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa16*pow2_spb79*spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mpppq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7)+(pow2_spa16*pow2_spb79)/(s789*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppqb2mq2pppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_235_6=spa24*spb62+spa34*spb63-spa45*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_4=spa56*(spa14*spb61+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_4=spa56*(spa14*spb61+spa24*spb62+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spab_5_67_4=(pow(spab_5_67_4,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57*spaa_5_67_1289_4)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_4)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spa18*pow3_spab_5_234_6*spab_4_235_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)+(pow2_spa18*pow3_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa14*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa18*pow2_spab_5_67_4*spa57)/(s4567*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb45)
-
(pow2_spa15*pow2_spab_5_67_9*spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppqb2mpq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_245_6=spa23*spb62-spa34*spb64-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_3=spa56*(spa13*spb61+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)+(pow2_spa18*pow3_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_234_6*spab_3_245_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa13*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppqb2mppq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_234_67_5*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_89_1))
-
(pow2_spa18*pow3_spab_5_234_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_8_79_6*spab_5_789_6)/(s6789*s789*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa57)/(spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpppqb2mq2ppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa34=(pow(spa34,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_3=spa45*(spa13*spb51+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa38*spb86+spa39*spb96)+spa47*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_56_1789_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spab_4_567_3=(pow(spab_4_567_3,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47*spaa_4_567_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spa13*spaa_4_123_56_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4*spaa_4_56_1789_3)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)
-
(pow2_spa14*pow2_spab_8_679_5*spa13*spab_4_123_5)/(s1234*s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow2_spab_4_567_3*spa47)/(s1289*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)+(pow2_spa18*pow3_spab_4_56_3)/(s3456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_567_9*spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpppqb2mpq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
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
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spaa_4_123_56_4)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_8_679_5*spab_4_123_5)/(s1234*s6789*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_4_567_9*spa47)/(spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qppppqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_456_3=-(spa13*(spa34*spb41+spa35*spb51+spa36*spb61))-spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_189_2=(pow(spab_3_189_2,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
((pow2_spa13*pow2_spaa_3_456_79_8*spaa_3_12_456_3)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)+(pow2_spa13*pow2_spaa_3_45_679_8*spaa_3_12_45_3)/(s6789*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa13*pow2_spab_8_123_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)+(pow2_spa18*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow2_spab_3_189_2*spa37)/(s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)
-
(pow2_spa18*pow3_spab_3_456_2)/(s1789*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_128_9*spa37)/(spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mq2pmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_34_5=-((spa36*spb53+spa46*spb54)*spb65)-(spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_234_5=-((spa26*spb52+spa36*spb53+spa46*spb54)*spb65)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_6_789_5=(pow(spab_6_789_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spab_6_789_5*pow2_spb79*spb15)/(s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)
-
(pow2_spb79*pow3_spab_6_234_5)/(s1789*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spb79*pow3_spab_6_34_5)/(s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)
-
(pow2_spab_4_567_9*pow2_spb57*spab_4_67_5)/(s4567*s567*spab_4_56_7*spb12*spb23*spb56*spb89*spbb_5_67_1289_3)
-
(pow2_spb57*pow2_spbb_5_34_128_9*spbb_5_67_34_5)/(s1289*spb12*spb34*spb45*spb56*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow2_spb57*spb15)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_5_1234_89_7*spbb_5_67_89_1)
-
(pow2_spb57*pow2_spbb_5_234_18_9*spbb_5_67_234_5)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mmq2pmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spab_3_128_9*pow2_spb47*spab_3_567_4)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_567_189_2))+(pow2_spab_8_123_4*pow2_spb47*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4)/(s3456*s456*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mmmq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3)/(s1789*s3456*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_89_7*spbb_3_456_789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpqb2mmmmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));

return(
(-((pow2_spa18*pow2_spb27)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56))
-
(pow2_spab_1_789_2*pow2_spb79)/(s1789*s789*spab_1_89_7*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmqb2mq2pmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_23_1789_5=spb42*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_789_5=spb41*(spa17*spb75+spa18*spb85+spa19*spb95)+spb42*(spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_46_5=-(spa34*spb54)+spa36*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))
-
(pow2_spab_5_123_4*pow2_spb79*spb14)/(s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)+(pow2_spb79*pow3_spab_5_23_4)/(s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)
-
(pow2_spab_3_128_9*pow2_spb47*spab_3_567_4*spb57)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow2_spb47*spb14*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spb57*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4*spab_3_46_5)/(s3456*s456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14*spbb_4_123_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4*spbb_4_23_1789_5)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmqb2mmq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_5=spb31*(spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_346_5=-(spa23*spb53)-spa24*spb54+spa26*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb57)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_346_5)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3)/(s2345*s345*spab_2_345_6*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmqb2mmmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_5=spa17*spb75+spa18*spb85+spa19*spb95;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb57)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_5)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79)/(s2345*s6789*spab_1_789_6*spb23*spb34*spb45*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmqb2mq2pmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
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
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_4=spb31*(spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spbb_3_12_6789_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_356_4=-(spa23*spb43)+spa25*spb54+spa26*spb64;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))
-
(pow2_spab_4_12_3*pow2_spb79*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb47)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow2_spb37*spb13*spb47)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_356_4)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_4)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13*spbb_3_12_6789_4)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmqb2mmq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_4=spa17*spb74+spa18*spb84+spa19*spb94;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_235_4=-(spa12*spb42)-spa13*spb43+spa15*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb47)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_4)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_235_4)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79)/(s1234*s234*spab_1_234_5*spb23*spb34*spb56*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q2Q3g2l_qpmmmqb2mq2pqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_3=spa17*spb73+spa18*spb83+spa19*spb93;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_245_3=-(spa12*spb32)+spa14*spb43+spa15*spb53;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb37)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_3)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_245_3)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79*spab_1_24_3)/(s1234*s234*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb89)
-
(pow2_spa13*pow2_spb79)/(s123*spa23*spab_1_23_4*spb45*spb56*spb67*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pqb2mpppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_456_3=-(spa13*(spa34*spb41+spa35*spb51+spa36*spb61))-spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_189_2=(pow(spab_3_189_2,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
((pow2_spa13*pow2_spaa_3_456_79_8*spaa_3_12_456_3)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)+(pow2_spa13*pow2_spaa_3_45_679_8*spaa_3_12_45_3)/(s6789*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa13*pow2_spab_8_123_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)+(pow2_spa18*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow2_spab_3_189_2*spa37)/(s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)
-
(pow2_spa18*pow3_spab_3_456_2)/(s1789*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_128_9*spa37)/(spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2ppqb2mppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
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
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spaa_4_123_56_4)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_8_679_5*spab_4_123_5)/(s1234*s6789*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_4_567_9*spa47)/(spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pppqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_234_67_5*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_89_1))
-
(pow2_spa18*pow3_spab_5_234_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_8_79_6*spab_5_789_6)/(s6789*s789*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa57)/(spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2ppppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7)+(pow2_spa16*pow2_spb79)/(s789*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpq2pqb2mppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa34=(pow(spa34,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_3=spa45*(spa13*spb51+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa38*spb86+spa39*spb96)+spa47*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_56_1789_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spab_4_567_3=(pow(spab_4_567_3,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47*spaa_4_567_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spa13*spaa_4_123_56_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4*spaa_4_56_1789_3)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)
-
(pow2_spa14*pow2_spab_8_679_5*spa13*spab_4_123_5)/(s1234*s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow2_spab_4_567_3*spa47)/(s1289*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)+(pow2_spa18*pow3_spab_4_56_3)/(s3456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_567_9*spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpq2ppqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_245_6=spa23*spb62-spa34*spb64-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_3=spa56*(spa13*spb61+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)+(pow2_spa18*pow3_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_234_6*spab_3_245_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa13*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpq2pppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_3_189_7=spa13*spb71+spa38*spb87+spa39*spb97;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_345_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7)+(pow2_spa18*pow2_spab_6_189_7*spab_3_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa16*pow2_spb79*spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmppq2pqb2mpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_235_6=spa24*spb62+spa34*spb63-spa45*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_4=spa56*(spa14*spb61+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_4=spa56*(spa14*spb61+spa24*spb62+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spab_5_67_4=(pow(spab_5_67_4,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57*spaa_5_67_1289_4)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_4)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spa18*pow3_spab_5_234_6*spab_4_235_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)+(pow2_spa18*pow3_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa14*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa18*pow2_spab_5_67_4*spa57)/(s4567*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb45)
-
(pow2_spa15*pow2_spab_5_67_9*spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmppq2ppqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_4_356_7=spa34*spb73-spa45*spb75-spa46*spb76;
const complex<T> spab_4_189_7=spa14*spb71+spa48*spb87+spa49*spb97;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
(-((pow2_spa18*pow2_spab_6_45_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7))+(pow2_spa18*pow2_spab_6_189_7*spab_4_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_4_356_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)+(pow2_spa16*pow2_spb79*spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpppq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_5_46_7=spa45*spb74-spa56*spb76;
const complex<T> spab_5_346_7=spa35*spb73+spa45*spb74-spa56*spb76;
const complex<T> spab_5_189_7=spa15*spb71+spa58*spb87+spa59*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7*spab_5_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_5_346_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)
-
(pow2_spa18*pow2_spab_6_45_7*spab_5_46_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7)+(pow2_spa18*pow2_spb57)/(s567*spa12*spa23*spa34*spa89*spab_4_56_7*spb56)+(pow2_spa16*pow2_spb79*spa15)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pqb2mmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_3=spa17*spb73+spa18*spb83+spa19*spb93;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_245_3=-(spa12*spb32)+spa14*spb43+spa15*spb53;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb37)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_3)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_245_3)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79*spab_1_24_3)/(s1234*s234*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb89)
-
(pow2_spa13*pow2_spb79)/(s123*spa23*spab_1_23_4*spb45*spb56*spb67*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pmqb2mmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_4=spa17*spb74+spa18*spb84+spa19*spb94;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_235_4=-(spa12*spb42)-spa13*spb43+spa15*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb47)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_4)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_235_4)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79)/(s1234*s234*spab_1_234_5*spb23*spb34*spb56*spb67*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pmmqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_5=spa17*spb75+spa18*spb85+spa19*spb95;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb57)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_5)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79)/(s2345*s6789*spab_1_789_6*spb23*spb34*spb45*spb67*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmq2pmmmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));

return(
(-((pow2_spa18*pow2_spb27)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56))
-
(pow2_spab_1_789_2*pow2_spb79)/(s1789*s789*spab_1_89_7*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmq2pqb2mmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
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
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_4=spb31*(spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spbb_3_12_6789_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_356_4=-(spa23*spb43)+spa25*spb54+spa26*spb64;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))
-
(pow2_spab_4_12_3*pow2_spb79*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb47)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow2_spb37*spb13*spb47)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_356_4)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_4)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13*spbb_3_12_6789_4)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmq2pmqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_5=spb31*(spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_346_5=-(spa23*spb53)-spa24*spb54+spa26*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb57)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_346_5)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3)/(s2345*s345*spab_2_345_6*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmq2pmmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3)/(s1789*s3456*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_89_7*spbb_3_456_789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmq2pqb2mmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_23_1789_5=spb42*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_789_5=spb41*(spa17*spb75+spa18*spb85+spa19*spb95)+spb42*(spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_46_5=-(spa34*spb54)+spa36*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))
-
(pow2_spab_5_123_4*pow2_spb79*spb14)/(s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)+(pow2_spb79*pow3_spab_5_23_4)/(s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)
-
(pow2_spab_3_128_9*pow2_spb47*spab_3_567_4*spb57)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow2_spb47*spb14*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spb57*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4*spab_3_46_5)/(s3456*s456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14*spbb_4_123_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4*spbb_4_23_1789_5)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmq2pmqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spab_3_128_9*pow2_spb47*spab_3_567_4)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_567_189_2))+(pow2_spab_8_123_4*pow2_spb47*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4)/(s3456*s456*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmmq2pqb2mqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_34_5=-((spa36*spb53+spa46*spb54)*spb65)-(spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_234_5=-((spa26*spb52+spa36*spb53+spa46*spb54)*spb65)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_6_789_5=(pow(spab_6_789_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spab_6_789_5*pow2_spb79*spb15)/(s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)
-
(pow2_spb79*pow3_spab_6_234_5)/(s1789*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spb79*pow3_spab_6_34_5)/(s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)
-
(pow2_spab_4_567_9*pow2_spb57*spab_4_67_5)/(s4567*s567*spab_4_56_7*spb12*spb23*spb56*spb89*spbb_5_67_1289_3)
-
(pow2_spb57*pow2_spbb_5_34_128_9*spbb_5_67_34_5)/(s1289*spb12*spb34*spb45*spb56*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow2_spb57*spb15)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_5_1234_89_7*spbb_5_67_89_1)
-
(pow2_spb57*pow2_spbb_5_234_18_9*spbb_5_67_234_5)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mq2ppppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_789_1=spa37*spb71+spa38*spb81+spa39*spb91;
const complex<T> spab_3_245_1=-(spa23*spb21)+spa34*spb41+spa35*spb51;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-((pow2_spaa_2_34_5679_8*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_3_245_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_3_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+pow2_spab_8_12_3/(s123*spa45*spa56*spa67*spa89*spab_4_23_1*spb23)+(pow2_spab_2_18_9*spa37)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mpq2pppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_789_1=spa47*spb71+spa48*spb81+spa49*spb91;
const complex<T> spab_4_235_1=-(spa24*spb21)-spa34*spb31+spa45*spb51;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-(pow2_spaa_2_34_5679_8/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_4_235_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_4_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mppq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_789_1=spa57*spb71+spa58*spb81+spa59*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));

return(
(pow2_spaa_2_345_679_8/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_5_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mpppq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));

return(
(pow2_spaa_2_3456_79_8/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_7_89_1)+pow2_spab_2_18_9/(s189*spa23*spa34*spa45*spa56*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpqb2mq2pppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_356_2=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_4=-(spa13*(spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_6789_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13*spaa_3_12_6789_4)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)+(pow2_spab_8_123_4*pow3_spa13)/(s1234*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)
-
(pow2_spa18*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa47)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_4_356_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_128_9*pow3_spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpqb2mpq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_346_2=-(spa35*spb32)-spa45*spb42+spa56*spb62;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_5=-(spa13*(spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa18*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa57)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_5_346_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpqb2mppq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_89_7*spaa_3_456_789_1))+(pow2_spa18*pow3_spab_3_189_2)/(s1289*s189*spa34*spa45*spa56*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmppqb2mq2ppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_5_46_3=-(spa45*spb43)+spa56*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_23_1789_5=-(spa24*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spa34*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_789_5=-(spa14*(spa57*spb71+spa58*spb81+spa59*spb91))-spa24*(spa57*spb72+spa58*spb82+spa59*spb92)-spa34*(spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14*spaa_4_123_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4*spaa_4_23_1789_5)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spab_8_679_5*pow3_spa14)/(s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)+(pow2_spa18*pow3_spab_4_567_3*spa57)/(s1289*s4567*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3*spab_5_46_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_567_9*pow3_spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmppqb2mpq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_567_3)/(s1289*s4567*spa12*spa45*spa56*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_7_456_3)
-
(pow2_spab_4_567_9*pow3_spa14)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmpppqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa15=(pow(spa15,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_67_4=(pow(spab_5_67_4,3));
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow3_spaa_5_34_67_5=(pow(spaa_5_34_67_5,3));
const complex<T> pow3_spaa_5_234_67_5=(pow(spaa_5_234_67_5,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spa18*pow3_spaa_5_34_67_5)/(s1289*spa12*spa34*spa45*spa56*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow3_spaa_5_234_67_5)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spab_8_79_6*pow3_spa15)/(s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)+(pow2_spa18*pow3_spab_5_234_6)/(s1789*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_34_6)/(s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)+(pow2_spa18*pow3_spab_5_67_4)/(s4567*s567*spa12*spa23*spa56*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb56)
-
(pow2_spab_5_67_9*pow3_spa15)/(spa12*spa23*spa34*spa45*spa56*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mq2pmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
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
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb37=(pow(spb37,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_456_12_3=(pow(spbb_3_456_12_3,3));
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))+(pow2_spb79*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow3_spb37)/(s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow3_spb37)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)+(pow2_spb79*pow3_spab_2_456_3)/(s1789*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow3_spbb_3_456_12_3)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)
-
(pow2_spb79*pow3_spab_2_45_3)/(s2345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow3_spbb_3_45_12_3)/(s6789*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mmq2pmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));

return(
((pow2_spb79*pow3_spab_5_123_4)/(s1234*s6789*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4)/(s234*s2345*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_8_123_4*pow3_spb47)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_89_1)+(pow2_spb79*pow3_spbb_4_56_123_4)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mmmq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5)/(s6789*s789*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5)/(s1789*s2345*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spab_8_679_5*pow3_spb57)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmqb2mmmmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-(pow2_spab_8_79_6/(s789*spa89*spab_7_89_1*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_2345_18_9/(s1789*s189*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmqb2mq2pmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_3=-(spb54*(spa15*spb31+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa68*spb83+spa69*spb93)-spb74*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))+(pow2_spb79*pow3_spab_5_123_4*spb13)/(s1234*s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_3_128_9*pow3_spb47)/(s1289*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow3_spb47*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47*spbb_4_567_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)
-
(pow2_spb79*pow3_spab_3_56_4)/(s3456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow3_spbb_4_56_123_4*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4*spbb_4_56_1789_3)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmqb2mmq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_3=-(spb65*(spa16*spb31+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_245_3=spa26*spb32-spa46*spb43-spa56*spb53;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5*spb13)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_245_3)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5)/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spbb_5_34_128_9*pow3_spb57)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmqb2mmmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_189_3=spa17*spb31+spa78*spb83+spa79*spb93;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb13)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_345_128_9/(s1289*s3456*spab_7_189_2*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_3)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmqb2mq2pmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_4=-(spb65*(spa16*spb41+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_4=-(spb65*(spa16*spb41+spa26*spb42+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_235_4=spa26*spb42+spa36*spb43-spa56*spb54;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa45*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spb79*pow3_spab_6_789_5*spb14)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_235_4)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spab_4_567_9*pow3_spb57)/(s4567*spab_4_56_7*spb12*spb23*spb56*spb67*spb89*spbb_5_67_1289_3)+(pow2_spbb_5_34_128_9*pow3_spb57*spbb_5_67_1289_4)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmqb2mmq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_356_4=spa37*spb43-spa57*spb54-spa67*spb64;
const complex<T> spab_7_189_4=spa17*spb41+spa78*spb84+spa79*spb94;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb14)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))+pow2_spbb_6_45_1238_9/(s456*s4567*spab_7_456_3*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_356_4)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_4)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qmmmmqb2mq2pqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_46_5=spa47*spb54-spa67*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_346_5=spa37*spb53+spa47*spb54-spa67*spb65;
const complex<T> spab_7_189_5=spa17*spb51+spa78*spb85+spa79*spb95;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spab_8_79_6*spb15)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spab_5_67_9/(s567*spa56*spab_7_56_4*spb12*spb23*spb34*spb89)+(pow2_spbb_6_45_1238_9*spab_7_46_5)/(s456*s4567*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_346_5)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_5)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pqb2mpppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa15=(pow(spa15,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_67_4=(pow(spab_5_67_4,3));
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow3_spaa_5_34_67_5=(pow(spaa_5_34_67_5,3));
const complex<T> pow3_spaa_5_234_67_5=(pow(spaa_5_234_67_5,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spa18*pow3_spaa_5_34_67_5)/(s1289*spa12*spa34*spa45*spa56*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow3_spaa_5_234_67_5)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spab_8_79_6*pow3_spa15)/(s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)+(pow2_spa18*pow3_spab_5_234_6)/(s1789*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_34_6)/(s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)+(pow2_spa18*pow3_spab_5_67_4)/(s4567*s567*spa12*spa23*spa56*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb56)
-
(pow2_spab_5_67_9*pow3_spa15)/(spa12*spa23*spa34*spa45*spa56*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2ppqb2mppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4)/(s189*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_567_3)/(s1289*s4567*spa12*spa45*spa56*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_7_456_3)
-
(pow2_spab_4_567_9*pow3_spa14)/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pppqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_89_7*spaa_3_456_789_1))+(pow2_spa18*pow3_spab_3_189_2)/(s1289*s189*spa34*spa45*spa56*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13)/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2ppppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));

return(
(pow2_spaa_2_3456_79_8/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_7_89_1)+pow2_spab_2_18_9/(s189*spa23*spa34*spa45*spa56*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppq2pqb2mppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa14=(pow(spa14,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_5_46_3=-(spa45*spb43)+spa56*spb63;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_23_1789_5=-(spa24*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spa34*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_789_5=-(spa14*(spa57*spb71+spa58*spb81+spa59*spb91))-spa24*(spa57*spb72+spa58*spb82+spa59*spb92)-spa34*(spa57*spb73+spa58*spb83+spa59*spb93);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_4_567_3=(pow(spab_4_567_3,3));
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_23_567_4=(pow(spaa_4_23_567_4,3));
const complex<T> pow3_spaa_4_23_56_4=(pow(spaa_4_23_56_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));

return(
(-((pow2_spa18*pow3_spaa_4_23_567_4*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))
-
(pow2_spaa_4_56_79_8*pow3_spa14*spaa_4_123_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_23_56_4*spaa_4_23_1789_5)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spab_8_679_5*pow3_spa14)/(s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow3_spab_4_23_5)/(s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)+(pow2_spa18*pow3_spab_4_567_3*spa57)/(s1289*s4567*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)
-
(pow2_spa18*pow3_spab_4_56_3*spab_5_46_3)/(s3456*s456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb45)
-
(pow2_spab_4_567_9*pow3_spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppq2ppqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_346_2=-(spa35*spb32)-spa45*spb42+spa56*spb62;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_5=-(spa13*(spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_5)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa18*pow3_spab_3_45_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa57)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_5_346_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)
-
(pow2_spab_3_128_9*pow3_spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppq2pppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_789_1=spa57*spb71+spa58*spb81+spa59*spb91;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));

return(
(pow2_spaa_2_345_679_8/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_5_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpppq2pqb2mpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_4_356_2=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spab_4_35_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_789_4=-(spa13*(spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_6789_4=-(spa13*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spa23*(spa46*spb62+spa47*spb72+spa48*spb82+spa49*spb92);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow3_spab_3_189_2=(pow(spab_3_189_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
(-((pow2_spaa_3_456_79_8*pow3_spa13*spaa_3_12_789_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1))
-
(pow2_spaa_3_45_679_8*pow3_spa13*spaa_3_12_6789_4)/(s6789*spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)+(pow2_spab_8_123_4*pow3_spa13)/(s1234*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)
-
(pow2_spa18*pow3_spab_3_45_2*spab_4_35_2)/(s2345*s345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow3_spab_3_189_2*spa47)/(s1289*s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)+(pow2_spa18*pow3_spab_3_456_2*spab_4_356_2)/(s1789*s3456*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb34)
-
(pow2_spab_3_128_9*pow3_spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpppq2ppqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_789_1=spa47*spb71+spa48*spb81+spa49*spb91;
const complex<T> spab_4_235_1=-(spa24*spb21)-spa34*spb31+spa45*spb51;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-(pow2_spaa_2_34_5679_8/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_4_235_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_4_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+(pow2_spab_2_18_9*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppppq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_3_789_1=spa37*spb71+spa38*spb81+spa39*spb91;
const complex<T> spab_3_245_1=-(spa23*spb21)+spa34*spb41+spa35*spb51;
const complex<T> spab_3_24_1=-(spa23*spb21)+spa34*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*spb63)-spa78*spb73+spa89*spb93)+spa24*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa25*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_5679_8=spa23*(-(spa58*spb53)-spa68*spb63-spa78*spb73+spa89*spb93)+spa24*(-(spa58*spb54)-spa68*spb64-spa78*spb74+spa89*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));
const complex<T> pow2_spaa_2_3456_79_8=(pow(spaa_2_3456_79_8,2));
const complex<T> pow2_spaa_2_345_679_8=(pow(spaa_2_345_679_8,2));
const complex<T> pow2_spaa_2_34_5679_8=(pow(spaa_2_34_5679_8,2));

return(
(-((pow2_spaa_2_34_5679_8*spab_3_24_1)/(s1234*s234*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1))+(pow2_spaa_2_345_679_8*spab_3_245_1)/(s2345*s6789*spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1)+(pow2_spaa_2_3456_79_8*spab_3_789_1)/(s1789*s789*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1)+pow2_spab_8_12_3/(s123*spa45*spa56*spa67*spa89*spab_4_23_1*spb23)+(pow2_spab_2_18_9*spa37)/(s189*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pqb2mmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_46_5=spa47*spb54-spa67*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_346_5=spa37*spb53+spa47*spb54-spa67*spb65;
const complex<T> spab_7_189_5=spa17*spb51+spa78*spb85+spa79*spb95;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));

return(
(-((pow2_spab_8_79_6*spb15)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spab_5_67_9/(s567*spa56*spab_7_56_4*spb12*spb23*spb34*spb89)+(pow2_spbb_6_45_1238_9*spab_7_46_5)/(s456*s4567*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_346_5)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_5)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pmqb2mmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_1238_9=spb64*(spa14*spb91+spa24*spb92+spa34*spb93-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92+spa35*spb93-spa58*spb98);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_356_4=spa37*spb43-spa57*spb54-spa67*spb64;
const complex<T> spab_7_189_4=spa17*spb41+spa78*spb84+spa79*spb94;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_45_1238_9=(pow(spbb_6_45_1238_9,2));
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb14)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))+pow2_spbb_6_45_1238_9/(s456*s4567*spab_7_456_3*spb12*spb23*spb45*spb56*spb89)
-
(pow2_spbb_6_345_128_9*spab_7_356_4)/(s1289*s3456*spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_4)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pmmqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+spa23*spb92-spa38*spb98)+spb64*(spa14*spb91+spa24*spb92-spa48*spb98)+spb65*(spa15*spb91+spa25*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_189_3=spa17*spb31+spa78*spb83+spa79*spb93;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spbb_6_345_128_9=(pow(spbb_6_345_128_9,2));
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-((pow2_spab_8_79_6*spb13)/(s789*spa89*spab_7_89_1*spb12*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_345_128_9/(s1289*s3456*spab_7_189_2*spb12*spb34*spb45*spb56*spb89)
-
(pow2_spbb_6_2345_18_9*spab_7_189_3)/(s1789*s189*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpq2pmmmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spbb_6_2345_18_9=(pow(spbb_6_2345_18_9,2));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));

return(
(-(pow2_spab_8_79_6/(s789*spa89*spab_7_89_1*spb23*spb34*spb45*spb56))
-
pow2_spbb_6_2345_18_9/(s1789*s189*spab_7_89_1*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmq2pqb2mmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_4=-(spb65*(spa16*spb41+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_4=-(spb65*(spa16*spb41+spa26*spb42+spa68*spb84+spa69*spb94))-spb75*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_235_4=spa26*spb42+spa36*spb43-spa56*spb54;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa45*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spb79*pow3_spab_6_789_5*spb14)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_235_4)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5*spab_6_35_4)/(s345*s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spab_4_567_9*pow3_spb57)/(s4567*spab_4_56_7*spb12*spb23*spb56*spb67*spb89*spbb_5_67_1289_3)+(pow2_spbb_5_34_128_9*pow3_spb57*spbb_5_67_1289_4)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmq2pmqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_3=-(spb65*(spa16*spb31+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_245_3=spa26*spb32-spa46*spb43-spa56*spb53;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5*spb13)/(s6789*s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5*spab_6_245_3)/(s1789*s2345*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)
-
(pow2_spb79*pow3_spab_6_34_5)/(s345*s3456*spab_6_345_2*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)+(pow2_spbb_5_34_128_9*pow3_spb57)/(s1289*spb12*spb34*spb45*spb56*spb67*spb89*spbb_5_34_1289_7*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow3_spb57*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57*spbb_5_67_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmq2pmmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb57=(pow(spb57,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_6_789_5=(pow(spab_6_789_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));

return(
((pow2_spb79*pow3_spab_6_789_5)/(s6789*s789*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)+(pow2_spb79*pow3_spab_6_234_5)/(s1789*s2345*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spab_8_679_5*pow3_spb57)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1)+(pow2_spbb_5_234_18_9*pow3_spb57)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_5_234_189_7*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmq2pqb2mmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_3=-(spb54*(spa15*spb31+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa68*spb83+spa69*spb93)-spb74*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_3=-(spb54*(spa15*spb31+spa57*spb73+spa58*spb83+spa59*spb93))-spb64*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa34*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))+(pow2_spb79*pow3_spab_5_123_4*spb13)/(s1234*s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4*spab_5_24_3)/(s234*s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_3_128_9*pow3_spb47)/(s1289*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow3_spb47*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47*spbb_4_567_189_3)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)
-
(pow2_spb79*pow3_spab_3_56_4)/(s3456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow3_spbb_4_56_123_4*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4*spbb_4_56_1789_3)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmq2pmqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb47=(pow(spb47,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_4_56_23_4=(pow(spbb_4_56_23_4,3));
const complex<T> pow3_spbb_4_56_123_4=(pow(spbb_4_56_123_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_5_123_4=(pow(spab_5_123_4,3));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));

return(
((pow2_spb79*pow3_spab_5_123_4)/(s1234*s6789*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_123_789_6)
-
(pow2_spb79*pow3_spab_5_23_4)/(s234*s2345*spab_5_234_1*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)+(pow2_spab_8_123_4*pow3_spb47)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)+(pow2_spbb_4_23_18_9*pow3_spb47)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_89_1)+(pow2_spb79*pow3_spbb_4_56_123_4)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_56_23_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmmq2pqb2mqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
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
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb37=(pow(spb37,3));
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spbb_3_456_12_3=(pow(spbb_3_456_12_3,3));
const complex<T> pow3_spbb_3_45_12_3=(pow(spbb_3_45_12_3,3));
const complex<T> pow3_spab_4_12_3=(pow(spab_4_12_3,3));
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa23*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))+(pow2_spb79*pow3_spab_4_12_3)/(s123*s1234*spab_4_23_1*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow3_spb37)/(s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow3_spb37)/(spa89*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)+(pow2_spb79*pow3_spab_2_456_3)/(s1789*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow3_spbb_3_456_12_3)/(s789*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)
-
(pow2_spb79*pow3_spab_2_45_3)/(s2345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow3_spbb_3_45_12_3)/(s6789*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mq2ppppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_5_46_7=spa45*spb74-spa56*spb76;
const complex<T> spab_5_346_7=spa35*spb73+spa45*spb74-spa56*spb76;
const complex<T> spab_5_189_7=spa15*spb71+spa58*spb87+spa59*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7*spab_5_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_5_346_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)
-
(pow2_spa18*pow2_spab_6_45_7*spab_5_46_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7)+(pow2_spa18*pow2_spb57)/(s567*spa12*spa23*spa34*spa89*spab_4_56_7*spb56)+(pow2_spa16*pow2_spb79*spa15)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mpq2pppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_4_356_7=spa34*spb73-spa45*spb75-spa46*spb76;
const complex<T> spab_4_189_7=spa14*spb71+spa48*spb87+spa49*spb97;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_45_7=(pow(spab_6_45_7,2));
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
(-((pow2_spa18*pow2_spab_6_45_7)/(s456*s4567*spa12*spa23*spa45*spa56*spa89*spab_3_456_7))+(pow2_spa18*pow2_spab_6_189_7*spab_4_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa18*pow2_spab_6_345_7*spab_4_356_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7)+(pow2_spa16*pow2_spb79*spa14)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mppq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_345_7=spa36*spb73+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_3_189_7=spa13*spb71+spa38*spb87+spa39*spb97;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow2_spab_6_345_7=(pow(spab_6_345_7,2));
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_345_7)/(s1289*s3456*spa12*spa34*spa45*spa56*spa89*spab_2_189_7)+(pow2_spa18*pow2_spab_6_189_7*spab_3_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7)+(pow2_spa16*pow2_spb79*spa13)/(s789*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mpppq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_189_7=spa16*spb71+spa68*spb87+spa69*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_6_189_7=(pow(spab_6_189_7,2));

return(
((pow2_spa18*pow2_spab_6_189_7)/(s1789*s189*spa23*spa34*spa45*spa56*spa89*spab_1_89_7)+(pow2_spa16*pow2_spb79)/(s789*spa23*spa34*spa45*spa56*spab_1_89_7*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppqb2mq2pppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_67_4=spa56*spb64+spa57*spb74;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_4_35_6=spa34*spb63-spa45*spb65;
const complex<T> spab_4_235_6=spa24*spb62+spa34*spb63-spa45*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_4=spa56*(spa14*spb61+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_67_1289_4=spa56*(spa14*spb61+spa24*spb62+spa48*spb86+spa49*spb96)+spa57*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spaa_5_67_1289_3=spa56*(spa13*spb61+spa23*spb62+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spab_5_67_4=(pow(spab_5_67_4,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57*spaa_5_67_1289_4)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_4)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
-
(pow2_spa18*pow3_spab_5_234_6*spab_4_235_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)+(pow2_spa18*pow3_spab_5_34_6*spab_4_35_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa14*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa18*pow2_spab_5_67_4*spa57)/(s4567*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4)+(pow2_spa18*pow3_spb46)/(s456*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*spb45)
-
(pow2_spa15*pow2_spab_5_67_9*spa14*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppqb2mpq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
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
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_3_245_6=spa23*spb62-spa34*spb64-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_67_189_3=spa56*(spa13*spb61+spa38*spb86+spa39*spb96)+spa57*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_5_67_189_2=spa56*(spa12*spb61+spa28*spb86+spa29*spb96)+spa57*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_5_34_67_5=-(spa35*(spa56*spb63+spa57*spb73))-spa45*(spa56*spb64+spa57*spb74);
const complex<T> spaa_5_34_1289_7=-(spa35*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93))-spa45*(spa17*spb41+spa27*spb42+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_5_34_6=(pow(spab_5_34_6,3));
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_34_67_5=(pow(spaa_5_34_67_5,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_34_67_5*spa57)/(s1289*spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_189_2))
-
(pow2_spa18*pow2_spaa_5_234_67_5*spa57*spaa_5_67_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)+(pow2_spa18*pow3_spab_5_34_6)/(s345*s3456*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6)
-
(pow2_spa18*pow3_spab_5_234_6*spab_3_245_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6)
-
(pow2_spa15*pow2_spab_8_79_6*spa13*spab_5_789_6)/(s6789*s789*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa13*spa57)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppqb2mppq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa15=(pow(spa15,2));
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_5_789_6=spa57*spb76+spa58*spb86+spa59*spb96;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> spab_5_234_6=spa25*spb62+spa35*spb63+spa45*spb64;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spaa_5_67_89_1=spa56*(spa18*spb86+spa19*spb96)+spa57*(spa18*spb87+spa19*spb97);
const complex<T> spaa_5_234_67_5=-(spa56*(spa25*spb62+spa35*spb63+spa45*spb64))-spa57*(spa25*spb72+spa35*spb73+spa45*spb74);
const complex<T> spaa_5_234_189_7=-(spa25*(spa17*spb21+spa78*spb82+spa79*spb92))-spa35*(spa17*spb31+spa78*spb83+spa79*spb93)-spa45*(spa17*spb41+spa78*spb84+spa79*spb94);
const complex<T> spaa_5_1234_89_7=-(spa78*(spa15*spb81+spa25*spb82+spa35*spb83+spa45*spb84))-spa79*(spa15*spb91+spa25*spb92+spa35*spb93+spa45*spb94);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow3_spab_5_234_6=(pow(spab_5_234_6,3));
const complex<T> pow2_spab_8_79_6=(pow(spab_8_79_6,2));
const complex<T> pow2_spab_5_67_9=(pow(spab_5_67_9,2));
const complex<T> pow2_spaa_5_234_67_5=(pow(spaa_5_234_67_5,2));

return(
(-((pow2_spa18*pow2_spaa_5_234_67_5*spa57)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_89_1))
-
(pow2_spa18*pow3_spab_5_234_6)/(s1789*s2345*spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_8_79_6*spab_5_789_6)/(s6789*s789*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6)
-
(pow2_spa15*pow2_spab_5_67_9*spa57)/(spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpppqb2mq2ppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow3_spa34=(pow(spa34,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_3_24_5=spa23*spb52-spa34*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_3=spa45*(spa13*spb51+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa38*spb86+spa39*spb96)+spa47*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_56_1789_3=spa45*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95)+spa46*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_56_3=(pow(spab_4_56_3,3));
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spab_4_567_3=(pow(spab_4_567_3,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47*spaa_4_567_189_3)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spa13*spaa_4_123_56_4)/(s789*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4*spaa_4_56_1789_3)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5*spab_3_24_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5)
-
(pow2_spa14*pow2_spab_8_679_5*spa13*spab_4_123_5)/(s1234*s6789*spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa18*pow2_spab_4_567_3*spa47)/(s1289*spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3)+(pow2_spa18*pow3_spab_4_56_3)/(s3456*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3)+(pow2_spa18*pow3_spb35)/(s345*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*spb34)
-
(pow2_spa14*pow2_spab_4_567_9*spa13*spa47)/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpppqb2mpq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
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
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa14=(pow(spa14,2));
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_56_23_4=spa45*(spa24*spb52+spa34*spb53)+spa46*(spa24*spb62+spa34*spb63);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_56_4=-(spa45*(spa14*spb51+spa24*spb52+spa34*spb53))-spa46*(spa14*spb61+spa24*spb62+spa34*spb63);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_4_23_5=(pow(spab_4_23_5,3));
const complex<T> pow3_spaa_4_56_23_4=(pow(spaa_4_56_23_4,3));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));
const complex<T> pow2_spaa_4_56_79_8=(pow(spaa_4_56_79_8,2));
const complex<T> pow2_spaa_4_23_567_4=(pow(spaa_4_23_567_4,2));

return(
(-((pow2_spa18*pow2_spaa_4_23_567_4*spa47)/(s189*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_89_1))+(pow2_spa14*pow2_spaa_4_56_79_8*spaa_4_123_56_4)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spaa_4_56_23_4)/(s1789*spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_789_1)+(pow2_spa18*pow3_spab_4_23_5)/(s234*s2345*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_8_679_5*spab_4_123_5)/(s1234*s6789*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5)
-
(pow2_spa14*pow2_spab_4_567_9*spa47)/(spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qppppqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb24=(pow(spb24,3));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+spa39*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> spaa_3_12_456_3=-(spa13*(spa34*spb41+spa35*spb51+spa36*spb61))-spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_12_45_3=-(spa13*(spa34*spb41+spa35*spb51))-spa23*(spa34*spb42+spa35*spb52);
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_3_456_2=(pow(spab_3_456_2,3));
const complex<T> pow3_spab_3_45_2=(pow(spab_3_45_2,3));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_189_2=(pow(spab_3_189_2,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));
const complex<T> pow2_spaa_3_456_79_8=(pow(spaa_3_456_79_8,2));
const complex<T> pow2_spaa_3_45_679_8=(pow(spaa_3_45_679_8,2));

return(
((pow2_spa13*pow2_spaa_3_456_79_8*spaa_3_12_456_3)/(s789*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)+(pow2_spa13*pow2_spaa_3_45_679_8*spaa_3_12_45_3)/(s6789*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1)
-
(pow2_spa13*pow2_spab_8_123_4*spab_3_12_4)/(s123*s1234*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4)+(pow2_spa18*pow3_spab_3_45_2)/(s2345*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2)+(pow2_spa18*pow2_spab_3_189_2*spa37)/(s189*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2)
-
(pow2_spa18*pow3_spab_3_456_2)/(s1789*spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2)+(pow2_spa18*pow3_spb24)/(s234*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*spb23)
-
(pow2_spa13*pow2_spab_3_128_9*spa37)/(spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*spb89))*complex<T>(0,-1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mq2pmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb15=SPB(i1,i5);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa46=(pow(spa46,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_34_5=-((spa36*spb53+spa46*spb54)*spb65)-(spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_234_5=-((spa26*spb52+spa36*spb53+spa46*spb54)*spb65)-(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*spb85+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s567=(spa56*spb65+spa57*spb75+spa67*spb76);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_6_34_5=(pow(spab_6_34_5,3));
const complex<T> pow3_spab_6_234_5=(pow(spab_6_234_5,3));
const complex<T> pow2_spbb_5_34_128_9=(pow(spbb_5_34_128_9,2));
const complex<T> pow2_spbb_5_234_18_9=(pow(spbb_5_234_18_9,2));
const complex<T> pow2_spab_8_679_5=(pow(spab_8_679_5,2));
const complex<T> pow2_spab_6_789_5=(pow(spab_6_789_5,2));
const complex<T> pow2_spab_4_567_9=(pow(spab_4_567_9,2));

return(
(-((pow2_spb79*pow3_spa46)/(s456*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb89))+(pow2_spab_6_789_5*pow2_spb79*spb15)/(s789*spab_6_789_1*spb12*spb23*spb34*spb45*spb89*spbb_5_1234_89_7)
-
(pow2_spb79*pow3_spab_6_234_5)/(s1789*spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*spb89*spbb_5_234_189_7)+(pow2_spb79*pow3_spab_6_34_5)/(s3456*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb89*spbb_5_34_1289_7)
-
(pow2_spab_4_567_9*pow2_spb57*spab_4_67_5)/(s4567*s567*spab_4_56_7*spb12*spb23*spb56*spb89*spbb_5_67_1289_3)
-
(pow2_spb57*pow2_spbb_5_34_128_9*spbb_5_67_34_5)/(s1289*spb12*spb34*spb45*spb56*spb89*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2)+(pow2_spab_8_679_5*pow2_spb57*spb15)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_5_1234_89_7*spbb_5_67_89_1)
-
(pow2_spb57*pow2_spbb_5_234_18_9*spbb_5_67_234_5)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mmq2pmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spab_3_128_9*pow2_spb47*spab_3_567_4)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_567_189_2))+(pow2_spab_8_123_4*pow2_spb47*spb14)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4)/(s3456*s456*spab_3_456_7*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mmmq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13)/(spa89*spb12*spb23*spb34*spb45*spb56*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3)/(s1789*s3456*spab_2_189_7*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_89_7*spbb_3_456_789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpqb2mmmmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));

return(
(-((pow2_spa18*pow2_spb27)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56))
-
(pow2_spab_1_789_2*pow2_spb79)/(s1789*s789*spab_1_89_7*spb23*spb34*spb45*spb56*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmqb2mq2pmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa35=(pow(spa35,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb47=(pow(spb47,2));
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_23_4=-(spb42*(spa25*spb54+spa26*spb64+spa27*spb74))-spb43*(spa35*spb54+spa36*spb64+spa37*spb74);
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+spa26*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_23_1789_5=spb42*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa13*spb51+spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_789_5=spb41*(spa17*spb75+spa18*spb85+spa19*spb95)+spb42*(spa27*spb75+spa28*spb85+spa29*spb95)+spb43*(spa37*spb75+spa38*spb85+spa39*spb95);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_46_5=-(spa34*spb54)+spa36*spb65;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s456=(spa45*spb54+spa46*spb64+spa56*spb65);
const complex<T> s4567=(spa45*spb54+spa46*spb64+spa56*spb65+spa47*spb74+spa57*spb75+spa67*spb76);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spbb_4_23_56_4=(pow(spbb_4_23_56_4,3));
const complex<T> pow3_spab_5_23_4=(pow(spab_5_23_4,3));
const complex<T> pow3_spab_3_56_4=(pow(spab_3_56_4,3));
const complex<T> pow2_spbb_4_56_123_4=(pow(spbb_4_56_123_4,2));
const complex<T> pow2_spbb_4_23_18_9=(pow(spbb_4_23_18_9,2));
const complex<T> pow2_spab_8_123_4=(pow(spab_8_123_4,2));
const complex<T> pow2_spab_5_123_4=(pow(spab_5_123_4,2));
const complex<T> pow2_spab_3_128_9=(pow(spab_3_128_9,2));

return(
(-((pow2_spb79*pow3_spa35)/(s345*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb89))
-
(pow2_spab_5_123_4*pow2_spb79*spb14)/(s6789*spab_5_234_1*spb12*spb23*spb34*spb67*spb89*spbb_4_123_789_6)+(pow2_spb79*pow3_spab_5_23_4)/(s2345*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb89*spbb_4_23_1789_6)
-
(pow2_spab_3_128_9*pow2_spb47*spab_3_567_4*spb57)/(s1289*s4567*spab_3_456_7*spb12*spb45*spb56*spb67*spb89*spbb_4_567_189_2)+(pow2_spab_8_123_4*pow2_spb47*spb14*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1)
-
(pow2_spb47*pow2_spbb_4_23_18_9*spb57*spbb_4_567_23_4)/(s189*spb23*spb34*spb45*spb56*spb67*spb89*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1)+(pow2_spb79*pow3_spab_3_56_4*spab_3_46_5)/(s3456*s456*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb89*spbb_4_56_1789_2)+(pow2_spb79*pow2_spbb_4_56_123_4*spb14*spbb_4_123_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1)
-
(pow2_spb79*pow3_spbb_4_23_56_4*spbb_4_23_1789_5)/(s1789*spb23*spb34*spb45*spb56*spb89*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmqb2mmq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_5=spb31*(spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_346_5=-(spa23*spb53)-spa24*spb54+spa26*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb57)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1))+(pow2_spab_8_12_3*pow2_spb37*spb13*spb57)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_346_5)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_5)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3)/(s2345*s345*spab_2_345_6*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmqb2mmmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_5=spa17*spb75+spa18*spb85+spa19*spb95;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb57)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_5)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79)/(s2345*s6789*spab_1_789_6*spb23*spb34*spb45*spb67*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmqb2mq2pmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
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
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa59=SPA(i5,i9);
const complex<T> spa58=SPA(i5,i8);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa49=SPA(i4,i9);
const complex<T> spa48=SPA(i4,i8);
const complex<T> spa47=SPA(i4,i7);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa29=SPA(i2,i9);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa24=(pow(spa24,3));
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb37=(pow(spb37,2));
const complex<T> spbb_3_456_789_1=-(spb43*(spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa57*spb71+spa58*spb81+spa59*spb91)-spb63*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_3_45_6789_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81+spa49*spb91))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91);
const complex<T> spbb_3_4567_89_1=(-(spa48*spb43)-spa58*spb53-spa68*spb63-spa78*spb73)*spb81+(-(spa49*spb43)-spa59*spb53-spa69*spb63-spa79*spb73)*spb91;
const complex<T> spbb_3_456_12_3=-(spb31*(spa14*spb43+spa15*spb53+spa16*spb63))-spb32*(spa24*spb43+spa25*spb53+spa26*spb63);
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(spa15*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_89_7=spb31*(spa18*spb87+spa19*spb97)+spb32*(spa28*spb87+spa29*spb97);
const complex<T> spbb_3_12_789_6=spb31*(spa17*spb76+spa18*spb86+spa19*spb96)+spb32*(spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spbb_3_12_789_4=spb31*(spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spbb_3_12_6789_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85+spa29*spb95);
const complex<T> spbb_3_12_6789_4=spb31*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spb32*(spa26*spb64+spa27*spb74+spa28*spb84+spa29*spb94);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_356_4=-(spa23*spb43)+spa25*spb54+spa26*spb64;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_2_189_3=spa12*spb31+spa28*spb83+spa29*spb93;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s345=(spa34*spb43+spa35*spb53+spa45*spb54);
const complex<T> s3456=(spa34*spb43+spa35*spb53+spa45*spb54+spa36*spb63+spa46*spb64+spa56*spb65);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1289=(spa12*spb21+spa18*spb81+spa28*spb82+spa19*spb91+spa29*spb92+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow3_spab_2_456_3=(pow(spab_2_456_3,3));
const complex<T> pow3_spab_2_45_3=(pow(spab_2_45_3,3));
const complex<T> pow2_spbb_3_456_12_3=(pow(spbb_3_456_12_3,2));
const complex<T> pow2_spbb_3_45_12_3=(pow(spbb_3_45_12_3,2));
const complex<T> pow2_spab_8_12_3=(pow(spab_8_12_3,2));
const complex<T> pow2_spab_4_12_3=(pow(spab_4_12_3,2));
const complex<T> pow2_spab_2_18_9=(pow(spab_2_18_9,2));

return(
(-((pow2_spb79*pow3_spa24)/(s234*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb89))
-
(pow2_spab_4_12_3*pow2_spb79*spb13)/(s1234*spab_4_23_1*spb12*spb23*spb56*spb67*spb89*spbb_3_12_6789_5)
-
(pow2_spab_2_18_9*pow2_spb37*spab_2_189_3*spb47)/(s1289*s189*spab_2_189_7*spb34*spb45*spb56*spb67*spb89*spbb_3_4567_89_1)+(pow2_spab_8_12_3*pow2_spb37*spb13*spb47)/(spa89*spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1)
-
(pow2_spb79*pow3_spab_2_456_3*spab_2_356_4)/(s1789*s3456*spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*spb89*spbb_3_456_789_1)+(pow2_spb79*pow2_spbb_3_456_12_3*spb13*spbb_3_12_789_4)/(s789*spb12*spb23*spb34*spb45*spb56*spb89*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1)+(pow2_spb79*pow3_spab_2_45_3*spab_2_35_4)/(s2345*s345*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb89*spbb_3_45_6789_1)+(pow2_spb79*pow2_spbb_3_45_12_3*spb13*spbb_3_12_6789_4)/(s6789*spb12*spb23*spb34*spb45*spb67*spb89*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmqb2mmq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb47=SPB(i4,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_4=spa17*spb74+spa18*spb84+spa19*spb94;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_235_4=-(spa12*spb42)-spa13*spb43+spa15*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb47)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_4)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_235_4)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79)/(s1234*s234*spab_1_234_5*spb23*spb34*spb56*spb67*spb89))*complex<T>(0,1)
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q2Q3g2l_qpmmmqb2mq2pqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb89=SPB(i8,i9);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb72=SPB(i7,i2);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb37=SPB(i3,i7);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb27=SPB(i2,i7);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa69=SPA(i6,i9);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=(pow(spb79,2));
const complex<T> pow2_spb27=(pow(spb27,2));
const complex<T> pow2_spa18=(pow(spa18,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_3=spa17*spb73+spa18*spb83+spa19*spb93;
const complex<T> spab_1_789_2=spa17*spb72+spa18*spb82+spa19*spb92;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_245_3=-(spa12*spb32)+spa14*spb43+spa15*spb53;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s789=(spa78*spb87+spa79*spb97+spa89*spb98);
const complex<T> s6789=(spa67*spb76+spa68*spb86+spa78*spb87+spa69*spb96+spa79*spb97+spa89*spb98);
const complex<T> s234=(spa23*spb32+spa24*spb42+spa34*spb43);
const complex<T> s2345=(spa23*spb32+spa24*spb42+spa34*spb43+spa25*spb52+spa35*spb53+spa45*spb54);
const complex<T> s189=(spa18*spb81+spa19*spb91+spa89*spb98);
const complex<T> s1789=(spa17*spb71+spa18*spb81+spa78*spb87+spa19*spb91+spa79*spb97+spa89*spb98);
const complex<T> s123=(spa12*spb21+spa13*spb31+spa23*spb32);
const complex<T> s1234=(spa12*spb21+spa13*spb31+spa23*spb32+spa14*spb41+spa24*spb42+spa34*spb43);
const complex<T> pow2_spab_1_789_2=(pow(spab_1_789_2,2));
const complex<T> pow2_spab_1_345_2=(pow(spab_1_345_2,2));
const complex<T> pow2_spab_1_34_2=(pow(spab_1_34_2,2));

return(
(-((pow2_spa18*pow2_spb27*spb37)/(s189*spa89*spab_1_89_7*spb23*spb34*spb45*spb56*spb67))
-
(pow2_spab_1_789_2*pow2_spb79*spab_1_789_3)/(s1789*s789*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb89)
-
(pow2_spab_1_345_2*pow2_spb79*spab_1_245_3)/(s2345*s6789*spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb89)+(pow2_spab_1_34_2*pow2_spb79*spab_1_24_3)/(s1234*s234*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb89)
-
(pow2_spa13*pow2_spb79)/(s123*spa23*spab_1_23_4*spb45*spb56*spb67*spb89))*complex<T>(0,1)
);
}




template <class T> complex<T> (*A2q2Q3g2l_Tree_Ptr_eval(long long hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 115423521:	 return &A2q2Q3g2l_qpmmmq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m m q2p qb2m qbm lp lbm
case 115424105:	 return &A2q2Q3g2l_qppppq2pqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p p q2p qb2m qbm lp lbm
case 115427105:	 return &A2q2Q3g2l_qpmmq2pmqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m q2p m qb2m qbm lp lbm
case 115427553:	 return &A2q2Q3g2l_qpmq2pmmqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p m m qb2m qbm lp lbm
case 115427609:	 return &A2q2Q3g2l_qpq2pmmmqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m m m qb2m qbm lp lbm
case 115431273:	 return &A2q2Q3g2l_qpppq2ppqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p q2p p qb2m qbm lp lbm
case 115432169:	 return &A2q2Q3g2l_qppq2pppqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p p p qb2m qbm lp lbm
case 115432281:	 return &A2q2Q3g2l_qpq2ppppqb2mqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p p p qb2m qbm lp lbm
case 115452193:	 return &A2q2Q3g2l_qpmmmqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m m qb2m q2p qbm lp lbm
case 115452777:	 return &A2q2Q3g2l_qppppqb2mq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p p qb2m q2p qbm lp lbm
case 115459361:	 return &A2q2Q3g2l_qpmmqb2mmq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m qb2m m q2p qbm lp lbm
case 115460257:	 return &A2q2Q3g2l_qpmqb2mmmq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m m m q2p qbm lp lbm
case 115460369:	 return &A2q2Q3g2l_qpqb2mmmmq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m m m q2p qbm lp lbm
case 115463529:	 return &A2q2Q3g2l_qpppqb2mpq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p qb2m p q2p qbm lp lbm
case 115464873:	 return &A2q2Q3g2l_qppqb2mppq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m p p q2p qbm lp lbm
case 115465041:	 return &A2q2Q3g2l_qpqb2mpppq2pqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p p p q2p qbm lp lbm
case 115484449:	 return &A2q2Q3g2l_qpmmq2pqb2mmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m q2p qb2m m qbm lp lbm
case 115484897:	 return &A2q2Q3g2l_qpmq2pmqb2mmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p m qb2m m qbm lp lbm
case 115484953:	 return &A2q2Q3g2l_qpq2pmmqb2mmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m m qb2m m qbm lp lbm
case 115488033:	 return &A2q2Q3g2l_qpmmqb2mq2pmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m qb2m q2p m qbm lp lbm
case 115488929:	 return &A2q2Q3g2l_qpmqb2mmq2pmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m m q2p m qbm lp lbm
case 115489041:	 return &A2q2Q3g2l_qpqb2mmmq2pmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m m q2p m qbm lp lbm
case 115492065:	 return &A2q2Q3g2l_qpmq2pqb2mmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p qb2m m m qbm lp lbm
case 115492121:	 return &A2q2Q3g2l_qpq2pmqb2mmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m qb2m m m qbm lp lbm
case 115492513:	 return &A2q2Q3g2l_qpmqb2mq2pmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m q2p m m qbm lp lbm
case 115492625:	 return &A2q2Q3g2l_qpqb2mmq2pmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m q2p m m qbm lp lbm
case 115493017:	 return &A2q2Q3g2l_qpq2pqb2mmmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p qb2m m m m qbm lp lbm
case 115493073:	 return &A2q2Q3g2l_qpqb2mq2pmmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m q2p m m m qbm lp lbm
case 115517289:	 return &A2q2Q3g2l_qpppq2pqb2mpqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p q2p qb2m p qbm lp lbm
case 115518185:	 return &A2q2Q3g2l_qppq2ppqb2mpqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p p qb2m p qbm lp lbm
case 115518297:	 return &A2q2Q3g2l_qpq2pppqb2mpqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p p qb2m p qbm lp lbm
case 115520873:	 return &A2q2Q3g2l_qpppqb2mq2ppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p qb2m q2p p qbm lp lbm
case 115522217:	 return &A2q2Q3g2l_qppqb2mpq2ppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m p q2p p qbm lp lbm
case 115522385:	 return &A2q2Q3g2l_qpqb2mppq2ppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p p q2p p qbm lp lbm
case 115528937:	 return &A2q2Q3g2l_qppq2pqb2mppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p qb2m p p qbm lp lbm
case 115529049:	 return &A2q2Q3g2l_qpq2ppqb2mppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p qb2m p p qbm lp lbm
case 115529385:	 return &A2q2Q3g2l_qppqb2mq2pppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m q2p p p qbm lp lbm
case 115529553:	 return &A2q2Q3g2l_qpqb2mpq2pppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p q2p p p qbm lp lbm
case 115530393:	 return &A2q2Q3g2l_qpq2pqb2mpppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp q2p qb2m p p p qbm lp lbm
case 115530449:	 return &A2q2Q3g2l_qpqb2mq2ppppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m q2p p p p qbm lp lbm
case 115685664:	 return &A2q2Q3g2l_qmmmmq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m m q2p qb2m qbp lp lbm
case 115686248:	 return &A2q2Q3g2l_qmpppq2pqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p p q2p qb2m qbp lp lbm
case 115689248:	 return &A2q2Q3g2l_qmmmq2pmqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m q2p m qb2m qbp lp lbm
case 115689696:	 return &A2q2Q3g2l_qmmq2pmmqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p m m qb2m qbp lp lbm
case 115689752:	 return &A2q2Q3g2l_qmq2pmmmqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m m m qb2m qbp lp lbm
case 115693416:	 return &A2q2Q3g2l_qmppq2ppqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p q2p p qb2m qbp lp lbm
case 115694312:	 return &A2q2Q3g2l_qmpq2pppqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p p p qb2m qbp lp lbm
case 115694424:	 return &A2q2Q3g2l_qmq2ppppqb2mqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p p p qb2m qbp lp lbm
case 115714336:	 return &A2q2Q3g2l_qmmmmqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m m qb2m q2p qbp lp lbm
case 115714920:	 return &A2q2Q3g2l_qmpppqb2mq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p p qb2m q2p qbp lp lbm
case 115721504:	 return &A2q2Q3g2l_qmmmqb2mmq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m qb2m m q2p qbp lp lbm
case 115722400:	 return &A2q2Q3g2l_qmmqb2mmmq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m m m q2p qbp lp lbm
case 115722512:	 return &A2q2Q3g2l_qmqb2mmmmq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m m m q2p qbp lp lbm
case 115725672:	 return &A2q2Q3g2l_qmppqb2mpq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p qb2m p q2p qbp lp lbm
case 115727016:	 return &A2q2Q3g2l_qmpqb2mppq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m p p q2p qbp lp lbm
case 115727184:	 return &A2q2Q3g2l_qmqb2mpppq2pqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p p p q2p qbp lp lbm
case 115746592:	 return &A2q2Q3g2l_qmmmq2pqb2mmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m q2p qb2m m qbp lp lbm
case 115747040:	 return &A2q2Q3g2l_qmmq2pmqb2mmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p m qb2m m qbp lp lbm
case 115747096:	 return &A2q2Q3g2l_qmq2pmmqb2mmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m m qb2m m qbp lp lbm
case 115750176:	 return &A2q2Q3g2l_qmmmqb2mq2pmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m qb2m q2p m qbp lp lbm
case 115751072:	 return &A2q2Q3g2l_qmmqb2mmq2pmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m m q2p m qbp lp lbm
case 115751184:	 return &A2q2Q3g2l_qmqb2mmmq2pmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m m q2p m qbp lp lbm
case 115754208:	 return &A2q2Q3g2l_qmmq2pqb2mmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p qb2m m m qbp lp lbm
case 115754264:	 return &A2q2Q3g2l_qmq2pmqb2mmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m qb2m m m qbp lp lbm
case 115754656:	 return &A2q2Q3g2l_qmmqb2mq2pmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m q2p m m qbp lp lbm
case 115754768:	 return &A2q2Q3g2l_qmqb2mmq2pmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m q2p m m qbp lp lbm
case 115755160:	 return &A2q2Q3g2l_qmq2pqb2mmmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p qb2m m m m qbp lp lbm
case 115755216:	 return &A2q2Q3g2l_qmqb2mq2pmmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m q2p m m m qbp lp lbm
case 115779432:	 return &A2q2Q3g2l_qmppq2pqb2mpqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p q2p qb2m p qbp lp lbm
case 115780328:	 return &A2q2Q3g2l_qmpq2ppqb2mpqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p p qb2m p qbp lp lbm
case 115780440:	 return &A2q2Q3g2l_qmq2pppqb2mpqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p p qb2m p qbp lp lbm
case 115783016:	 return &A2q2Q3g2l_qmppqb2mq2ppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p qb2m q2p p qbp lp lbm
case 115784360:	 return &A2q2Q3g2l_qmpqb2mpq2ppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m p q2p p qbp lp lbm
case 115784528:	 return &A2q2Q3g2l_qmqb2mppq2ppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p p q2p p qbp lp lbm
case 115791080:	 return &A2q2Q3g2l_qmpq2pqb2mppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p qb2m p p qbp lp lbm
case 115791192:	 return &A2q2Q3g2l_qmq2ppqb2mppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p qb2m p p qbp lp lbm
case 115791528:	 return &A2q2Q3g2l_qmpqb2mq2pppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m q2p p p qbp lp lbm
case 115791696:	 return &A2q2Q3g2l_qmqb2mpq2pppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p q2p p p qbp lp lbm
case 115792536:	 return &A2q2Q3g2l_qmq2pqb2mpppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm q2p qb2m p p p qbp lp lbm
case 115792592:	 return &A2q2Q3g2l_qmqb2mq2ppppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m q2p p p p qbp lp lbm
case 130103585:	 return &A2q2Q3g2l_qpmmmq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m m q2p qb2m qbm lm lbp
case 130104169:	 return &A2q2Q3g2l_qppppq2pqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p p q2p qb2m qbm lm lbp
case 130107169:	 return &A2q2Q3g2l_qpmmq2pmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m q2p m qb2m qbm lm lbp
case 130107617:	 return &A2q2Q3g2l_qpmq2pmmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p m m qb2m qbm lm lbp
case 130107673:	 return &A2q2Q3g2l_qpq2pmmmqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m m m qb2m qbm lm lbp
case 130111337:	 return &A2q2Q3g2l_qpppq2ppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p q2p p qb2m qbm lm lbp
case 130112233:	 return &A2q2Q3g2l_qppq2pppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p p p qb2m qbm lm lbp
case 130112345:	 return &A2q2Q3g2l_qpq2ppppqb2mqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p p p qb2m qbm lm lbp
case 130132257:	 return &A2q2Q3g2l_qpmmmqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m m qb2m q2p qbm lm lbp
case 130132841:	 return &A2q2Q3g2l_qppppqb2mq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p p qb2m q2p qbm lm lbp
case 130139425:	 return &A2q2Q3g2l_qpmmqb2mmq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m qb2m m q2p qbm lm lbp
case 130140321:	 return &A2q2Q3g2l_qpmqb2mmmq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m m m q2p qbm lm lbp
case 130140433:	 return &A2q2Q3g2l_qpqb2mmmmq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m m m q2p qbm lm lbp
case 130143593:	 return &A2q2Q3g2l_qpppqb2mpq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p qb2m p q2p qbm lm lbp
case 130144937:	 return &A2q2Q3g2l_qppqb2mppq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m p p q2p qbm lm lbp
case 130145105:	 return &A2q2Q3g2l_qpqb2mpppq2pqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p p p q2p qbm lm lbp
case 130164513:	 return &A2q2Q3g2l_qpmmq2pqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m q2p qb2m m qbm lm lbp
case 130164961:	 return &A2q2Q3g2l_qpmq2pmqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p m qb2m m qbm lm lbp
case 130165017:	 return &A2q2Q3g2l_qpq2pmmqb2mmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m m qb2m m qbm lm lbp
case 130168097:	 return &A2q2Q3g2l_qpmmqb2mq2pmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m qb2m q2p m qbm lm lbp
case 130168993:	 return &A2q2Q3g2l_qpmqb2mmq2pmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m m q2p m qbm lm lbp
case 130169105:	 return &A2q2Q3g2l_qpqb2mmmq2pmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m m q2p m qbm lm lbp
case 130172129:	 return &A2q2Q3g2l_qpmq2pqb2mmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m q2p qb2m m m qbm lm lbp
case 130172185:	 return &A2q2Q3g2l_qpq2pmqb2mmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p m qb2m m m qbm lm lbp
case 130172577:	 return &A2q2Q3g2l_qpmqb2mq2pmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m qb2m q2p m m qbm lm lbp
case 130172689:	 return &A2q2Q3g2l_qpqb2mmq2pmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m m q2p m m qbm lm lbp
case 130173081:	 return &A2q2Q3g2l_qpq2pqb2mmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p qb2m m m m qbm lm lbp
case 130173137:	 return &A2q2Q3g2l_qpqb2mq2pmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m q2p m m m qbm lm lbp
case 130197353:	 return &A2q2Q3g2l_qpppq2pqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p q2p qb2m p qbm lm lbp
case 130198249:	 return &A2q2Q3g2l_qppq2ppqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p p qb2m p qbm lm lbp
case 130198361:	 return &A2q2Q3g2l_qpq2pppqb2mpqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p p qb2m p qbm lm lbp
case 130200937:	 return &A2q2Q3g2l_qpppqb2mq2ppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p qb2m q2p p qbm lm lbp
case 130202281:	 return &A2q2Q3g2l_qppqb2mpq2ppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m p q2p p qbm lm lbp
case 130202449:	 return &A2q2Q3g2l_qpqb2mppq2ppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p p q2p p qbm lm lbp
case 130209001:	 return &A2q2Q3g2l_qppq2pqb2mppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p q2p qb2m p p qbm lm lbp
case 130209113:	 return &A2q2Q3g2l_qpq2ppqb2mppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p p qb2m p p qbm lm lbp
case 130209449:	 return &A2q2Q3g2l_qppqb2mq2pppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p qb2m q2p p p qbm lm lbp
case 130209617:	 return &A2q2Q3g2l_qpqb2mpq2pppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m p q2p p p qbm lm lbp
case 130210457:	 return &A2q2Q3g2l_qpq2pqb2mpppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp q2p qb2m p p p qbm lm lbp
case 130210513:	 return &A2q2Q3g2l_qpqb2mq2ppppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp qb2m q2p p p p qbm lm lbp
case 130365728:	 return &A2q2Q3g2l_qmmmmq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m m q2p qb2m qbp lm lbp
case 130366312:	 return &A2q2Q3g2l_qmpppq2pqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p p q2p qb2m qbp lm lbp
case 130369312:	 return &A2q2Q3g2l_qmmmq2pmqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m q2p m qb2m qbp lm lbp
case 130369760:	 return &A2q2Q3g2l_qmmq2pmmqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p m m qb2m qbp lm lbp
case 130369816:	 return &A2q2Q3g2l_qmq2pmmmqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m m m qb2m qbp lm lbp
case 130373480:	 return &A2q2Q3g2l_qmppq2ppqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p q2p p qb2m qbp lm lbp
case 130374376:	 return &A2q2Q3g2l_qmpq2pppqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p p p qb2m qbp lm lbp
case 130374488:	 return &A2q2Q3g2l_qmq2ppppqb2mqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p p p qb2m qbp lm lbp
case 130394400:	 return &A2q2Q3g2l_qmmmmqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m m qb2m q2p qbp lm lbp
case 130394984:	 return &A2q2Q3g2l_qmpppqb2mq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p p qb2m q2p qbp lm lbp
case 130401568:	 return &A2q2Q3g2l_qmmmqb2mmq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m qb2m m q2p qbp lm lbp
case 130402464:	 return &A2q2Q3g2l_qmmqb2mmmq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m m m q2p qbp lm lbp
case 130402576:	 return &A2q2Q3g2l_qmqb2mmmmq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m m m q2p qbp lm lbp
case 130405736:	 return &A2q2Q3g2l_qmppqb2mpq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p qb2m p q2p qbp lm lbp
case 130407080:	 return &A2q2Q3g2l_qmpqb2mppq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m p p q2p qbp lm lbp
case 130407248:	 return &A2q2Q3g2l_qmqb2mpppq2pqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p p p q2p qbp lm lbp
case 130426656:	 return &A2q2Q3g2l_qmmmq2pqb2mmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m q2p qb2m m qbp lm lbp
case 130427104:	 return &A2q2Q3g2l_qmmq2pmqb2mmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p m qb2m m qbp lm lbp
case 130427160:	 return &A2q2Q3g2l_qmq2pmmqb2mmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m m qb2m m qbp lm lbp
case 130430240:	 return &A2q2Q3g2l_qmmmqb2mq2pmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m qb2m q2p m qbp lm lbp
case 130431136:	 return &A2q2Q3g2l_qmmqb2mmq2pmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m m q2p m qbp lm lbp
case 130431248:	 return &A2q2Q3g2l_qmqb2mmmq2pmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m m q2p m qbp lm lbp
case 130434272:	 return &A2q2Q3g2l_qmmq2pqb2mmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m q2p qb2m m m qbp lm lbp
case 130434328:	 return &A2q2Q3g2l_qmq2pmqb2mmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p m qb2m m m qbp lm lbp
case 130434720:	 return &A2q2Q3g2l_qmmqb2mq2pmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m qb2m q2p m m qbp lm lbp
case 130434832:	 return &A2q2Q3g2l_qmqb2mmq2pmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m m q2p m m qbp lm lbp
case 130435224:	 return &A2q2Q3g2l_qmq2pqb2mmmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p qb2m m m m qbp lm lbp
case 130435280:	 return &A2q2Q3g2l_qmqb2mq2pmmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m q2p m m m qbp lm lbp
case 130459496:	 return &A2q2Q3g2l_qmppq2pqb2mpqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p q2p qb2m p qbp lm lbp
case 130460392:	 return &A2q2Q3g2l_qmpq2ppqb2mpqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p p qb2m p qbp lm lbp
case 130460504:	 return &A2q2Q3g2l_qmq2pppqb2mpqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p p qb2m p qbp lm lbp
case 130463080:	 return &A2q2Q3g2l_qmppqb2mq2ppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p qb2m q2p p qbp lm lbp
case 130464424:	 return &A2q2Q3g2l_qmpqb2mpq2ppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m p q2p p qbp lm lbp
case 130464592:	 return &A2q2Q3g2l_qmqb2mppq2ppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p p q2p p qbp lm lbp
case 130471144:	 return &A2q2Q3g2l_qmpq2pqb2mppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p q2p qb2m p p qbp lm lbp
case 130471256:	 return &A2q2Q3g2l_qmq2ppqb2mppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p p qb2m p p qbp lm lbp
case 130471592:	 return &A2q2Q3g2l_qmpqb2mq2pppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p qb2m q2p p p qbp lm lbp
case 130471760:	 return &A2q2Q3g2l_qmqb2mpq2pppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m p q2p p p qbp lm lbp
case 130472600:	 return &A2q2Q3g2l_qmq2pqb2mpppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm q2p qb2m p p p qbp lm lbp
case 130472656:	 return &A2q2Q3g2l_qmqb2mq2ppppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm qb2m q2p p p p qbp lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q2Q3g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q2Q3g2l_Tree_Ptr_eval(long long hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q2Q3g2l_Tree_Ptr_eval(long long hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q2Q3g2l_Tree_Ptr_eval(long long hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A2q2Q3g2l_Tree_Ptr_eval(long long hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}
