/*
 * A0_2q4g2l_eval.cpp
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
 * The 2 quarks 4g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmmmmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(SPA(i1,i2)*spb52)-spa13*spb53-spa14*spb54)*spb81)-(SPA(i2,i7)*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=SPA(i2,i6)*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,1)*(-((spa46*pow(spab_4_56_8,2))/((s45+s46+s56)*spa45*spa56*spab_6_45_3*spb12*spb23*spb78)
)+
(spab_6_234_5*pow(spbb_5_234_17_8,2))/(spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78*(s78+spb71*SPA(i1,i7
)+
spb81*SPA(i1,i8))*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(spab_6_34_5*pow(spbb_5_34_127_8,2))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow(spab_7_68_5,2)*SPB(i1,i5))/(spa78*spab_6_78_1*spb12*spb23*spb34*spb45*(s78+spa68*spb86+spa67*SPB(i7,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmmmmpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(SPA(i1,i2)*spb52)-spa13*spb53-spa14*spb54)*spb81)-(SPA(i2,i7)*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=SPA(i2,i6)*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,-1)*(-((spa46*pow(spab_4_56_8,2))/((s45+s46+s56)*spa45*spa56*spab_6_45_3*spb12*spb23*spb78)
)+
(spab_6_234_5*pow(spbb_5_234_17_8,2))/(spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78*(s78+spb71*SPA(i1,i7
)+
spb81*SPA(i1,i8))*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(spab_6_34_5*pow(spbb_5_34_127_8,2))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow(spab_7_68_5,2)*SPB(i1,i5))/(spa78*spab_6_78_1*spb12*spb23*spb34*spb45*(s78+spa68*spb86+spa67*SPB(i7,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmmmpmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=pow(spb46,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+SPA(i2,i6)*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_8=-(spa34*SPB(i8,i4))-spa35*SPB(i8,i5)-spa36*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
)-
(pow3_spb46*spb14*pow(spab_7_123_4,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1
)-
(pow3_spb46*spbb_4_23_56_4*pow(spbb_4_23_17_8,2))/(spb23*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1
)-
(pow2_spb68*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_4_123_78_6*(s23+s24+s34+spa12*spb21+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow2_spb68*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb78*spbb_4_23_178_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb46*spab_3_56_4*pow(spab_3_456_8,2))/((s45+s46+s56)*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmmmpmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=pow(spb46,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+SPA(i2,i6)*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_8=-(spa34*SPB(i8,i4))-spa35*SPB(i8,i5)-spa36*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
)-
(pow3_spb46*spb14*pow(spab_7_123_4,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1
)-
(pow3_spb46*spbb_4_23_56_4*pow(spbb_4_23_17_8,2))/(spb23*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1
)-
(pow2_spb68*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_4_123_78_6*(s23+s24+s34+spa12*spb21+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow2_spb68*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb78*spbb_4_23_178_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb46*spab_3_56_4*pow(spab_3_456_8,2))/((s45+s46+s56)*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmmpmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=pow(spb36,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(SPA(i1,i5)*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)
)-
(pow3_spb36*spb13*pow(spab_7_12_3,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1
)-
(pow2_spb68*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1
)+
(pow3_spb36*spab_2_456_3*pow(spab_2_17_8,2))/(spab_2_345_6*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_3_456_78_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow2_spb68*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb68*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmmpmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=pow(spb36,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(SPA(i1,i5)*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)
)-
(pow3_spb36*spb13*pow(spab_7_12_3,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1
)-
(pow2_spb68*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1
)+
(pow3_spb36*spab_2_456_3*pow(spab_2_17_8,2))/(spab_2_345_6*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_3_456_78_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow2_spb68*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb68*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmmpppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*SPB(i6,i5))-spa78*(-(spa23*spb83)-spa24*spb84-spa25*SPB(i8,i5));
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((spab_2_34_1*pow(spaa_2_34_568_7,2))/((s23+s24+s34)*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)-
(pow(spab_7_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa78*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(pow(spab_2_17_8,2)*SPA(i2,i6))/(spa23*spa34*spa45*spa56*spab_6_78_1*(s78+spa17*spb71+spb81*SPA(i1,i8))*SPB(i7,i8)
)+
(spab_2_345_1*pow(spaa_2_345_68_7,2))/(spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1*(s23+s24+s34+spa35*spb53+spa45*spb54+spa25*SPB(i5,i2))*(s78+spa67*SPB(i7,i6
)+
spa68*SPB(i8,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmmpppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*SPB(i6,i5))-spa78*(-(spa23*spb83)-spa24*spb84-spa25*SPB(i8,i5));
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((spab_2_34_1*pow(spaa_2_34_568_7,2))/((s23+s24+s34)*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)-
(pow(spab_7_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa78*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(pow(spab_2_17_8,2)*SPA(i2,i6))/(spa23*spa34*spa45*spa56*spab_6_78_1*(s78+spa17*spb71+spb81*SPA(i1,i8))*SPB(i7,i8)
)+
(spab_2_345_1*pow(spaa_2_345_68_7,2))/(spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1*(s23+s24+s34+spa35*spb53+spa45*spb54+spa25*SPB(i5,i2))*(s78+spa67*SPB(i7,i6
)+
spa68*SPB(i8,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmpmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb78)
)+
(pow2_spb68*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_78_6*spb23*spb34*spb45*spb78*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))
)-
(pow2_spb68*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb78*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow(spa17,2)*pow(SPB(i2,i6),3))/(spa78*spab_1_78_6*spb23*spb34*spb45*spb56*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmpmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb78)
)+
(pow2_spb68*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_78_6*spb23*spb34*spb45*spb78*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))
)-
(pow2_spb68*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb78*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow(spa17,2)*pow(SPB(i2,i6),3))/(spa78*spab_1_78_6*spb23*spb34*spb45*spb56*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmpmppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_123_4=spa17*spb41+SPA(i2,i7)*spb42+SPA(i3,i7)*spb43;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*SPB(i5,i1)+spa23*spb52);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_68_7,2))/(spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1*(s78+spa67*spb76+spa68*spb86))
)-
(pow3_spa13*spab_3_12_4*pow(spab_7_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa17*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa17*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa78*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa17*spa36*pow(spab_3_456_2,3))/(spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_345_2*(s78+spa17*spb71+spa18*spb81)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa36*pow(spab_3_456_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*SPB(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmpmppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_123_4=spa17*spb41+SPA(i2,i7)*spb42+SPA(i3,i7)*spb43;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*SPB(i5,i1)+spa23*spb52);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_68_7,2))/(spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1*(s78+spa67*spb76+spa68*spb86))
)-
(pow3_spa13*spab_3_12_4*pow(spab_7_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa17*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa17*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa78*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa17*spa36*pow(spab_3_456_2,3))/(spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_345_2*(s78+spa17*spb71+spa18*spb81)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa36*pow(spab_3_456_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*SPB(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmppmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*SPB(i6,i2)))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,1)*(-((pow2_spa17*spa46*pow(spaa_4_23_56_4,3))/(spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1*(s78+spa17*spb71+spa18*spb81))
)+
(pow2_spa17*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa17*spa46*pow(spab_4_56_3,3))/((s45+s46+s56)*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_7_68_5,2))/(spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_234_5*(s78+spa67*spb76+spa68*spb86)*(s23+s24+s34+spa12*spb21+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa17*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa78*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)-
(pow3_spa14*spa46*pow(spab_4_56_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*SPB(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmppmpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*SPB(i6,i2)))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,-1)*(-((pow2_spa17*spa46*pow(spaa_4_23_56_4,3))/(spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1*(s78+spa17*spb71+spa18*spb81))
)+
(pow2_spa17*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa17*spa46*pow(spab_4_56_3,3))/((s45+s46+s56)*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_7_68_5,2))/(spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_234_5*(s78+spa67*spb76+spa68*spb86)*(s23+s24+s34+spa12*spb21+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa17*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa78*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)-
(pow3_spa14*spa46*pow(spab_4_56_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*SPB(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmpppmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spa78=SPA(i7,i8);
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
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*SPB(i6,i2)+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,1)*(-((pow2_spa17*pow(spab_5_34_6,3))/((s34+s35+s45)*spa12*spa34*spa45*spa78*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6)))
)-
(pow2_spa17*pow(SPB(i4,i6),3))/((s45+s46+s56)*spa12*spa23*spa78*spab_3_45_6*SPB(i4,i5)*SPB(i5,i6)
)+
(pow(SPA(i1,i5),3)*pow(SPB(i6,i8),2))/(spa12*spa23*spa34*spa45*spab_1_78_6*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))*SPB(i7,i8)
)+
(pow2_spa17*pow(spab_5_234_6,3))/(spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_345_6*(s34+s35+s45+spa23*SPB(i3,i2
)+
spa24*SPB(i4,i2
)+
spa25*SPB(i5,i2))*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmpppmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spa78=SPA(i7,i8);
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
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*SPB(i6,i2)+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,-1)*(-((pow2_spa17*pow(spab_5_34_6,3))/((s34+s35+s45)*spa12*spa34*spa45*spa78*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6)))
)-
(pow2_spa17*pow(SPB(i4,i6),3))/((s45+s46+s56)*spa12*spa23*spa78*spab_3_45_6*SPB(i4,i5)*SPB(i5,i6)
)+
(pow(SPA(i1,i5),3)*pow(SPB(i6,i8),2))/(spa12*spa23*spa34*spa45*spab_1_78_6*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))*SPB(i7,i8)
)+
(pow2_spa17*pow(spab_5_234_6,3))/(spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_345_6*(s34+s35+s45+spa23*SPB(i3,i2
)+
spa24*SPB(i4,i2
)+
spa25*SPB(i5,i2))*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpmmmpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb78)
)+
(pow2_spb68*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_78_6*spb23*spb34*spb45*spb78*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))
)-
(pow2_spb68*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb78*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow(spa17,2)*pow(SPB(i2,i6),3))/(spa78*spab_1_78_6*spb23*spb34*spb45*spb56*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpmmmpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb76=SPB(i7,i6);
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
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb78)
)+
(pow2_spb68*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_78_6*spb23*spb34*spb45*spb78*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))
)-
(pow2_spb68*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb78*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow(spa17,2)*pow(SPB(i2,i6),3))/(spa78*spab_1_78_6*spb23*spb34*spb45*spb56*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpmmpmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=pow(spb36,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(SPA(i1,i5)*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)
)-
(pow3_spb36*spb13*pow(spab_7_12_3,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1
)-
(pow2_spb68*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1
)+
(pow3_spb36*spab_2_456_3*pow(spab_2_17_8,2))/(spab_2_345_6*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_3_456_78_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow2_spb68*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb68*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpmmpmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb36=pow(spb36,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_3_45_678_1=-(spb43*(spa46*spb61+spa47*spb71+spa48*spb81))-spb53*(spa56*spb61+spa57*spb71+spa58*spb81);
const complex<T> spbb_3_456_78_1=(-(spa47*spb43)-spa57*spb53-spa67*spb63)*spb71+(-(spa48*spb43)-spa58*spb53-spa68*spb63)*spb81;
const complex<T> spbb_3_45_12_3=-((spa14*spb31+spa24*spb32)*spb43)-(SPA(i1,i5)*spb31+spa25*spb32)*spb53;
const complex<T> spbb_3_12_78_6=spb31*(spa17*spb76+spa18*spb86)+spb32*(spa27*spb76+spa28*spb86);
const complex<T> spbb_3_12_678_5=spb31*(spa16*spb65+spa17*spb75+spa18*spb85)+spb32*(spa26*spb65+spa27*spb75+spa28*spb85);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=spa14*spb31+spa24*spb32;
const complex<T> spab_2_456_3=spa24*spb43+spa25*spb53+spa26*spb63;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb78)
)-
(pow3_spb36*spb13*pow(spab_7_12_3,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_3_12_78_6*spbb_3_456_78_1
)-
(pow2_spb68*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_3_12_678_5*spbb_3_12_78_6*spbb_3_45_678_1
)+
(pow3_spb36*spab_2_456_3*pow(spab_2_17_8,2))/(spab_2_345_6*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_3_456_78_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow2_spb68*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb78*spbb_3_12_678_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb68*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb78*spbb_3_45_678_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpmpmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=pow(spb46,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+SPA(i2,i6)*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_8=-(spa34*SPB(i8,i4))-spa35*SPB(i8,i5)-spa36*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,-1)*(-((pow2_spb68*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
)-
(pow3_spb46*spb14*pow(spab_7_123_4,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1
)-
(pow3_spb46*spbb_4_23_56_4*pow(spbb_4_23_17_8,2))/(spb23*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1
)-
(pow2_spb68*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_4_123_78_6*(s23+s24+s34+spa12*spb21+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow2_spb68*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb78*spbb_4_23_178_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb46*spab_3_56_4*pow(spab_3_456_8,2))/((s45+s46+s56)*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpmpmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=pow(spb46,3);
const complex<T> pow2_spb68=pow(spb68,2);
const complex<T> spbb_4_56_78_1=-(spb54*(spa57*spb71+spa58*spb81))-spb64*(spa67*spb71+spa68*spb81);
const complex<T> spbb_4_56_178_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82);
const complex<T> spbb_4_23_56_4=spb42*(spa25*spb54+SPA(i2,i6)*spb64)+spb43*(spa35*spb54+spa36*spb64);
const complex<T> spbb_4_23_17_8=spb42*(spa12*spb81-spa27*spb87)+spb43*(spa13*spb81-spa37*spb87);
const complex<T> spbb_4_23_178_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86);
const complex<T> spbb_4_123_78_6=(spa17*spb41+spa27*spb42+spa37*spb43)*spb76+(spa18*spb41+spa28*spb42+spa38*spb43)*spb86;
const complex<T> spab_7_123_4=spa17*spb41+spa27*spb42+spa37*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_8=-(spa34*SPB(i8,i4))-spa35*SPB(i8,i5)-spa36*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,1)*(-((pow2_spb68*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb78)
)-
(pow3_spb46*spb14*pow(spab_7_123_4,2))/(spa78*spb12*spb23*spb34*spb45*spb56*spbb_4_123_78_6*spbb_4_56_78_1
)-
(pow3_spb46*spbb_4_23_56_4*pow(spbb_4_23_17_8,2))/(spb23*spb34*spb45*spb56*spb78*(s78+spa17*spb71+spa18*spb81)*spbb_4_23_178_6*spbb_4_56_178_2*spbb_4_56_78_1
)-
(pow2_spb68*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb78*(s78+spa67*spb76+spa68*spb86)*spbb_4_123_78_6*(s23+s24+s34+spa12*spb21+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow2_spb68*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb78*spbb_4_23_178_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb46*spab_3_56_4*pow(spab_3_456_8,2))/((s45+s46+s56)*spab_3_45_6*spb12*spb45*spb56*spb78*spbb_4_56_178_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpmpppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spa78=SPA(i7,i8);
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
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*SPB(i6,i2)+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,-1)*(-((pow2_spa17*pow(spab_5_34_6,3))/((s34+s35+s45)*spa12*spa34*spa45*spa78*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6)))
)-
(pow2_spa17*pow(SPB(i4,i6),3))/((s45+s46+s56)*spa12*spa23*spa78*spab_3_45_6*SPB(i4,i5)*SPB(i5,i6)
)+
(pow(SPA(i1,i5),3)*pow(SPB(i6,i8),2))/(spa12*spa23*spa34*spa45*spab_1_78_6*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))*SPB(i7,i8)
)+
(pow2_spa17*pow(spab_5_234_6,3))/(spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_345_6*(s34+s35+s45+spa23*SPB(i3,i2
)+
spa24*SPB(i4,i2
)+
spa25*SPB(i5,i2))*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpmpppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spa78=SPA(i7,i8);
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
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_234_6=spa25*SPB(i6,i2)+spa35*spb63+spa45*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_2_345_6=-(spa23*spb63)-spa24*spb64-spa25*spb65;
const complex<T> spab_1_78_6=spa17*spb76+spa18*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,1)*(-((pow2_spa17*pow(spab_5_34_6,3))/((s34+s35+s45)*spa12*spa34*spa45*spa78*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6)))
)-
(pow2_spa17*pow(SPB(i4,i6),3))/((s45+s46+s56)*spa12*spa23*spa78*spab_3_45_6*SPB(i4,i5)*SPB(i5,i6)
)+
(pow(SPA(i1,i5),3)*pow(SPB(i6,i8),2))/(spa12*spa23*spa34*spa45*spab_1_78_6*(s78+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8))*SPB(i7,i8)
)+
(pow2_spa17*pow(spab_5_234_6,3))/(spa23*spa34*spa45*spa78*spab_1_78_6*spab_2_345_6*(s34+s35+s45+spa23*SPB(i3,i2
)+
spa24*SPB(i4,i2
)+
spa25*SPB(i5,i2))*(s78+spa17*SPB(i7,i1
)+
spa18*SPB(i8,i1))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qppmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(SPA(i1,i2)*spb52)-spa13*spb53-spa14*spb54)*spb81)-(SPA(i2,i7)*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=SPA(i2,i6)*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,-1)*(-((spa46*pow(spab_4_56_8,2))/((s45+s46+s56)*spa45*spa56*spab_6_45_3*spb12*spb23*spb78)
)+
(spab_6_234_5*pow(spbb_5_234_17_8,2))/(spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78*(s78+spb71*SPA(i1,i7
)+
spb81*SPA(i1,i8))*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(spab_6_34_5*pow(spbb_5_34_127_8,2))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow(spab_7_68_5,2)*SPB(i1,i5))/(spa78*spab_6_78_1*spb12*spb23*spb34*spb45*(s78+spa68*spb86+spa67*SPB(i7,i6))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qppmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb85=SPB(i8,i5);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb78=SPB(i7,i8);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spbb_5_34_127_8=spb53*(spa13*spb81+spa23*spb82-spa37*spb87)+spb54*(spa14*spb81+spa24*spb82-spa47*spb87);
const complex<T> spbb_5_234_17_8=-((-(SPA(i1,i2)*spb52)-spa13*spb53-spa14*spb54)*spb81)-(SPA(i2,i7)*spb52+spa37*spb53+spa47*spb54)*spb87;
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=SPA(i2,i6)*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;

return(
complex<T>(0,1)*(-((spa46*pow(spab_4_56_8,2))/((s45+s46+s56)*spa45*spa56*spab_6_45_3*spb12*spb23*spb78)
)+
(spab_6_234_5*pow(spbb_5_234_17_8,2))/(spab_6_345_2*spab_6_78_1*spb23*spb34*spb45*spb78*(s78+spb71*SPA(i1,i7
)+
spb81*SPA(i1,i8))*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(spab_6_34_5*pow(spbb_5_34_127_8,2))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb78*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow(spab_7_68_5,2)*SPB(i1,i5))/(spa78*spab_6_78_1*spb12*spb23*spb34*spb45*(s78+spa68*spb86+spa67*SPB(i7,i6))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qppmppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*SPB(i6,i2)))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,-1)*(-((pow2_spa17*spa46*pow(spaa_4_23_56_4,3))/(spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1*(s78+spa17*spb71+spa18*spb81))
)+
(pow2_spa17*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa17*spa46*pow(spab_4_56_3,3))/((s45+s46+s56)*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_7_68_5,2))/(spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_234_5*(s78+spa67*spb76+spa68*spb86)*(s23+s24+s34+spa12*spb21+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa17*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa78*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)-
(pow3_spa14*spa46*pow(spab_4_56_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qppmppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_68_5=-(spa67*spb65)+spa78*spb85;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_56_8=-(spa45*spb85)-spa46*spb86;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_78_1=spa45*(spa17*spb75+spa18*spb85)+spa46*(spa17*spb76+spa18*spb86);
const complex<T> spaa_4_56_178_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*SPB(i6,i2)))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_178_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83);
const complex<T> spaa_4_123_78_6=-(spa67*(spa14*spb71+spa24*spb72+spa34*spb73))-spa68*(spa14*spb81+spa24*spb82+spa34*spb83);
const complex<T> s78=spa78*spb87;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;

return(
complex<T>(0,1)*(-((pow2_spa17*spa46*pow(spaa_4_23_56_4,3))/(spa23*spa34*spa45*spa56*spa78*spaa_4_23_178_6*spaa_4_56_178_2*spaa_4_56_78_1*(s78+spa17*spb71+spa18*spb81))
)+
(pow2_spa17*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa78*spaa_4_23_178_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa17*spa46*pow(spab_4_56_3,3))/((s45+s46+s56)*spa12*spa45*spa56*spa78*spaa_4_56_178_2*spab_6_45_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_7_68_5,2))/(spa12*spa23*spa34*spa78*spaa_4_123_78_6*spab_1_234_5*(s78+spa67*spb76+spa68*spb86)*(s23+s24+s34+spa12*spb21+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa17*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa78*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)-
(pow3_spa14*spa46*pow(spab_4_56_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_4_123_78_6*spaa_4_56_78_1*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpppmpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_123_4=spa17*spb41+SPA(i2,i7)*spb42+SPA(i3,i7)*spb43;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*SPB(i5,i1)+spa23*spb52);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_68_7,2))/(spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1*(s78+spa67*spb76+spa68*spb86))
)-
(pow3_spa13*spab_3_12_4*pow(spab_7_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa17*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa17*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa78*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa17*spa36*pow(spab_3_456_2,3))/(spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_345_2*(s78+spa17*spb71+spa18*spb81)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa36*pow(spab_3_456_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpppmpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa17=pow(spa17,2);
const complex<T> spab_7_123_4=spa17*spb41+SPA(i2,i7)*spb42+SPA(i3,i7)*spb43;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_8=-(spa34*spb84)-spa35*spb85-spa36*spb86;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_68_7=spa34*(-(spa67*spb64)+spa78*spb84)+spa35*(-(spa67*spb65)+spa78*spb85);
const complex<T> spaa_3_45_678_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85);
const complex<T> spaa_3_456_78_1=-(spa17*(-(spa34*spb74)-spa35*spb75-spa36*spb76))-spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*SPB(i5,i1)+spa23*spb52);
const complex<T> spaa_3_12_78_6=-(spa13*(spa67*spb71+spa68*spb81))-spa23*(spa67*spb72+spa68*spb82);
const complex<T> spaa_3_12_678_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82);
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_68_7,2))/(spa12*spa23*spa34*spa45*spa78*spaa_3_12_678_5*spaa_3_12_78_6*spaa_3_45_678_1*(s78+spa67*spb76+spa68*spb86))
)-
(pow3_spa13*spab_3_12_4*pow(spab_7_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa78*spaa_3_12_678_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa17*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa78*spaa_3_45_678_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa17*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa78*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa17*spa36*pow(spab_3_456_2,3))/(spa34*spa45*spa56*spa78*spaa_3_456_78_1*spab_6_345_2*(s78+spa17*spb71+spa18*spb81)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa36*pow(spab_3_456_8,2))/(spa12*spa23*spa34*spa45*spa56*spaa_3_12_78_6*spaa_3_456_78_1*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qppppmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*SPB(i6,i5))-spa78*(-(spa23*spb83)-spa24*spb84-spa25*SPB(i8,i5));
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((spab_2_34_1*pow(spaa_2_34_568_7,2))/((s23+s24+s34)*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)-
(pow(spab_7_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa78*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(pow(spab_2_17_8,2)*SPA(i2,i6))/(spa23*spa34*spa45*spa56*spab_6_78_1*(s78+spa17*spb71+spb81*SPA(i1,i8))*SPB(i7,i8)
)+
(spab_2_345_1*pow(spaa_2_345_68_7,2))/(spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1*(s23+s24+s34+spa35*spb53+spa45*spb54+spa25*SPB(i5,i2))*(s78+spa67*SPB(i7,i6
)+
spa68*SPB(i8,i6))))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qppppmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb84=SPB(i8,i4);
const complex<T> spb83=SPB(i8,i3);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa68=SPA(i6,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_7_12_3=spa17*spb31+spa27*spb32;
const complex<T> spab_6_78_1=spa67*spb71+spa68*spb81;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_17_8=spa12*spb81-spa27*spb87;
const complex<T> spaa_2_345_68_7=spa67*(-(spa23*spb63)-spa24*spb64-spa25*SPB(i6,i5))-spa78*(-(spa23*spb83)-spa24*spb84-spa25*SPB(i8,i5));
const complex<T> spaa_2_34_568_7=spa23*(-(spa57*spb53)-spa67*spb63+spa78*spb83)+spa24*(-(spa57*spb54)-spa67*spb64+spa78*spb84);
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((spab_2_34_1*pow(spaa_2_34_568_7,2))/((s23+s24+s34)*spa23*spa34*spa56*spa78*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)-
(pow(spab_7_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa78*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(pow(spab_2_17_8,2)*SPA(i2,i6))/(spa23*spa34*spa45*spa56*spab_6_78_1*(s78+spa17*spb71+spb81*SPA(i1,i8))*SPB(i7,i8)
)+
(spab_2_345_1*pow(spaa_2_345_68_7,2))/(spa23*spa34*spa45*spa78*spab_5_234_1*spab_6_78_1*(s23+s24+s34+spa35*spb53+spa45*spb54+spa25*SPB(i5,i2))*(s78+spa67*SPB(i7,i6
)+
spa68*SPB(i8,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmppppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i7),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, class T> complex<T>  A2q4g2l_qmmmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPB(i6,i8),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpppppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i7),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmppppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i7),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i7, int i8, class T> complex<T>  A2q4g2l_qpmmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPB(i6,i8),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i7,i8)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i8, int i7, class T> complex<T>  A2q4g2l_qmmmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPB(i6,i8),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpppppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i7),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i7,i8)))
);
}


template <int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i7, class T> complex<T>  A2q4g2l_qpmmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPB(i6,i8),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i7,i8)))
);
}




template <class T> complex<T> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 1360802:	 return &A2q4g2l_qpmmmmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m m m m qbm lp lbm
case 1360820:	 return &A2q4g2l_qppmmmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p m m m qbm lp lbm
case 1360910:	 return &A2q4g2l_qpmpmmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m p m m qbm lp lbm
case 1361450:	 return &A2q4g2l_qpmmpmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m m p m qbm lp lbm
case 1361576:	 return &A2q4g2l_qppppmqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p p p m qbm lp lbm
case 1364690:	 return &A2q4g2l_qpmmmpqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m m m p qbm lp lbm
case 1364816:	 return &A2q4g2l_qpppmpqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p p m p qbm lp lbm
case 1365356:	 return &A2q4g2l_qppmppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p m p p qbm lp lbm
case 1365446:	 return &A2q4g2l_qpmpppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp m p p p qbm lp lbm
case 1365464:	 return &A2q4g2l_qpppppqbmlplbm_eval<0,1,2,3,4,5,6,7>;//qp p p p p qbm lp lbm
case 1368577:	 return &A2q4g2l_qmmmmmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m m m m qbp lp lbm
case 1368595:	 return &A2q4g2l_qmpmmmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p m m m qbp lp lbm
case 1368685:	 return &A2q4g2l_qmmpmmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m p m m qbp lp lbm
case 1369225:	 return &A2q4g2l_qmmmpmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m m p m qbp lp lbm
case 1369351:	 return &A2q4g2l_qmpppmqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p p p m qbp lp lbm
case 1372465:	 return &A2q4g2l_qmmmmpqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m m m p qbp lp lbm
case 1372591:	 return &A2q4g2l_qmppmpqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p p m p qbp lp lbm
case 1373131:	 return &A2q4g2l_qmpmppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p m p p qbp lp lbm
case 1373221:	 return &A2q4g2l_qmmpppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm m p p p qbp lp lbm
case 1373239:	 return &A2q4g2l_qmppppqbplplbm_eval<0,1,2,3,4,5,6,7>;//qm p p p p qbp lp lbm
case 1594082:	 return &A2q4g2l_qpmmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m m m m qbm lm lbp
case 1594100:	 return &A2q4g2l_qppmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p m m m qbm lm lbp
case 1594190:	 return &A2q4g2l_qpmpmmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m p m m qbm lm lbp
case 1594730:	 return &A2q4g2l_qpmmpmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m m p m qbm lm lbp
case 1594856:	 return &A2q4g2l_qppppmqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p p p m qbm lm lbp
case 1597970:	 return &A2q4g2l_qpmmmpqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m m m p qbm lm lbp
case 1598096:	 return &A2q4g2l_qpppmpqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p p m p qbm lm lbp
case 1598636:	 return &A2q4g2l_qppmppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p m p p qbm lm lbp
case 1598726:	 return &A2q4g2l_qpmpppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp m p p p qbm lm lbp
case 1598744:	 return &A2q4g2l_qpppppqbmlmlbp_eval<0,1,2,3,4,5,6,7>;//qp p p p p qbm lm lbp
case 1601857:	 return &A2q4g2l_qmmmmmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m m m m qbp lm lbp
case 1601875:	 return &A2q4g2l_qmpmmmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p m m m qbp lm lbp
case 1601965:	 return &A2q4g2l_qmmpmmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m p m m qbp lm lbp
case 1602505:	 return &A2q4g2l_qmmmpmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m m p m qbp lm lbp
case 1602631:	 return &A2q4g2l_qmpppmqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p p p m qbp lm lbp
case 1605745:	 return &A2q4g2l_qmmmmpqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m m m p qbp lm lbp
case 1605871:	 return &A2q4g2l_qmppmpqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p p m p qbp lm lbp
case 1606411:	 return &A2q4g2l_qmpmppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p m p p qbp lm lbp
case 1606501:	 return &A2q4g2l_qmmpppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm m p p p qbp lm lbp
case 1606519:	 return &A2q4g2l_qmppppqbplmlbp_eval<0,1,2,3,4,5,6,7>;//qm p p p p qbp lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q4g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2q4g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}
