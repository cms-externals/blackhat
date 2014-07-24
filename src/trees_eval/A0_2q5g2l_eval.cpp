/*
 * A0_2q5g2l_eval.cpp
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
 * The 2 quarks 5g 2 lepton amplitudes
 *
 *
 */


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmmmmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
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
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_4567_9=spb64*(-(spa45*SPB(i9,i5))-spa46*spb96-spa47*spb97)+spb65*(spa45*SPB(i9,i4)-spa56*spb96-spa57*spb97);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+SPA(i2,i3)*spb92-spa38*spb98)+spb64*(spa14*spb91+SPA(i2,i4)*spb92-spa48*spb98)+spb65*(spa15*spb91+SPA(i2,i5)*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_45_6=spa47*spb64+spa57*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_345_6=SPA(i3,i7)*spb63+spa47*spb64+spa57*spb65;
const complex<T> spab_7_189_6=spa17*SPB(i6,i1)+spa78*SPB(i8,i6)+spa79*spb96;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((spa57*pow(spab_5_67_9,2))/((s56+s57+s67)*spa56*spa67*spab_7_56_4*spb12*spb23*spb34*spb98)
)+
(spab_7_345_6*pow(spbb_6_345_128_9,2))/(spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb98*(s18+s19+s89+spa12*spb21+spa28*spb82+spb92*SPA(i2,i9))*(s45+s46+s56+spb43*SPA(i3,i4
)+
spb53*SPA(i3,i5
)+
spb63*SPA(i3,i6))
)+
(pow(spab_8_79_6,2)*SPB(i1,i6))/((s78+s79+s89)*spab_7_89_1*spb12*spb23*spb34*spb45*spb56*SPA(i9,i8)
)+
(spab_7_189_6*pow(spbb_6_2345_18_9,2))/((s18+s19+s89)*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))
)-
(spab_7_45_6*pow(spbb_6_45_4567_9,2))/((s45+s46+s56)*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb98*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmmmmpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
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
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_4567_9=spb64*(-(spa45*SPB(i9,i5))-spa46*spb96-spa47*spb97)+spb65*(spa45*SPB(i9,i4)-spa56*spb96-spa57*spb97);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+SPA(i2,i3)*spb92-spa38*spb98)+spb64*(spa14*spb91+SPA(i2,i4)*spb92-spa48*spb98)+spb65*(spa15*spb91+SPA(i2,i5)*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_45_6=spa47*spb64+spa57*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_345_6=SPA(i3,i7)*spb63+spa47*spb64+spa57*spb65;
const complex<T> spab_7_189_6=spa17*SPB(i6,i1)+spa78*SPB(i8,i6)+spa79*spb96;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((spa57*pow(spab_5_67_9,2))/((s56+s57+s67)*spa56*spa67*spab_7_56_4*spb12*spb23*spb34*spb98)
)+
(spab_7_345_6*pow(spbb_6_345_128_9,2))/(spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb98*(s18+s19+s89+spa12*spb21+spa28*spb82+spb92*SPA(i2,i9))*(s45+s46+s56+spb43*SPA(i3,i4
)+
spb53*SPA(i3,i5
)+
spb63*SPA(i3,i6))
)+
(pow(spab_8_79_6,2)*SPB(i1,i6))/((s78+s79+s89)*spab_7_89_1*spb12*spb23*spb34*spb45*spb56*SPA(i9,i8)
)+
(spab_7_189_6*pow(spbb_6_2345_18_9,2))/((s18+s19+s89)*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))
)-
(spab_7_45_6*pow(spbb_6_45_4567_9,2))/((s45+s46+s56)*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb98*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmmmpmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
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
const complex<T> pow3_spb57=pow(spb57,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_67_5=spb53*(spa36*spb65+spa37*spb75)+spb54*(spa46*spb65+spa47*spb75);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_67_5=(spa26*spb52+spa36*spb53+spa46*spb54)*spb65+(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*SPB(i8,i5)+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa46,4))/((s45+s46+s56)*spa45*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb98)
)-
(pow3_spb57*spbb_5_234_67_5*pow(spbb_5_234_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1
)-
(pow3_spb57*spbb_5_34_67_5*pow(spbb_5_34_128_9,2))/(spb12*spb34*spb45*spb56*spb67*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)*spb98*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2
)-
(pow2_spb79*pow(spab_6_234_5,4))/(spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_5_234_189_7*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(pow3_spb57*spb15*pow(spab_8_679_5,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_6_34_5,4))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb98*spbb_5_34_1289_7*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow3_spb57*spab_4_67_5*pow(spab_4_567_9,2))/((s56+s57+s67)*spab_4_56_7*spb12*spb23*spb56*spb67*spb98*spbb_5_67_1289_3*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))
)-
(pow2_spb79*spb15*pow(spab_6_789_5,3))/((s78+s79+s89)*spab_6_789_1*spb12*spb23*spb34*spb45*spb98*spbb_5_1234_89_7*(s67+s78+s79+s89+spa69*spb96+spa68*SPB(i8,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmmmpmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
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
const complex<T> pow3_spb57=pow(spb57,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_67_5=spb53*(spa36*spb65+spa37*spb75)+spb54*(spa46*spb65+spa47*spb75);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_67_5=(spa26*spb52+spa36*spb53+spa46*spb54)*spb65+(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*SPB(i8,i5)+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa46,4))/((s45+s46+s56)*spa45*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb98)
)-
(pow3_spb57*spbb_5_234_67_5*pow(spbb_5_234_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1
)-
(pow3_spb57*spbb_5_34_67_5*pow(spbb_5_34_128_9,2))/(spb12*spb34*spb45*spb56*spb67*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)*spb98*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2
)-
(pow2_spb79*pow(spab_6_234_5,4))/(spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_5_234_189_7*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(pow3_spb57*spb15*pow(spab_8_679_5,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_6_34_5,4))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb98*spbb_5_34_1289_7*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow3_spb57*spab_4_67_5*pow(spab_4_567_9,2))/((s56+s57+s67)*spab_4_56_7*spb12*spb23*spb56*spb67*spb98*spbb_5_67_1289_3*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))
)-
(pow2_spb79*spb15*pow(spab_6_789_5,3))/((s78+s79+s89)*spab_6_789_1*spb12*spb23*spb34*spb45*spb98*spbb_5_1234_89_7*(s67+s78+s79+s89+spa69*spb96+spa68*SPB(i8,i6))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmmpmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> pow3_spb47=pow(spb47,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_189_4=spb42*(spa12*spb41+spa28*spb84+spa29*spb94)+spb43*(spa13*spb41+spa38*spb84+spa39*spb94);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb98)
)+
(pow3_spb47*spbb_4_23_189_4*pow(spbb_4_23_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1
)-
(pow2_spb79*spb14*pow(spbb_4_56_123_4,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1
)-
(pow2_spb79*pow(spbb_4_56_23_4,4))/(spb23*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1
)-
(pow2_spb79*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb67*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*spb98*spbb_4_123_789_6*(s12+s23+s24+s34+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow3_spb47*spab_3_567_4*pow(spab_3_128_9,2))/(spab_3_456_7*spb12*spb45*spb56*spb67*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*spb98*spbb_4_567_189_2*(s45+s46+s56+s67+spa57*spb75+spb74*SPA(i4,i7))
)-
(pow3_spb47*spb14*pow(spab_8_123_4,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb98*spbb_4_23_1789_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)-
(pow2_spb79*pow(spab_3_56_4,4))/((s45+s46+s56)*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb98*spbb_4_56_1789_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmmpmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> pow3_spb47=pow(spb47,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_189_4=spb42*(spa12*spb41+spa28*spb84+spa29*spb94)+spb43*(spa13*spb41+spa38*spb84+spa39*spb94);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb98)
)+
(pow3_spb47*spbb_4_23_189_4*pow(spbb_4_23_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1
)-
(pow2_spb79*spb14*pow(spbb_4_56_123_4,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1
)-
(pow2_spb79*pow(spbb_4_56_23_4,4))/(spb23*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1
)-
(pow2_spb79*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb67*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*spb98*spbb_4_123_789_6*(s12+s23+s24+s34+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow3_spb47*spab_3_567_4*pow(spab_3_128_9,2))/(spab_3_456_7*spb12*spb45*spb56*spb67*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*spb98*spbb_4_567_189_2*(s45+s46+s56+s67+spa57*spb75+spb74*SPA(i4,i7))
)-
(pow3_spb47*spb14*pow(spab_8_123_4,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb98*spbb_4_23_1789_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)-
(pow2_spb79*pow(spab_3_56_4,4))/((s45+s46+s56)*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb98*spbb_4_56_1789_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmpmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> pow3_spb37=pow(spb37,3);
const complex<T> pow2_spb79=pow(spb79,2);
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
const complex<T> spab_2_189_3=spa12*spb31+spa28*SPB(i8,i3)+spa29*SPB(i9,i3);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb98)
)-
(pow2_spb79*spb13*pow(spbb_3_456_12_3,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1
)-
(pow2_spb79*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb67*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)*spb98*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1
)+
(pow2_spb79*pow(spab_2_456_3,4))/(spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_3_456_789_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow3_spb37*spb13*pow(spab_8_12_3,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1*SPA(i9,i8)
)-
(pow2_spb79*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb67*spb98*spbb_3_12_6789_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb79*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb98*spbb_3_45_6789_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb37*spab_2_189_3*pow(spab_2_18_9,2))/((s18+s19+s89)*spab_2_189_7*spb34*spb45*spb56*spb67*spb98*spbb_3_4567_89_1*(s12+s18+s19+s89+spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmpmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> pow3_spb37=pow(spb37,3);
const complex<T> pow2_spb79=pow(spb79,2);
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
const complex<T> spab_2_189_3=spa12*spb31+spa28*SPB(i8,i3)+spa29*SPB(i9,i3);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb98)
)-
(pow2_spb79*spb13*pow(spbb_3_456_12_3,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1
)-
(pow2_spb79*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb67*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)*spb98*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1
)+
(pow2_spb79*pow(spab_2_456_3,4))/(spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_3_456_789_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow3_spb37*spb13*pow(spab_8_12_3,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1*SPA(i9,i8)
)-
(pow2_spb79*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb67*spb98*spbb_3_12_6789_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb79*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb98*spbb_3_45_6789_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb37*spab_2_189_3*pow(spab_2_18_9,2))/((s18+s19+s89)*spab_2_189_7*spb34*spb45*spb56*spb67*spb98*spbb_3_4567_89_1*(s12+s18+s19+s89+spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmppppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
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
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_789_1=spa27*spb71+spa28*spb81+SPA(i2,i9)*spb91;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*SPB(i6,i3))-spa78*spb73+spa89*spb93)+spa24*(-(spa68*SPB(i6,i4))-spa78*spb74+spa89*spb94)+spa25*(-(spa68*SPB(i6,i5))-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_1234_8=spa24*(spa18*spb41+spa28*spb42+SPA(i3,i8)*spb43)+spa23*(spa18*spb31+spa28*spb32-SPA(i4,i8)*spb43);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((spab_2_34_1*pow(spaa_2_34_1234_8,2))/((s23+s24+s34)*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)+
(spab_2_789_1*pow(spaa_2_3456_79_8,2))/((s78+s79+s89)*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow(spab_8_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa67*spa89*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(spab_2_345_1*pow(spaa_2_345_679_8,2))/(spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1*(s23+s24+s34+spa25*SPB(i5,i2
)+
spa35*SPB(i5,i3
)+
spa45*SPB(i5,i4))*(s78+s79+s89+spa67*spb76+spa69*spb96+spa68*SPB(i8,i6))
)-
(spa27*pow(spab_2_18_9,2)*SPB(i1,i9))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb91*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmppppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
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
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_789_1=spa27*spb71+spa28*spb81+SPA(i2,i9)*spb91;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*SPB(i6,i3))-spa78*spb73+spa89*spb93)+spa24*(-(spa68*SPB(i6,i4))-spa78*spb74+spa89*spb94)+spa25*(-(spa68*SPB(i6,i5))-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_1234_8=spa24*(spa18*spb41+spa28*spb42+SPA(i3,i8)*spb43)+spa23*(spa18*spb31+spa28*spb32-SPA(i4,i8)*spb43);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((spab_2_34_1*pow(spaa_2_34_1234_8,2))/((s23+s24+s34)*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)+
(spab_2_789_1*pow(spaa_2_3456_79_8,2))/((s78+s79+s89)*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow(spab_8_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa67*spa89*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(spab_2_345_1*pow(spaa_2_345_679_8,2))/(spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1*(s23+s24+s34+spa25*SPB(i5,i2
)+
spa35*SPB(i5,i3
)+
spa45*SPB(i5,i4))*(s78+s79+s89+spa67*spb76+spa69*spb96+spa68*SPB(i8,i6))
)-
(spa27*pow(spab_2_18_9,2)*SPB(i1,i9))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb91*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmpmmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb67=SPB(i6,i7);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_2=spa17*SPB(i7,i2)+spa18*SPB(i8,i2)+spa19*SPB(i9,i2);
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb67*spb98)
)+
(pow2_spb79*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb98*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+s79+s89+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)+
(pow(spa18,2)*pow(SPB(i2,i7),3))/((s18+s19+s89)*spab_1_89_7*spb23*spb34*spb45*spb56*spb67*SPA(i9,i8)
)-
(pow2_spb79*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb98*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow2_spb79*pow(spab_1_789_2,3))/((s78+s79+s89)*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmpmmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb67=SPB(i6,i7);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_2=spa17*SPB(i7,i2)+spa18*SPB(i8,i2)+spa19*SPB(i9,i2);
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb67*spb98)
)+
(pow2_spb79*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb98*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+s79+s89+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)+
(pow(spa18,2)*pow(SPB(i2,i7),3))/((s18+s19+s89)*spab_1_89_7*spb23*spb34*spb45*spb56*spb67*SPA(i9,i8)
)-
(pow2_spb79*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb98*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow2_spb79*pow(spab_1_789_2,3))/((s78+s79+s89)*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmpmpppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+SPA(i3,i9)*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_456_12_3=spa13*(spa34*spb41+spa35*spb51+spa36*spb61)+spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*spb51+spa23*spb52);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow3_spa13*spaa_3_456_12_3*pow(spaa_3_456_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)
)-
(pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_679_8,2))/(spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)
)-
(pow3_spa13*spab_3_12_4*pow(spab_8_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa18*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa18*spa37*pow(spab_3_189_2,3))/((s18+s19+s89)*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2*(s12+s18+s19+s89+spa28*spb82+spb92*SPA(i2,i9))
)-
(pow2_spa18*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa18*pow(spab_3_456_2,4))/(spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa37*pow(spab_3_128_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmpmpppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+SPA(i3,i9)*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_456_12_3=spa13*(spa34*spb41+spa35*spb51+spa36*spb61)+spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*spb51+spa23*spb52);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow3_spa13*spaa_3_456_12_3*pow(spaa_3_456_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)
)-
(pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_679_8,2))/(spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)
)-
(pow3_spa13*spab_3_12_4*pow(spab_8_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa18*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa18*spa37*pow(spab_3_189_2,3))/((s18+s19+s89)*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2*(s12+s18+s19+s89+spa28*spb82+spb92*SPA(i2,i9))
)-
(pow2_spa18*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa18*pow(spab_3_456_2,4))/(spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa37*pow(spab_3_128_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmppmppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_4=spa45*(spa47*spb75+spa48*spb85+spa49*spb95)+spa46*(spa47*spb76+spa48*spb86+spa49*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spa18*spa47*pow(spaa_4_23_567_4,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1)
)-
(pow2_spa18*pow(spaa_4_23_56_4,4))/(spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1*(s18+s19+s78+s79+s89+spa17*spb71)
)+
(pow3_spa14*spaa_4_56_789_4*pow(spaa_4_56_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1
)+
(pow2_spa18*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa18*pow(spab_4_56_3,4))/((s45+s46+s56)*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_8_679_5,2))/(spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*(s12+s23+s24+s34+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa18*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)+
(pow2_spa18*spa47*pow(spab_4_567_3,3))/(spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*(s45+s46+s56+s67+spa57*spb75+spa47*SPB(i7,i4))
)-
(pow3_spa14*spa47*pow(spab_4_567_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmppmppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_4=spa45*(spa47*spb75+spa48*spb85+spa49*spb95)+spa46*(spa47*spb76+spa48*spb86+spa49*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spa18*spa47*pow(spaa_4_23_567_4,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1)
)-
(pow2_spa18*pow(spaa_4_23_56_4,4))/(spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1*(s18+s19+s78+s79+s89+spa17*spb71)
)+
(pow3_spa14*spaa_4_56_789_4*pow(spaa_4_56_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1
)+
(pow2_spa18*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa18*pow(spab_4_56_3,4))/((s45+s46+s56)*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_8_679_5,2))/(spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*(s12+s23+s24+s34+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa18*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)+
(pow2_spa18*spa47*pow(spab_4_567_3,3))/(spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*(s45+s46+s56+s67+spa57*spb75+spa47*SPB(i7,i4))
)-
(pow3_spa14*spa47*pow(spab_4_567_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmpppmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
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
const complex<T> pow3_spa15=pow(spa15,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+SPA(i5,i8)*spb86+SPA(i5,i9)*spb96;
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
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((pow2_spa18*spa57*pow(spaa_5_234_67_5,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
)-
(pow2_spa18*spa57*pow(spaa_5_34_67_5,3))/(spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)
)+
(pow2_spa18*pow(spab_5_34_6,4))/((s34+s35+s45)*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)+
(pow2_spa18*spa57*pow(spab_5_67_4,3))/((s56+s57+s67)*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow3_spa15*spab_5_789_6*pow(spab_8_79_6,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6*(s67+s78+s79+s89+spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)-
(pow2_spa18*pow(spab_5_234_6,4))/(spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa23*spb32+spa24*spb42+spa25*SPB(i5,i2))
)-
(pow2_spa18*pow(SPB(i4,i6),4))/((s45+s46+s56)*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*SPB(i4,i5)*SPB(i5,i6)
)-
(pow3_spa15*spa57*pow(spab_5_67_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmpppmpqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
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
const complex<T> pow3_spa15=pow(spa15,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+SPA(i5,i8)*spb86+SPA(i5,i9)*spb96;
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
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((pow2_spa18*spa57*pow(spaa_5_234_67_5,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
)-
(pow2_spa18*spa57*pow(spaa_5_34_67_5,3))/(spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)
)+
(pow2_spa18*pow(spab_5_34_6,4))/((s34+s35+s45)*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)+
(pow2_spa18*spa57*pow(spab_5_67_4,3))/((s56+s57+s67)*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow3_spa15*spab_5_789_6*pow(spab_8_79_6,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6*(s67+s78+s79+s89+spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)-
(pow2_spa18*pow(spab_5_234_6,4))/(spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa23*spb32+spa24*spb42+spa25*SPB(i5,i2))
)-
(pow2_spa18*pow(SPB(i4,i6),4))/((s45+s46+s56)*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*SPB(i4,i5)*SPB(i5,i6)
)-
(pow3_spa15*spa57*pow(spab_5_67_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmppppmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
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
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*SPB(i7,i3)+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+SPA(i6,i8)*spb87+SPA(i6,i9)*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*((pow2_spa18*pow(spab_6_189_7,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow2_spa18*pow(spab_6_45_7,3))/((s45+s46+s56)*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow2_spa18*pow(SPB(i5,i7),3))/((s56+s57+s67)*spa12*spa23*spa34*spa89*spab_4_56_7*SPB(i5,i6)*SPB(i6,i7)
)+
(pow(spa16,3)*pow(SPB(i7,i9),2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*SPB(i8,i9)
)+
(pow2_spa18*pow(spab_6_345_7,3))/(spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7*(s45+s46+s56+spa34*SPB(i4,i3
)+
spa35*SPB(i5,i3
)+
spa36*SPB(i6,i3))*(s18+s19+s89+spa12*SPB(i2,i1
)+
spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmppppmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
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
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*SPB(i7,i3)+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+SPA(i6,i8)*spb87+SPA(i6,i9)*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*((pow2_spa18*pow(spab_6_189_7,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow2_spa18*pow(spab_6_45_7,3))/((s45+s46+s56)*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow2_spa18*pow(SPB(i5,i7),3))/((s56+s57+s67)*spa12*spa23*spa34*spa89*spab_4_56_7*SPB(i5,i6)*SPB(i6,i7)
)+
(pow(spa16,3)*pow(SPB(i7,i9),2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*SPB(i8,i9)
)+
(pow2_spa18*pow(spab_6_345_7,3))/(spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7*(s45+s46+s56+spa34*SPB(i4,i3
)+
spa35*SPB(i5,i3
)+
spa36*SPB(i6,i3))*(s18+s19+s89+spa12*SPB(i2,i1
)+
spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmmmmpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb67=SPB(i6,i7);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_2=spa17*SPB(i7,i2)+spa18*SPB(i8,i2)+spa19*SPB(i9,i2);
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb67*spb98)
)+
(pow2_spb79*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb98*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+s79+s89+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)+
(pow(spa18,2)*pow(SPB(i2,i7),3))/((s18+s19+s89)*spab_1_89_7*spb23*spb34*spb45*spb56*spb67*SPA(i9,i8)
)-
(pow2_spb79*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb98*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow2_spb79*pow(spab_1_789_2,3))/((s78+s79+s89)*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmmmmpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb79=SPB(i7,i9);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb67=SPB(i6,i7);
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
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> spab_1_789_6=spa17*spb76+spa18*spb86+spa19*spb96;
const complex<T> spab_1_789_2=spa17*SPB(i7,i2)+spa18*SPB(i8,i2)+spa19*SPB(i9,i2);
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+SPA(i1,i5)*spb52;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa13,3))/((s12+s13+s23)*spa12*spa23*spab_1_23_4*spb45*spb56*spb67*spb98)
)+
(pow2_spb79*pow(spab_1_345_2,3))/(spab_1_234_5*spab_1_789_6*spb23*spb34*spb45*spb67*spb98*(s23+s24+s34+spb52*SPA(i2,i5
)+
spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*(s78+s79+s89+spb76*SPA(i6,i7
)+
spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)+
(pow(spa18,2)*pow(SPB(i2,i7),3))/((s18+s19+s89)*spab_1_89_7*spb23*spb34*spb45*spb56*spb67*SPA(i9,i8)
)-
(pow2_spb79*pow(spab_1_34_2,3))/((s23+s24+s34)*spab_1_234_5*spab_1_23_4*spb23*spb34*spb56*spb67*spb98*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)+
(pow2_spb79*pow(spab_1_789_2,3))/((s78+s79+s89)*spab_1_789_6*spab_1_89_7*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmmmpmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> pow3_spb37=pow(spb37,3);
const complex<T> pow2_spb79=pow(spb79,2);
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
const complex<T> spab_2_189_3=spa12*spb31+spa28*SPB(i8,i3)+spa29*SPB(i9,i3);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb98)
)-
(pow2_spb79*spb13*pow(spbb_3_456_12_3,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1
)-
(pow2_spb79*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb67*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)*spb98*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1
)+
(pow2_spb79*pow(spab_2_456_3,4))/(spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_3_456_789_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow3_spb37*spb13*pow(spab_8_12_3,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1*SPA(i9,i8)
)-
(pow2_spb79*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb67*spb98*spbb_3_12_6789_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb79*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb98*spbb_3_45_6789_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb37*spab_2_189_3*pow(spab_2_18_9,2))/((s18+s19+s89)*spab_2_189_7*spb34*spb45*spb56*spb67*spb98*spbb_3_4567_89_1*(s12+s18+s19+s89+spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmmmpmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> pow3_spb37=pow(spb37,3);
const complex<T> pow2_spb79=pow(spb79,2);
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
const complex<T> spab_2_189_3=spa12*spb31+spa28*SPB(i8,i3)+spa29*SPB(i9,i3);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa24,4))/((s23+s24+s34)*spa23*spa34*spab_2_34_5*spab_4_23_1*spb56*spb67*spb98)
)-
(pow2_spb79*spb13*pow(spbb_3_456_12_3,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_3_12_789_6*spbb_3_12_89_7*spbb_3_456_789_1
)-
(pow2_spb79*spb13*pow(spbb_3_45_12_3,3))/(spb12*spb23*spb34*spb45*spb67*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)*spb98*spbb_3_12_6789_5*spbb_3_12_789_6*spbb_3_45_6789_1
)+
(pow2_spb79*pow(spab_2_456_3,4))/(spab_2_189_7*spab_2_345_6*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_3_456_789_1*(s34+s35+s45+spa46*spb64+spa56*spb65+spb63*SPA(i3,i6))
)-
(pow3_spb37*spb13*pow(spab_8_12_3,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_3_12_89_7*spbb_3_4567_89_1*SPA(i9,i8)
)-
(pow2_spb79*spb13*pow(spab_4_12_3,3))/((s12+s13+s23)*spab_4_23_1*spb12*spb23*spb56*spb67*spb98*spbb_3_12_6789_5*(s12+s13+s23+s24+s34+spa14*SPB(i4,i1))
)-
(pow2_spb79*pow(spab_2_45_3,4))/((s34+s35+s45)*spab_2_345_6*spab_2_34_5*spb34*spb45*spb67*spb98*spbb_3_45_6789_1*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)+
(pow3_spb37*spab_2_189_3*pow(spab_2_18_9,2))/((s18+s19+s89)*spab_2_189_7*spb34*spb45*spb56*spb67*spb98*spbb_3_4567_89_1*(s12+s18+s19+s89+spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmmpmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> pow3_spb47=pow(spb47,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_189_4=spb42*(spa12*spb41+spa28*spb84+spa29*spb94)+spb43*(spa13*spb41+spa38*spb84+spa39*spb94);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb98)
)+
(pow3_spb47*spbb_4_23_189_4*pow(spbb_4_23_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1
)-
(pow2_spb79*spb14*pow(spbb_4_56_123_4,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1
)-
(pow2_spb79*pow(spbb_4_56_23_4,4))/(spb23*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1
)-
(pow2_spb79*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb67*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*spb98*spbb_4_123_789_6*(s12+s23+s24+s34+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow3_spb47*spab_3_567_4*pow(spab_3_128_9,2))/(spab_3_456_7*spb12*spb45*spb56*spb67*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*spb98*spbb_4_567_189_2*(s45+s46+s56+s67+spa57*spb75+spb74*SPA(i4,i7))
)-
(pow3_spb47*spb14*pow(spab_8_123_4,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb98*spbb_4_23_1789_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)-
(pow2_spb79*pow(spab_3_56_4,4))/((s45+s46+s56)*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb98*spbb_4_56_1789_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmmpmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb86=SPB(i8,i6);
const complex<T> spb84=SPB(i8,i4);
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
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
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
const complex<T> pow3_spb47=pow(spb47,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_4_56_789_1=-(spb54*(spa57*spb71+spa58*spb81+spa59*spb91))-spb64*(spa67*spb71+spa68*spb81+spa69*spb91);
const complex<T> spbb_4_567_89_1=(-(spa58*spb54)-spa68*spb64-spa78*spb74)*spb81+(-(spa59*spb54)-spa69*spb64-spa79*spb74)*spb91;
const complex<T> spbb_4_567_189_2=-(spb54*(spa15*spb21+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa68*spb82+spa69*spb92)-spb74*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_4_56_23_4=-((spa25*spb42+spa35*spb43)*spb54)-(spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_56_1789_2=-(spb54*(spa15*spb21+spa57*spb72+spa58*spb82+spa59*spb92))-spb64*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spbb_4_56_123_4=-((spa15*spb41+spa25*spb42+spa35*spb43)*spb54)-(spa16*spb41+spa26*spb42+spa36*spb43)*spb64;
const complex<T> spbb_4_23_18_9=spb42*(spa12*spb91-spa28*spb98)+spb43*(spa13*spb91-spa38*spb98);
const complex<T> spbb_4_23_189_7=spb42*(spa12*spb71+spa28*spb87+spa29*spb97)+spb43*(spa13*spb71+spa38*spb87+spa39*spb97);
const complex<T> spbb_4_23_189_4=spb42*(spa12*spb41+spa28*spb84+spa29*spb94)+spb43*(spa13*spb41+spa38*spb84+spa39*spb94);
const complex<T> spbb_4_23_1789_6=spb42*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa13*spb61+spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spbb_4_123_89_7=(spa18*spb41+spa28*spb42+spa38*spb43)*spb87+(spa19*spb41+spa29*spb42+spa39*spb43)*spb97;
const complex<T> spbb_4_123_789_6=spb41*(spa17*spb76+spa18*spb86+spa19*spb96)+spb42*(spa27*spb76+spa28*spb86+spa29*spb96)+spb43*(spa37*spb76+spa38*spb86+spa39*spb96);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_5_123_4=spa15*spb41+spa25*spb42+spa35*spb43;
const complex<T> spab_3_567_4=spa35*spb54+spa36*spb64+spa37*spb74;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_45_6=-(spa34*spb64)-spa35*spb65;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa35,4))/((s34+s35+s45)*spa34*spa45*spab_3_45_6*spab_5_34_2*spb12*spb67*spb98)
)+
(pow3_spb47*spbb_4_23_189_4*pow(spbb_4_23_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_4_23_189_7*spbb_4_567_189_2*spbb_4_567_89_1
)-
(pow2_spb79*spb14*pow(spbb_4_56_123_4,3))/((s78+s79+s89)*spb12*spb23*spb34*spb45*spb56*spb98*spbb_4_123_789_6*spbb_4_123_89_7*spbb_4_56_789_1
)-
(pow2_spb79*pow(spbb_4_56_23_4,4))/(spb23*spb34*spb45*spb56*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_4_23_1789_6*spbb_4_23_189_7*spbb_4_56_1789_2*spbb_4_56_789_1
)-
(pow2_spb79*spb14*pow(spab_5_123_4,3))/(spab_5_234_1*spb12*spb23*spb34*spb67*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*spb98*spbb_4_123_789_6*(s12+s23+s24+s34+spa13*spb31+spb41*SPA(i1,i4))
)+
(pow3_spb47*spab_3_567_4*pow(spab_3_128_9,2))/(spab_3_456_7*spb12*spb45*spb56*spb67*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*spb98*spbb_4_567_189_2*(s45+s46+s56+s67+spa57*spb75+spb74*SPA(i4,i7))
)-
(pow3_spb47*spb14*pow(spab_8_123_4,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_4_123_89_7*spbb_4_567_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_5_23_4,4))/((s23+s24+s34)*spab_5_234_1*spab_5_34_2*spb23*spb34*spb67*spb98*spbb_4_23_1789_6*(s23+s24+s34+s35+s45+spa25*SPB(i5,i2))
)-
(pow2_spb79*pow(spab_3_56_4,4))/((s45+s46+s56)*spab_3_456_7*spab_3_45_6*spb12*spb45*spb56*spb98*spbb_4_56_1789_2*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmpmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
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
const complex<T> pow3_spb57=pow(spb57,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_67_5=spb53*(spa36*spb65+spa37*spb75)+spb54*(spa46*spb65+spa47*spb75);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_67_5=(spa26*spb52+spa36*spb53+spa46*spb54)*spb65+(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*SPB(i8,i5)+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((pow2_spb79*pow(spa46,4))/((s45+s46+s56)*spa45*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb98)
)-
(pow3_spb57*spbb_5_234_67_5*pow(spbb_5_234_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1
)-
(pow3_spb57*spbb_5_34_67_5*pow(spbb_5_34_128_9,2))/(spb12*spb34*spb45*spb56*spb67*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)*spb98*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2
)-
(pow2_spb79*pow(spab_6_234_5,4))/(spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_5_234_189_7*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(pow3_spb57*spb15*pow(spab_8_679_5,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_6_34_5,4))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb98*spbb_5_34_1289_7*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow3_spb57*spab_4_67_5*pow(spab_4_567_9,2))/((s56+s57+s67)*spab_4_56_7*spb12*spb23*spb56*spb67*spb98*spbb_5_67_1289_3*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))
)-
(pow2_spb79*spb15*pow(spab_6_789_5,3))/((s78+s79+s89)*spab_6_789_1*spb12*spb23*spb34*spb45*spb98*spbb_5_1234_89_7*(s67+s78+s79+s89+spa69*spb96+spa68*SPB(i8,i6))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmpmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
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
const complex<T> pow3_spb57=pow(spb57,3);
const complex<T> pow2_spb79=pow(spb79,2);
const complex<T> spbb_5_67_89_1=-(spb65*(spa68*spb81+spa69*spb91))-spb75*(spa78*spb81+spa79*spb91);
const complex<T> spbb_5_67_189_2=-(spb65*(spa16*spb21+spa68*spb82+spa69*spb92))-spb75*(spa17*spb21+spa78*spb82+spa79*spb92);
const complex<T> spbb_5_67_1289_3=-(spb65*(spa16*spb31+spa26*spb32+spa68*spb83+spa69*spb93))-spb75*(spa17*spb31+spa27*spb32+spa78*spb83+spa79*spb93);
const complex<T> spbb_5_34_67_5=spb53*(spa36*spb65+spa37*spb75)+spb54*(spa46*spb65+spa47*spb75);
const complex<T> spbb_5_34_128_9=spb53*(spa13*spb91+spa23*spb92-spa38*spb98)+spb54*(spa14*spb91+spa24*spb92-spa48*spb98);
const complex<T> spbb_5_34_1289_7=spb53*(spa13*spb71+spa23*spb72+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa24*spb72+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_234_67_5=(spa26*spb52+spa36*spb53+spa46*spb54)*spb65+(spa27*spb52+spa37*spb53+spa47*spb54)*spb75;
const complex<T> spbb_5_234_18_9=-((-(spa12*spb52)-spa13*spb53-spa14*spb54)*spb91)-(spa28*spb52+spa38*spb53+spa48*spb54)*spb98;
const complex<T> spbb_5_234_189_7=spb52*(spa12*spb71+spa28*spb87+spa29*spb97)+spb53*(spa13*spb71+spa38*spb87+spa39*spb97)+spb54*(spa14*spb71+spa48*spb87+spa49*spb97);
const complex<T> spbb_5_1234_89_7=(spa18*spb51+spa28*spb52+spa38*spb53+spa48*spb54)*spb87+(spa19*spb51+spa29*spb52+spa39*spb53+spa49*spb54)*spb97;
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_6_789_5=spa67*spb75+spa68*SPB(i8,i5)+spa69*spb95;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_34_5=spa36*spb53+spa46*spb54;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_6_234_5=spa26*spb52+spa36*spb53+spa46*spb54;
const complex<T> spab_4_67_5=spa46*spb65+spa47*spb75;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((pow2_spb79*pow(spa46,4))/((s45+s46+s56)*spa45*spa56*spab_4_56_7*spab_6_45_3*spb12*spb23*spb98)
)-
(pow3_spb57*spbb_5_234_67_5*pow(spbb_5_234_18_9,2))/((s18+s19+s89)*spb23*spb34*spb45*spb56*spb67*spb98*spbb_5_234_189_7*spbb_5_67_189_2*spbb_5_67_89_1
)-
(pow3_spb57*spbb_5_34_67_5*pow(spbb_5_34_128_9,2))/(spb12*spb34*spb45*spb56*spb67*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)*spb98*spbb_5_34_1289_7*spbb_5_67_1289_3*spbb_5_67_189_2
)-
(pow2_spb79*pow(spab_6_234_5,4))/(spab_6_345_2*spab_6_789_1*spb23*spb34*spb45*(s18+s19+s78+s79+s89+spa17*spb71)*spb98*spbb_5_234_189_7*(s34+s35+s45+spa23*spb32+spa24*spb42+spb52*SPA(i2,i5))
)-
(pow3_spb57*spb15*pow(spab_8_679_5,2))/(spb12*spb23*spb34*spb45*spb56*spb67*spbb_5_1234_89_7*spbb_5_67_89_1*SPA(i9,i8)
)+
(pow2_spb79*pow(spab_6_34_5,4))/((s34+s35+s45)*spab_6_345_2*spab_6_45_3*spb12*spb34*spb45*spb98*spbb_5_34_1289_7*(s34+s35+s45+s46+s56+spa36*SPB(i6,i3))
)+
(pow3_spb57*spab_4_67_5*pow(spab_4_567_9,2))/((s56+s57+s67)*spab_4_56_7*spb12*spb23*spb56*spb67*spb98*spbb_5_67_1289_3*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))
)-
(pow2_spb79*spb15*pow(spab_6_789_5,3))/((s78+s79+s89)*spab_6_789_1*spb12*spb23*spb34*spb45*spb98*spbb_5_1234_89_7*(s67+s78+s79+s89+spa69*spb96+spa68*SPB(i8,i6))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmppppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
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
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*SPB(i7,i3)+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+SPA(i6,i8)*spb87+SPA(i6,i9)*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*((pow2_spa18*pow(spab_6_189_7,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow2_spa18*pow(spab_6_45_7,3))/((s45+s46+s56)*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow2_spa18*pow(SPB(i5,i7),3))/((s56+s57+s67)*spa12*spa23*spa34*spa89*spab_4_56_7*SPB(i5,i6)*SPB(i6,i7)
)+
(pow(spa16,3)*pow(SPB(i7,i9),2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*SPB(i8,i9)
)+
(pow2_spa18*pow(spab_6_345_7,3))/(spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7*(s45+s46+s56+spa34*SPB(i4,i3
)+
spa35*SPB(i5,i3
)+
spa36*SPB(i6,i3))*(s18+s19+s89+spa12*SPB(i2,i1
)+
spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmppppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spa89=SPA(i8,i9);
const complex<T> spa79=SPA(i7,i9);
const complex<T> spa78=SPA(i7,i8);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
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
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_6_45_7=spa46*spb74+spa56*spb75;
const complex<T> spab_6_345_7=spa36*SPB(i7,i3)+spa46*spb74+spa56*spb75;
const complex<T> spab_6_189_7=spa16*spb71+SPA(i6,i8)*spb87+SPA(i6,i9)*spb97;
const complex<T> spab_4_56_7=-(spa45*spb75)-spa46*spb76;
const complex<T> spab_3_456_7=-(spa34*spb74)-spa35*spb75-spa36*spb76;
const complex<T> spab_2_189_7=spa12*spb71+spa28*spb87+spa29*spb97;
const complex<T> spab_1_89_7=spa18*spb87+spa19*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*((pow2_spa18*pow(spab_6_189_7,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa89*spab_1_89_7*spab_2_189_7*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow2_spa18*pow(spab_6_45_7,3))/((s45+s46+s56)*spa12*spa23*spa45*spa56*spa89*spab_3_456_7*spab_4_56_7*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow2_spa18*pow(SPB(i5,i7),3))/((s56+s57+s67)*spa12*spa23*spa34*spa89*spab_4_56_7*SPB(i5,i6)*SPB(i6,i7)
)+
(pow(spa16,3)*pow(SPB(i7,i9),2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spab_1_89_7*SPB(i8,i9)
)+
(pow2_spa18*pow(spab_6_345_7,3))/(spa12*spa34*spa45*spa56*spa89*spab_2_189_7*spab_3_456_7*(s45+s46+s56+spa34*SPB(i4,i3
)+
spa35*SPB(i5,i3
)+
spa36*SPB(i6,i3))*(s18+s19+s89+spa12*SPB(i2,i1
)+
spa28*SPB(i8,i2
)+
spa29*SPB(i9,i2))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qppmmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
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
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_4567_9=spb64*(-(spa45*SPB(i9,i5))-spa46*spb96-spa47*spb97)+spb65*(spa45*SPB(i9,i4)-spa56*spb96-spa57*spb97);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+SPA(i2,i3)*spb92-spa38*spb98)+spb64*(spa14*spb91+SPA(i2,i4)*spb92-spa48*spb98)+spb65*(spa15*spb91+SPA(i2,i5)*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_45_6=spa47*spb64+spa57*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_345_6=SPA(i3,i7)*spb63+spa47*spb64+spa57*spb65;
const complex<T> spab_7_189_6=spa17*SPB(i6,i1)+spa78*SPB(i8,i6)+spa79*spb96;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((spa57*pow(spab_5_67_9,2))/((s56+s57+s67)*spa56*spa67*spab_7_56_4*spb12*spb23*spb34*spb98)
)+
(spab_7_345_6*pow(spbb_6_345_128_9,2))/(spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb98*(s18+s19+s89+spa12*spb21+spa28*spb82+spb92*SPA(i2,i9))*(s45+s46+s56+spb43*SPA(i3,i4
)+
spb53*SPA(i3,i5
)+
spb63*SPA(i3,i6))
)+
(pow(spab_8_79_6,2)*SPB(i1,i6))/((s78+s79+s89)*spab_7_89_1*spb12*spb23*spb34*spb45*spb56*SPA(i9,i8)
)+
(spab_7_189_6*pow(spbb_6_2345_18_9,2))/((s18+s19+s89)*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))
)-
(spab_7_45_6*pow(spbb_6_45_4567_9,2))/((s45+s46+s56)*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb98*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qppmmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb82=SPB(i8,i2);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
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
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
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
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spbb_6_45_4567_9=spb64*(-(spa45*SPB(i9,i5))-spa46*spb96-spa47*spb97)+spb65*(spa45*SPB(i9,i4)-spa56*spb96-spa57*spb97);
const complex<T> spbb_6_345_128_9=spb63*(spa13*spb91+SPA(i2,i3)*spb92-spa38*spb98)+spb64*(spa14*spb91+SPA(i2,i4)*spb92-spa48*spb98)+spb65*(spa15*spb91+SPA(i2,i5)*spb92-spa58*spb98);
const complex<T> spbb_6_2345_18_9=-((-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65)*spb91)-(spa28*spb62+spa38*spb63+spa48*spb64+spa58*spb65)*spb98;
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_7_45_6=spa47*spb64+spa57*spb65;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_7_345_6=SPA(i3,i7)*spb63+spa47*spb64+spa57*spb65;
const complex<T> spab_7_189_6=spa17*SPB(i6,i1)+spa78*SPB(i8,i6)+spa79*spb96;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_5_67_9=-(spa56*spb96)-spa57*spb97;
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((spa57*pow(spab_5_67_9,2))/((s56+s57+s67)*spa56*spa67*spab_7_56_4*spb12*spb23*spb34*spb98)
)+
(spab_7_345_6*pow(spbb_6_345_128_9,2))/(spab_7_189_2*spab_7_456_3*spb12*spb34*spb45*spb56*spb98*(s18+s19+s89+spa12*spb21+spa28*spb82+spb92*SPA(i2,i9))*(s45+s46+s56+spb43*SPA(i3,i4
)+
spb53*SPA(i3,i5
)+
spb63*SPA(i3,i6))
)+
(pow(spab_8_79_6,2)*SPB(i1,i6))/((s78+s79+s89)*spab_7_89_1*spb12*spb23*spb34*spb45*spb56*SPA(i9,i8)
)+
(spab_7_189_6*pow(spbb_6_2345_18_9,2))/((s18+s19+s89)*spab_7_189_2*spab_7_89_1*spb23*spb34*spb45*spb56*spb98*(s18+s19+s78+s79+s89+spa17*SPB(i7,i1))
)-
(spab_7_45_6*pow(spbb_6_45_4567_9,2))/((s45+s46+s56)*spab_7_456_3*spab_7_56_4*spb12*spb23*spb45*spb56*spb98*(s45+s46+s56+s57+s67+spa47*SPB(i7,i4))))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qppmpppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
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
const complex<T> pow3_spa15=pow(spa15,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+SPA(i5,i8)*spb86+SPA(i5,i9)*spb96;
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
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,1)*(-((pow2_spa18*spa57*pow(spaa_5_234_67_5,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
)-
(pow2_spa18*spa57*pow(spaa_5_34_67_5,3))/(spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)
)+
(pow2_spa18*pow(spab_5_34_6,4))/((s34+s35+s45)*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)+
(pow2_spa18*spa57*pow(spab_5_67_4,3))/((s56+s57+s67)*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow3_spa15*spab_5_789_6*pow(spab_8_79_6,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6*(s67+s78+s79+s89+spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)-
(pow2_spa18*pow(spab_5_234_6,4))/(spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa23*spb32+spa24*spb42+spa25*SPB(i5,i2))
)-
(pow2_spa18*pow(SPB(i4,i6),4))/((s45+s46+s56)*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*SPB(i4,i5)*SPB(i5,i6)
)-
(pow3_spa15*spa57*pow(spab_5_67_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qppmpppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa39=SPA(i3,i9);
const complex<T> spa38=SPA(i3,i8);
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
const complex<T> pow3_spa15=pow(spa15,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_79_6=-(spa78*spb76)+spa89*spb96;
const complex<T> spab_7_56_4=-(spa57*spb54)-spa67*spb64;
const complex<T> spab_5_789_6=spa57*spb76+SPA(i5,i8)*spb86+SPA(i5,i9)*spb96;
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
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s57=spa57*spb75;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;

return(
complex<T>(0,-1)*(-((pow2_spa18*spa57*pow(spaa_5_234_67_5,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_5_234_189_7*spaa_5_67_189_2*spaa_5_67_89_1)
)-
(pow2_spa18*spa57*pow(spaa_5_34_67_5,3))/(spa12*spa34*spa45*spa56*spa67*spa89*spaa_5_34_1289_7*spaa_5_67_1289_3*spaa_5_67_189_2*(s18+s19+s89+spa12*spb21+spa28*spb82+spa29*spb92)
)+
(pow2_spa18*pow(spab_5_34_6,4))/((s34+s35+s45)*spa12*spa34*spa45*spa89*spaa_5_34_1289_7*spab_2_345_6*spab_3_45_6*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)+
(pow2_spa18*spa57*pow(spab_5_67_4,3))/((s56+s57+s67)*spa12*spa23*spa56*spa67*spa89*spaa_5_67_1289_3*spab_7_56_4*(s45+s46+s56+s57+s67+spb74*SPA(i4,i7))
)-
(pow3_spa15*spab_5_789_6*pow(spab_8_79_6,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa89*spaa_5_1234_89_7*spab_1_789_6*(s67+s78+s79+s89+spb86*SPA(i6,i8
)+
spb96*SPA(i6,i9))
)-
(pow2_spa18*pow(spab_5_234_6,4))/(spa23*spa34*spa45*spa89*spaa_5_234_189_7*spab_1_789_6*spab_2_345_6*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa23*spb32+spa24*spb42+spa25*SPB(i5,i2))
)-
(pow2_spa18*pow(SPB(i4,i6),4))/((s45+s46+s56)*spa12*spa23*spa89*spab_3_45_6*spab_7_56_4*SPB(i4,i5)*SPB(i5,i6)
)-
(pow3_spa15*spa57*pow(spab_5_67_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_5_1234_89_7*spaa_5_67_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpppmppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_4=spa45*(spa47*spb75+spa48*spb85+spa49*spb95)+spa46*(spa47*spb76+spa48*spb86+spa49*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow2_spa18*spa47*pow(spaa_4_23_567_4,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1)
)-
(pow2_spa18*pow(spaa_4_23_56_4,4))/(spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1*(s18+s19+s78+s79+s89+spa17*spb71)
)+
(pow3_spa14*spaa_4_56_789_4*pow(spaa_4_56_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1
)+
(pow2_spa18*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa18*pow(spab_4_56_3,4))/((s45+s46+s56)*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_8_679_5,2))/(spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*(s12+s23+s24+s34+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa18*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)+
(pow2_spa18*spa47*pow(spab_4_567_3,3))/(spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*(s45+s46+s56+s67+spa57*spb75+spa47*SPB(i7,i4))
)-
(pow3_spa14*spa47*pow(spab_4_567_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpppmppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa14=pow(spa14,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_679_5=-(spa68*spb65)-spa78*spb75+spa89*spb95;
const complex<T> spab_7_456_3=-(spa47*spb43)-spa57*spb53-spa67*spb63;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_4_567_9=-(spa45*spb95)-spa46*spb96-spa47*spb97;
const complex<T> spab_4_567_3=spa45*spb53+spa46*spb63+spa47*spb73;
const complex<T> spab_4_56_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_23_5=spa24*spb52+spa34*spb53;
const complex<T> spab_4_123_5=spa14*spb51+spa24*spb52+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spaa_4_56_79_8=spa45*(-(spa78*spb75)+spa89*spb95)+spa46*(-(spa78*spb76)+spa89*spb96);
const complex<T> spaa_4_56_789_4=spa45*(spa47*spb75+spa48*spb85+spa49*spb95)+spa46*(spa47*spb76+spa48*spb86+spa49*spb96);
const complex<T> spaa_4_56_789_1=spa45*(spa17*spb75+spa18*spb85+spa19*spb95)+spa46*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_4_567_89_1=-(spa18*(-(spa45*spb85)-spa46*spb86-spa47*spb87))-spa19*(-(spa45*spb95)-spa46*spb96-spa47*spb97);
const complex<T> spaa_4_567_189_2=spa45*(spa12*spb51+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa28*spb86+spa29*spb96)+spa47*(spa12*spb71+spa28*spb87+spa29*spb97);
const complex<T> spaa_4_56_1789_2=spa45*(spa12*spb51+spa27*spb75+spa28*spb85+spa29*spb95)+spa46*(spa12*spb61+spa27*spb76+spa28*spb86+spa29*spb96);
const complex<T> spaa_4_23_567_4=-(spa24*(spa45*spb52+spa46*spb62+spa47*spb72))-spa34*(spa45*spb53+spa46*spb63+spa47*spb73);
const complex<T> spaa_4_23_56_4=-(spa24*(spa45*spb52+spa46*spb62))-spa34*(spa45*spb53+spa46*spb63);
const complex<T> spaa_4_23_189_7=-(spa24*(spa17*spb21+spa78*spb82+spa79*spb92))-spa34*(spa17*spb31+spa78*spb83+spa79*spb93);
const complex<T> spaa_4_23_1789_6=-(spa24*(spa16*spb21+spa67*spb72+spa68*spb82+spa69*spb92))-spa34*(spa16*spb31+spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> spaa_4_123_89_7=-(spa78*(spa14*spb81+spa24*spb82+spa34*spb83))-spa79*(spa14*spb91+spa24*spb92+spa34*spb93);
const complex<T> spaa_4_123_789_6=-(spa14*(spa67*spb71+spa68*spb81+spa69*spb91))-spa24*(spa67*spb72+spa68*spb82+spa69*spb92)-spa34*(spa67*spb73+spa68*spb83+spa69*spb93);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s67=spa67*spb76;
const complex<T> s56=spa56*spb65;
const complex<T> s46=spa46*spb64;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow2_spa18*spa47*pow(spaa_4_23_567_4,3))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spa89*spaa_4_23_189_7*spaa_4_567_189_2*spaa_4_567_89_1)
)-
(pow2_spa18*pow(spaa_4_23_56_4,4))/(spa23*spa34*spa45*spa56*spa89*spaa_4_23_1789_6*spaa_4_23_189_7*spaa_4_56_1789_2*spaa_4_56_789_1*(s18+s19+s78+s79+s89+spa17*spb71)
)+
(pow3_spa14*spaa_4_56_789_4*pow(spaa_4_56_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_4_123_789_6*spaa_4_123_89_7*spaa_4_56_789_1
)+
(pow2_spa18*pow(spab_4_23_5,4))/((s23+s24+s34)*spa23*spa34*spa67*spa89*spaa_4_23_1789_6*spab_1_234_5*spab_2_34_5*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)-
(pow2_spa18*pow(spab_4_56_3,4))/((s45+s46+s56)*spa12*spa45*spa56*spa89*spaa_4_56_1789_2*spab_6_45_3*spab_7_456_3*(s34+s35+s45+s46+s56+spb63*SPA(i3,i6))
)-
(pow3_spa14*spab_4_123_5*pow(spab_8_679_5,2))/(spa12*spa23*spa34*spa67*spa89*spaa_4_123_789_6*spab_1_234_5*(s67+s78+s79+s89+spa68*spb86+spa69*spb96)*(s12+s23+s24+s34+spa13*spb31+spa14*SPB(i4,i1))
)-
(pow2_spa18*pow(SPB(i3,i5),4))/((s34+s35+s45)*spa12*spa67*spa89*spab_2_34_5*spab_6_45_3*SPB(i3,i4)*SPB(i4,i5)
)+
(pow2_spa18*spa47*pow(spab_4_567_3,3))/(spa12*spa45*spa56*spa67*spa89*spaa_4_567_189_2*spab_7_456_3*(s12+s18+s19+s89+spa28*spb82+spa29*spb92)*(s45+s46+s56+s67+spa57*spb75+spa47*SPB(i7,i4))
)-
(pow3_spa14*spa47*pow(spab_4_567_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_4_123_89_7*spaa_4_567_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qppppmpqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+SPA(i3,i9)*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_456_12_3=spa13*(spa34*spb41+spa35*spb51+spa36*spb61)+spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*spb51+spa23*spb52);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((pow3_spa13*spaa_3_456_12_3*pow(spaa_3_456_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)
)-
(pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_679_8,2))/(spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)
)-
(pow3_spa13*spab_3_12_4*pow(spab_8_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa18*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa18*spa37*pow(spab_3_189_2,3))/((s18+s19+s89)*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2*(s12+s18+s19+s89+spa28*spb82+spb92*SPA(i2,i9))
)-
(pow2_spa18*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa18*pow(spab_3_456_2,4))/(spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa37*pow(spab_3_128_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qppppmpqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb92=SPB(i9,i2);
const complex<T> spb91=SPB(i9,i1);
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
const complex<T> spa38=SPA(i3,i8);
const complex<T> spa37=SPA(i3,i7);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=pow(spa13,3);
const complex<T> pow2_spa18=pow(spa18,2);
const complex<T> spab_8_123_4=spa18*spb41+spa28*spb42+spa38*spb43;
const complex<T> spab_7_189_2=spa17*spb21+spa78*spb82+spa79*spb92;
const complex<T> spab_6_345_2=-(spa36*spb32)-spa46*spb42-spa56*spb52;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_45_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_189_2=spa13*spb21+spa38*spb82+SPA(i3,i9)*spb92;
const complex<T> spab_3_128_9=spa13*spb91+spa23*spb92-spa38*spb98;
const complex<T> spab_3_12_4=spa13*spb41+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_456_79_8=spa78*(-(spa34*spb74)-spa35*spb75-spa36*spb76)-spa89*(-(spa34*spb94)-spa35*spb95-spa36*spb96);
const complex<T> spaa_3_45_679_8=spa34*(-(spa68*spb64)-spa78*spb74+spa89*spb94)+spa35*(-(spa68*spb65)-spa78*spb75+spa89*spb95);
const complex<T> spaa_3_456_789_1=spa34*(spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa17*spb75+spa18*spb85+spa19*spb95)+spa36*(spa17*spb76+spa18*spb86+spa19*spb96);
const complex<T> spaa_3_45_6789_1=spa34*(spa16*spb64+spa17*spb74+spa18*spb84+spa19*spb94)+spa35*(spa16*spb65+spa17*spb75+spa18*spb85+spa19*spb95);
const complex<T> spaa_3_4567_89_1=-(spa18*(-(spa34*spb84)-spa35*spb85-spa36*spb86-spa37*spb87))-spa19*(-(spa34*spb94)-spa35*spb95-spa36*spb96-spa37*spb97);
const complex<T> spaa_3_456_12_3=spa13*(spa34*spb41+spa35*spb51+spa36*spb61)+spa23*(spa34*spb42+spa35*spb52+spa36*spb62);
const complex<T> spaa_3_45_12_3=spa34*(spa13*spb41+spa23*spb42)+spa35*(spa13*spb51+spa23*spb52);
const complex<T> spaa_3_12_89_7=-(spa13*(spa78*spb81+spa79*spb91))-spa23*(spa78*spb82+spa79*spb92);
const complex<T> spaa_3_12_789_6=-(spa13*(spa67*spb71+spa68*spb81+spa69*spb91))-spa23*(spa67*spb72+spa68*spb82+spa69*spb92);
const complex<T> spaa_3_12_6789_5=-(spa13*(spa56*spb61+spa57*spb71+spa58*spb81+spa59*spb91))-spa23*(spa56*spb62+spa57*spb72+spa58*spb82+spa59*spb92);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s45=spa45*spb54;
const complex<T> s35=spa35*spb53;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((pow3_spa13*spaa_3_456_12_3*pow(spaa_3_456_79_8,2))/((s78+s79+s89)*spa12*spa23*spa34*spa45*spa56*spa89*spaa_3_12_789_6*spaa_3_12_89_7*spaa_3_456_789_1)
)-
(pow3_spa13*spaa_3_45_12_3*pow(spaa_3_45_679_8,2))/(spa12*spa23*spa34*spa45*spa67*spa89*spaa_3_12_6789_5*spaa_3_12_789_6*spaa_3_45_6789_1*(s78+s79+s89+spa67*spb76+spa68*spb86+spa69*spb96)
)-
(pow3_spa13*spab_3_12_4*pow(spab_8_123_4,2))/((s12+s13+s23)*spa12*spa23*spa56*spa67*spa89*spaa_3_12_6789_5*spab_1_23_4*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4))
)-
(pow2_spa18*pow(spab_3_45_2,4))/((s34+s35+s45)*spa34*spa45*spa67*spa89*spaa_3_45_6789_1*spab_5_34_2*spab_6_345_2*(s23+s24+s34+s35+s45+spb52*SPA(i2,i5))
)+
(pow2_spa18*spa37*pow(spab_3_189_2,3))/((s18+s19+s89)*spa34*spa45*spa56*spa67*spa89*spaa_3_4567_89_1*spab_7_189_2*(s12+s18+s19+s89+spa28*spb82+spb92*SPA(i2,i9))
)-
(pow2_spa18*pow(SPB(i2,i4),4))/((s23+s24+s34)*spa56*spa67*spa89*spab_1_23_4*spab_5_34_2*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa18*pow(spab_3_456_2,4))/(spa34*spa45*spa56*spa89*spaa_3_456_789_1*spab_6_345_2*spab_7_189_2*(s18+s19+s78+s79+s89+spa17*spb71)*(s34+s35+s45+spa46*spb64+spa56*spb65+spa36*SPB(i6,i3))
)-
(pow3_spa13*spa37*pow(spab_3_128_9,2))/(spa12*spa23*spa34*spa45*spa56*spa67*spaa_3_12_89_7*spaa_3_4567_89_1*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpppppmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
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
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_789_1=spa27*spb71+spa28*spb81+SPA(i2,i9)*spb91;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*SPB(i6,i3))-spa78*spb73+spa89*spb93)+spa24*(-(spa68*SPB(i6,i4))-spa78*spb74+spa89*spb94)+spa25*(-(spa68*SPB(i6,i5))-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_1234_8=spa24*(spa18*spb41+spa28*spb42+SPA(i3,i8)*spb43)+spa23*(spa18*spb31+spa28*spb32-SPA(i4,i8)*spb43);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,1)*(-((spab_2_34_1*pow(spaa_2_34_1234_8,2))/((s23+s24+s34)*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)+
(spab_2_789_1*pow(spaa_2_3456_79_8,2))/((s78+s79+s89)*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow(spab_8_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa67*spa89*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(spab_2_345_1*pow(spaa_2_345_679_8,2))/(spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1*(s23+s24+s34+spa25*SPB(i5,i2
)+
spa35*SPB(i5,i3
)+
spa45*SPB(i5,i4))*(s78+s79+s89+spa67*spb76+spa69*spb96+spa68*SPB(i8,i6))
)-
(spa27*pow(spab_2_18_9,2)*SPB(i1,i9))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb91*SPB(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpppppmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb98=SPB(i9,i8);
const complex<T> spb97=SPB(i9,i7);
const complex<T> spb96=SPB(i9,i6);
const complex<T> spb95=SPB(i9,i5);
const complex<T> spb94=SPB(i9,i4);
const complex<T> spb93=SPB(i9,i3);
const complex<T> spb91=SPB(i9,i1);
const complex<T> spb87=SPB(i8,i7);
const complex<T> spb81=SPB(i8,i1);
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb73=SPB(i7,i3);
const complex<T> spb71=SPB(i7,i1);
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
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa28=SPA(i2,i8);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa19=SPA(i1,i9);
const complex<T> spa18=SPA(i1,i8);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_8_12_3=spa18*spb31+spa28*spb32;
const complex<T> spab_7_89_1=spa78*spb81+spa79*spb91;
const complex<T> spab_6_789_1=spa67*spb71+spa68*spb81+spa69*spb91;
const complex<T> spab_5_234_1=-(spa25*spb21)-spa35*spb31-spa45*spb41;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_789_1=spa27*spb71+spa28*spb81+SPA(i2,i9)*spb91;
const complex<T> spab_2_345_1=spa23*spb31+spa24*spb41+spa25*SPB(i5,i1);
const complex<T> spab_2_34_1=spa23*spb31+spa24*spb41;
const complex<T> spab_2_18_9=spa12*spb91-spa28*spb98;
const complex<T> spaa_2_3456_79_8=spa78*(-(spa23*spb73)-spa24*spb74-spa25*spb75-spa26*spb76)-spa89*(-(spa23*spb93)-spa24*spb94-spa25*spb95-spa26*spb96);
const complex<T> spaa_2_345_679_8=spa23*(-(spa68*SPB(i6,i3))-spa78*spb73+spa89*spb93)+spa24*(-(spa68*SPB(i6,i4))-spa78*spb74+spa89*spb94)+spa25*(-(spa68*SPB(i6,i5))-spa78*spb75+spa89*spb95);
const complex<T> spaa_2_34_1234_8=spa24*(spa18*spb41+spa28*spb42+SPA(i3,i8)*spb43)+spa23*(spa18*spb31+spa28*spb32-SPA(i4,i8)*spb43);
const complex<T> s89=spa89*spb98;
const complex<T> s79=spa79*spb97;
const complex<T> s78=spa78*spb87;
const complex<T> s34=spa34*spb43;
const complex<T> s24=spa24*spb42;
const complex<T> s23=spa23*spb32;
const complex<T> s19=spa19*spb91;
const complex<T> s18=spa18*spb81;
const complex<T> s13=spa13*spb31;
const complex<T> s12=spa12*spb21;

return(
complex<T>(0,-1)*(-((spab_2_34_1*pow(spaa_2_34_1234_8,2))/((s23+s24+s34)*spa23*spa34*spa56*spa67*spa89*spab_4_23_1*spab_5_234_1*(s12+s13+s23+s24+s34+spb41*SPA(i1,i4)))
)+
(spab_2_789_1*pow(spaa_2_3456_79_8,2))/((s78+s79+s89)*spa23*spa34*spa45*spa56*spa89*spab_6_789_1*spab_7_89_1*(s18+s19+s78+s79+s89+spb71*SPA(i1,i7))
)-
(pow(spab_8_12_3,2)*SPB(i1,i3))/((s12+s13+s23)*spa45*spa56*spa67*spa89*spab_4_23_1*SPB(i1,i2)*SPB(i2,i3)
)+
(spab_2_345_1*pow(spaa_2_345_679_8,2))/(spa23*spa34*spa45*spa67*spa89*spab_5_234_1*spab_6_789_1*(s23+s24+s34+spa25*SPB(i5,i2
)+
spa35*SPB(i5,i3
)+
spa45*SPB(i5,i4))*(s78+s79+s89+spa67*spb76+spa69*spb96+spa68*SPB(i8,i6))
)-
(spa27*pow(spab_2_18_9,2)*SPB(i1,i9))/((s18+s19+s89)*spa23*spa34*spa45*spa56*spa67*spab_7_89_1*spb91*SPB(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmpppppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i8),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i6,i7)*SPA(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8, int i9, class T> complex<T>  A2q5g2l_qmmmmmmqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
(complex<T>(0,1)*pow(SPB(i7,i9),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i6,i7)*SPB(i8,i9))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qppppppqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i8),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i6,i7)*SPA(i8,i9)))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmpppppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i8),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i6,i7)*SPA(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i8, int i9, class T> complex<T>  A2q5g2l_qpmmmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
(complex<T>(0,1)*pow(SPB(i7,i9),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i6,i7)*SPB(i8,i9))
);
}


template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i9, int i8, class T> complex<T>  A2q5g2l_qmmmmmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
(complex<T>(0,-1)*pow(SPB(i7,i9),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i6,i7)*SPB(i8,i9))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qppppppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i8),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i5,i6)*SPA(i6,i7)*SPA(i8,i9)))
);
}


template <int i7, int i6, int i5, int i4, int i3, int i2, int i1, int i9, int i8, class T> complex<T>  A2q5g2l_qpmmmmmqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
(complex<T>(0,-1)*pow(SPB(i7,i9),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i5,i6)*SPB(i6,i7)*SPB(i8,i9))
);
}




template <class T> complex<T> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&)
{
	switch (hc) {

case 8164802:	 return &A2q5g2l_qpmmmmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m m m m qbm lp lbm
case 8164820:	 return &A2q5g2l_qppmmmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p m m m m qbm lp lbm
case 8164910:	 return &A2q5g2l_qpmpmmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m p m m m qbm lp lbm
case 8165450:	 return &A2q5g2l_qpmmpmmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m p m m qbm lp lbm
case 8168690:	 return &A2q5g2l_qpmmmpmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m m p m qbm lp lbm
case 8169464:	 return &A2q5g2l_qpppppmqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p p p m qbm lp lbm
case 8188130:	 return &A2q5g2l_qpmmmmpqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m m m m p qbm lp lbm
case 8188904:	 return &A2q5g2l_qppppmpqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p p m p qbm lp lbm
case 8192144:	 return &A2q5g2l_qpppmppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p m p p qbm lp lbm
case 8192684:	 return &A2q5g2l_qppmpppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p m p p p qbm lp lbm
case 8192774:	 return &A2q5g2l_qpmppppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp m p p p p qbm lp lbm
case 8192792:	 return &A2q5g2l_qppppppqbmlplbm_eval<0,1,2,3,4,5,6,7,8>;//qp p p p p p qbm lp lbm
case 8211457:	 return &A2q5g2l_qmmmmmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m m m m qbp lp lbm
case 8211475:	 return &A2q5g2l_qmpmmmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p m m m m qbp lp lbm
case 8211565:	 return &A2q5g2l_qmmpmmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m p m m m qbp lp lbm
case 8212105:	 return &A2q5g2l_qmmmpmmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m p m m qbp lp lbm
case 8215345:	 return &A2q5g2l_qmmmmpmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m m p m qbp lp lbm
case 8216119:	 return &A2q5g2l_qmppppmqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p p p m qbp lp lbm
case 8234785:	 return &A2q5g2l_qmmmmmpqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m m m m p qbp lp lbm
case 8235559:	 return &A2q5g2l_qmpppmpqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p p m p qbp lp lbm
case 8238799:	 return &A2q5g2l_qmppmppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p m p p qbp lp lbm
case 8239339:	 return &A2q5g2l_qmpmpppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p m p p p qbp lp lbm
case 8239429:	 return &A2q5g2l_qmmppppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm m p p p p qbp lp lbm
case 8239447:	 return &A2q5g2l_qmpppppqbplplbm_eval<0,1,2,3,4,5,6,7,8>;//qm p p p p p qbp lp lbm
case 9564482:	 return &A2q5g2l_qpmmmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m m m m qbm lm lbp
case 9564500:	 return &A2q5g2l_qppmmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p m m m m qbm lm lbp
case 9564590:	 return &A2q5g2l_qpmpmmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m p m m m qbm lm lbp
case 9565130:	 return &A2q5g2l_qpmmpmmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m p m m qbm lm lbp
case 9568370:	 return &A2q5g2l_qpmmmpmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m m p m qbm lm lbp
case 9569144:	 return &A2q5g2l_qpppppmqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p p p m qbm lm lbp
case 9587810:	 return &A2q5g2l_qpmmmmpqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m m m m p qbm lm lbp
case 9588584:	 return &A2q5g2l_qppppmpqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p p m p qbm lm lbp
case 9591824:	 return &A2q5g2l_qpppmppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p m p p qbm lm lbp
case 9592364:	 return &A2q5g2l_qppmpppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p m p p p qbm lm lbp
case 9592454:	 return &A2q5g2l_qpmppppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp m p p p p qbm lm lbp
case 9592472:	 return &A2q5g2l_qppppppqbmlmlbp_eval<0,1,2,3,4,5,6,7,8>;//qp p p p p p qbm lm lbp
case 9611137:	 return &A2q5g2l_qmmmmmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m m m m qbp lm lbp
case 9611155:	 return &A2q5g2l_qmpmmmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p m m m m qbp lm lbp
case 9611245:	 return &A2q5g2l_qmmpmmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m p m m m qbp lm lbp
case 9611785:	 return &A2q5g2l_qmmmpmmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m p m m qbp lm lbp
case 9615025:	 return &A2q5g2l_qmmmmpmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m m p m qbp lm lbp
case 9615799:	 return &A2q5g2l_qmppppmqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p p p m qbp lm lbp
case 9634465:	 return &A2q5g2l_qmmmmmpqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m m m m p qbp lm lbp
case 9635239:	 return &A2q5g2l_qmpppmpqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p p m p qbp lm lbp
case 9638479:	 return &A2q5g2l_qmppmppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p m p p qbp lm lbp
case 9639019:	 return &A2q5g2l_qmpmpppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p m p p p qbp lm lbp
case 9639109:	 return &A2q5g2l_qmmppppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm m p p p p qbp lm lbp
case 9639127:	 return &A2q5g2l_qmpppppqbplmlbp_eval<0,1,2,3,4,5,6,7,8>;//qm p p p p p qbp lm lbp
	
	default: // We return the zero pointer for all other helicity combinations
		//cout << " A2q5g2l_Tree_Ptr_eval : Missing entry for helcode=" << hex << hc << dec << endl;
		return 0;
	}
}

template complex<R> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&);
template complex<RHP> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&);
template complex<RVHP> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&);

#if BH_USE_GMP

template complex<RGMP> (*A2q5g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&);
#endif

}
