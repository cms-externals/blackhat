/*
 * A0_2q3g2l_eval.cpp
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

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q3g2l_qmmppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_5_67_1=spa56*spb61+spa57*spb71;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_2_34_1=spa23*spb31+spa24*SPB(i4,i1);
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;
const complex<T> spaa_2_34_57_6=spa23*(-(spa56*SPB(i5,i3))+spa67*SPB(i7,i3))+spa24*(-(spa56*SPB(i5,i4))+spa67*SPB(i7,i4));

return(
complex<T>(0,1)*(-((pow(spab_6_12_3,2)*SPB(i1,i3))/(spa45*spa67*spab_4_23_1*(spa12*spb21+spa23*spb32+spb31*SPA(i1,i3))*SPB(i1,i2)*SPB(i2,i3))
)+
(pow(spab_2_16_7,2)*SPA(i2,i5))/(spa23*spa34*spa45*spab_5_67_1*(spa16*spb61+spa67*spb76+spb71*SPA(i1,i7))*SPB(i6,i7)
)+
(spab_2_34_1*pow(spaa_2_34_57_6,2))/(spa23*spa34*spa67*spab_4_23_1*spab_5_67_1*(spa23*spb32+spa24*SPB(i4,i2
)+
spa34*SPB(i4,i3))*(spa67*spb76+spa56*SPB(i6,i5
)+
spa57*SPB(i7,i5))))
);
}
template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q3g2l_qpmpmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
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
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa27=SPA(i2,i7);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb35=(pow(spb35,3));
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spbb_3_45_67_1=-(spb43*(spa46*spb61+spa47*spb71))-spb53*(spa56*spb61+spa57*spb71);
const complex<T> spbb_3_12_67_5=spb31*(spa16*spb65+spa17*spb75)+spb32*(spa26*spb65+spa27*spb75);
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_4_23_1=-(spa24*spb21)-spa34*spb31;
const complex<T> spab_4_12_3=SPA(i1,i4)*spb31+spa24*spb32;
const complex<T> spab_2_45_3=spa24*spb43+SPA(i2,i5)*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_2_16_7=spa12*spb71-spa26*spb76;

return(
complex<T>(0,1)*((pow3_spb35*spb13*pow(spab_6_12_3,2))/(spa67*spb12*spb23*spb34*spb45*spbb_3_12_67_5*spbb_3_45_67_1
)+
(pow2_spb57*spb13*pow(spab_4_12_3,3))/(spab_4_23_1*spb12*spb23*spb67*(spa56*spb65+spa57*spb75+spa67*spb76)*spbb_3_12_67_5*(spa12*spb21+spa23*spb32+spb31*SPA(i1,i3))
)-
(pow3_spb35*spab_2_45_3*pow(spab_2_16_7,2))/(spab_2_34_5*spb34*spb45*spb67*(spa16*spb61+spa17*spb71+spa67*spb76)*spbb_3_45_67_1*(spa34*spb43+spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))
)+
(pow2_spb57*pow(spa24,4))/(spa23*spa34*spab_2_34_5*spab_4_23_1*spb67*(spa23*spb32+spa34*spb43+spa24*SPB(i4,i2))))
);
}
template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q3g2l_qmpppqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,1)*pow(SPA(i1,i6),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i6,i7)))
);
}

}
