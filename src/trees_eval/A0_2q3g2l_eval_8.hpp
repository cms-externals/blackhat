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

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T>  A2q3g2l_qmpmpqbplmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
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
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;
const complex<T> spab_3_45_2=spa34*spb42+spa35*SPB(i5,i2);
const complex<T> spab_3_12_4=spa13*SPB(i4,i1)+spa23*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spaa_3_45_67_1=spa34*(spa16*spb64+spa17*spb74)+spa35*(spa16*spb65+spa17*spb75);
const complex<T> spaa_3_12_67_5=-(spa13*(spa56*spb61+spa57*spb71))-spa23*(spa56*spb62+spa57*spb72);

return(
complex<T>(0,1)*(-((pow3_spa13*spab_3_12_4*pow(spab_6_57_4,2))/(spa12*spa23*spa67*spaa_3_12_67_5*spab_1_23_4*(spa56*spb65+spa57*spb75+spa67*spb76)*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1)))
)-
(pow2_spa16*pow(SPB(i2,i4),4))/(spa67*spab_1_23_4*spab_5_34_2*(spa23*spb32+spa34*spb43+spb42*SPA(i2,i4))*SPB(i2,i3)*SPB(i3,i4)
)+
(pow2_spa16*spa35*pow(spab_3_45_2,3))/(spa34*spa45*spa67*spaa_3_45_67_1*spab_5_34_2*(spa16*spb61+spa17*spb71+spa67*spb76)*(spa34*spb43+spa45*spb54+spa35*SPB(i5,i3))
)-
(pow3_spa13*spa35*pow(spab_3_45_7,2))/(spa12*spa23*spa34*spa45*spaa_3_12_67_5*spaa_3_45_67_1*SPB(i6,i7)))
);
}
template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q3g2l_qppmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb74=SPB(i7,i4);
const complex<T> spb71=SPB(i7,i1);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa57=SPA(i5,i7);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spbb_4_23_16_7=spb42*(SPA(i1,i2)*spb71-SPA(i2,i6)*spb76)+spb43*(SPA(i1,i3)*spb71-SPA(i3,i6)*spb76);
const complex<T> spab_6_57_4=-(spa56*spb54)+spa67*spb74;
const complex<T> spab_5_67_1=spa56*spb61+spa57*spb71;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_23_4=SPA(i2,i5)*spb42+spa35*spb43;
const complex<T> spab_3_45_7=-(spa34*spb74)-spa35*spb75;

return(
complex<T>(0,1)*(-((spab_5_23_4*pow(spbb_4_23_16_7,2))/(spab_5_34_2*spab_5_67_1*spb23*spb34*spb67*(spa67*spb76+spb61*SPA(i1,i6
)+
spb71*SPA(i1,i7))*(spa34*spb43+spb32*SPA(i2,i3
)+
spb42*SPA(i2,i4)))
)+
(spa35*pow(spab_3_45_7,2))/(spa34*spa45*spab_5_34_2*spb12*spb67*(spa34*spb43+spa45*spb54+spa35*SPB(i5,i3))
)-
(pow(spab_6_57_4,2)*SPB(i1,i4))/(spa67*spab_5_67_1*spb12*spb23*spb34*(spa57*spb75+spa67*spb76+spa56*SPB(i6,i5))))
);
}
template <int i5, int i4, int i3, int i2, int i1, int i6, int i7, class T> complex<T>  A2q3g2l_qpmmmqbmlmlbp_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
(complex<T>(0,1)*pow(SPB(i5,i7),2))/(SPB(i1,i2)*SPB(i2,i3)*SPB(i3,i4)*SPB(i4,i5)*SPB(i6,i7))
);
}

}
