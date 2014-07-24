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

template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q3g2l_qmpmmqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb67=SPB(i6,i7);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb57=SPB(i5,i7);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spb57=(pow(spb57,2));
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;
const complex<T> spab_1_34_2=spa13*spb32+SPA(i1,i4)*spb42;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;

return(
complex<T>(0,-1)*(-((pow2_spb57*pow(spab_1_34_2,3))/(spab_1_23_4*spab_1_67_5*spb23*spb34*spb67*(spa23*spb32+spb42*SPA(i2,i4
)+
spb43*SPA(i3,i4))*(spa67*spb76+spb65*SPA(i5,i6
)+
spb75*SPA(i5,i7)))
)+
(pow2_spb57*pow(spa13,3))/(spa12*spa23*spab_1_23_4*spb45*spb67*(spa23*spb32+spa12*SPB(i2,i1
)+
spa13*SPB(i3,i1))
)-
(pow(spa16,2)*pow(SPB(i2,i5),3))/(spa67*spab_1_67_5*spb23*spb34*spb45*(spa67*spb76+spa16*SPB(i6,i1
)+
spa17*SPB(i7,i1))))
);
}
template <int i5, int i4, int i3, int i2, int i1, int i7, int i6, class T> complex<T>  A2q3g2l_qpmppqbmlplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{
const complex<T> spb76=SPB(i7,i6);
const complex<T> spb75=SPB(i7,i5);
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spa67=SPA(i6,i7);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa17=SPA(i1,i7);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow2_spa16=(pow(spa16,2));
const complex<T> spab_4_23_5=spa24*SPB(i5,i2)+spa34*spb53;
const complex<T> spab_2_34_5=-(spa23*spb53)-spa24*spb54;
const complex<T> spab_1_67_5=spa16*spb65+spa17*spb75;

return(
complex<T>(0,-1)*(-((pow2_spa16*pow(SPB(i3,i5),3))/(spa12*spa67*spab_2_34_5*(spa34*spb43+spb53*SPA(i3,i5
)+
spb54*SPA(i4,i5))*SPB(i3,i4)*SPB(i4,i5))
)+
(pow(SPA(i1,i4),3)*pow(SPB(i5,i7),2))/(spa12*spa23*spa34*spab_1_67_5*(spa67*spb76+spb65*SPA(i5,i6
)+
spb75*SPA(i5,i7))*SPB(i6,i7)
)+
(pow2_spa16*pow(spab_4_23_5,3))/(spa23*spa34*spa67*spab_1_67_5*spab_2_34_5*(spa34*spb43+spa23*SPB(i3,i2
)+
spa24*SPB(i4,i2))*(spa67*spb76+spa16*SPB(i6,i1
)+
spa17*SPB(i7,i1))))
);
}
template <int i1, int i2, int i3, int i4, int i5, int i7, int i6, class T> complex<T>  A2q3g2l_qmpppqbplplbm_eval(const eval_param<T>& ep, const mass_param_coll& masses)
{

return(
-((complex<T>(0,-1)*pow(SPA(i1,i6),2))/(SPA(i1,i2)*SPA(i2,i3)*SPA(i3,i4)*SPA(i4,i5)*SPA(i6,i7)))
);
}

}
