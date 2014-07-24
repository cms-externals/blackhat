
#include "R_6g_eval.h"
#include "eval_param.h"


using namespace std;
namespace BH  {


#define _ONLY_X_PART 1
#define X 2
//#define _VERBOSE 1

#define SPA(i,j) ep.spa(i,j)
#define SPB(i,j) ep.spb(i,j)

/*

 mpmppp, mppmpp rational terms from Lance using Haralds scripts. Can be further cleaned up
 complex conjugation used explicity for some rational parts; we would thus get wrong values 
 for complex momenta.
	
*/

int helcode_g(const process& pro);

// added by Daniel
template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6g1m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R6 : -+++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
+complex<T>(0,1)/complex<T>(3,0)*
		((pow(ep.spa(i1,i4),3)
		*ep.spa(i3,i5)*(ep.spab(i1,i2,i4) + ep.spab(i1,i3,i4)))/(ep.spa(i1,i2)*ep.spa(i1,i6)*ep.spa(i2,i3)*pow(ep.spa(i3,i4),2)*pow(ep.spa(i4,i5),2)*ep.spa(i5,i6)) + pow(ep.spab(i1,i2,i6) + 	ep.spab(i1,i3,i6),3)/(ep.s(i1,i2,i3)*ep.spa(i1,i2)*ep.spa(i2,i3)*pow(ep.spa(i4,i5),2)*(ep.spab(i3,i1,i6) +ep.spab(i3,i2,i6))) - pow(ep.spab(i1,i3,i2) + ep.spab(i1,i4,i2),3)/(ep.s(i2,i3,i4)*ep.spa(i1,i6)*pow(ep.spa(i3,i4),2)*ep.spa(i5,i6)*(ep.spab(i5,i3,i2) +ep.spab(i5,i4,i2))) - (pow(ep.spa(i1,i3),3)*ep.spa(i2,i4)*ep.spb(i3,i2))/(ep.spa(i1,i6)*pow(ep.spa(i2,i3),2)*pow(ep.spa(i3,i4),2)*ep.spa(i4,i5)*ep.spa(i5,i6)) - (pow(ep.spa(i1,i5),3)*ep.spa(i4,i6)*ep.spb(i6,i5))/(ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*pow(ep.spa(i4,i5),2)*pow(ep.spa(i5,i6),2)) + (pow(ep.spb(i6,i2),3)*((ep.spb(i3,i2)*ep.spb(i4,i3))/(ep.spa(i4,i5)*(ep.spab(i5,i3,i2) + ep.spab(i5,i4,i2))) - ep.spb(i5,i3)/(ep.spa(i3,i4)*ep.spa(i4,i5)) - (ep.spb(i5,i4)*ep.spb(i6,i5))/(ep.spa(i3,i4)*(ep.spab(i3,i1,i6) + ep.spab(i3,i2,i6)))))/(ep.s(i3,i4,i5)*ep.spb(i2,i1)*ep.spb(i6,i1))
		);

		#endif
;
}

// added by Daniel
template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6g1p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R6 : +-----");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
-std::conj(R6g1m<i1,i2,i3,i4,i5,i6>(ep,mpc));

		#endif
;
}

// added by Daniel
template <class T> complex<T> R6g6p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R6 : ++++++");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
+(T(1)/complex<T>(0,-3)*(
		- ep.spa(0,1)*ep.spa(2,3)*ep.spb(2,1)*ep.spb(3,0)
		- ep.spa(0,1)*ep.spa(2,4)*ep.spb(2,1)*ep.spb(4,0)
		- ep.spa(0,1)*ep.spa(3,4)*ep.spb(3,1)*ep.spb(4,0)
		- ep.spa(0,2)*ep.spa(3,4)*ep.spb(3,2)*ep.spb(4,0)
		- ep.spa(1,2)*ep.spa(3,4)*ep.spb(3,2)*ep.spb(4,1)
		- ep.spa(0,1)*ep.spa(2,5)*ep.spb(2,1)*ep.spb(5,0)
		- ep.spa(0,1)*ep.spa(3,5)*ep.spb(3,1)*ep.spb(5,0)
		- ep.spa(0,2)*ep.spa(3,5)*ep.spb(3,2)*ep.spb(5,0)
		- ep.spa(0,1)*ep.spa(4,5)*ep.spb(4,1)*ep.spb(5,0)
		- ep.spa(0,2)*ep.spa(4,5)*ep.spb(4,2)*ep.spb(5,0)
		- ep.spa(0,3)*ep.spa(4,5)*ep.spb(4,3)*ep.spb(5,0)
		- ep.spa(1,2)*ep.spa(3,5)*ep.spb(3,2)*ep.spb(5,1)
		- ep.spa(1,2)*ep.spa(4,5)*ep.spb(4,2)*ep.spb(5,1)
		- ep.spa(1,3)*ep.spa(4,5)*ep.spb(4,3)*ep.spb(5,1)
		- ep.spa(2,3)*ep.spa(4,5)*ep.spb(4,3)*ep.spb(5,2)
		))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)
				*ep.spa(3,4)*ep.spa(4,5))
#endif
	;

}

template <class T> complex<T> R6g6m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
#if _VERBOSE
	_MESSAGE("R6 : ------");
#endif
	return complex<T>(0,0)
	#if _ONLY_X_PART
+(T(1)/complex<T>(0,-3)*(
		- ep.spb(1,0)*ep.spb(3,2)*ep.spa(1,2)*ep.spa(0,3)
		- ep.spb(1,0)*ep.spb(4,2)*ep.spa(1,2)*ep.spa(0,4)
		- ep.spb(1,0)*ep.spb(4,3)*ep.spa(1,3)*ep.spa(0,4)
		- ep.spb(2,0)*ep.spb(4,3)*ep.spa(2,3)*ep.spa(0,4)
		- ep.spb(2,1)*ep.spb(4,3)*ep.spa(2,3)*ep.spa(1,4)
		- ep.spb(1,0)*ep.spb(5,2)*ep.spa(1,2)*ep.spa(0,5)
		- ep.spb(1,0)*ep.spb(5,3)*ep.spa(1,3)*ep.spa(0,5)
		- ep.spb(2,0)*ep.spb(5,3)*ep.spa(2,3)*ep.spa(0,5)
		- ep.spb(1,0)*ep.spb(5,4)*ep.spa(1,4)*ep.spa(0,5)
		- ep.spb(2,0)*ep.spb(5,4)*ep.spa(2,4)*ep.spa(0,5)
		- ep.spb(3,0)*ep.spb(5,4)*ep.spa(3,4)*ep.spa(0,5)
		- ep.spb(2,1)*ep.spb(5,3)*ep.spa(2,3)*ep.spa(1,5)
		- ep.spb(2,1)*ep.spb(5,4)*ep.spa(2,4)*ep.spa(1,5)
		- ep.spb(3,1)*ep.spb(5,4)*ep.spa(3,4)*ep.spa(1,5)
		- ep.spb(3,2)*ep.spb(5,4)*ep.spa(3,4)*ep.spa(2,5)
		))/(ep.spb(1,0)*ep.spb(5,0)*ep.spb(2,1)*ep.spb(3,2)
				*ep.spb(4,3)*ep.spb(5,4))
#endif
	;

}



template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6g2m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
  return  complex<T>(0,1)*((complex<T>(-X,0)/complex<T>(9,0)*pow(ep.spa(i1,i2),3))/(ep.spa(i1,i6)*ep.spa(i2,i3)*
ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i6))
#if _ONLY_X_PART
+
(complex<T>(1,0)/complex<T>(3,0)*pow(-(ep.spa(i1,i2)*ep.spb(i3,i2))+ep.spa(i1,i4)*ep.spb(i4,i3),3)*
ep.spa(i3,i5))/(ep.s(i2,i3,i4)*ep.spa(i1,i6)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i6)*
ep.spb(i3,i2)*(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*ep.spa(i1,i2)*(-(ep.spa(i1,i2)*ep.spb(i3,i2))+ep.spa(i1,i4)*ep.spb(i4,i3))*
(ep.spa(i1,i4)*ep.spb(i4,i3)+complex<T>(2,0)*(-(ep.spa(i1,i2)*ep.spb(i3,i2)) +
   ep.spa(i1,i4)*ep.spb(i4,i3))))/(ep.s(i2,i3,i4)*ep.spa(i1,i6)*ep.spa(i3,i4)*ep.spa(i4,i5)*
ep.spa(i5,i6)*ep.spb(i3,i2))+(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spb(i6,i3),3))/
 (pow(ep.spa(i4,i5),2)*ep.spb(i2,i1)*ep.spb(i3,i2)*ep.spb(i6,i1)) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spa(i1,i2),3)*pow(ep.spb(i6,i4),2)*(ep.s(i4,i5)+ep.s(i5,i6)))/
 (ep.s(i1,i2,i3)*ep.spa(i2,i3)*ep.spa(i4,i5)*ep.spa(i5,i6)*(-(ep.spa(i1,i2)*ep.spb(i4,i2)) -
 ep.spa(i1,i3)*ep.spb(i4,i3))*(ep.spa(i1,i3)*ep.spb(i6,i1)+ep.spa(i2,i3)*ep.spb(i6,i2))) +
(complex<T>(-1,0)/complex<T>(3,0)*ep.s(i3,i5)*(ep.spa(i1,i4)*ep.spb(i3,i1)+ep.spa(i2,i4)*ep.spb(i3,i2))*
(ep.spa(i1,i4)*ep.spb(i6,i1)+ep.spa(i2,i4)*ep.spb(i6,i2))*(ep.spa(i1,i5)*ep.spb(i6,i1) +
 ep.spa(i2,i5)*ep.spb(i6,i2)))/(pow(ep.spa(i3,i4),2)*pow(ep.spa(i4,i5),2)*ep.spb(i2,i1)*
(ep.spa(i1,i6)*ep.spb(i3,i1)+ep.spa(i2,i6)*ep.spb(i3,i2))*(-(ep.spa(i3,i5)*ep.spb(i3,i2)) -
 ep.spa(i4,i5)*ep.spb(i4,i2))*ep.spb(i6,i1)) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spa(i1,i4)*ep.spb(i6,i1)+ep.spa(i2,i4)*ep.spb(i6,i2),2)*ep.spa(i3,i5)*
ep.spb(i6,i3))/(pow(ep.spa(i3,i4),2)*pow(ep.spa(i4,i5),2)*ep.spb(i2,i1)*
(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*ep.spb(i6,i1)) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spa(i1,i5),2)*pow(ep.spb(i4,i3),2)*
(-(ep.spa(i1,i6)*ep.spa(i4,i5)*ep.spb(i4,i3))-ep.spa(i5,i6)*
  (-(ep.spa(i1,i2)*ep.spb(i3,i2))+ep.spa(i1,i4)*ep.spb(i4,i3)))*ep.spb(i6,i5))/
 (pow(ep.spa(i5,i6),2)*ep.s(i2,i3,i4)*ep.spa(i4,i5)*ep.spb(i3,i2)*
(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*
(-(ep.spa(i1,i2)*ep.spb(i4,i2))-ep.spa(i1,i3)*ep.spb(i4,i3))) +
(complex<T>(1,0)/complex<T>(6,0)*pow(ep.spa(i1,i5)*ep.spb(i6,i1)+ep.spa(i2,i5)*ep.spb(i6,i2),2)*
(ep.s(i1,i2)*ep.spa(i4,i5)+complex<T>(-2,0)*(ep.s(i1,i5)*ep.spa(i4,i5) +
   ep.s(i2,i5)*ep.spa(i4,i5)-ep.spa(i1,i5)*ep.spa(i3,i4)*ep.spb(i3,i1) -
   ep.spa(i2,i5)*ep.spa(i3,i4)*ep.spb(i3,i2)))*ep.spb(i6,i5))/
 (pow(ep.spa(i4,i5),2)*ep.spa(i3,i4)*ep.spa(i5,i6)*ep.spb(i2,i1)*
(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*ep.spb(i6,i1)*
(ep.spa(i1,i3)*ep.spb(i6,i1)+ep.spa(i2,i3)*ep.spb(i6,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*pow(ep.spa(i1,i2),3)*ep.spa(i3,i5)*ep.spb(i6,i4)*ep.spb(i6,i5))/
 (ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i6)*(-(ep.spa(i1,i2)*ep.spb(i4,i2)) -
 ep.spa(i1,i3)*ep.spb(i4,i3))*(ep.spa(i1,i3)*ep.spb(i6,i1)+ep.spa(i2,i3)*ep.spb(i6,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*ep.spa(i1,i2)*ep.spa(i1,i5)*ep.spb(i4,i3)*(ep.s(i3,i5)*ep.spa(i1,i5) +
 ep.s(i4,i5)*ep.spa(i1,i5)+ep.spa(i1,i6)*ep.spa(i3,i5)*ep.spb(i6,i3) +
 ep.spa(i1,i6)*ep.spa(i4,i5)*ep.spb(i6,i4))*ep.spb(i6,i5))/(ep.s(i2,i3,i4)*ep.spa(i3,i4)*
ep.spa(i4,i5)*ep.spa(i5,i6)*(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*
(-(ep.spa(i1,i2)*ep.spb(i4,i2))-ep.spa(i1,i3)*ep.spb(i4,i3))) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spa(i1,i2),3)*pow(ep.spb(i5,i3),2)*(-ep.s(i3,i4)-ep.s(i4,i5)))/
 (ep.s(i3,i4,i5)*ep.spa(i1,i6)*ep.spa(i3,i4)*ep.spa(i4,i5)*(ep.spa(i1,i6)*ep.spb(i3,i1) +
 ep.spa(i2,i6)*ep.spb(i3,i2))*(ep.spa(i1,i2)*ep.spb(i5,i1)+ep.spa(i2,i6)*ep.spb(i6,i5))) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spa(i1,i2),2)*pow(ep.spb(i5,i3),2)*(-(ep.s(i3,i5)*ep.spa(i1,i5)) +
 ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spb(i3,i2)+ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spb(i4,i2)))/
 (ep.spa(i1,i6)*ep.spa(i3,i4)*ep.spa(i4,i5)*(ep.spa(i1,i6)*ep.spb(i3,i1) +
 ep.spa(i2,i6)*ep.spb(i3,i2))*(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*
(ep.spa(i1,i2)*ep.spb(i5,i1)+ep.spa(i2,i6)*ep.spb(i6,i5))) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spb(i5,i3),2)*ep.spa(i1,i2)*ep.spa(i2,i4)*ep.spa(i3,i5)*
(ep.spa(i1,i5)*ep.spb(i6,i1)+ep.spa(i2,i5)*ep.spb(i6,i2))*ep.spb(i6,i5))/
 (pow(ep.spa(i3,i4),2)*ep.spa(i4,i5)*(ep.spa(i1,i6)*ep.spb(i3,i1) +
 ep.spa(i2,i6)*ep.spb(i3,i2))*(-(ep.spa(i3,i5)*ep.spb(i3,i2))-ep.spa(i4,i5)*ep.spb(i4,i2))*
ep.spb(i6,i1)*(ep.spa(i1,i2)*ep.spb(i5,i1)+ep.spa(i2,i6)*ep.spb(i6,i5))) +
(complex<T>(1,0)/complex<T>(6,0)*pow(ep.spa(i1,i2),2)*(-((ep.spa(i1,i4)*ep.spb(i4,i3))/ep.spb(i3,i2)) -
 (ep.spa(i2,i5)*ep.spb(i6,i5))/ep.spb(i6,i1)))/(ep.spa(i1,i6)*ep.spa(i2,i3)*ep.spa(i3,i4)*
ep.spa(i4,i5)*ep.spa(i5,i6)) +
(complex<T>(-1,0)/complex<T>(6,0)*(-(((-pow(ep.s(i2,i3,i4),2)+pow(ep.spa(i2,i3),2)*
pow(ep.spb(i3,i2),2))*ep.spa(i1,i4)*ep.spa(i2,i4)*ep.spb(i4,i3)*
(-(ep.spa(i1,i2)*ep.spb(i4,i2))-ep.spa(i1,i3)*ep.spb(i4,i3))*
(-(ep.spa(i1,i4)*ep.spa(i2,i3)*ep.spb(i4,i3))+ep.spa(i2,i4)*
(-(ep.spa(i1,i2)*ep.spb(i4,i2))-ep.spa(i1,i3)*ep.spb(i4,i3))))/
   (pow(-ep.s(i2,i4)-ep.s(i3,i4),3)*ep.s(i2,i3,i4)*ep.spb(i3,i2))) -
 (((pow(ep.spa(i1,i6),2)*pow(ep.spb(i6,i1),2))/ep.s(i2,i3,i4)-ep.s(i2,i3,i4))*
   ep.spa(i1,i5)*ep.spa(i2,i5)*ep.spb(i6,i5)*(ep.spa(i1,i2)*ep.spb(i5,i1) +
ep.spa(i2,i6)*ep.spb(i6,i5))*(ep.spa(i1,i6)*ep.spa(i2,i5)*ep.spb(i6,i5) +
ep.spa(i1,i5)*(ep.spa(i1,i2)*ep.spb(i5,i1)+ep.spa(i2,i6)*ep.spb(i6,i5))))/
  (pow(ep.s(i1,i6)-ep.s(i2,i3,i4),3)*ep.spb(i6,i1))))/(ep.spa(i1,i6)*ep.spa(i2,i3)*
ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i6))
#endif
);


}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gmppmpp(const eval_param<T>& ep,const mass_param_coll& mpc){


const complex<T> INT1(1,0);
const complex<T> INT2(2,0);
const complex<T> INT3(3,0);
const complex<T> INT4(4,0);
const complex<T> INT5(5,0);
const complex<T> INT6(6,0);
const complex<T> INT8(8,0);
const complex<T> INT9(9,0);

{
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb16=SPB(i1,i6);
const complex<T> spb14=SPB(i1,i4);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spba_6_1_4=spa14*spb61;
const complex<T> spba_6_12_4=spa14*spb61+spa24*spb62;
const complex<T> spba_6_123_4=spa14*spb61+spa24*spb62+spa34*spb63;
const complex<T> spba_5_16_4=spa14*spb51+spa46*spb65;
const complex<T> spba_5_126_4=spa14*spb51+spa24*spb52+spa46*spb65;
const complex<T> spba_5_1236_4=spa14*spb51+spa24*spb52+spa34*spb53+spa46*spb65;
const complex<T> spba_3_456_4=spa45*spb53+spa46*spb63;
const complex<T> spba_3_45_4=spa45*spb53;
const complex<T> spba_2_3456_4=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spba_2_345_4=-(spa34*spb32)+spa45*spb52;
const complex<T> spba_2_34_4=-(spa34*spb32);
const complex<T> spab_6_35_4=spa36*spb43-spa56*spb54;
const complex<T> spab_6_14_6=spa16*spb61+spa46*spb64;
const complex<T> spab_6_14_2=spa16*spb21-spa46*spb42;
const complex<T> spab_5_36_4=spa35*spb43+spa56*spb64;
const complex<T> spab_5_23_6=spa25*spb62+spa35*spb63;
const complex<T> spab_5_23_4=spa25*spb42+spa35*spb43;
const complex<T> spab_5_14_6=spa15*spb61+spa45*spb64;
const complex<T> spab_5_14_2=spa15*spb21-spa45*spb42;
const complex<T> spab_4_456_3=spa45*spb53+spa46*spb63;
const complex<T> spab_4_45_3=spa45*spb53;
const complex<T> spab_4_3456_2=-(spa34*spb32)+spa45*spb52+spa46*spb62;
const complex<T> spab_4_345_2=-(spa34*spb32)+spa45*spb52;
const complex<T> spab_4_34_2=-(spa34*spb32);
const complex<T> spab_4_1_6=spa14*spb61;
const complex<T> spab_4_16_5=spa14*spb51+spa46*spb65;
const complex<T> spab_4_12_6=spa14*spb61+spa24*spb62;
const complex<T> spab_4_126_5=spa14*spb51+spa24*spb52+spa46*spb65;
const complex<T> spab_4_123_6=spa14*spb61+spa24*spb62+spa34*spb63;
const complex<T> spab_4_1236_5=spa14*spb51+spa24*spb52+spa34*spb53+spa46*spb65;
const complex<T> spab_3_56_4=spa35*spb54+spa36*spb64;
const complex<T> spab_3_56_2=spa35*spb52+spa36*spb62;
const complex<T> spab_3_25_4=spa23*spb42+spa35*spb54;
const complex<T> spab_3_14_6=spa13*spb61-spa34*spb64;
const complex<T> spab_3_14_2=spa13*spb21+spa34*spb42;
const complex<T> spab_2_35_4=-(spa23*spb43)+spa25*spb54;
const complex<T> spab_2_14_6=spa12*spb61-spa24*spb64;
const complex<T> spab_2_14_2=spa12*spb21+spa24*spb42;
const complex<T> spab_1_56_4=spa15*spb54+spa16*spb64;
const complex<T> spab_1_5_3=spa15*spb53;
const complex<T> spab_1_46_5=-(spa14*spb54)+spa16*spb65;
const complex<T> spab_1_4_5=-(spa14*spb54);
const complex<T> spab_1_45_6=-(spa14*spb64)-spa15*spb65;
const complex<T> spab_1_456_6=-(spa14*spb64)-spa15*spb65;
const complex<T> spab_1_456_4=spa15*spb54+spa16*spb64;
const complex<T> spab_1_456_3=spa14*spb43+spa15*spb53+spa16*spb63;
const complex<T> spab_1_45_5=-(spa14*spb54);
const complex<T> spab_1_45_4=spa15*spb54;
const complex<T> spab_1_45_3=spa14*spb43+spa15*spb53;
const complex<T> spab_1_4_3=spa14*spb43;
const complex<T> spab_1_3_5=-(spa13*spb53);
const complex<T> spab_1_35_4=-(spa13*spb43)+spa15*spb54;
const complex<T> spab_1_34_5=-(spa13*spb53)-spa14*spb54;
const complex<T> spab_1_345_6=-(spa13*spb63)-spa14*spb64-spa15*spb65;
const complex<T> spab_1_3456_6=-(spa13*spb63)-spa14*spb64-spa15*spb65;
const complex<T> spab_1_3456_4=-(spa13*spb43)+spa15*spb54+spa16*spb64;
const complex<T> spab_1_3456_3=spa14*spb43+spa15*spb53+spa16*spb63;
const complex<T> spab_1_3456_2=spa13*spb32+spa14*spb42+spa15*spb52+spa16*spb62;
const complex<T> spab_1_345_5=-(spa13*spb53)-spa14*spb54;
const complex<T> spab_1_345_4=-(spa13*spb43)+spa15*spb54;
const complex<T> spab_1_345_3=spa14*spb43+spa15*spb53;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_4=-(spa13*spb43);
const complex<T> spab_1_34_3=spa14*spb43;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spab_1_2345_6=-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65;
const complex<T> spab_1_2345_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spab_1_2345_4=-(spa12*spb42)-spa13*spb43+spa15*spb54;
const complex<T> spab_1_2345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_234_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_16_5=spa16*spb65;
const complex<T> spab_1_12_6=-(spa12*spb62);
const complex<T> spab_1_126_5=-(spa12*spb52)+spa16*spb65;
const complex<T> spab_1_123_6=-(spa12*spb62)-spa13*spb63;
const complex<T> spab_1_1236_5=-(spa12*spb52)-spa13*spb53+spa16*spb65;
const complex<T> spaa_1_6_5_3=-(spa16*spa35*spb65);
const complex<T> spaa_1_6_3_5=spa16*spa35*spb63;
const complex<T> spaa_1_6_25_3=spa16*(spa23*spb62-spa35*spb65);
const complex<T> spaa_1_6_2_3=spa16*spa23*spb62;
const complex<T> spaa_1_6_23_5=spa16*(spa25*spb62+spa35*spb63);
const complex<T> spaa_1_56_2_3=spa23*(spa15*spb52+spa16*spb62);
const complex<T> spaa_1_2_6_5=spa12*spa56*spb62;
const complex<T> spaa_1_26_5_3=-(spa35*(-(spa12*spb52)+spa16*spb65));
const complex<T> spaa_1_26_3_5=spa35*(-(spa12*spb32)+spa16*spb63);
const complex<T> spaa_1_26_35_6=spa12*(-(spa36*spb32)-spa56*spb52)+spa16*(spa36*spb63+spa56*spb65);
const complex<T> spaa_1_26_35_2=spa12*(spa23*spb32+spa25*spb52)+spa16*(-(spa23*spb63)-spa25*spb65);
const complex<T> spaa_1_2_56_3=spa12*(spa35*spb52+spa36*spb62);
const complex<T> spaa_1_2_5_3=spa12*spa35*spb52;
const complex<T> spaa_1_23_6_5=-(spa56*(-(spa12*spb62)-spa13*spb63));
const complex<T> spaa_1_2_36_5=spa12*(-(spa35*spb32)+spa56*spb62);
const complex<T> spaa_1_2_3_5=-(spa12*spa35*spb32);
const complex<T> s56=(spa56*spb65);
const complex<T> s46=(spa46*spb64);
const complex<T> s45=(spa45*spb54);
const complex<T> s35=(spa35*spb53);
const complex<T> s34=(spa34*spb43);
const complex<T> s26=(spa26*spb62);
const complex<T> s24=(spa24*spb42);
const complex<T> s23=(spa23*spb32);
const complex<T> s16=(spa16*spb61);
const complex<T> s15=(spa15*spb51);
const complex<T> s14=(spa14*spb41);
const complex<T> s13=(spa13*spb31);
const complex<T> s12=(spa12*spb21);

return(
complex<T>(0,1)*(-((INT1*(s12+s14+s24)*spa13*spa14*spab_3_14_2*spb23)/(INT6*s12*spa23*spa36*spa56*spab_3_56_4*spab_5_14_2))-(INT1*(s14+s16+s46)*spa14*spa15*spab_5_14_6*spb56)/(INT6*s16*spa23*spa25*spa56*spab_3_14_6*spab_5_23_4)+(INT1*spa14*spa15*(spa14/spa23+(INT3*spa15*spab_1_24_3)/(spa25*spab_1_23_4))*spb23*spb56)/(INT6*(s23+s24+s34)*spa25*spa56*spab_1_34_2)+(INT1*spa13*spa14*(-(spa14/spa56)-(INT3*spa13*spab_1_46_5)/(spa36*spab_1_56_4))*spb23*spb56)/(INT6*(s45+s46+s56)*spa23*spa36*spab_1_45_6)+(INT1*(s23+s24+INT2*s34)*spa12*spa24*spab_1_34_2*spab_4_34_2*(-(spa24*spab_1_34_2)+spa12*spba_2_34_4))/(INT6*pow(-s23-s24,2)*s34*(s23+s24+s34)*spa16*spa23*spa25*spa34*spa56)+(INT1*(s12+s45+s46+s56)*spa13*spa34*spab_1_456_3*spab_4_456_3*(-(spa34*spab_1_456_3)+spa13*spba_3_456_4))/(INT6*s12*(s45+s46+s56)*pow(-s12+s45+s46+s56,2)*spa12*spa23*spa36*spa45*spa56)+(INT1*(s15+INT2*s16+s56)*spa15*spa45*spab_1_16_5*spab_4_16_5*(spa45*spab_1_16_5+spa15*spba_5_16_4))/(INT6*s16*pow(-s15-s56,2)*(s15+s16+s56)*spa16*spa23*spa25*spa34*spa56)+(INT1*(s12+s13+s23+s45)*spa16*spa46*spab_1_123_6*spab_4_123_6*(spa46*spab_1_123_6+spa16*spba_6_123_4))/(INT6*(s12+s13+s23)*pow(s12+s13+s23-s45,2)*s45*spa12*spa23*spa36*spa45*spa56)-(INT1*spab_3_56_2*spb23*((spa13*spa16)/spa36+(spa12*spb24)/spb14)*pow(s12+s14+s24,2))/(INT3*spa12*spa23*spa56*spab_3_56_4*spab_5_14_2*spab_6_35_4*spb12)+(INT1*spab_5_23_6*(-((spa12*spa15)/spa25)-(spa16*spb46)/spb14)*spb56*pow(s14+s16+s46,2))/(INT3*spa16*spa23*spa56*spab_2_35_4*spab_3_14_6*spab_5_23_4*spb16)-(INT1*spab_1_34_5*((spa16*spa35*spaa_1_26_35_6)/spa36+spa15*spa36*spab_1_45_3)*spb35*pow(spa13,2))/(INT2*(s34+s35+s45)*spa12*spa16*spa23*spa35*spa36*spa56*spab_1_35_4*spb45)+(INT1*(s34+s35+INT2*s45)*spa34*spab_1_45_3*spab_4_45_3*(-(spa34*spab_1_45_3)+spa13*spba_3_45_4)*pow(spa13,2))/(INT6*pow(-s34-s35,2)*s45*(s34+s35+s45)*spa12*spa16*spa23*spa35*spa36*spa45)-(INT1*spa12*(INT1+(INT3*spa15*spa23)/(spa13*spa25))*spb23*pow(spa14,2))/(INT6*spa16*spa23*spa25*spa34*spa56*spb34)-(INT1*spa16*(INT1+(INT3*spa13*spa56)/(spa15*spa36))*spb56*pow(spa14,2))/(INT6*spa12*spa23*spa36*spa45*spa56*spb45)+(INT1*(INT1+(INT3*spa16*spa35)/(spa15*spa36))*spb35*pow(spa13,2)*pow(spa14,2))/(INT6*spa12*spa16*spa23*spa35*spa36*spa45*spb45)-(INT8*spa13*pow(spa14,3))/(INT9*spa12*spa16*spa23*spa34*spa35*spa56)-(INT8*spa15*pow(spa14,3))/(INT9*spa12*spa16*spa23*spa35*spa45*spa56)+(INT2*pow(spa14,4))/(INT3*spa12*spa16*spa23*spa34*spa45*spa56)-(INT1*((spa12*spa35*spaa_1_26_35_2)/spa25-spa13*spa25*spab_1_34_5)*spab_1_45_3*spb35*pow(spa15,2))/(INT2*(s34+s35+s45)*spa12*spa16*spa23*spa25*spa35*spa56*spab_1_35_4*spb34)+(INT1*(s12+s16+s26+s34)*spa45*spab_1_126_5*spab_4_126_5*(spa45*spab_1_126_5+spa15*spba_5_126_4)*pow(spa15,2))/(INT6*(s12+s16+s26)*pow(s12+s16+s26-s34,2)*s34*spa12*spa16*spa25*spa34*spa35*spa56)+(INT1*(INT1+(INT3*spa12*spa35)/(spa13*spa25))*spb35*pow(spa14,2)*pow(spa15,2))/(INT6*spa12*spa16*spa25*spa34*spa35*spa56*spb34)+(INT1*(s16+s34+s35+s45)*spa12*spab_1_345_2*spab_4_345_2*(-(spa24*spab_1_345_2)+spa12*spba_2_345_4)*pow(spa24,2))/(INT6*s16*(s34+s35+s45)*pow(-s16+s34+s35+s45,2)*spa16*spa23*spa25*spa26*spa34*spa45)-(INT1*spa12*spa14*spa15*(-(spa12*spab_1_2345_2)-spa15*spab_1_2345_5))/(INT2*spa13*spa16*spa23*spa56*spab_1_2345_4*pow(spa25,2))+(INT1*spa12*spa14*spa15*(-(spa12*spab_1_234_2)+spa15*spab_1_234_5))/(INT2*spa13*spa16*spa23*spa56*spab_1_234_4*pow(spa25,2))-(INT1*spa12*spa14*spa15*(-(spa12*spab_1_345_2)+spa15*spab_1_345_5))/(INT2*spa13*spa16*spa23*spa56*spab_1_345_4*pow(spa25,2))+(INT1*spa12*spa14*spa15*(-(spa12*spab_1_34_2)-spa15*spab_1_34_5))/(INT2*spa13*spa16*spa23*spa56*spab_1_34_4*pow(spa25,2))+(INT1*(s15+INT2*s16+s56)*spa24*spa45*spab_1_16_5*spab_4_16_5*pow(spa15,2))/(INT2*s16*(s15+s56)*(s15+s16+s56)*spa16*spa23*spa34*spa56*pow(spa25,2))+(INT1*(s14+s16+s46)*spab_2_14_6*spab_5_14_6*spb56*pow(spa15,2))/(INT2*s16*spa23*spa56*spab_3_14_6*spab_5_23_4*spb46*pow(spa25,2))+(INT1*spa15*spb25*((-(INT3*spa15*spa23)+spa13*spa25)*(spaa_1_56_2_3+spaa_1_6_25_3)-INT3*spa12*spa15*spa23*spa35*spb25)*pow(spa14,3))/(INT6*spa13*spa23*spa34*spa56*spab_1_34_2*spab_3_14_6*pow(spa16,2)*pow(spa25,2))-(INT1*(s23+s24+INT2*s34)*spa12*spa15*spab_1_34_2*spab_4_34_2*pow(spa24,2))/(INT2*(s23+s24)*s34*(s23+s24+s34)*spa16*spa23*spa34*spa56*pow(spa25,2))-(INT1*(s16+s34+s35+s45)*spa24*spab_1_345_2*spab_4_345_2*pow(spa14,2)*(-((INT2*spa12*spa15*spa24*spa45)/(pow(spa14,2)*pow(spa25,2)))+(INT2*spa12*spa16*spa24*spa46)/(pow(spa14,2)*pow(spa26,2))))/(INT4*s16*(s16-s34-s35-s45)*(s34+s35+s45)*spa16*spa23*spa34*spa45*spa56)+(INT1*spa12*spa14*spa16*(-(spa12*spab_1_2345_2)+spa16*spab_1_2345_6))/(INT2*spa13*spa15*spa23*spa56*spab_1_2345_4*pow(spa26,2))-(INT1*spa12*spa14*spa16*(-(spa12*spab_1_3456_2)+spa16*spab_1_3456_6))/(INT2*spa13*spa15*spa23*spa56*spab_1_3456_4*pow(spa26,2))+(INT1*spa12*spa14*spa16*(-(spa12*spab_1_345_2)-spa16*spab_1_345_6))/(INT2*spa13*spa15*spa23*spa56*spab_1_345_4*pow(spa26,2))-(INT2*pow(spa14,4)*(-((INT1*spa12*spa16*spa24*spa46)/(INT2*pow(spa14,2)*pow(spa26,2)))+(-((INT1*spa16*spa24*(-(spa24*spab_1_3456_2)+spa12*spba_2_3456_4))/(INT6*spa26*pow(spa14,2)))-(spa16*spa46*(spa46*spab_1_3456_2+spa16*spab_4_3456_2)*pow(spa12,2)*pow(spa24,2))/(pow(spa14,4)*pow(spa26,3)))/s12))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa14,4)*(-((INT1*spa12*spa16*spa24*spa46)/(INT2*pow(spa14,2)*pow(spa26,2)))+(-((INT1*spa16*spa24*(-(spa24*spab_1_3456_2)+spa12*spba_2_3456_4))/(INT6*spa26*pow(spa14,2)))-(spa16*spa46*(spa46*spab_1_3456_2+spa16*spab_4_3456_2)*pow(spa12,2)*pow(spa24,2))/(pow(spa14,4)*pow(spa26,3)))/s12))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT1*(s12+s16+s26+s34)*spa15*spab_1_126_5*spab_4_126_5*pow(spa14,2)*(-((INT2*spa12*spa15*spa24*spa45)/(pow(spa14,2)*pow(spa25,2)))+(INT2*spa13*spa15*spa34*spa45)/(pow(spa14,2)*pow(spa35,2))))/(INT4*(s12+s16+s26)*s34*(-s12-s16-s26+s34)*spa12*spa16*spa23*spa34*spa56)-(INT1*spa13*spa14*spa15*(-(spa13*spab_1_345_3)-spa15*spab_1_345_5))/(INT2*spa12*spa16*spa23*spa56*spab_1_345_4*pow(spa35,2))+(INT1*spa13*spa14*spa15*(-(spa13*spab_1_34_3)+spa15*spab_1_34_5))/(INT2*spa12*spa16*spa23*spa56*spab_1_34_4*pow(spa35,2))-(INT1*spa13*spa14*spa15*(-(spa13*spab_1_45_3)+spa15*spab_1_45_5))/(INT2*spa12*spa16*spa23*spa56*spab_1_45_4*pow(spa35,2))+(INT1*spa36*(INT2*spa13*spa26*(spaa_1_26_5_3+spaa_1_2_56_3)-INT3*spa16*spa23*(INT2*spaa_1_2_56_3+spaa_1_6_5_3))*spb26*pow(spa14,3))/(INT6*spa12*spa13*spa34*spa56*spab_1_34_5*spab_3_14_2*pow(spa26,2)*pow(spa35,2))+(INT1*spa25*(-(INT2*spa15*spa26*(spaa_1_26_3_5+spaa_1_6_23_5))+INT3*spa12*spa56*(spaa_1_2_3_5+INT2*spaa_1_6_23_5))*spb26*pow(spa14,3))/(INT6*spa15*spa16*spa23*spa45*spab_1_45_3*spab_5_14_6*pow(spa26,2)*pow(spa35,2))-(INT2*pow(spa14,4)*(-((INT1*spa13*spa15*spa34*spa45)/(INT2*pow(spa14,2)*pow(spa35,2)))+((INT1*spa13*spa34*spa45*spab_1_4_3)/(INT6*spa35*pow(spa14,2))+(spa15*pow(spa45,2)*spab_1_4_3*pow(spa13,2)*pow(spa34,2))/(pow(spa14,4)*pow(spa35,3)))/s34))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa14,4)*(-((INT1*spa13*spa15*spa34*spa45)/(INT2*pow(spa14,2)*pow(spa35,2)))+((INT1*spa13*spa34*spa45*spab_1_4_3)/(INT6*spa35*pow(spa14,2))+(spa15*pow(spa45,2)*spab_1_4_3*pow(spa13,2)*pow(spa34,2))/(pow(spa14,4)*pow(spa35,3)))/s34))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT1*(INT2*s12+s16+s26)*spa46*spab_1_12_6*spab_4_12_6*pow(spa14,2)*(-((INT2*spa12*spa16*spa24*spa46)/(pow(spa14,2)*pow(spa26,2)))+(INT2*spa13*spa16*spa34*spa46)/(pow(spa14,2)*pow(spa36,2))))/(INT4*s12*(s16+s26)*(s12+s16+s26)*spa12*spa23*spa34*spa45*spa56)-(INT1*(s34+s35+INT2*s45)*spa13*spab_1_45_3*spab_4_45_3*pow(spa14,2)*(-((INT2*spa13*spa15*spa34*spa45)/(pow(spa14,2)*pow(spa35,2)))+(INT2*spa13*spa16*spa34*spa46)/(pow(spa14,2)*pow(spa36,2))))/(INT4*(s34+s35)*s45*(s34+s35+s45)*spa12*spa16*spa23*spa45*spa56)-(INT1*spa13*spa14*spa16*(-(spa13*spab_1_3456_3)-spa16*spab_1_3456_6))/(INT2*spa12*spa15*spa23*spa56*spab_1_3456_4*pow(spa36,2))+(INT1*spa13*spa14*spa16*(-(spa13*spab_1_345_3)+spa16*spab_1_345_6))/(INT2*spa12*spa15*spa23*spa56*spab_1_345_4*pow(spa36,2))-(INT1*spa13*spa14*spa16*(-(spa13*spab_1_456_3)+spa16*spab_1_456_6))/(INT2*spa12*spa15*spa23*spa56*spab_1_456_4*pow(spa36,2))+(INT1*spa13*spa14*spa16*(-(spa13*spab_1_45_3)-spa16*spab_1_45_6))/(INT2*spa12*spa15*spa23*spa56*spab_1_45_4*pow(spa36,2))+(INT1*(s12+s45+s46+s56)*spa34*spa46*spab_1_456_3*spab_4_456_3*pow(spa13,2))/(INT2*s12*(s12-s45-s46-s56)*(s45+s46+s56)*spa12*spa23*spa45*spa56*pow(spa36,2))+(INT1*(s12+s14+s24)*spab_3_14_2*spab_6_14_2*spb23*pow(spa13,2))/(INT2*s12*spa23*spa56*spab_3_56_4*spab_5_14_2*spb24*pow(spa36,2))+(INT1*spa13*spb36*((-(spa15*spa36)+INT3*spa13*spa56)*(spaa_1_23_6_5+spaa_1_2_36_5)+INT3*spa13*spa16*spa35*spa56*spb36)*pow(spa14,3))/(INT6*spa15*spa23*spa45*spa56*spab_1_45_6*spab_5_14_2*pow(spa12,2)*pow(spa36,2))-(INT2*pow(spa14,4)*(-((INT1*spa13*spa15*spa34*spa45)/(INT2*pow(spa14,2)*pow(spa35,2)))+((INT1*spa15*spa34*(spa45*spab_1_1236_5+spa15*spba_5_1236_4))/(INT6*spa35*pow(spa14,2))-(spa13*spa34*(-(spa34*spab_1_1236_5)+spa13*spab_4_1236_5)*pow(spa15,2)*pow(spa45,2))/(pow(spa14,4)*pow(spa35,3)))/s45))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa14,4)*(-((INT1*spa13*spa15*spa34*spa45)/(INT2*pow(spa14,2)*pow(spa35,2)))+((INT1*spa15*spa34*(spa45*spab_1_1236_5+spa15*spba_5_1236_4))/(INT6*spa35*pow(spa14,2))-(spa13*spa34*(-(spa34*spab_1_1236_5)+spa13*spab_4_1236_5)*pow(spa15,2)*pow(spa45,2))/(pow(spa14,4)*pow(spa35,3)))/s45))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT1*(INT2*s12+s16+s26)*spa16*spab_1_12_6*spab_4_12_6*(spa46*spab_1_12_6+spa16*spba_6_12_4)*pow(spa46,2))/(INT6*s12*pow(-s16-s26,2)*(s12+s16+s26)*spa12*spa26*spa34*spa36*spa45*spa56)-(INT1*(s12+s13+s23+s45)*spa13*spa16*spab_1_123_6*spab_4_123_6*pow(spa46,2))/(INT2*(s12+s13+s23)*s45*(-s12-s13-s23+s45)*spa12*spa23*spa45*spa56*pow(spa36,2))-(INT2*pow(spa14,4)*(-((INT1*spa12*spa16*spa24*spa46)/(INT2*pow(spa14,2)*pow(spa26,2)))+((INT1*spa12*spa16*spa46*spba_6_1_4)/(INT6*spa26*pow(spa14,2))+(pow(spa12,2)*spa24*spab_4_1_6*pow(spa16,2)*pow(spa46,2))/(pow(spa14,4)*pow(spa26,3)))/s16))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa14,4)*(-((INT1*spa12*spa16*spa24*spa46)/(INT2*pow(spa14,2)*pow(spa26,2)))+((INT1*spa12*spa16*spa46*spba_6_1_4)/(INT6*spa26*pow(spa14,2))+(pow(spa12,2)*spa24*spab_4_1_6*pow(spa16,2)*pow(spa46,2))/(pow(spa14,4)*pow(spa26,3)))/s16))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT1*spa12*spa14*spb26*spb35*pow(spaa_1_26_35_2,2))/(INT6*(s34+s35+s45)*spa23*spa25*spa26*spa35*spab_1_34_5*spab_1_45_3*spab_2_35_4)-(INT1*spa12*(-((spa16*spa25*spaa_1_26_35_6)/spa26)-(spa15*spa26*spa35*spab_1_45_3)/spa25)*spb26*spb35*pow(spaa_1_26_35_2,2))/(INT2*(s34+s35+s45)*spa23*spa25*spa26*spa35*spa56*spab_1_34_5*spab_1_35_4*spab_1_45_3*spab_2_35_4)+(INT1*spa14*spa16*spb26*spb35*pow(spaa_1_26_35_6,2))/(INT6*(s34+s35+s45)*spa26*spa35*spa36*spa56*spab_1_34_5*spab_1_45_3*spab_6_35_4)+(INT1*spa16*(-((spa12*spa36*spaa_1_26_35_2)/spa26)+(spa13*spa26*spa35*spab_1_34_5)/spa36)*spb26*spb35*pow(spaa_1_26_35_6,2))/(INT2*(s34+s35+s45)*spa23*spa26*spa35*spa36*spa56*spab_1_34_5*spab_1_35_4*spab_1_45_3*spab_6_35_4)+(INT1*spa12*spa15*((spa12*pow(spab_1_2345_2,2))/spb24-(spa15*pow(spab_1_2345_5,2))/spb45))/(INT2*spa13*spa16*spa23*spa56*spab_1_2345_4*pow(spa25,2))-(INT1*spa12*spa16*((spa12*pow(spab_1_2345_2,2))/spb24+(spa16*pow(spab_1_2345_6,2))/spb46))/(INT2*spa13*spa15*spa23*spa56*spab_1_2345_4*pow(spa26,2))-(INT1*spa12*spa15*((spa12*pow(spab_1_234_2,2))/spb24+(spa15*pow(spab_1_234_5,2))/spb45))/(INT2*spa13*spa16*spa23*spa56*spab_1_234_4*pow(spa25,2))+(INT1*spa12*spa15*spb23*pow(spab_1_24_3,2))/(INT2*(s23+s24+s34)*spa16*spa56*spab_1_23_4*spb34*pow(spa25,2))+(INT1*spa13*spa16*((spa13*pow(spab_1_3456_3,2))/spb34-(spa16*pow(spab_1_3456_6,2))/spb46))/(INT2*spa12*spa15*spa23*spa56*spab_1_3456_4*pow(spa36,2))+(INT1*spa12*spa16*((spa12*pow(spab_1_3456_2,2))/spb24+(spa16*pow(spab_1_3456_6,2))/spb46))/(INT2*spa13*spa15*spa23*spa56*spab_1_3456_4*pow(spa26,2))+(INT1*spa13*spa15*((spa13*pow(spab_1_345_3,2))/spb34-(spa15*pow(spab_1_345_5,2))/spb45))/(INT2*spa12*spa16*spa23*spa56*spab_1_345_4*pow(spa35,2))+(INT1*spa12*spa15*((spa12*pow(spab_1_345_2,2))/spb24+(spa15*pow(spab_1_345_5,2))/spb45))/(INT2*spa13*spa16*spa23*spa56*spab_1_345_4*pow(spa25,2))-(INT1*spa12*spa16*((spa12*pow(spab_1_345_2,2))/spb24-(spa16*pow(spab_1_345_6,2))/spb46))/(INT2*spa13*spa15*spa23*spa56*spab_1_345_4*pow(spa26,2))-(INT1*spa13*spa16*((spa13*pow(spab_1_345_3,2))/spb34+(spa16*pow(spab_1_345_6,2))/spb46))/(INT2*spa12*spa15*spa23*spa56*spab_1_345_4*pow(spa36,2))-(INT1*spa12*spa15*((spa12*pow(spab_1_34_2,2))/spb24-(spa15*pow(spab_1_34_5,2))/spb45))/(INT2*spa13*spa16*spa23*spa56*spab_1_34_4*pow(spa25,2))-(INT1*spa13*spa15*((spa13*pow(spab_1_34_3,2))/spb34+(spa15*pow(spab_1_34_5,2))/spb45))/(INT2*spa12*spa16*spa23*spa56*spab_1_34_4*pow(spa35,2))+(INT1*spa13*spa16*((spa13*pow(spab_1_456_3,2))/spb34+(spa16*pow(spab_1_456_6,2))/spb46))/(INT2*spa12*spa15*spa23*spa56*spab_1_456_4*pow(spa36,2))+(INT1*spa13*spa15*((spa13*pow(spab_1_45_3,2))/spb34+(spa15*pow(spab_1_45_5,2))/spb45))/(INT2*spa12*spa16*spa23*spa56*spab_1_45_4*pow(spa35,2))-(INT1*spa13*spa16*((spa13*pow(spab_1_45_3,2))/spb34-(spa16*pow(spab_1_45_6,2))/spb46))/(INT2*spa12*spa15*spa23*spa56*spab_1_45_4*pow(spa36,2))-(INT1*spa13*spa16*spb56*pow(spab_1_46_5,2))/(INT2*(s45+s46+s56)*spa12*spa23*spab_1_56_4*spb45*pow(spa36,2))+(INT1*(s14+s16+s46)*spa12*spb26*((spa15*spa26*spab_5_14_6)/spa25-(spa16*spa25*spab_6_14_6)/spa26+(INT1*spa14*spa56*spb46)/INT3)*pow(spab_2_14_6,2))/(INT2*s16*spa23*spa25*spa26*spa56*spab_2_35_4*spab_3_14_6*spab_5_14_6*spb46)+(INT1*spab_3_25_4*spb23*spb46*pow(spab_2_14_6,3))/(INT3*spab_5_14_6*spab_5_23_4*spb14*spb16*pow(spa23,2)*pow(spab_2_35_4,2))+(INT1*(s12+s14+s24)*spa16*(-((spa12*spa36*spab_2_14_2)/spa26)+(spa13*spa26*spab_3_14_2)/spa36+(INT1*spa14*spa23*spb24)/INT3)*spb26*pow(spab_6_14_2,2))/(INT2*s12*spa23*spa26*spa36*spa56*spab_3_14_2*spab_5_14_2*spab_6_35_4*spb24)+(INT1*spab_5_36_4*spb24*spb56*pow(spab_6_14_2,3))/(INT3*spab_3_14_2*spab_3_56_4*spb12*spb14*pow(spa56,2)*pow(spab_6_35_4,2))-(INT1*(INT3*spa14+(INT3*spa15*spa23*spab_1_24_3)/(spa25*spab_1_23_4)-(INT2*spa15*spa26*spb23)/(spa56*spb34))*spb56*pow(spa15,2)*pow(spb23,2))/(INT6*(s23+s24+s34)*spa25*spa56*spab_1_34_2*spab_5_23_4)-(INT1*spab_6_14_2*pow(s12+s14+s24,2)*((spab_6_14_2*spb14*spb26*pow(spa16,2))/spab_3_56_2-(spa12*spa26*spa36*spb23*pow(spb24,2))/spab_3_56_4))/(INT3*s12*spa26*spa36*spa56*spab_5_14_2*spb14*pow(spab_6_35_4,2))+(pow(spa12,2)*pow(spa15,2)*pow(spb25,2))/(spa13*spa16*spa23*spa56*spb24*spb45*pow(spa25,2))+(pow(spa12,2)*pow(spa16,2)*pow(spb26,2))/(spa13*spa15*spa23*spa56*spb24*spb46*pow(spa26,2))+(spa12*spa23*pow(spa14,3)*((INT1*spa13*(spaa_1_26_5_3+spaa_1_6_25_3)*(spa15*spa26*spa35*spb25+spa16*spa25*spa36*spb26))/(INT3*spa12*spa23*spa25*spa26)+((spaa_1_26_5_3*spaa_1_2_5_3+spaa_1_2_5_3*spaa_1_6_2_3+spaa_1_6_25_3*spaa_1_6_2_3)*pow(spa13,2))/(pow(spa12,2)*pow(spa23,2))-(INT1*spa16*spa35*spb25*spb26*((spa13*spa16)/(spa12*spa26)+(spa13*spa35)/(spa23*spa25)+(INT1*pow(spa13,2))/(INT3*spa12*spa23)+(spa15*spa35)/pow(spa25,2)+(spa16*spa36)/pow(spa26,2)))/INT2+((spa13*spa26+spa12*spa36)*spaa_1_6_25_3*spb26*pow(spa16,2))/(pow(spa12,2)*pow(spa26,2))+((spa15*spa23+spa13*spa25)*spaa_1_26_5_3*spb25*pow(spa35,2))/(pow(spa23,2)*pow(spa25,2))-(INT1*spa13*pow(spa15,2)*pow(spa35,2)*pow(spb25,3))/(INT3*spa23*spa25*spab_1_34_2)+(INT1*pow(spa13,2)*pow(spa16,2)*pow(spa35,2)*pow(spb25,2)*pow(spb26,2))/(INT3*spa12*spa23*spab_1_34_2*spab_3_14_2)-(INT1*spa13*pow(spa16,2)*pow(spa36,2)*pow(spb26,3))/(INT3*spa12*spa26*spab_3_14_2)))/(spa13*spa34*spa56*spab_1_34_5*spab_3_14_6*pow(spa16,2)*pow(spa35,2))+(INT1*spa12*spab_1_24_3*spb23*(spab_1_24_3+spa12*spb23))/(INT6*(s23+s24+s34)*spa16*spa23*spa25*spa56*pow(spb34,2))-(INT1*spab_1_45_3*(spab_1_4_3+INT2*spab_1_5_3)*spb35*pow(spa15,2))/(INT6*(s34+s35+s45)*spa12*spa16*spa25*spa35*spa56*pow(spb34,2))+(INT1*spb24*pow(spab_1_24_3,3))/(INT3*(s23+s24+s34)*spa16*spa23*spa56*spab_5_23_4*pow(spb34,2))+(INT1*pow(spa13,2)*pow(spa15,2)*pow(spb35,2))/(spa12*spa16*spa23*spa56*spb34*spb45*pow(spa35,2))-(INT1*spb26*pow(spa12,2)*pow(spaa_1_26_35_2,2)*pow(spb35,2))/(INT3*(s34+s35+s45)*spa23*spa25*spa26*spab_1_34_5*spab_1_45_3*pow(spab_2_35_4,2))-(INT1*spb26*pow(spa16,2)*pow(spaa_1_26_35_6,2)*pow(spb35,2))/(INT3*(s34+s35+s45)*spa26*spa36*spa56*spab_1_34_5*spab_1_45_3*pow(spab_6_35_4,2))-(INT1*pow(spb24,3)*pow(spb35,3))/(INT3*spb12*spb14*spb34*spb45*pow(spab_6_35_4,2))+(INT1*spa26*spb26*pow(spab_1_35_4,4)*pow(spb35,4))/(INT3*(s34+s35+s45)*spab_1_34_5*spab_1_45_3*spb34*spb45*pow(spab_2_35_4,2)*pow(spab_6_35_4,2))+(INT1*pow(spa13,2)*pow(spa16,2)*pow(spb36,2))/(spa12*spa15*spa23*spa56*spb34*spb46*pow(spa36,2))-(spa16*spa56*pow(spa14,3)*(-((INT1*spa15*(spaa_1_26_3_5+spaa_1_2_36_5)*(-(spa12*spa25*spa36*spb26)-spa13*spa26*spa35*spb36))/(INT3*spa16*spa26*spa36*spa56))-((-(spa16*spa25)-spa15*spa26)*spaa_1_2_36_5*spb26*pow(spa12,2))/(pow(spa16,2)*pow(spa26,2))+(INT1*spa12*spa35*spb26*spb36*(-((spa12*spa15)/(spa16*spa26))-(spa15*spa35)/(spa36*spa56)-(INT1*pow(spa15,2))/(INT3*spa16*spa56)-(spa12*spa25)/pow(spa26,2)-(spa13*spa35)/pow(spa36,2)))/INT2+((spaa_1_2_36_5*spaa_1_2_6_5+spaa_1_26_3_5*spaa_1_6_3_5+spaa_1_2_6_5*spaa_1_6_3_5)*pow(spa15,2))/(pow(spa16,2)*pow(spa56,2))-((-(spa15*spa36)-spa13*spa56)*spaa_1_26_3_5*spb36*pow(spa35,2))/(pow(spa36,2)*pow(spa56,2))-(INT1*spa15*pow(spa12,2)*pow(spa25,2)*pow(spb26,3))/(INT3*spa16*spa26*spab_5_14_6)-(INT1*pow(spa12,2)*pow(spa15,2)*pow(spa35,2)*pow(spb26,2)*pow(spb36,2))/(INT3*spa16*spa56*spab_1_45_6*spab_5_14_6)+(INT1*spa15*pow(spa13,2)*pow(spa35,2)*pow(spb36,3))/(INT3*spa36*spa56*spab_1_45_6)))/(spa15*spa23*spa45*spab_1_45_3*spab_5_14_2*pow(spa12,2)*pow(spa35,2))+(INT1*spa16*spab_1_46_5*spb56*(spab_1_46_5-spa16*spb56))/(INT6*(s45+s46+s56)*spa12*spa23*spa36*spa56*pow(spb45,2))-(INT1*spab_1_34_5*(INT2*spab_1_3_5+spab_1_4_5)*spb35*pow(spa13,2))/(INT6*(s34+s35+s45)*spa12*spa16*spa23*spa35*spa36*pow(spb45,2))-(INT1*spb46*pow(spab_1_46_5,3))/(INT3*(s45+s46+s56)*spa12*spa23*spa56*spab_3_56_4*pow(spb45,2))-(INT1*pow(spab_1_35_4,4)*pow(spb35,3))/(INT3*(s34+s35+s45)*spa12*spa16*spa35*spab_2_35_4*spab_6_35_4*pow(spb34,2)*pow(spb45,2))+(INT1*spab_2_14_6*pow(s14+s16+s46,2)*(-((spab_2_14_6*spb14*spb26*pow(spa12,2))/spab_5_23_6)+(spa16*spa25*spa26*spb56*pow(spb46,2))/spab_5_23_4))/(INT3*s16*spa23*spa25*spa26*spab_3_14_6*spb14*pow(spab_2_35_4,2))-(INT1*pow(spb35,3)*pow(spb46,3))/(INT3*spb14*spb16*spb34*spb45*pow(spab_2_35_4,2))+(INT1*spb23*(INT3*spa14+(INT3*spa13*spa56*spab_1_46_5)/(spa36*spab_1_56_4)-(INT2*spa13*spa26*spb56)/(spa23*spb45))*pow(spa13,2)*pow(spb56,2))/(INT6*(s45+s46+s56)*spa23*spa36*spab_1_45_6*spab_3_56_4))
);
}


}

/* testing only
template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gmpmpmp(const eval_param<T>& ep,const mass_param_coll& mpc){
	return complex<T>(0,0);
}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gmmpmpp(const eval_param<T>& ep,const mass_param_coll& mpc){
	return complex<T>(0,0);
}
*/


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gmpmppp(const eval_param<T>& ep,const mass_param_coll& mpc){

const complex<T> INT1(1,0);
const complex<T> INT2(2,0);
const complex<T> INT3(3,0);
const complex<T> INT4(4,0);
const complex<T> INT5(5,0);
const complex<T> INT6(6,0);
const complex<T> INT8(8,0);
const complex<T> INT9(9,0);
{
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb52=SPB(i5,i2);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb36=SPB(i3,i6);
const complex<T> spb35=SPB(i3,i5);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb31=SPB(i3,i1);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb16=SPB(i1,i6);
const complex<T> spb13=SPB(i1,i3);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa46=SPA(i4,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> spba_6_1_3=spa13*spb61;
const complex<T> spba_6_12_3=spa13*spb61+spa23*spb62;
const complex<T> spba_5_16_3=spa13*spb51+spa36*spb65;
const complex<T> spba_5_126_3=spa13*spb51+spa23*spb52+spa36*spb65;
const complex<T> spba_4_156_3=spa13*spb41+spa35*spb54+spa36*spb64;
const complex<T> spba_4_1256_3=spa13*spb41+spa23*spb42+spa35*spb54+spa36*spb64;
const complex<T> spba_2_3456_3=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spba_2_345_3=spa34*spb42+spa35*spb52;
const complex<T> spba_2_34_3=spa34*spb42;
const complex<T> spab_6_45_3=-(spa46*spb43)-spa56*spb53;
const complex<T> spab_6_4_3=-(spa46*spb43);
const complex<T> spab_6_2_3=spa26*spb32;
const complex<T> spab_6_13_6=spa16*spb61+spa36*spb63;
const complex<T> spab_6_13_2=spa16*spb21-spa36*spb32;
const complex<T> spab_6_12_3=spa16*spb31+spa26*spb32;
const complex<T> spab_5_24_3=spa25*spb32-spa45*spb43;
const complex<T> spab_5_13_6=spa15*spb61+spa35*spb63;
const complex<T> spab_5_13_2=spa15*spb21-spa35*spb32;
const complex<T> spab_4_13_6=spa14*spb61+spa34*spb63;
const complex<T> spab_4_13_2=spa14*spb21-spa34*spb32;
const complex<T> spab_3_3456_2=spa34*spb42+spa35*spb52+spa36*spb62;
const complex<T> spab_3_345_2=spa34*spb42+spa35*spb52;
const complex<T> spab_3_34_2=spa34*spb42;
const complex<T> spab_3_1_6=spa13*spb61;
const complex<T> spab_3_16_5=spa13*spb51+spa36*spb65;
const complex<T> spab_3_156_4=spa13*spb41+spa35*spb54+spa36*spb64;
const complex<T> spab_3_12_6=spa13*spb61+spa23*spb62;
const complex<T> spab_3_126_5=spa13*spb51+spa23*spb52+spa36*spb65;
const complex<T> spab_3_1256_4=spa13*spb41+spa23*spb42+spa35*spb54+spa36*spb64;
const complex<T> spab_2_45_6=-(spa24*spb64)-spa25*spb65;
const complex<T> spab_2_45_3=spa24*spb43+spa25*spb53;
const complex<T> spab_2_13_6=spa12*spb61-spa23*spb63;
const complex<T> spab_1_5_4=spa15*spb54;
const complex<T> spab_1_45_3=spa14*spb43+spa15*spb53;
const complex<T> spab_1_4_2=spa14*spb42;
const complex<T> spab_1_35_4=-(spa13*spb43)+spa15*spb54;
const complex<T> spab_1_3_4=-(spa13*spb43);
const complex<T> spab_1_34_5=-(spa13*spb53)-spa14*spb54;
const complex<T> spab_1_345_6=-(spa13*spb63)-spa14*spb64-spa15*spb65;
const complex<T> spab_1_3456_6=-(spa13*spb63)-spa14*spb64-spa15*spb65;
const complex<T> spab_1_3456_3=spa14*spb43+spa15*spb53+spa16*spb63;
const complex<T> spab_1_3456_2=spa13*spb32+spa14*spb42+spa15*spb52+spa16*spb62;
const complex<T> spab_1_345_5=-(spa13*spb53)-spa14*spb54;
const complex<T> spab_1_345_3=spa14*spb43+spa15*spb53;
const complex<T> spab_1_345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_34_4=-(spa13*spb43);
const complex<T> spab_1_34_3=spa14*spb43;
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_3_2=spa13*spb32;
const complex<T> spab_1_25_4=-(spa12*spb42)+spa15*spb54;
const complex<T> spab_1_2_4=-(spa12*spb42);
const complex<T> spab_1_24_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spab_1_2345_6=-(spa12*spb62)-spa13*spb63-spa14*spb64-spa15*spb65;
const complex<T> spab_1_2345_5=-(spa12*spb52)-spa13*spb53-spa14*spb54;
const complex<T> spab_1_2345_3=-(spa12*spb32)+spa14*spb43+spa15*spb53;
const complex<T> spab_1_2345_2=spa13*spb32+spa14*spb42+spa15*spb52;
const complex<T> spab_1_234_4=-(spa12*spb42)-spa13*spb43;
const complex<T> spab_1_234_3=-(spa12*spb32)+spa14*spb43;
const complex<T> spab_1_234_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_23_3=-(spa12*spb32);
const complex<T> spab_1_23_2=spa13*spb32;
const complex<T> spab_1_16_5=spa16*spb65;
const complex<T> spab_1_156_4=spa15*spb54+spa16*spb64;
const complex<T> spab_1_12_6=-(spa12*spb62);
const complex<T> spab_1_126_5=-(spa12*spb52)+spa16*spb65;
const complex<T> spab_1_1256_4=-(spa12*spb42)+spa15*spb54+spa16*spb64;
const complex<T> spaa_5_46_13_6=-(spa45*(spa16*spb41+spa36*spb43))+spa56*(spa16*spb61+spa36*spb63);
const complex<T> spaa_5_13_6_2=-(spa26*(spa15*spb61+spa35*spb63));
const complex<T> spaa_1_6_5_4=-(spa16*spa45*spb65);
const complex<T> spaa_1_6_25_4=spa16*(spa24*spb62-spa45*spb65);
const complex<T> spaa_1_56_24_5=spa15*(spa25*spb52+spa45*spb54)+spa16*(spa25*spb62+spa45*spb64);
const complex<T> spaa_1_5_2_4=spa15*spa24*spb52;
const complex<T> spaa_1_36_13_2=spa12*spa13*spb31+spa16*(spa12*spb61-spa23*spb63);
const complex<T> spaa_1_3_24_5=spa13*(spa25*spb32-spa45*spb43);
const complex<T> spaa_1_26_4_5=spa45*(-(spa12*spb42)+spa16*spb64);
const complex<T> spaa_1_26_45_6=spa12*(-(spa46*spb42)-spa56*spb52)+spa16*(spa46*spb64+spa56*spb65);
const complex<T> spaa_1_26_45_2=spa12*(spa24*spb42+spa25*spb52)+spa16*(-(spa24*spb64)-spa25*spb65);
const complex<T> spaa_1_2_56_4=spa12*(spa45*spb52+spa46*spb62);
const complex<T> spaa_1_2_5_4=spa12*spa45*spb52;
const complex<T> s56=(spa56*spb65);
const complex<T> s45=(spa45*spb54);
const complex<T> s35=(spa35*spb53);
const complex<T> s34=(spa34*spb43);
const complex<T> s26=(spa26*spb62);
const complex<T> s25=(spa25*spb52);
const complex<T> s24=(spa24*spb42);
const complex<T> s23=(spa23*spb32);
const complex<T> s16=(spa16*spb61);
const complex<T> s15=(spa15*spb51);
const complex<T> s13=(spa13*spb31);
const complex<T> s12=(spa12*spb21);

return(
complex<T>(0,1)*((INT1*(s12+s13+s23)*spa13*spa16*spab_6_13_2*spb26)/(INT6*spa12*spa26*spa45*spa56*spab_4_13_2*spab_6_12_3*spb12)-(INT1*(s24+s25+s45)*spa12*spa13*spab_2_13_6*spb26)/(INT6*spa16*spa25*spa26*spa45*spab_2_45_3*spab_4_13_6*spb16)+(INT2*(s24+s25+s45)*spa13*spa15*spab_5_13_6*spb25*spb36)/(INT3*spa16*spa45*spa56*spab_2_45_3*spab_4_13_6*spab_5_24_3*spb16)+(INT1*(s16+s34+s35+s45)*spa12*spa23*spab_1_345_2*spab_3_345_2*(-(spa23*spab_1_345_2)+spa12*spba_2_345_3))/(INT6*s16*(s34+s35+s45)*pow(-s16+s34+s35+s45,2)*spa16*spa25*spa26*spa34*spa45)+(INT1*(s23+s24+INT2*s34)*spa12*spa23*spab_1_34_2*spab_3_34_2*(-(spa23*spab_1_34_2)+spa12*spba_2_34_3))/(INT6*pow(-s23-s24,2)*s34*(s23+s24+s34)*spa16*spa24*spa25*spa34*spa56)+(INT1*(s15+s16+s23+s56)*spa14*spa34*spab_1_156_4*spab_3_156_4*(spa34*spab_1_156_4+spa14*spba_4_156_3))/(INT6*s23*(s15+s16+s56)*pow(s15+s16-s23+s56,2)*spa16*spa23*spa24*spa45*spa56)+(INT1*(INT2*s12+s16+s26)*spa16*spa36*spab_1_12_6*spab_3_12_6*(spa36*spab_1_12_6+spa16*spba_6_12_3))/(INT6*s12*pow(-s16-s26,2)*(s12+s16+s26)*spa12*spa26*spa34*spa45*spa56)+(INT1*spb26*spb36*pow(s24+s25+s45,2)*(spa56*spaa_1_36_13_2+spa16*spaa_5_13_6_2-(spa56*spab_5_13_6*spb13*pow(spa12,2))/(spa25*spb36)))/(INT3*spa16*spa26*spa45*spa56*spab_2_45_3*spab_4_13_6*spab_5_24_3*spb13*spb16)+(INT1*spa14*spb24*pow(spa13,2))/(INT6*spa16*spa23*spa24*spa45*spa56*spb23)+(INT1*spa12*spb24*pow(spa13,2))/(INT6*spa16*spa24*spa25*spa34*spa56*spb34)+(INT1*spab_5_13_6*spb25*spb36*pow(spa13,2))/(INT3*spa16*spa45*spa56*spab_2_45_3*spab_4_13_6*spb16)-(INT1*spa15*spb24*(spaa_1_3_24_5-INT4*spa15*spa24*spb24)*spb56*pow(spa13,2))/(INT6*(s23+s24+s34)*spa24*spa25*spa45*spa56*spab_1_23_4*spab_1_34_2)-(INT2*pow(spa13,3))/(INT9*spa16*spa23*spa24*spa45*spa56)-(INT2*spa14*pow(spa13,3))/(INT9*spa12*spa16*spa24*spa34*spa45*spa56)-(INT1*spa14*spb25*spb26*pow(spa13,3))/(INT6*spa16*spa34*spa45*spa56*spab_1_34_5*spab_4_13_6)+(INT1*(spa45*spab_1_34_5-spa16*spab_4_13_6)*((spa15*spa26*spb25)/spa16+(spa25*spa46*spb26)/spa45)*pow(spa13,3))/(INT3*spa16*spa25*spa26*spa34*spa45*spa56*spab_1_34_5*spab_4_13_6)+(INT1*spa15*spa25*spb45*spb56*pow(spa13,3))/(INT6*spa16*spa23*spa24*spa45*spa56*spab_1_23_4*spab_2_13_6)+(INT1*spab_1_23_4*spab_1_34_2*spb24*pow(spa14,2))/(INT2*(s23+s24+s34)*spa16*spa24*spa45*spa56*spab_1_24_3*spb23)+(INT1*(s12+s16+s26+s34)*spa35*spab_1_126_5*spab_3_126_5*(spa35*spab_1_126_5+spa15*spba_5_126_3)*pow(spa15,2))/(INT6*(s12+s16+s26)*pow(s12+s16+s26-s34,2)*s34*spa12*spa16*spa25*spa34*spa45*spa56)-(INT1*spb45*pow(spa13,2)*pow(spa15,2))/(INT6*spa12*spa16*spa25*spa34*spa45*spa56*spb34)+(spb25*pow(spa13,3)*(-((INT1*(spa16*spa24*spa25-spa12*spa26*spa45)*spb26)/(INT2*spa26))+(spa12*spa24*spa45*spb25*pow(spa15,2))/(spa14*spa16*spa25)))/(spa16*spa25*spa34*spa45*spa56*spab_1_34_5*spab_4_13_6)-(INT1*spa15*(spa24*spab_1_34_2+spa16*spab_4_13_6)*spb25*pow(spa13,3))/(INT6*spa24*spa25*spa34*spa56*spab_1_34_2*spab_4_13_6*pow(spa16,2))-(INT1*spa13*spa14*(-(spa12*spab_1_234_2)-spa14*spab_1_234_4))/(INT2*spa16*spa45*spa56*spab_1_234_3*pow(spa24,2))+(INT1*spa13*spa14*(-(spa12*spab_1_23_2)+spa14*spab_1_23_4))/(INT2*spa16*spa45*spa56*spab_1_23_3*pow(spa24,2))-(INT1*spa13*spa14*(-(spa12*spab_1_34_2)+spa14*spab_1_34_4))/(INT2*spa16*spa45*spa56*spab_1_34_3*pow(spa24,2))+(INT1*(s15+s16+s23+s56)*spa34*spab_1_156_4*spab_3_156_4*pow(spa14,2))/(INT2*s23*(-s15-s16+s23-s56)*(s15+s16+s56)*spa16*spa45*spa56*pow(spa24,2))-(INT2*pow(spa13,4)*(-((INT1*spa12*spa14*spa23*spa34)/(INT2*pow(spa13,2)*pow(spa24,2)))+((INT1*spa12*spa23*spa34*spab_1_3_2)/(INT6*spa24*pow(spa13,2))+(spa14*pow(spa34,2)*spab_1_3_2*pow(spa12,2)*pow(spa23,2))/(pow(spa13,4)*pow(spa24,3)))/s23))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa13,4)*(-((INT1*spa12*spa14*spa23*spa34)/(INT2*pow(spa13,2)*pow(spa24,2)))+((INT1*spa12*spa23*spa34*spab_1_3_2)/(INT6*spa24*pow(spa13,2))+(spa14*pow(spa34,2)*spab_1_3_2*pow(spa12,2)*pow(spa23,2))/(pow(spa13,4)*pow(spa24,3)))/s23))/(spa12*spa16*spa23*spa34*spa45*spa56)+(pow(spa13,3)*((spa46*spb26)/(spa26*spa45)-(spa15*spa24*(-(spa15*spb25)+spa16*spb26))/(spa14*spa16*pow(spa25,2))))/(spa34*spa45*spa56*spab_1_34_5)-(INT1*(s23+s24+INT2*s34)*spab_1_34_2*spab_3_34_2*pow(spa13,2)*(-((INT2*spa12*spa14*spa23*spa34)/(pow(spa13,2)*pow(spa24,2)))+(INT2*spa12*spa15*spa23*spa35)/(pow(spa13,2)*pow(spa25,2))))/(INT4*(s23+s24)*s34*(s23+s24+s34)*spa16*spa34*spa45*spa56)+(INT1*spa12*spab_1_23_4*(spa15*spa24*spab_1_23_4-spa12*spa25*spab_1_34_2)*spb24)/(INT2*(s23+s24+s34)*spa16*spa24*spa56*spab_1_24_3*spb34*pow(spa25,2))+(INT1*spa15*spab_1_25_4*pow(spa13,2))/(INT2*spa14*spa16*spa34*spa56*spb34*pow(spa25,2))-(INT1*spa15*spa45*(spaa_1_5_2_4+INT2*spaa_1_6_25_4)*spb56*pow(spa13,3))/(INT2*spa14*spa16*spa24*spa34*spa56*spab_1_34_2*spab_4_13_6*pow(spa25,2))-(spa15*(-((spa15*spa24*spb25)/(spa16*spab_4_13_6))+(spa16*spa45*spb56)/(spa24*spab_1_34_2)+(spa15*spa45*spb25*spb56)/(spab_1_34_2*spab_4_13_6))*pow(spa13,3))/(spa14*spa16*spa34*spa56*pow(spa25,2))-(INT1*spa13*(-(spa12*spab_1_2345_2)-spa15*spab_1_2345_5)*pow(spa15,2))/(INT2*spa14*spa16*spa45*spa56*spab_1_2345_3*pow(spa25,2))+(INT1*spa13*(-(spa12*spab_1_234_2)+spa15*spab_1_234_5)*pow(spa15,2))/(INT2*spa14*spa16*spa45*spa56*spab_1_234_3*pow(spa25,2))-(INT1*spa13*(-(spa12*spab_1_345_2)+spa15*spab_1_345_5)*pow(spa15,2))/(INT2*spa14*spa16*spa45*spa56*spab_1_345_3*pow(spa25,2))+(INT1*spa13*(-(spa12*spab_1_34_2)-spa15*spab_1_34_5)*pow(spa15,2))/(INT2*spa14*spa16*spa45*spa56*spab_1_34_3*pow(spa25,2))-(INT1*spa12*spaa_1_26_45_2*spb26*spb45*(-((spa15*spaa_1_26_4_5)/pow(spa25,2))+(spa16*spaa_1_26_45_6)/pow(spa26,2)))/(INT2*(s34+s35+s45)*spa45*spa56*spab_1_34_5*spab_1_45_3*spab_2_45_3)-(INT1*(s24+s25+s45)*spa12*spa15*spab_2_13_6*spb26*(spab_5_13_6/(spa16*pow(spa25,2))-spab_6_13_6/(spa15*pow(spa26,2))))/(INT2*spa45*spa56*spab_2_45_3*spab_4_13_6*spb16*spb36)-(INT1*(s16+s34+s35+s45)*spab_1_345_2*spab_3_345_2*pow(spa13,2)*(-((INT2*spa12*spa15*spa23*spa35)/(pow(spa13,2)*pow(spa25,2)))+(INT2*spa12*spa16*spa23*spa36)/(pow(spa13,2)*pow(spa26,2))))/(INT4*s16*(s16-s34-s35-s45)*(s34+s35+s45)*spa16*spa34*spa45*spa56)+(INT1*spa13*spa16*(-(spa12*spab_1_2345_2)+spa16*spab_1_2345_6))/(INT2*spa14*spa45*spa56*spab_1_2345_3*pow(spa26,2))-(INT1*spa13*spa16*(-(spa12*spab_1_3456_2)+spa16*spab_1_3456_6))/(INT2*spa14*spa45*spa56*spab_1_3456_3*pow(spa26,2))+(INT1*spa13*spa16*(-(spa12*spab_1_345_2)-spa16*spab_1_345_6))/(INT2*spa14*spa45*spa56*spab_1_345_3*pow(spa26,2))-(INT2*pow(spa13,4)*(-((INT1*spa12*spa16*spa23*spa36)/(INT2*pow(spa13,2)*pow(spa26,2)))+(-((INT1*spa16*spa23*(-(spa23*spab_1_3456_2)+spa12*spba_2_3456_3))/(INT6*spa26*pow(spa13,2)))-(spa16*spa36*(spa36*spab_1_3456_2+spa16*spab_3_3456_2)*pow(spa12,2)*pow(spa23,2))/(pow(spa13,4)*pow(spa26,3)))/s12))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa13,4)*(-((INT1*spa12*spa16*spa23*spa36)/(INT2*pow(spa13,2)*pow(spa26,2)))+(-((INT1*spa16*spa23*(-(spa23*spab_1_3456_2)+spa12*spba_2_3456_3))/(INT6*spa26*pow(spa13,2)))-(spa16*spa36*(spa36*spab_1_3456_2+spa16*spab_3_3456_2)*pow(spa12,2)*pow(spa23,2))/(pow(spa13,4)*pow(spa26,3)))/s12))/(spa12*spa16*spa23*spa34*spa45*spa56)-(INT2*pow(spa13,4)*(-((INT1*spa12*spa14*spa23*spa34)/(INT2*pow(spa13,2)*pow(spa24,2)))+((INT1*spa14*spa23*(spa34*spab_1_1256_4+spa14*spba_4_1256_3))/(INT6*spa24*pow(spa13,2))-(spa12*spa23*(-(spa23*spab_1_1256_4)+spa12*spab_3_1256_4)*pow(spa14,2)*pow(spa34,2))/(pow(spa13,4)*pow(spa24,3)))/s34))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa13,4)*(-((INT1*spa12*spa14*spa23*spa34)/(INT2*pow(spa13,2)*pow(spa24,2)))+((INT1*spa14*spa23*(spa34*spab_1_1256_4+spa14*spba_4_1256_3))/(INT6*spa24*pow(spa13,2))-(spa12*spa23*(-(spa23*spab_1_1256_4)+spa12*spab_3_1256_4)*pow(spa14,2)*pow(spa34,2))/(pow(spa13,4)*pow(spa24,3)))/s34))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT1*(s15+INT2*s16+s56)*spa15*spab_1_16_5*spab_3_16_5*(spa35*spab_1_16_5+spa15*spba_5_16_3)*pow(spa35,2))/(INT6*s16*pow(-s15-s56,2)*(s15+s16+s56)*spa16*spa23*spa25*spa34*spa45*spa56)-(INT1*(s12+s16+s26+s34)*spab_1_126_5*spab_3_126_5*pow(spa15,2)*pow(spa35,2))/(INT2*(s12+s16+s26)*s34*(-s12-s16-s26+s34)*spa16*spa34*spa45*spa56*pow(spa25,2))+(INT1*(s15+INT2*s16+s56)*spab_1_16_5*spab_3_16_5*pow(spa15,2)*pow(spa35,2))/(INT2*s16*(s15+s56)*(s15+s16+s56)*spa16*spa34*spa45*spa56*pow(spa25,2))-(INT1*(INT2*s12+s16+s26)*spa16*spab_1_12_6*spab_3_12_6*pow(spa36,2))/(INT2*s12*(s16+s26)*(s12+s16+s26)*spa34*spa45*spa56*pow(spa26,2))-(INT2*pow(spa13,4)*(-((INT1*spa12*spa16*spa23*spa36)/(INT2*pow(spa13,2)*pow(spa26,2)))+((INT1*spa12*spa16*spa36*spba_6_1_3)/(INT6*spa26*pow(spa13,2))+(pow(spa12,2)*spa23*spab_3_1_6*pow(spa16,2)*pow(spa36,2))/(pow(spa13,4)*pow(spa26,3)))/s16))/(spa12*spa16*spa23*spa34*spa45*spa56)+(INT2*pow(spa13,4)*(-((INT1*spa12*spa16*spa23*spa36)/(INT2*pow(spa13,2)*pow(spa26,2)))+((INT1*spa12*spa16*spa36*spba_6_1_3)/(INT6*spa26*pow(spa13,2))+(pow(spa12,2)*spa23*spab_3_1_6*pow(spa16,2)*pow(spa36,2))/(pow(spa13,4)*pow(spa26,3)))/s16))/(spa12*spa16*spa23*spa34*spa45*spa56)-(INT1*spa46*(spa45*spab_1_34_5-spa12*spab_4_13_2)*spb26*pow(spa13,3))/(INT6*spa12*spa26*spa34*spa56*spab_1_34_5*spab_4_13_2*pow(spa45,2))-(INT1*spa24*(spaa_1_2_5_4+INT2*spaa_1_6_25_4)*spb26*pow(spa13,3)*((spa15*spa45)/pow(spa25,2)-(spa16*spa46)/pow(spa26,2)))/(INT2*spa14*spa16*spa34*spa56*spab_1_34_5*spab_4_13_6*pow(spa45,2))+(INT1*(INT2*spaa_1_2_56_4+spaa_1_6_5_4)*spb26*pow(spa13,3)*pow(spa46,2))/(INT2*spa14*spa34*spa56*spab_1_34_5*spab_4_13_2*pow(spa26,2)*pow(spa45,2))+(INT1*spa12*spa16*spb26*spb45*pow(spaa_1_26_45_6,2))/(INT2*(s34+s35+s45)*spa45*spa56*spab_1_34_5*spab_1_45_3*spab_6_45_3*pow(spa26,2))-(INT1*spb24*spb56*pow(spa15,2)*pow(spaa_1_56_24_5,2))/(INT2*(s23+s24+s34)*spa45*spa56*spab_1_24_3*spab_1_34_2*spab_5_24_3*pow(spa25,2))+(INT1*pow(spa15,2)*((spa12*pow(spab_1_2345_2,2))/spb23-(spa15*pow(spab_1_2345_5,2))/spb35))/(INT2*spa14*spa16*spa45*spa56*spab_1_2345_3*pow(spa25,2))-(INT1*spa16*((spa12*pow(spab_1_2345_2,2))/spb23+(spa16*pow(spab_1_2345_6,2))/spb36))/(INT2*spa14*spa45*spa56*spab_1_2345_3*pow(spa26,2))+(INT1*spa14*((spa12*pow(spab_1_234_2,2))/spb23-(spa14*pow(spab_1_234_4,2))/spb34))/(INT2*spa16*spa45*spa56*spab_1_234_3*pow(spa24,2))-(INT1*pow(spa15,2)*((spa12*pow(spab_1_234_2,2))/spb23+(spa15*pow(spab_1_234_5,2))/spb35))/(INT2*spa14*spa16*spa45*spa56*spab_1_234_3*pow(spa25,2))-(INT1*spa14*((spa12*pow(spab_1_23_2,2))/spb23+(spa14*pow(spab_1_23_4,2))/spb34))/(INT2*spa16*spa45*spa56*spab_1_23_3*pow(spa24,2))+(INT1*spa16*((spa12*pow(spab_1_3456_2,2))/spb23+(spa16*pow(spab_1_3456_6,2))/spb36))/(INT2*spa14*spa45*spa56*spab_1_3456_3*pow(spa26,2))+(INT1*pow(spa15,2)*((spa12*pow(spab_1_345_2,2))/spb23+(spa15*pow(spab_1_345_5,2))/spb35))/(INT2*spa14*spa16*spa45*spa56*spab_1_345_3*pow(spa25,2))-(INT1*spa16*((spa12*pow(spab_1_345_2,2))/spb23-(spa16*pow(spab_1_345_6,2))/spb36))/(INT2*spa14*spa45*spa56*spab_1_345_3*pow(spa26,2))+(INT1*spa14*((spa12*pow(spab_1_34_2,2))/spb23+(spa14*pow(spab_1_34_4,2))/spb34))/(INT2*spa16*spa45*spa56*spab_1_34_3*pow(spa24,2))-(INT1*pow(spa15,2)*((spa12*pow(spab_1_34_2,2))/spb23-(spa15*pow(spab_1_34_5,2))/spb35))/(INT2*spa14*spa16*spa45*spa56*spab_1_34_3*pow(spa25,2))-(INT1*spb45*pow(spa15,2)*pow(spab_1_35_4,2))/(INT2*(s34+s35+s45)*spa16*spa56*spab_1_45_3*spb34*pow(spa25,2))+(INT1*spa12*spaa_1_26_45_2*spb26*spb45*(spaa_1_26_45_2+spa12*spa45*spb45))/(INT6*(s34+s35+s45)*spa25*spa26*spa45*spab_1_34_5*pow(spab_2_45_3,2))-(INT1*spa46*spb45*pow(spab_5_13_2,3))/(INT3*(s12+s13+s23)*spab_6_12_3*spb12*spb13*pow(spa45,2)*pow(spa56,2))+(INT1*(s24+s25+s45)*spa13*spa15*spb56*pow(spab_5_13_6,2))/(INT6*spa16*spa25*spa45*spa56*spab_2_45_6*spab_4_13_6*spab_5_24_3*spb16)-(INT1*(s24+s25+s45)*spb56*pow(spa15,2)*pow(spab_5_13_6,2))/(INT2*spa16*spa45*spa56*spab_4_13_6*spab_5_24_3*spb16*spb36*pow(spa25,2))-(INT1*(s24+s25+s45)*spa12*spb56*pow(spab_5_13_6,3))/(INT3*spa25*spa45*spa56*spab_2_45_3*spab_2_45_6*spab_4_13_6*spab_5_24_3*spb16)-(INT1*spa46*spb36*spb45*pow(spab_5_13_6,3))/(INT3*spa56*spab_2_45_3*spab_2_45_6*spab_5_24_3*spb13*spb16*pow(spa45,2))+(INT1*spa15*spab_1_34_2*spab_5_13_6*spb36*pow(s24+s25+s45,2))/(INT3*spa16*spa45*spa56*spab_2_45_3*spab_4_13_6*spb16*pow(spab_5_24_3,2))-(INT1*(s12+s13+s23)*spa16*spb26*pow(spab_6_13_2,2))/(INT2*spa45*spa56*spab_4_13_2*spab_6_12_3*spb12*spb23*pow(spa26,2))+(INT1*spaa_5_46_13_6*spb26*pow(spab_6_13_2,2))/(INT3*spa26*spa45*spab_4_13_2*spab_6_12_3*spb12*spb13*pow(spa56,2))-(INT1*spa16*spaa_1_26_45_6*spb26*spb45*(spaa_1_26_45_6+spa16*spa45*spb45))/(INT6*(s34+s35+s45)*spa26*spa45*spa56*spab_1_34_5*pow(spab_6_45_3,2))-(INT1*spa14*spab_1_34_2*(spab_1_3_2+INT2*spab_1_4_2)*spb24)/(INT6*(s23+s24+s34)*spa16*spa24*spa45*spa56*pow(spb23,2))+(spa12*pow(spa14,2)*pow(spb24,2))/(spa16*spa45*spa56*spb23*spb34*pow(spa24,2))-(INT5*spa13*spa24*spb56*pow(spa15,3)*pow(spb24,3))/(INT6*(s23+s24+s34)*spa25*spa45*spa56*spab_1_23_4*spab_1_34_2*spab_5_24_3)-(INT1*pow(spb24,3))/(INT3*spb12*spb13*spb34*pow(spa56,2))-(INT1*(spa12*spa45*spab_6_2_3+spa14*spa25*spab_6_4_3)*spb56*pow(spa15,3)*pow(spb24,4))/(INT3*(s23+s24+s34)*spa25*spa45*spab_1_23_4*spab_1_34_2*spab_5_24_3*spb23*spb34*pow(spa56,2))-(INT1*pow(spa13,3)*pow(spa15,2)*pow(spb25,2))/(INT2*spa25*spa34*spa56*spab_1_34_2*spab_4_13_6*pow(spa16,2))+(spa12*pow(spa15,3)*pow(spb25,2))/(spa14*spa16*spa45*spa56*spb23*spb35*pow(spa25,2))-(INT1*spa12*pow(spa13,3)*pow(spa15,2)*pow(spb25,3))/(INT3*spa25*spa34*spa56*spab_1_34_2*spab_1_34_5*spab_4_13_6*pow(spa16,2))+(spa12*pow(spa16,2)*pow(spb26,2))/(spa14*spa45*spa56*spb23*spb36*pow(spa26,2))+(INT1*pow(spa13,3)*pow(spa46,2)*pow(spb26,2))/(INT2*spa26*spa34*spa56*spab_1_34_5*spab_4_13_2*pow(spa45,2))+(INT1*spa14*pow(spa13,3)*pow(spb25,2)*pow(spb26,2))/(INT3*spa34*spa56*spab_1_34_2*spab_1_34_5*spab_4_13_2*spab_4_13_6)-(INT1*spa24*pow(spa13,3)*pow(spa46,2)*pow(spb26,3))/(INT3*spa26*spa34*spa56*spab_1_34_5*spab_4_13_2*spab_4_13_6*pow(spa45,2))-(INT1*spa12*spab_1_23_4*(INT2*spab_1_2_4+spab_1_3_4)*spb24)/(INT6*(s23+s24+s34)*spa16*spa24*spa25*spa56*pow(spb34,2))+(INT1*spab_1_35_4*(spab_1_3_4+INT2*spab_1_5_4)*spb45*pow(spa15,2))/(INT6*(s34+s35+s45)*spa12*spa16*spa25*spa45*spa56*pow(spb34,2))+(INT1*spab_1_45_3*spb35*pow(spab_1_35_4,3))/(INT3*(s34+s35+s45)*spa12*spa16*spa45*spab_2_45_3*spab_6_45_3*pow(spb34,2))-(INT1*pow(spab_1_24_3,3)*pow(spb24,3))/(INT3*(s23+s24+s34)*spa16*spa24*spa56*spab_5_24_3*pow(spb23,2)*pow(spb34,2))-(INT1*pow(spb24,3)*pow(spb36,3))/(INT3*spb13*spb16*spb23*spb34*pow(spab_5_24_3,2))-(INT1*s26*pow(spab_1_45_3,3)*pow(spb45,3))/(INT3*(s34+s35+s45)*spab_1_34_5*spb34*pow(spab_2_45_3,2)*pow(spab_6_45_3,2))+(INT1*(s45+s56)*pow(spa13,3)*pow(SPB(i4,i6),2))/(INT3*(s12+s13+s23)*spa23*spa45*spa56*spab_1_23_4*spab_2_13_6))
);
}


}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gpmmpmm(const eval_param<T>& ep,const mass_param_coll& mpc){

return (-std::conj(R6gmppmpp<i1,i2,i3,i4,i5,i6>(ep,mpc)));

}

template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gpmpmmm(const eval_param<T>& ep,const mass_param_coll& mpc){

return (-std::conj(R6gmpmppp<i1,i2,i3,i4,i5,i6>(ep,mpc)));

}



template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6gmmmppp(const eval_param<T>& ep,const mass_param_coll& mpc){


const complex<T> INT1(1,0);
const complex<T> INT2(2,0);
const complex<T> INT3(3,0);
const complex<T> INT4(4,0);
const complex<T> INT5(5,0);
const complex<T> INT6(6,0);
const complex<T> INT8(8,0);
const complex<T> INT9(9,0);

{
const complex<T> spb65=SPB(i6,i5);
const complex<T> spb64=SPB(i6,i4);
const complex<T> spb63=SPB(i6,i3);
const complex<T> spb62=SPB(i6,i2);
const complex<T> spb61=SPB(i6,i1);
const complex<T> spb56=SPB(i5,i6);
const complex<T> spb54=SPB(i5,i4);
const complex<T> spb53=SPB(i5,i3);
const complex<T> spb51=SPB(i5,i1);
const complex<T> spb46=SPB(i4,i6);
const complex<T> spb45=SPB(i4,i5);
const complex<T> spb43=SPB(i4,i3);
const complex<T> spb42=SPB(i4,i2);
const complex<T> spb41=SPB(i4,i1);
const complex<T> spb34=SPB(i3,i4);
const complex<T> spb32=SPB(i3,i2);
const complex<T> spb26=SPB(i2,i6);
const complex<T> spb25=SPB(i2,i5);
const complex<T> spb24=SPB(i2,i4);
const complex<T> spb23=SPB(i2,i3);
const complex<T> spb21=SPB(i2,i1);
const complex<T> spb16=SPB(i1,i6);
const complex<T> spb12=SPB(i1,i2);
const complex<T> spa56=SPA(i5,i6);
const complex<T> spa45=SPA(i4,i5);
const complex<T> spa36=SPA(i3,i6);
const complex<T> spa35=SPA(i3,i5);
const complex<T> spa34=SPA(i3,i4);
const complex<T> spa26=SPA(i2,i6);
const complex<T> spa25=SPA(i2,i5);
const complex<T> spa24=SPA(i2,i4);
const complex<T> spa23=SPA(i2,i3);
const complex<T> spa16=SPA(i1,i6);
const complex<T> spa15=SPA(i1,i5);
const complex<T> spa14=SPA(i1,i4);
const complex<T> spa13=SPA(i1,i3);
const complex<T> spa12=SPA(i1,i2);
const complex<T> pow3_spb46=(pow(spb46,3));
const complex<T> pow3_spa13=(pow(spa13,3));
const complex<T> pow2_spb46=(pow(spb46,2));
const complex<T> pow2_spa13=(pow(spa13,2));
const complex<T> spab_5_34_6=spa35*spb63+spa45*spb64;
const complex<T> spab_5_34_2=-(spa35*spb32)-spa45*spb42;
const complex<T> spab_5_16_4=spa15*spb41+spa56*spb64;
const complex<T> spab_5_16_2=spa15*spb21+spa56*spb62;
const complex<T> spab_3_5_6=-(spa35*spb65);
const complex<T> spab_3_4_6=-(spa34*spb64);
const complex<T> spab_3_2_6=spa23*spb62;
const complex<T> spab_3_1_6=spa13*spb61;
const complex<T> spab_3_16_2=spa13*spb21+spa36*spb62;
const complex<T> spab_3_12_6=spa13*spb61+spa23*spb62;
const complex<T> spab_1_6_4=spa16*spb64;
const complex<T> spab_1_5_4=spa15*spb54;
const complex<T> spab_1_3_4=-(spa13*spb43);
const complex<T> spab_1_34_2=spa13*spb32+spa14*spb42;
const complex<T> spab_1_2_4=-(spa12*spb42);
const complex<T> spab_1_23_4=-(spa12*spb42)-spa13*spb43;
const complex<T> s56=(spa56*spb65);
const complex<T> s45=(spa45*spb54);
const complex<T> s35=(spa35*spb53);
const complex<T> s34=(spa34*spb43);
const complex<T> s26=(spa26*spb62);
const complex<T> s24=(spa24*spb42);
const complex<T> s23=(spa23*spb32);
const complex<T> s16=(spa16*spb61);
const complex<T> s15=(spa15*spb51);
const complex<T> s12=(spa12*spb21);
const complex<T> pow3_spab_3_12_6=(pow(spab_3_12_6,3));
const complex<T> pow3_spab_1_23_4=(pow(spab_1_23_4,3));

return(
complex<T>(0,1)*(-((INT4*pow3_spab_3_12_6)/(INT9*(s34+s35+s45)*spa34*spa45*spab_5_34_2*spb12*spb16))-
((s12+INT2*s16+s26)*spa12*spa23*spab_3_16_2*spb26*(spa12*spab_3_16_2*spb16-(s12+s16+s26)*spa23*spb26))/(INT6*s16*pow(-s12-s26,2)*(s12+s16+s26)*spa34*spa45*spab_5_16_2*spb12)+
(INT2*(pow3_spab_3_12_6/((s34+s35+s45)*spa34*spa45*spab_5_34_2*spb12*spb16)-
pow3_spab_1_23_4/((s23+s24+s34)*spa16*spa56*spab_5_34_2*spb23*spb34)))/INT3+(INT4*pow3_spab_1_23_4)/(INT9*(s23+s24+s34)*spa16*spa56*spab_5_34_2*spb23*spb34)-
((s23+s24+INT2*s34)*spa12*spa23*spab_1_34_2*spb24*((s23+s24+s34)*spa12*spb24-spa23*spab_1_34_2*spb34))/(INT6*pow(-s23-s24,2)*s34*(s23+s24+s34)*spa16*spa56*spab_5_34_2*spb23)-
((INT2*s34+s35+s45)*spa35*spab_5_34_6*spb45*spb56*(-(spa34*spab_5_34_6*spb45)+
(s34+s35+s45)*spa35*spb56))/(INT6*s34*pow(-s35-s45,2)*(s34+s35+s45)*spa45*spab_5_34_2*spb12*spb16)-
((s15+INT2*s16+s56)*spa15*spab_5_16_4*spb45*spb56*(-((s15+s16+s56)*spa15*spb45)+
spa16*spab_5_16_4*spb56))/(INT6*s16*pow(-s15-s56,2)*(s15+s16+s56)*spa56*spab_5_16_2*spb23*spb34)+
(((pow2_spa13*(INT3*spab_1_2_4+spab_1_3_4))/(spa16*spa34)+
(pow3_spa13*spa25*spb23)/(spa16*spa34*spa45)-
(pow2_spb46*(INT3*spab_1_5_4+spab_1_6_4))/(spb16*spb34)+
(pow3_spb46*spa56*spb25)/(spb12*spb16*spb34)-
((spa13/spa34+(spab_1_2_4-spab_1_5_4)/(s23+s24+s34)+
spb46/spb16)*pow(spab_1_23_4,2))/(spa16*spb34)))/(INT6*spa56*spab_5_34_2*spb23)+
(((pow2_spa13*(spab_3_1_6+INT3*spab_3_2_6))/(spa16*spa34)-
(pow3_spa13*spa25*spb12)/(spa16*spa34*spa56)-
(pow2_spb46*(spab_3_4_6+INT3*spab_3_5_6))/(spb16*spb34)-
(pow3_spb46*spa45*spb25)/(spb16*spb23*spb34)-
((-(spa13/spa16)+
(spab_3_2_6-spab_3_5_6)/(s12+s16+s26)-
spb46/spb34)*pow(spab_3_12_6,2))/(spa34*spb16)))/(INT6*spa45*spab_5_16_2*spb12))
);
}

}


template <int i1, int i2, int i3, int i4, int i5, int i6, class T> complex<T> R6g2p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
  return  complex<T>(0,1)*((complex<T>(-X,0)/complex<T>(9,0)*pow(ep.spb(i1,i2),3))/(ep.spb(i1,i6)*ep.spb(i2,i3)*
ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i6))
#if _ONLY_X_PART
+
(complex<T>(1,0)/complex<T>(3,0)*pow(-(ep.spb(i1,i2)*ep.spa(i3,i2))+ep.spb(i1,i4)*ep.spa(i4,i3),3)*
ep.spb(i3,i5))/(ep.s(i2,i3,i4)*ep.spb(i1,i6)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i6)*
ep.spa(i3,i2)*(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*ep.spb(i1,i2)*(-(ep.spb(i1,i2)*ep.spa(i3,i2))+ep.spb(i1,i4)*ep.spa(i4,i3))*
(ep.spb(i1,i4)*ep.spa(i4,i3)+complex<T>(2,0)*(-(ep.spb(i1,i2)*ep.spa(i3,i2)) +
   ep.spb(i1,i4)*ep.spa(i4,i3))))/(ep.s(i2,i3,i4)*ep.spb(i1,i6)*ep.spb(i3,i4)*ep.spb(i4,i5)*
ep.spb(i5,i6)*ep.spa(i3,i2))+(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spa(i6,i3),3))/
 (pow(ep.spb(i4,i5),2)*ep.spa(i2,i1)*ep.spa(i3,i2)*ep.spa(i6,i1)) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spb(i1,i2),3)*pow(ep.spa(i6,i4),2)*(ep.s(i4,i5)+ep.s(i5,i6)))/
 (ep.s(i1,i2,i3)*ep.spb(i2,i3)*ep.spb(i4,i5)*ep.spb(i5,i6)*(-(ep.spb(i1,i2)*ep.spa(i4,i2)) -
 ep.spb(i1,i3)*ep.spa(i4,i3))*(ep.spb(i1,i3)*ep.spa(i6,i1)+ep.spb(i2,i3)*ep.spa(i6,i2))) +
(complex<T>(-1,0)/complex<T>(3,0)*ep.s(i3,i5)*(ep.spb(i1,i4)*ep.spa(i3,i1)+ep.spb(i2,i4)*ep.spa(i3,i2))*
(ep.spb(i1,i4)*ep.spa(i6,i1)+ep.spb(i2,i4)*ep.spa(i6,i2))*(ep.spb(i1,i5)*ep.spa(i6,i1) +
 ep.spb(i2,i5)*ep.spa(i6,i2)))/(pow(ep.spb(i3,i4),2)*pow(ep.spb(i4,i5),2)*ep.spa(i2,i1)*
(ep.spb(i1,i6)*ep.spa(i3,i1)+ep.spb(i2,i6)*ep.spa(i3,i2))*(-(ep.spb(i3,i5)*ep.spa(i3,i2)) -
 ep.spb(i4,i5)*ep.spa(i4,i2))*ep.spa(i6,i1)) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spb(i1,i4)*ep.spa(i6,i1)+ep.spb(i2,i4)*ep.spa(i6,i2),2)*ep.spb(i3,i5)*
ep.spa(i6,i3))/(pow(ep.spb(i3,i4),2)*pow(ep.spb(i4,i5),2)*ep.spa(i2,i1)*
(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*ep.spa(i6,i1)) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spb(i1,i5),2)*pow(ep.spa(i4,i3),2)*
(-(ep.spb(i1,i6)*ep.spb(i4,i5)*ep.spa(i4,i3))-ep.spb(i5,i6)*
  (-(ep.spb(i1,i2)*ep.spa(i3,i2))+ep.spb(i1,i4)*ep.spa(i4,i3)))*ep.spa(i6,i5))/
 (pow(ep.spb(i5,i6),2)*ep.s(i2,i3,i4)*ep.spb(i4,i5)*ep.spa(i3,i2)*
(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*
(-(ep.spb(i1,i2)*ep.spa(i4,i2))-ep.spb(i1,i3)*ep.spa(i4,i3))) +
(complex<T>(1,0)/complex<T>(6,0)*pow(ep.spb(i1,i5)*ep.spa(i6,i1)+ep.spb(i2,i5)*ep.spa(i6,i2),2)*
(ep.s(i1,i2)*ep.spb(i4,i5)+complex<T>(-2,0)*(ep.s(i1,i5)*ep.spb(i4,i5) +
   ep.s(i2,i5)*ep.spb(i4,i5)-ep.spb(i1,i5)*ep.spb(i3,i4)*ep.spa(i3,i1) -
   ep.spb(i2,i5)*ep.spb(i3,i4)*ep.spa(i3,i2)))*ep.spa(i6,i5))/
 (pow(ep.spb(i4,i5),2)*ep.spb(i3,i4)*ep.spb(i5,i6)*ep.spa(i2,i1)*
(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*ep.spa(i6,i1)*
(ep.spb(i1,i3)*ep.spa(i6,i1)+ep.spb(i2,i3)*ep.spa(i6,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*pow(ep.spb(i1,i2),3)*ep.spb(i3,i5)*ep.spa(i6,i4)*ep.spa(i6,i5))/
 (ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i6)*(-(ep.spb(i1,i2)*ep.spa(i4,i2)) -
 ep.spb(i1,i3)*ep.spa(i4,i3))*(ep.spb(i1,i3)*ep.spa(i6,i1)+ep.spb(i2,i3)*ep.spa(i6,i2))) +
(complex<T>(-1,0)/complex<T>(6,0)*ep.spb(i1,i2)*ep.spb(i1,i5)*ep.spa(i4,i3)*(ep.s(i3,i5)*ep.spb(i1,i5) +
 ep.s(i4,i5)*ep.spb(i1,i5)+ep.spb(i1,i6)*ep.spb(i3,i5)*ep.spa(i6,i3) +
 ep.spb(i1,i6)*ep.spb(i4,i5)*ep.spa(i6,i4))*ep.spa(i6,i5))/(ep.s(i2,i3,i4)*ep.spb(i3,i4)*
ep.spb(i4,i5)*ep.spb(i5,i6)*(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*
(-(ep.spb(i1,i2)*ep.spa(i4,i2))-ep.spb(i1,i3)*ep.spa(i4,i3))) +
(complex<T>(1,0)/complex<T>(3,0)*pow(ep.spb(i1,i2),3)*pow(ep.spa(i5,i3),2)*(-ep.s(i3,i4)-ep.s(i4,i5)))/
 (ep.s(i3,i4,i5)*ep.spb(i1,i6)*ep.spb(i3,i4)*ep.spb(i4,i5)*(ep.spb(i1,i6)*ep.spa(i3,i1) +
 ep.spb(i2,i6)*ep.spa(i3,i2))*(ep.spb(i1,i2)*ep.spa(i5,i1)+ep.spb(i2,i6)*ep.spa(i6,i5))) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spb(i1,i2),2)*pow(ep.spa(i5,i3),2)*(-(ep.s(i3,i5)*ep.spb(i1,i5)) +
 ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spa(i3,i2)+ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spa(i4,i2)))/
 (ep.spb(i1,i6)*ep.spb(i3,i4)*ep.spb(i4,i5)*(ep.spb(i1,i6)*ep.spa(i3,i1) +
 ep.spb(i2,i6)*ep.spa(i3,i2))*(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*
(ep.spb(i1,i2)*ep.spa(i5,i1)+ep.spb(i2,i6)*ep.spa(i6,i5))) +
(complex<T>(-1,0)/complex<T>(3,0)*pow(ep.spa(i5,i3),2)*ep.spb(i1,i2)*ep.spb(i2,i4)*ep.spb(i3,i5)*
(ep.spb(i1,i5)*ep.spa(i6,i1)+ep.spb(i2,i5)*ep.spa(i6,i2))*ep.spa(i6,i5))/
 (pow(ep.spb(i3,i4),2)*ep.spb(i4,i5)*(ep.spb(i1,i6)*ep.spa(i3,i1) +
 ep.spb(i2,i6)*ep.spa(i3,i2))*(-(ep.spb(i3,i5)*ep.spa(i3,i2))-ep.spb(i4,i5)*ep.spa(i4,i2))*
ep.spa(i6,i1)*(ep.spb(i1,i2)*ep.spa(i5,i1)+ep.spb(i2,i6)*ep.spa(i6,i5))) +
(complex<T>(1,0)/complex<T>(6,0)*pow(ep.spb(i1,i2),2)*(-((ep.spb(i1,i4)*ep.spa(i4,i3))/ep.spa(i3,i2)) -
 (ep.spb(i2,i5)*ep.spa(i6,i5))/ep.spa(i6,i1)))/(ep.spb(i1,i6)*ep.spb(i2,i3)*ep.spb(i3,i4)*
ep.spb(i4,i5)*ep.spb(i5,i6)) +
(-complex<T>(1,0)/complex<T>(6,0)*(-(((-pow(ep.s(i2,i3,i4),2)+pow(ep.spb(i2,i3),2)*
pow(ep.spa(i3,i2),2))*ep.spb(i1,i4)*ep.spb(i2,i4)*ep.spa(i4,i3)*
(-(ep.spb(i1,i2)*ep.spa(i4,i2))-ep.spb(i1,i3)*ep.spa(i4,i3))*
(-(ep.spb(i1,i4)*ep.spb(i2,i3)*ep.spa(i4,i3))+ep.spb(i2,i4)*
(-(ep.spb(i1,i2)*ep.spa(i4,i2))-ep.spb(i1,i3)*ep.spa(i4,i3))))/
   (pow(-ep.s(i2,i4)-ep.s(i3,i4),3)*ep.s(i2,i3,i4)*ep.spa(i3,i2))) -
 (((pow(ep.spb(i1,i6),2)*pow(ep.spa(i6,i1),2))/ep.s(i2,i3,i4)-ep.s(i2,i3,i4))*
   ep.spb(i1,i5)*ep.spb(i2,i5)*ep.spa(i6,i5)*(ep.spb(i1,i2)*ep.spa(i5,i1) +
ep.spb(i2,i6)*ep.spa(i6,i5))*(ep.spb(i1,i6)*ep.spb(i2,i5)*ep.spa(i6,i5) +
ep.spb(i1,i5)*(ep.spb(i1,i2)*ep.spa(i5,i1)+ep.spb(i2,i6)*ep.spa(i6,i5))))/
  (pow(ep.s(i1,i6)-ep.s(i2,i3,i4),3)*ep.spa(i6,i1))))/(ep.spb(i1,i6)*ep.spb(i2,i3)*
ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i6))
#endif
);
}





//R6g6p and R6g6m already in rat_ampl;

template <class T> complex<T> R6g3(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<0,1,2,3,4,5>(ep,mpc);}  //++----
template <class T> complex<T> R6g6(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<1,2,3,4,5,0>(ep,mpc);}  //-++---
template <class T> complex<T> R6g12(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<2,3,4,5,0,1>(ep,mpc);} //--++--
template <class T> complex<T> R6g24(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<3,4,5,0,1,2>(ep,mpc);} //---++-
template <class T> complex<T> R6g48(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<4,5,0,1,2,3>(ep,mpc);} //----++
template <class T> complex<T> R6g33(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2p<5,0,1,2,3,4>(ep,mpc);} //+----+
//
template <class T> complex<T> R6g30(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<5,0,1,2,3,4>(ep,mpc);} //-++++-
template <class T> complex<T> R6g15(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<4,5,0,1,2,3>(ep,mpc);} //++++--
template <class T> complex<T> R6g39(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<3,4,5,0,1,2>(ep,mpc);} //+++--+
template <class T> complex<T> R6g51(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<2,3,4,5,0,1>(ep,mpc);} //++--++
template <class T> complex<T> R6g57(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<1,2,3,4,5,0>(ep,mpc);} //+--+++
template <class T> complex<T> R6g60(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g2m<0,1,2,3,4,5>(ep,mpc);} //--++++


template <class T> complex<T> R6g58(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<0,1,2,3,4,5>(ep,mpc);} //-+-+++
template <class T> complex<T> R6g53(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<1,2,3,4,5,0>(ep,mpc);} //+-+-++ 
template <class T> complex<T> R6g43(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<2,3,4,5,0,1>(ep,mpc);} //++-+-+
template <class T> complex<T> R6g23(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<3,4,5,0,1,2>(ep,mpc);} //+++-+-
template <class T> complex<T> R6g46(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<4,5,0,1,2,3>(ep,mpc);} //-+++-+
template <class T> complex<T> R6g29(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmppp<5,0,1,2,3,4>(ep,mpc);} //+-+++-

template <class T> complex<T> R6g54(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmppmpp<0,1,2,3,4,5>(ep,mpc);} //-++-++
template <class T> complex<T> R6g45(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmppmpp<1,2,3,4,5,0>(ep,mpc);} //+-++-+
template <class T> complex<T> R6g27(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmppmpp<2,3,4,5,0,1>(ep,mpc);} //++-++-

template <class T> complex<T> R6g5(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<0,1,2,3,4,5>(em,pmc);} //+-+---
template <class T> complex<T> R6g10(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<1,2,3,4,5,0>(em,pmc);} //-+-+-- 
template <class T> complex<T> R6g20(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<2,3,4,5,0,1>(em,pmc);} //--+-+-
template <class T> complex<T> R6g40(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<3,4,5,0,1,2>(em,pmc);} //---+-+
template <class T> complex<T> R6g17(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<4,5,0,1,2,3>(em,pmc);} //+---+-
template <class T> complex<T> R6g34(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmpmmm<5,0,1,2,3,4>(em,pmc);} //-+---+

template <class T> complex<T> R6g9(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmmpmm<0,1,2,3,4,5>(em,pmc);} //+--+--
template <class T> complex<T> R6g18(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmmpmm<1,2,3,4,5,0>(em,pmc);} //-+--+- 
template <class T> complex<T> R6g36(const eval_param<T>& em,const mass_param_coll& pmc){return R6gpmmpmm<2,3,4,5,0,1>(em,pmc);} //--+--+

template <class T> complex<T> R6g56(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<0,1,2,3,4,5>(ep,mpc);} //---+++
template <class T> complex<T> R6g49(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<1,2,3,4,5,0>(ep,mpc);} //+---++
template <class T> complex<T> R6g35(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<2,3,4,5,0,1>(ep,mpc);} //++---+
template <class T> complex<T> R6g7(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<3,4,5,0,1,2>(ep,mpc);} //+++---
template <class T> complex<T> R6g14(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<4,5,0,1,2,3>(ep,mpc);} //-+++--
template <class T> complex<T> R6g28(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmmppp<5,0,1,2,3,4>(ep,mpc);} //--+++-

template <class T> complex<T> R6g63(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g6p(ep,mpc);} //++++++
template <class T> complex<T> R6g62(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<0,1,2,3,4,5>(ep,mpc);} // -+++++
template <class T> complex<T> R6g61(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<1,2,3,4,5,0>(ep,mpc);} // +-++++
template <class T> complex<T> R6g59(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<2,3,4,5,0,1>(ep,mpc);} // ++-+++
template <class T> complex<T> R6g55(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<3,4,5,0,1,2>(ep,mpc);} // +++-++
template <class T> complex<T> R6g47(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<4,5,0,1,2,3>(ep,mpc);} // ++++-+
template <class T> complex<T> R6g31(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1m<5,0,1,2,3,4>(ep,mpc);} // +++++-

template <class T> complex<T> R6g0(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g6m(ep,mpc);} //------
template <class T> complex<T> R6g1(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<0,1,2,3,4,5>(ep,mpc);} // +-----
template <class T> complex<T> R6g2(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<1,2,3,4,5,0>(ep,mpc);} // -+----
template <class T> complex<T> R6g4(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<2,3,4,5,0,1>(ep,mpc);} // --+---
template <class T> complex<T> R6g8(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<3,4,5,0,1,2>(ep,mpc);} // ---+--
template <class T> complex<T> R6g16(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<4,5,0,1,2,3>(ep,mpc);} // ----+-
template <class T> complex<T> R6g32(const eval_param<T>& ep,const mass_param_coll& mpc){return R6g1p<5,0,1,2,3,4>(ep,mpc);} // -----+

/* for testing only
template <class T> complex<T> R6g42(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmpmpmp<0,1,2,3,4,5>(ep,mpc);} //-+-+-+ 

template <class T> complex<T> R6g52(const eval_param<T>& ep,const mass_param_coll& mpc){return R6gmmpmpp<0,1,2,3,4,5>(ep,mpc);} //--+-++ 
*/

	//------------------  NF   -----------------------------
#if 1
template <class T> complex<T> R6g3_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<0,1,2,3,4,5>(ep,mpc);}  //++----
template <class T> complex<T> R6g6_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<1,2,3,4,5,0>(ep,mpc);}  //-++---
template <class T> complex<T> R6g12_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<2,3,4,5,0,1>(ep,mpc);} //--++--
template <class T> complex<T> R6g24_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<3,4,5,0,1,2>(ep,mpc);} //---++-
template <class T> complex<T> R6g48_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<4,5,0,1,2,3>(ep,mpc);} //----++
template <class T> complex<T> R6g33_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2p<5,0,1,2,3,4>(ep,mpc);} //+----+
//
template <class T> complex<T> R6g30_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<5,0,1,2,3,4>(ep,mpc);} //-++++-
template <class T> complex<T> R6g15_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<4,5,0,1,2,3>(ep,mpc);} //++++--
template <class T> complex<T> R6g39_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<3,4,5,0,1,2>(ep,mpc);} //+++--+
template <class T> complex<T> R6g51_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<2,3,4,5,0,1>(ep,mpc);} //++--++
template <class T> complex<T> R6g57_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<1,2,3,4,5,0>(ep,mpc);} //+--+++
template <class T> complex<T> R6g60_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g2m<0,1,2,3,4,5>(ep,mpc);} //--++++


template <class T> complex<T> R6g58_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<0,1,2,3,4,5>(ep,mpc);} //-+-+++
template <class T> complex<T> R6g53_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<1,2,3,4,5,0>(ep,mpc);} //+-+-++ 
template <class T> complex<T> R6g43_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<2,3,4,5,0,1>(ep,mpc);} //++-+-+
template <class T> complex<T> R6g23_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<3,4,5,0,1,2>(ep,mpc);} //+++-+-
template <class T> complex<T> R6g46_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<4,5,0,1,2,3>(ep,mpc);} //-+++-+
template <class T> complex<T> R6g29_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmpmppp<5,0,1,2,3,4>(ep,mpc);} //+-+++-

template <class T> complex<T> R6g54_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmppmpp<0,1,2,3,4,5>(ep,mpc);} //-++-++
template <class T> complex<T> R6g45_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmppmpp<1,2,3,4,5,0>(ep,mpc);} //+-++-+
template <class T> complex<T> R6g27_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmppmpp<2,3,4,5,0,1>(ep,mpc);} //++-++-

template <class T> complex<T> R6g5_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<0,1,2,3,4,5>(em,pmc);} //+-+---
template <class T> complex<T> R6g10_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<1,2,3,4,5,0>(em,pmc);} //-+-+-- 
template <class T> complex<T> R6g20_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<2,3,4,5,0,1>(em,pmc);} //--+-+-
template <class T> complex<T> R6g40_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<3,4,5,0,1,2>(em,pmc);} //---+-+
template <class T> complex<T> R6g17_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<4,5,0,1,2,3>(em,pmc);} //+---+-
template <class T> complex<T> R6g34_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmpmmm<5,0,1,2,3,4>(em,pmc);} //-+---+

template <class T> complex<T> R6g9_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmmpmm<0,1,2,3,4,5>(em,pmc);} //+--+--
template <class T> complex<T> R6g18_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmmpmm<1,2,3,4,5,0>(em,pmc);} //-+--+- 
template <class T> complex<T> R6g36_nf(const eval_param<T>& em,const mass_param_coll& pmc){return -R6gpmmpmm<2,3,4,5,0,1>(em,pmc);} //--+--+

template <class T> complex<T> R6g56_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<0,1,2,3,4,5>(ep,mpc);} //---+++
template <class T> complex<T> R6g49_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<1,2,3,4,5,0>(ep,mpc);} //+---++
template <class T> complex<T> R6g35_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<2,3,4,5,0,1>(ep,mpc);} //++---+
template <class T> complex<T> R6g7_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<3,4,5,0,1,2>(ep,mpc);} //+++---
template <class T> complex<T> R6g14_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<4,5,0,1,2,3>(ep,mpc);} //-+++--
template <class T> complex<T> R6g28_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6gmmmppp<5,0,1,2,3,4>(ep,mpc);} //--+++-

template <class T> complex<T> R6g63_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g6p(ep,mpc);} //++++++
template <class T> complex<T> R6g62_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<0,1,2,3,4,5>(ep,mpc);} // -+++++
template <class T> complex<T> R6g61_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<1,2,3,4,5,0>(ep,mpc);} // +-++++
template <class T> complex<T> R6g59_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<2,3,4,5,0,1>(ep,mpc);} // ++-+++
template <class T> complex<T> R6g55_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<3,4,5,0,1,2>(ep,mpc);} // +++-++
template <class T> complex<T> R6g47_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<4,5,0,1,2,3>(ep,mpc);} // ++++-+
template <class T> complex<T> R6g31_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1m<5,0,1,2,3,4>(ep,mpc);} // +++++-

template <class T> complex<T> R6g0_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g6m(ep,mpc);} //------
template <class T> complex<T> R6g1_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<0,1,2,3,4,5>(ep,mpc);} // +-----
template <class T> complex<T> R6g2_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<1,2,3,4,5,0>(ep,mpc);} // -+----
template <class T> complex<T> R6g4_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<2,3,4,5,0,1>(ep,mpc);} // --+---
template <class T> complex<T> R6g8_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<3,4,5,0,1,2>(ep,mpc);} // ---+--
template <class T> complex<T> R6g16_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<4,5,0,1,2,3>(ep,mpc);} // ----+-
template <class T> complex<T> R6g32_nf(const eval_param<T>& ep,const mass_param_coll& mpc){return -R6g1p<5,0,1,2,3,4>(ep,mpc);} // -----+
#endif

template <class T> complex<T> (*R6g_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {

switch (hc) {
case 3 : return &R6g3;
case 6 : return &R6g6;
case 12: return &R6g12;
case 24: return &R6g24;
case 48: return &R6g48;
case 33: return &R6g33;
//
case 30: return &R6g30;
case 15: return &R6g15;
case 39: return &R6g39;
case 51: return &R6g51;
case 57: return &R6g57;
case 60: return &R6g60;
//
case 58: return &R6g58;
case 53: return &R6g53;
case 43: return &R6g43;
case 23: return &R6g23;
case 46: return &R6g46;
case 29: return &R6g29;
//
case 5: return &R6g5;
case 10: return &R6g10;
case 20: return &R6g20;
case 40: return &R6g40;
case 17: return &R6g17;
case 34: return &R6g34;
//
case 54: return &R6g54;
case 45: return &R6g45;
case 27: return &R6g27;
//
case 9: return &R6g9;
case 18: return &R6g18;
case 36: return &R6g36;
//
case 56: return &R6g56;
case 49: return &R6g49;
case 35: return &R6g35;
case 7: return &R6g7;
case 14: return &R6g14;
case 28: return &R6g28;
//
//
case 63: return &R6g63;
//
case 62: return &R6g62;
case 61: return &R6g61;
case 59: return &R6g59;
case 55: return &R6g55;
case 47: return &R6g47;
case 31: return &R6g31;
//
case 0: return &R6g0;
//
case 1: return &R6g1;
case 2: return &R6g2;
case 4: return &R6g4;
case 8: return &R6g8;
case 16: return &R6g16;
case 32: return &R6g32;
/*testing only
case 42: return &R6g42;
case 52: return &R6g52;
*/
default: return 0;
	}
}


			//------------------  NF   -----------------------------

template <class T> complex<T> (*R6g_Ptr_nf_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {

switch (hc) {
case 3 : return &R6g3_nf;
case 6 : return &R6g6_nf;
case 12: return &R6g12_nf;
case 24: return &R6g24_nf;
case 48: return &R6g48_nf;
case 33: return &R6g33_nf;
//
case 30: return &R6g30_nf;
case 15: return &R6g15_nf;
case 39: return &R6g39_nf;
case 51: return &R6g51_nf;
case 57: return &R6g57_nf;
case 60: return &R6g60_nf;
//
case 58: return &R6g58_nf;
case 53: return &R6g53_nf;
case 43: return &R6g43_nf;
case 23: return &R6g23_nf;
case 46: return &R6g46_nf;
case 29: return &R6g29_nf;
//
case 5: return &R6g5_nf;
case 10: return &R6g10_nf;
case 20: return &R6g20_nf;
case 40: return &R6g40_nf;
case 17: return &R6g17_nf;
case 34: return &R6g34_nf;
//
case 54: return &R6g54_nf;
case 45: return &R6g45_nf;
case 27: return &R6g27_nf;
//
case 9: return &R6g9_nf;
case 18: return &R6g18_nf;
case 36: return &R6g36_nf;
//
case 56: return &R6g56_nf;
case 49: return &R6g49_nf;
case 35: return &R6g35_nf;
case 7: return &R6g7_nf;
case 14: return &R6g14_nf;
case 28: return &R6g28_nf;
//
//
case 63: return &R6g63_nf;
//
case 62: return &R6g62_nf;
case 61: return &R6g61_nf;
case 59: return &R6g59_nf;
case 55: return &R6g55_nf;
case 47: return &R6g47_nf;
case 31: return &R6g31_nf;
//
case 0: return &R6g0_nf;
//
case 1: return &R6g1_nf;
case 2: return &R6g2_nf;
case 4: return &R6g4_nf;
case 8: return &R6g8_nf;
case 16: return &R6g16_nf;
case 32: return &R6g32_nf;
//
	default: return 0;
	}
}

#if 1
template <class T> complex<T> R6g(const process& p,const eval_param<T>& ep,const mass_param_coll& mpc){

switch (helcode_g(p)) {
case 3 : return R6g3(ep,mpc);
case 6 : return R6g6(ep,mpc);
case 12: return R6g12(ep,mpc);
case 24: return R6g24(ep,mpc);
case 48: return R6g48(ep,mpc);
case 33: return R6g33(ep,mpc);
//
case 30: return R6g30(ep,mpc);
case 15: return R6g15(ep,mpc);
case 39: return R6g39(ep,mpc);
case 51: return R6g51(ep,mpc);
case 57: return R6g57(ep,mpc);
case 60: return R6g60(ep,mpc);
//
case 58: return R6g58(ep,mpc);
case 53: return R6g53(ep,mpc);
case 43: return R6g43(ep,mpc);
case 23: return R6g23(ep,mpc);
case 46: return R6g46(ep,mpc);
case 29: return R6g29(ep,mpc);
//
case 5: return R6g5(ep,mpc);
case 10: return R6g10(ep,mpc);
case 20: return R6g20(ep,mpc);
case 40: return R6g40(ep,mpc);
case 17: return R6g17(ep,mpc);
case 34: return R6g34(ep,mpc);
//
case 54: return R6g54(ep,mpc);
case 45: return R6g45(ep,mpc);
case 27: return R6g27(ep,mpc);
//
case 9: return R6g9(ep,mpc);
case 18: return R6g18(ep,mpc);
case 36: return R6g36(ep,mpc);
//
case 56: return R6g56(ep,mpc);
case 49: return R6g49(ep,mpc);
case 35: return R6g35(ep,mpc);
case 7: return R6g7(ep,mpc);
case 14: return R6g14(ep,mpc);
case 28: return R6g28(ep,mpc);
//
//
case 63: return R6g63(ep,mpc);
//
case 62: return R6g62(ep,mpc);
case 61: return R6g61(ep,mpc);
case 59: return R6g59(ep,mpc);
case 55: return R6g55(ep,mpc);
case 47: return R6g47(ep,mpc);
case 31: return R6g31(ep,mpc);
//
case 0: return R6g0(ep,mpc);
//
case 1: return R6g1(ep,mpc);
case 2: return R6g2(ep,mpc);
case 4: return R6g4(ep,mpc);
case 8: return R6g8(ep,mpc);
case 16: return R6g16(ep,mpc);
case 32: return R6g32(ep,mpc);

default: {_WARNING3("using unknown known rational term for ",p," returned 0;"); return complex<T>(0,0);};

}
}



template complex<R> R6g(const process& p,const eval_param<R>& ep,const mass_param_coll& mpc);
template complex<RHP> R6g(const process& p,const eval_param<RHP>& ep,const mass_param_coll& mpc);
template complex<RVHP> R6g(const process& p,const eval_param<RVHP>& ep,const mass_param_coll& mpc);
#endif


template complex<R> (*R6g_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R6g_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R6g_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R6g_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif



template complex<R> (*R6g_Ptr_nf_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R6g_Ptr_nf_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R6g_Ptr_nf_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R6g_Ptr_nf_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif


}

