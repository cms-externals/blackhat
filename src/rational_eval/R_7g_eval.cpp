

#include "R_7g_eval.h"
#include "eval_param.h" 
#include "R_alln_eval.h"
 using namespace std;
namespace BH   {


#define _ONLY_X_PART 1
#define X 2
//#define _VERBOSE 1

int helcode_g(const process& pro);

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T> R7gap(const eval_param<T>& ep,const mass_param_coll& mpc)
{
  return complex<T>(0,0);

}

//#if _FAST_COMPILE == 0

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T> R7g1m(const eval_param<T>& ep,const mass_param_coll& mpc)
{
  return complex<T>(0,1)/T(3)*(((ep.spa(i1,i2)*ep.spa(i1,i3)*(ep.spa(i1,i2)*ep.spa(i3,i1)*ep.spb(i2,i3) +
          ep.spa(i1,i2)*ep.spa(i4,i1)*ep.spb(i2,i4) + ep.spa(i1,i2)*ep.spa(i5,i1)*ep.spb(i2,i5) +
          ep.spa(i1,i2)*ep.spa(i6,i1)*ep.spb(i2,i6) + ep.spa(i1,i2)*ep.spa(i7,i1)*ep.spb(i2,i7) +
          ep.spa(i1,i3)*ep.spa(i4,i1)*ep.spb(i3,i4) + ep.spa(i1,i3)*ep.spa(i5,i1)*ep.spb(i3,i5) +
          ep.spa(i1,i3)*ep.spa(i6,i1)*ep.spb(i3,i6) + ep.spa(i1,i3)*ep.spa(i7,i1)*ep.spb(i3,i7)))/
      ep.spa(i2,i3) + (ep.spa(i1,i3)*ep.spa(i1,i4)*
        (ep.spa(i1,i3)*ep.spa(i4,i1)*ep.spb(i3,i4) + ep.spa(i1,i3)*ep.spa(i5,i1)*ep.spb(i3,i5) +
          ep.spa(i1,i3)*ep.spa(i6,i1)*ep.spb(i3,i6) + ep.spa(i1,i3)*ep.spa(i7,i1)*ep.spb(i3,i7) +
          ep.spa(i1,i4)*ep.spa(i5,i1)*ep.spb(i4,i5) + ep.spa(i1,i4)*ep.spa(i6,i1)*ep.spb(i4,i6) +
          ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7)))/ep.spa(i3,i4) +
     (ep.spa(i1,i4)*ep.spa(i1,i5)*(ep.spa(i1,i4)*ep.spa(i5,i1)*ep.spb(i4,i5) +
          ep.spa(i1,i4)*ep.spa(i6,i1)*ep.spb(i4,i6) + ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7) +
          ep.spa(i1,i5)*ep.spa(i6,i1)*ep.spb(i5,i6) + ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7)))/
      ep.spa(i4,i5) + (pow(ep.spa(i1,i6),2)*ep.spa(i1,i7)*ep.spa(i7,i1)*ep.spb(i6,i7))/ep.spa(i6,i7) +
     (ep.spa(i1,i5)*ep.spa(i1,i6)*(ep.spa(i1,i5)*ep.spa(i6,i1)*ep.spb(i5,i6) +
          ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7) + ep.spa(i1,i6)*ep.spa(i7,i1)*ep.spb(i6,i7)))/
      ep.spa(i5,i6) + (ep.spa(i2,i3)*ep.spa(i4,i5)*
        pow(ep.spa(i1,i3)*ep.spa(i5,i1)*ep.spb(i3,i5) + ep.spa(i1,i3)*ep.spa(i6,i1)*ep.spb(i3,i6) +
          ep.spa(i1,i3)*ep.spa(i7,i1)*ep.spb(i3,i7) + ep.spa(i1,i4)*ep.spa(i5,i1)*ep.spb(i4,i5) +
          ep.spa(i1,i4)*ep.spa(i6,i1)*ep.spb(i4,i6) + ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i5,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i5) +
          ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i6,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i6) +
          ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i4)*ep.spb(i2,i4))*
        (ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spb(i2,i4))*
        (ep.spa(i1,i5)*ep.spa(i3,i2)*ep.spb(i5,i3) + ep.spa(i1,i5)*ep.spa(i4,i2)*ep.spb(i5,i4) +
          ep.spa(i1,i6)*ep.spa(i3,i2)*ep.spb(i6,i3) + ep.spa(i1,i6)*ep.spa(i4,i2)*ep.spb(i6,i4) +
          ep.spa(i1,i7)*ep.spa(i3,i2)*ep.spb(i7,i3) + ep.spa(i1,i7)*ep.spa(i4,i2)*ep.spb(i7,i4))*
        (ep.spa(i1,i5)*ep.spa(i3,i3)*ep.spb(i5,i3) + ep.spa(i1,i5)*ep.spa(i4,i3)*ep.spb(i5,i4) +
          ep.spa(i1,i6)*ep.spa(i3,i3)*ep.spb(i6,i3) + ep.spa(i1,i6)*ep.spa(i4,i3)*ep.spb(i6,i4) +
          ep.spa(i1,i7)*ep.spa(i3,i3)*ep.spb(i7,i3) + ep.spa(i1,i7)*ep.spa(i4,i3)*ep.spb(i7,i4))*
        ep.s(i3,i4)) + (ep.spa(i3,i4)*ep.spa(i5,i6)*
        pow(ep.spa(i1,i4)*ep.spa(i6,i1)*ep.spb(i4,i6) + ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7) +
          ep.spa(i1,i5)*ep.spa(i6,i1)*ep.spb(i5,i6) + ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i6,i1)*ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i5,i6) +
          ep.spa(i1,i3)*pow(ep.spa(i4,i5),2)*ep.spa(i6,i1)*ep.spb(i3,i4)*ep.spb(i5,i4)*ep.spb(i5,i6) +
          ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i5,i7) +
          ep.spa(i1,i3)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*ep.spb(i3,i4)*ep.spb(i5,i4)*ep.spb(i5,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spb(i2,i4) + ep.spa(i1,i2)*ep.spa(i5,i5)*ep.spb(i2,i5) +
          ep.spa(i1,i3)*ep.spa(i4,i5)*ep.spb(i3,i4) + ep.spa(i1,i3)*ep.spa(i5,i5)*ep.spb(i3,i5))*
        (ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spb(i2,i4) + ep.spa(i1,i2)*ep.spa(i5,i6)*ep.spb(i2,i5) +
          ep.spa(i1,i3)*ep.spa(i4,i6)*ep.spb(i3,i4) + ep.spa(i1,i3)*ep.spa(i5,i6)*ep.spb(i3,i5))*
        (ep.spa(i1,i6)*ep.spa(i4,i3)*ep.spb(i6,i4) + ep.spa(i1,i6)*ep.spa(i5,i3)*ep.spb(i6,i5) +
          ep.spa(i1,i7)*ep.spa(i4,i3)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i3)*ep.spb(i7,i5))*
        (ep.spa(i1,i6)*ep.spa(i4,i4)*ep.spb(i6,i4) + ep.spa(i1,i6)*ep.spa(i5,i4)*ep.spb(i6,i5) +
          ep.spa(i1,i7)*ep.spa(i4,i4)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i4)*ep.spb(i7,i5))*
        ep.s(i4,i5)) + (ep.spa(i4,i5)*ep.spa(i6,i7)*
        pow(ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7) + ep.spa(i1,i6)*ep.spa(i7,i1)*ep.spb(i6,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i5)*ep.spb(i6,i5)*ep.spb(i6,i7) +
          ep.spa(i1,i3)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i3,i5)*ep.spb(i6,i5)*ep.spb(i6,i7) +
          ep.spa(i1,i4)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i4,i5)*ep.spb(i6,i5)*ep.spb(i6,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i5,i6)*ep.spb(i2,i5) + ep.spa(i1,i2)*ep.spa(i6,i6)*ep.spb(i2,i6) +
          ep.spa(i1,i3)*ep.spa(i5,i6)*ep.spb(i3,i5) + ep.spa(i1,i3)*ep.spa(i6,i6)*ep.spb(i3,i6) +
          ep.spa(i1,i4)*ep.spa(i5,i6)*ep.spb(i4,i5) + ep.spa(i1,i4)*ep.spa(i6,i6)*ep.spb(i4,i6))*
        (ep.spa(i1,i2)*ep.spa(i5,i7)*ep.spb(i2,i5) + ep.spa(i1,i2)*ep.spa(i6,i7)*ep.spb(i2,i6) +
          ep.spa(i1,i3)*ep.spa(i5,i7)*ep.spb(i3,i5) + ep.spa(i1,i3)*ep.spa(i6,i7)*ep.spb(i3,i6) +
          ep.spa(i1,i4)*ep.spa(i5,i7)*ep.spb(i4,i5) + ep.spa(i1,i4)*ep.spa(i6,i7)*ep.spb(i4,i6))*
        (ep.spa(i1,i7)*ep.spa(i5,i4)*ep.spb(i7,i5) + ep.spa(i1,i7)*ep.spa(i6,i4)*ep.spb(i7,i6))*
        (ep.spa(i1,i7)*ep.spa(i5,i5)*ep.spb(i7,i5) + ep.spa(i1,i7)*ep.spa(i6,i5)*ep.spb(i7,i6))*
        ep.s(i5,i6)) + (ep.spa(i2,i3)*ep.spa(i5,i6)*
        pow(ep.spa(i1,i3)*ep.spa(i6,i1)*ep.spb(i3,i6) + ep.spa(i1,i3)*ep.spa(i7,i1)*ep.spb(i3,i7) +
          ep.spa(i1,i4)*ep.spa(i6,i1)*ep.spb(i4,i6) + ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7) +
          ep.spa(i1,i5)*ep.spa(i6,i1)*ep.spb(i5,i6) + ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i6,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i6) +
          ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i6,i1)*ep.spb(i2,i3)*ep.spb(i4,i6)*
           ep.spb(i5,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i6,i1)*ep.spb(i2,i4)*
           ep.spb(i4,i6)*ep.spb(i5,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i7,i1)*
           ep.spb(i2,i3)*ep.spb(i4,i7)*ep.spb(i5,i3) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i4,i7)*
           ep.spb(i5,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i6,i1)*ep.spb(i2,i3)*
           ep.spb(i4,i3)*ep.spb(i5,i6) + ep.spa(i1,i2)*pow(ep.spa(i3,i5),2)*ep.spa(i6,i1)*
           ep.spb(i2,i3)*ep.spb(i5,i3)*ep.spb(i5,i6) +
          ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i6,i1)*ep.spb(i2,i4)*ep.spb(i5,i3)*
           ep.spb(i5,i6) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i6,i1)*ep.spb(i2,i3)*
           ep.spb(i5,i4)*ep.spb(i5,i6) + ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i6,i1)*
           ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i5,i6) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*
           ep.spb(i5,i7) + ep.spa(i1,i2)*pow(ep.spa(i3,i5),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*
           ep.spb(i5,i3)*ep.spb(i5,i7) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i7,i1)*
           ep.spb(i2,i4)*ep.spb(i5,i3)*ep.spb(i5,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i5,i4)*
           ep.spb(i5,i7) + ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i5,i4)*ep.spb(i5,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spb(i2,i4) +
          ep.spa(i1,i2)*ep.spa(i5,i5)*ep.spb(i2,i5))*
        (ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spb(i2,i4) +
          ep.spa(i1,i2)*ep.spa(i5,i6)*ep.spb(i2,i5))*
        (ep.spa(i1,i6)*ep.spa(i3,i2)*ep.spb(i6,i3) + ep.spa(i1,i6)*ep.spa(i4,i2)*ep.spb(i6,i4) +
          ep.spa(i1,i6)*ep.spa(i5,i2)*ep.spb(i6,i5) + ep.spa(i1,i7)*ep.spa(i3,i2)*ep.spb(i7,i3) +
          ep.spa(i1,i7)*ep.spa(i4,i2)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i2)*ep.spb(i7,i5))*
        (ep.spa(i1,i6)*ep.spa(i3,i3)*ep.spb(i6,i3) + ep.spa(i1,i6)*ep.spa(i4,i3)*ep.spb(i6,i4) +
          ep.spa(i1,i6)*ep.spa(i5,i3)*ep.spb(i6,i5) + ep.spa(i1,i7)*ep.spa(i3,i3)*ep.spb(i7,i3) +
          ep.spa(i1,i7)*ep.spa(i4,i3)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i3)*ep.spb(i7,i5))*
        ep.s(i3,i4,i5)) + (ep.spa(i3,i4)*ep.spa(i6,i7)*
        pow(ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7) + ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7) +
          ep.spa(i1,i6)*ep.spa(i7,i1)*ep.spb(i6,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i5,i7) +
          ep.spa(i1,i3)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*ep.spb(i3,i4)*ep.spb(i5,i4)*ep.spb(i5,i7) +
          ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i5,i7)*
           ep.spb(i6,i4) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i5)*
           ep.spb(i5,i7)*ep.spb(i6,i4) + ep.spa(i1,i3)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*
           ep.spb(i3,i4)*ep.spb(i5,i7)*ep.spb(i6,i4) +
          ep.spa(i1,i3)*ep.spa(i4,i5)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i3,i5)*ep.spb(i5,i7)*
           ep.spb(i6,i4) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i5,i4)*ep.spb(i6,i7) + ep.spa(i1,i3)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*
           ep.spb(i3,i4)*ep.spb(i5,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*pow(ep.spa(i4,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i6,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i5)*ep.spb(i6,i4)*
           ep.spb(i6,i7) + ep.spa(i1,i3)*pow(ep.spa(i4,i6),2)*ep.spa(i7,i1)*ep.spb(i3,i4)*
           ep.spb(i6,i4)*ep.spb(i6,i7) + ep.spa(i1,i3)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i3,i5)*ep.spb(i6,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i6,i5)*
           ep.spb(i6,i7) + ep.spa(i1,i2)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i5)*
           ep.spb(i6,i5)*ep.spb(i6,i7) + ep.spa(i1,i3)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i3,i4)*ep.spb(i6,i5)*ep.spb(i6,i7) +
          ep.spa(i1,i3)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i3,i5)*ep.spb(i6,i5)*ep.spb(i6,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spb(i2,i4) + ep.spa(i1,i2)*ep.spa(i5,i6)*ep.spb(i2,i5) +
          ep.spa(i1,i2)*ep.spa(i6,i6)*ep.spb(i2,i6) + ep.spa(i1,i3)*ep.spa(i4,i6)*ep.spb(i3,i4) +
          ep.spa(i1,i3)*ep.spa(i5,i6)*ep.spb(i3,i5) + ep.spa(i1,i3)*ep.spa(i6,i6)*ep.spb(i3,i6))*
        (ep.spa(i1,i2)*ep.spa(i4,i7)*ep.spb(i2,i4) + ep.spa(i1,i2)*ep.spa(i5,i7)*ep.spb(i2,i5) +
          ep.spa(i1,i2)*ep.spa(i6,i7)*ep.spb(i2,i6) + ep.spa(i1,i3)*ep.spa(i4,i7)*ep.spb(i3,i4) +
          ep.spa(i1,i3)*ep.spa(i5,i7)*ep.spb(i3,i5) + ep.spa(i1,i3)*ep.spa(i6,i7)*ep.spb(i3,i6))*
        (ep.spa(i1,i7)*ep.spa(i4,i3)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i3)*ep.spb(i7,i5) +
          ep.spa(i1,i7)*ep.spa(i6,i3)*ep.spb(i7,i6))*
        (ep.spa(i1,i7)*ep.spa(i4,i4)*ep.spb(i7,i4) + ep.spa(i1,i7)*ep.spa(i5,i4)*ep.spb(i7,i5) +
          ep.spa(i1,i7)*ep.spa(i6,i4)*ep.spb(i7,i6))*ep.s(i4,i5,i6)) +
     (ep.spa(i2,i3)*ep.spa(i6,i7)*pow(ep.spa(i1,i3)*ep.spa(i7,i1)*ep.spb(i3,i7) +
          ep.spa(i1,i4)*ep.spa(i7,i1)*ep.spb(i4,i7) + ep.spa(i1,i5)*ep.spa(i7,i1)*ep.spb(i5,i7) +
          ep.spa(i1,i6)*ep.spa(i7,i1)*ep.spb(i6,i7),3)*
        (ep.spa(i1,i2)*pow(ep.spa(i3,i4),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i4,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i7)*
           ep.spb(i5,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i4,i7)*ep.spb(i5,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i5)*ep.spa(i7,i1)*
           ep.spb(i2,i3)*ep.spb(i4,i3)*ep.spb(i5,i7) +
          ep.spa(i1,i2)*pow(ep.spa(i3,i5),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i5,i3)*ep.spb(i5,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i5,i3)*
           ep.spb(i5,i7) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i3)*
           ep.spb(i5,i4)*ep.spb(i5,i7) + ep.spa(i1,i2)*pow(ep.spa(i4,i5),2)*ep.spa(i7,i1)*
           ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i5,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i7)*
           ep.spb(i6,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i4,i7)*ep.spb(i6,i3) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i5)*ep.spb(i4,i7)*ep.spb(i6,i3) +
          ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i3,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i5,i7)*
           ep.spb(i6,i3) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i5,i7)*ep.spb(i6,i3) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i5)*ep.spb(i5,i7)*ep.spb(i6,i3) +
          ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i5,i7)*
           ep.spb(i6,i4) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*
           ep.spb(i5,i7)*ep.spb(i6,i4) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i5)*ep.spb(i5,i7)*ep.spb(i6,i4) +
          ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i3,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i4,i3)*
           ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i4)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*
           ep.spb(i4,i5)*ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i3,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i3)*ep.spb(i5,i3)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i4,i5)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i5,i3)*
           ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*
           ep.spb(i5,i4)*ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i4,i5)*ep.spa(i4,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i4)*ep.spb(i5,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*pow(ep.spa(i3,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i3)*ep.spb(i6,i3)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i4,i6)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i6,i3)*
           ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i5)*
           ep.spb(i6,i3)*ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i4,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i3)*ep.spb(i6,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*pow(ep.spa(i4,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i4)*ep.spb(i6,i4)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i5)*ep.spb(i6,i4)*
           ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*ep.spb(i2,i3)*
           ep.spb(i6,i5)*ep.spb(i6,i7) + ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spa(i5,i6)*ep.spa(i7,i1)*
           ep.spb(i2,i4)*ep.spb(i6,i5)*ep.spb(i6,i7) +
          ep.spa(i1,i2)*pow(ep.spa(i5,i6),2)*ep.spa(i7,i1)*ep.spb(i2,i5)*ep.spb(i6,i5)*ep.spb(i6,i7)))/
      ((ep.spa(i1,i2)*ep.spa(i3,i6)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i6)*ep.spb(i2,i4) +
          ep.spa(i1,i2)*ep.spa(i5,i6)*ep.spb(i2,i5) + ep.spa(i1,i2)*ep.spa(i6,i6)*ep.spb(i2,i6))*
        (ep.spa(i1,i2)*ep.spa(i3,i7)*ep.spb(i2,i3) + ep.spa(i1,i2)*ep.spa(i4,i7)*ep.spb(i2,i4) +
          ep.spa(i1,i2)*ep.spa(i5,i7)*ep.spb(i2,i5) + ep.spa(i1,i2)*ep.spa(i6,i7)*ep.spb(i2,i6))*
        (ep.spa(i1,i7)*ep.spa(i3,i2)*ep.spb(i7,i3) + ep.spa(i1,i7)*ep.spa(i4,i2)*ep.spb(i7,i4) +
          ep.spa(i1,i7)*ep.spa(i5,i2)*ep.spb(i7,i5) + ep.spa(i1,i7)*ep.spa(i6,i2)*ep.spb(i7,i6))*
        (ep.spa(i1,i7)*ep.spa(i3,i3)*ep.spb(i7,i3) + ep.spa(i1,i7)*ep.spa(i4,i3)*ep.spb(i7,i4) +
          ep.spa(i1,i7)*ep.spa(i5,i3)*ep.spb(i7,i5) + ep.spa(i1,i7)*ep.spa(i6,i3)*ep.spb(i7,i6))*
        ep.s(i3,i4,i5,i6)))/
   (ep.spa(i1,i2)*ep.spa(i2,i3)*ep.spa(i3,i4)*ep.spa(i4,i5)*ep.spa(i5,i6)*ep.spa(i6,i7)*
     ep.spa(i7,i1)));
}

template <int i1, int i2, int i3, int i4, int i5, int i6, int i7, class T> complex<T> R7g1p(const eval_param<T>& ep,const mass_param_coll& mpc)
{
  return complex<T>(0,1)/T(3)*(((ep.spb(i1,i2)*ep.spb(i1,i3)*(ep.spb(i1,i2)*ep.spb(i3,i1)*ep.spa(i2,i3) +
          ep.spb(i1,i2)*ep.spb(i4,i1)*ep.spa(i2,i4) + ep.spb(i1,i2)*ep.spb(i5,i1)*ep.spa(i2,i5) +
          ep.spb(i1,i2)*ep.spb(i6,i1)*ep.spa(i2,i6) + ep.spb(i1,i2)*ep.spb(i7,i1)*ep.spa(i2,i7) +
          ep.spb(i1,i3)*ep.spb(i4,i1)*ep.spa(i3,i4) + ep.spb(i1,i3)*ep.spb(i5,i1)*ep.spa(i3,i5) +
          ep.spb(i1,i3)*ep.spb(i6,i1)*ep.spa(i3,i6) + ep.spb(i1,i3)*ep.spb(i7,i1)*ep.spa(i3,i7)))/
      ep.spb(i2,i3) + (ep.spb(i1,i3)*ep.spb(i1,i4)*
        (ep.spb(i1,i3)*ep.spb(i4,i1)*ep.spa(i3,i4) + ep.spb(i1,i3)*ep.spb(i5,i1)*ep.spa(i3,i5) +
          ep.spb(i1,i3)*ep.spb(i6,i1)*ep.spa(i3,i6) + ep.spb(i1,i3)*ep.spb(i7,i1)*ep.spa(i3,i7) +
          ep.spb(i1,i4)*ep.spb(i5,i1)*ep.spa(i4,i5) + ep.spb(i1,i4)*ep.spb(i6,i1)*ep.spa(i4,i6) +
          ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7)))/ep.spb(i3,i4) +
     (ep.spb(i1,i4)*ep.spb(i1,i5)*(ep.spb(i1,i4)*ep.spb(i5,i1)*ep.spa(i4,i5) +
          ep.spb(i1,i4)*ep.spb(i6,i1)*ep.spa(i4,i6) + ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7) +
          ep.spb(i1,i5)*ep.spb(i6,i1)*ep.spa(i5,i6) + ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7)))/
      ep.spb(i4,i5) + (pow(ep.spb(i1,i6),2)*ep.spb(i1,i7)*ep.spb(i7,i1)*ep.spa(i6,i7))/ep.spb(i6,i7) +
     (ep.spb(i1,i5)*ep.spb(i1,i6)*(ep.spb(i1,i5)*ep.spb(i6,i1)*ep.spa(i5,i6) +
          ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7) + ep.spb(i1,i6)*ep.spb(i7,i1)*ep.spa(i6,i7)))/
      ep.spb(i5,i6) + (ep.spb(i2,i3)*ep.spb(i4,i5)*
        pow(ep.spb(i1,i3)*ep.spb(i5,i1)*ep.spa(i3,i5) + ep.spb(i1,i3)*ep.spb(i6,i1)*ep.spa(i3,i6) +
          ep.spb(i1,i3)*ep.spb(i7,i1)*ep.spa(i3,i7) + ep.spb(i1,i4)*ep.spb(i5,i1)*ep.spa(i4,i5) +
          ep.spb(i1,i4)*ep.spb(i6,i1)*ep.spa(i4,i6) + ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i5,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i5) +
          ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i6,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i6) +
          ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i4)*ep.spa(i2,i4))*
        (ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spa(i2,i4))*
        (ep.spb(i1,i5)*ep.spb(i3,i2)*ep.spa(i5,i3) + ep.spb(i1,i5)*ep.spb(i4,i2)*ep.spa(i5,i4) +
          ep.spb(i1,i6)*ep.spb(i3,i2)*ep.spa(i6,i3) + ep.spb(i1,i6)*ep.spb(i4,i2)*ep.spa(i6,i4) +
          ep.spb(i1,i7)*ep.spb(i3,i2)*ep.spa(i7,i3) + ep.spb(i1,i7)*ep.spb(i4,i2)*ep.spa(i7,i4))*
        (ep.spb(i1,i5)*ep.spb(i3,i3)*ep.spa(i5,i3) + ep.spb(i1,i5)*ep.spb(i4,i3)*ep.spa(i5,i4) +
          ep.spb(i1,i6)*ep.spb(i3,i3)*ep.spa(i6,i3) + ep.spb(i1,i6)*ep.spb(i4,i3)*ep.spa(i6,i4) +
          ep.spb(i1,i7)*ep.spb(i3,i3)*ep.spa(i7,i3) + ep.spb(i1,i7)*ep.spb(i4,i3)*ep.spa(i7,i4))*
        ep.s(i3,i4)) + (ep.spb(i3,i4)*ep.spb(i5,i6)*
        pow(ep.spb(i1,i4)*ep.spb(i6,i1)*ep.spa(i4,i6) + ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7) +
          ep.spb(i1,i5)*ep.spb(i6,i1)*ep.spa(i5,i6) + ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i6,i1)*ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i5,i6) +
          ep.spb(i1,i3)*pow(ep.spb(i4,i5),2)*ep.spb(i6,i1)*ep.spa(i3,i4)*ep.spa(i5,i4)*ep.spa(i5,i6) +
          ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i5,i7) +
          ep.spb(i1,i3)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*ep.spa(i3,i4)*ep.spa(i5,i4)*ep.spa(i5,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spa(i2,i4) + ep.spb(i1,i2)*ep.spb(i5,i5)*ep.spa(i2,i5) +
          ep.spb(i1,i3)*ep.spb(i4,i5)*ep.spa(i3,i4) + ep.spb(i1,i3)*ep.spb(i5,i5)*ep.spa(i3,i5))*
        (ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spa(i2,i4) + ep.spb(i1,i2)*ep.spb(i5,i6)*ep.spa(i2,i5) +
          ep.spb(i1,i3)*ep.spb(i4,i6)*ep.spa(i3,i4) + ep.spb(i1,i3)*ep.spb(i5,i6)*ep.spa(i3,i5))*
        (ep.spb(i1,i6)*ep.spb(i4,i3)*ep.spa(i6,i4) + ep.spb(i1,i6)*ep.spb(i5,i3)*ep.spa(i6,i5) +
          ep.spb(i1,i7)*ep.spb(i4,i3)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i3)*ep.spa(i7,i5))*
        (ep.spb(i1,i6)*ep.spb(i4,i4)*ep.spa(i6,i4) + ep.spb(i1,i6)*ep.spb(i5,i4)*ep.spa(i6,i5) +
          ep.spb(i1,i7)*ep.spb(i4,i4)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i4)*ep.spa(i7,i5))*
        ep.s(i4,i5)) + (ep.spb(i4,i5)*ep.spb(i6,i7)*
        pow(ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7) + ep.spb(i1,i6)*ep.spb(i7,i1)*ep.spa(i6,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i5)*ep.spa(i6,i5)*ep.spa(i6,i7) +
          ep.spb(i1,i3)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i3,i5)*ep.spa(i6,i5)*ep.spa(i6,i7) +
          ep.spb(i1,i4)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i4,i5)*ep.spa(i6,i5)*ep.spa(i6,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i5,i6)*ep.spa(i2,i5) + ep.spb(i1,i2)*ep.spb(i6,i6)*ep.spa(i2,i6) +
          ep.spb(i1,i3)*ep.spb(i5,i6)*ep.spa(i3,i5) + ep.spb(i1,i3)*ep.spb(i6,i6)*ep.spa(i3,i6) +
          ep.spb(i1,i4)*ep.spb(i5,i6)*ep.spa(i4,i5) + ep.spb(i1,i4)*ep.spb(i6,i6)*ep.spa(i4,i6))*
        (ep.spb(i1,i2)*ep.spb(i5,i7)*ep.spa(i2,i5) + ep.spb(i1,i2)*ep.spb(i6,i7)*ep.spa(i2,i6) +
          ep.spb(i1,i3)*ep.spb(i5,i7)*ep.spa(i3,i5) + ep.spb(i1,i3)*ep.spb(i6,i7)*ep.spa(i3,i6) +
          ep.spb(i1,i4)*ep.spb(i5,i7)*ep.spa(i4,i5) + ep.spb(i1,i4)*ep.spb(i6,i7)*ep.spa(i4,i6))*
        (ep.spb(i1,i7)*ep.spb(i5,i4)*ep.spa(i7,i5) + ep.spb(i1,i7)*ep.spb(i6,i4)*ep.spa(i7,i6))*
        (ep.spb(i1,i7)*ep.spb(i5,i5)*ep.spa(i7,i5) + ep.spb(i1,i7)*ep.spb(i6,i5)*ep.spa(i7,i6))*
        ep.s(i5,i6)) + (ep.spb(i2,i3)*ep.spb(i5,i6)*
        pow(ep.spb(i1,i3)*ep.spb(i6,i1)*ep.spa(i3,i6) + ep.spb(i1,i3)*ep.spb(i7,i1)*ep.spa(i3,i7) +
          ep.spb(i1,i4)*ep.spb(i6,i1)*ep.spa(i4,i6) + ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7) +
          ep.spb(i1,i5)*ep.spb(i6,i1)*ep.spa(i5,i6) + ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i6,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i6) +
          ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i6,i1)*ep.spa(i2,i3)*ep.spa(i4,i6)*
           ep.spa(i5,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i6,i1)*ep.spa(i2,i4)*
           ep.spa(i4,i6)*ep.spa(i5,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i7,i1)*
           ep.spa(i2,i3)*ep.spa(i4,i7)*ep.spa(i5,i3) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i4,i7)*
           ep.spa(i5,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i6,i1)*ep.spa(i2,i3)*
           ep.spa(i4,i3)*ep.spa(i5,i6) + ep.spb(i1,i2)*pow(ep.spb(i3,i5),2)*ep.spb(i6,i1)*
           ep.spa(i2,i3)*ep.spa(i5,i3)*ep.spa(i5,i6) +
          ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i6,i1)*ep.spa(i2,i4)*ep.spa(i5,i3)*
           ep.spa(i5,i6) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i6,i1)*ep.spa(i2,i3)*
           ep.spa(i5,i4)*ep.spa(i5,i6) + ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i6,i1)*
           ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i5,i6) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*
           ep.spa(i5,i7) + ep.spb(i1,i2)*pow(ep.spb(i3,i5),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*
           ep.spa(i5,i3)*ep.spa(i5,i7) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i7,i1)*
           ep.spa(i2,i4)*ep.spa(i5,i3)*ep.spa(i5,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i5,i4)*
           ep.spa(i5,i7) + ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i5,i4)*ep.spa(i5,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spa(i2,i4) +
          ep.spb(i1,i2)*ep.spb(i5,i5)*ep.spa(i2,i5))*
        (ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spa(i2,i4) +
          ep.spb(i1,i2)*ep.spb(i5,i6)*ep.spa(i2,i5))*
        (ep.spb(i1,i6)*ep.spb(i3,i2)*ep.spa(i6,i3) + ep.spb(i1,i6)*ep.spb(i4,i2)*ep.spa(i6,i4) +
          ep.spb(i1,i6)*ep.spb(i5,i2)*ep.spa(i6,i5) + ep.spb(i1,i7)*ep.spb(i3,i2)*ep.spa(i7,i3) +
          ep.spb(i1,i7)*ep.spb(i4,i2)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i2)*ep.spa(i7,i5))*
        (ep.spb(i1,i6)*ep.spb(i3,i3)*ep.spa(i6,i3) + ep.spb(i1,i6)*ep.spb(i4,i3)*ep.spa(i6,i4) +
          ep.spb(i1,i6)*ep.spb(i5,i3)*ep.spa(i6,i5) + ep.spb(i1,i7)*ep.spb(i3,i3)*ep.spa(i7,i3) +
          ep.spb(i1,i7)*ep.spb(i4,i3)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i3)*ep.spa(i7,i5))*
        ep.s(i3,i4,i5)) + (ep.spb(i3,i4)*ep.spb(i6,i7)*
        pow(ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7) + ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7) +
          ep.spb(i1,i6)*ep.spb(i7,i1)*ep.spa(i6,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i5,i7) +
          ep.spb(i1,i3)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*ep.spa(i3,i4)*ep.spa(i5,i4)*ep.spa(i5,i7) +
          ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i5,i7)*
           ep.spa(i6,i4) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i5)*
           ep.spa(i5,i7)*ep.spa(i6,i4) + ep.spb(i1,i3)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*
           ep.spa(i3,i4)*ep.spa(i5,i7)*ep.spa(i6,i4) +
          ep.spb(i1,i3)*ep.spb(i4,i5)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i3,i5)*ep.spa(i5,i7)*
           ep.spa(i6,i4) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i5,i4)*ep.spa(i6,i7) + ep.spb(i1,i3)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*
           ep.spa(i3,i4)*ep.spa(i5,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*pow(ep.spb(i4,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i6,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i5)*ep.spa(i6,i4)*
           ep.spa(i6,i7) + ep.spb(i1,i3)*pow(ep.spb(i4,i6),2)*ep.spb(i7,i1)*ep.spa(i3,i4)*
           ep.spa(i6,i4)*ep.spa(i6,i7) + ep.spb(i1,i3)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i3,i5)*ep.spa(i6,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i6,i5)*
           ep.spa(i6,i7) + ep.spb(i1,i2)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i5)*
           ep.spa(i6,i5)*ep.spa(i6,i7) + ep.spb(i1,i3)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i3,i4)*ep.spa(i6,i5)*ep.spa(i6,i7) +
          ep.spb(i1,i3)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i3,i5)*ep.spa(i6,i5)*ep.spa(i6,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spa(i2,i4) + ep.spb(i1,i2)*ep.spb(i5,i6)*ep.spa(i2,i5) +
          ep.spb(i1,i2)*ep.spb(i6,i6)*ep.spa(i2,i6) + ep.spb(i1,i3)*ep.spb(i4,i6)*ep.spa(i3,i4) +
          ep.spb(i1,i3)*ep.spb(i5,i6)*ep.spa(i3,i5) + ep.spb(i1,i3)*ep.spb(i6,i6)*ep.spa(i3,i6))*
        (ep.spb(i1,i2)*ep.spb(i4,i7)*ep.spa(i2,i4) + ep.spb(i1,i2)*ep.spb(i5,i7)*ep.spa(i2,i5) +
          ep.spb(i1,i2)*ep.spb(i6,i7)*ep.spa(i2,i6) + ep.spb(i1,i3)*ep.spb(i4,i7)*ep.spa(i3,i4) +
          ep.spb(i1,i3)*ep.spb(i5,i7)*ep.spa(i3,i5) + ep.spb(i1,i3)*ep.spb(i6,i7)*ep.spa(i3,i6))*
        (ep.spb(i1,i7)*ep.spb(i4,i3)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i3)*ep.spa(i7,i5) +
          ep.spb(i1,i7)*ep.spb(i6,i3)*ep.spa(i7,i6))*
        (ep.spb(i1,i7)*ep.spb(i4,i4)*ep.spa(i7,i4) + ep.spb(i1,i7)*ep.spb(i5,i4)*ep.spa(i7,i5) +
          ep.spb(i1,i7)*ep.spb(i6,i4)*ep.spa(i7,i6))*ep.s(i4,i5,i6)) +
     (ep.spb(i2,i3)*ep.spb(i6,i7)*pow(ep.spb(i1,i3)*ep.spb(i7,i1)*ep.spa(i3,i7) +
          ep.spb(i1,i4)*ep.spb(i7,i1)*ep.spa(i4,i7) + ep.spb(i1,i5)*ep.spb(i7,i1)*ep.spa(i5,i7) +
          ep.spb(i1,i6)*ep.spb(i7,i1)*ep.spa(i6,i7),3)*
        (ep.spb(i1,i2)*pow(ep.spb(i3,i4),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i4,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i7)*
           ep.spa(i5,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i4,i7)*ep.spa(i5,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i5)*ep.spb(i7,i1)*
           ep.spa(i2,i3)*ep.spa(i4,i3)*ep.spa(i5,i7) +
          ep.spb(i1,i2)*pow(ep.spb(i3,i5),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i5,i3)*ep.spa(i5,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i5,i3)*
           ep.spa(i5,i7) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i3)*
           ep.spa(i5,i4)*ep.spa(i5,i7) + ep.spb(i1,i2)*pow(ep.spb(i4,i5),2)*ep.spb(i7,i1)*
           ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i5,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i7)*
           ep.spa(i6,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i4,i7)*ep.spa(i6,i3) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i5)*ep.spa(i4,i7)*ep.spa(i6,i3) +
          ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i3,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i5,i7)*
           ep.spa(i6,i3) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i5,i7)*ep.spa(i6,i3) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i5)*ep.spa(i5,i7)*ep.spa(i6,i3) +
          ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i5,i7)*
           ep.spa(i6,i4) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*
           ep.spa(i5,i7)*ep.spa(i6,i4) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i5)*ep.spa(i5,i7)*ep.spa(i6,i4) +
          ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i3,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i4,i3)*
           ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i4)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*
           ep.spa(i4,i5)*ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i3,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i3)*ep.spa(i5,i3)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i4,i5)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i5,i3)*
           ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*
           ep.spa(i5,i4)*ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i4,i5)*ep.spb(i4,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i4)*ep.spa(i5,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*pow(ep.spb(i3,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i3)*ep.spa(i6,i3)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i4,i6)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i6,i3)*
           ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i5)*
           ep.spa(i6,i3)*ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i4,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i3)*ep.spa(i6,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*pow(ep.spb(i4,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i4)*ep.spa(i6,i4)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i5)*ep.spa(i6,i4)*
           ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*ep.spa(i2,i3)*
           ep.spa(i6,i5)*ep.spa(i6,i7) + ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spb(i5,i6)*ep.spb(i7,i1)*
           ep.spa(i2,i4)*ep.spa(i6,i5)*ep.spa(i6,i7) +
          ep.spb(i1,i2)*pow(ep.spb(i5,i6),2)*ep.spb(i7,i1)*ep.spa(i2,i5)*ep.spa(i6,i5)*ep.spa(i6,i7)))/
      ((ep.spb(i1,i2)*ep.spb(i3,i6)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i6)*ep.spa(i2,i4) +
          ep.spb(i1,i2)*ep.spb(i5,i6)*ep.spa(i2,i5) + ep.spb(i1,i2)*ep.spb(i6,i6)*ep.spa(i2,i6))*
        (ep.spb(i1,i2)*ep.spb(i3,i7)*ep.spa(i2,i3) + ep.spb(i1,i2)*ep.spb(i4,i7)*ep.spa(i2,i4) +
          ep.spb(i1,i2)*ep.spb(i5,i7)*ep.spa(i2,i5) + ep.spb(i1,i2)*ep.spb(i6,i7)*ep.spa(i2,i6))*
        (ep.spb(i1,i7)*ep.spb(i3,i2)*ep.spa(i7,i3) + ep.spb(i1,i7)*ep.spb(i4,i2)*ep.spa(i7,i4) +
          ep.spb(i1,i7)*ep.spb(i5,i2)*ep.spa(i7,i5) + ep.spb(i1,i7)*ep.spb(i6,i2)*ep.spa(i7,i6))*
        (ep.spb(i1,i7)*ep.spb(i3,i3)*ep.spa(i7,i3) + ep.spb(i1,i7)*ep.spb(i4,i3)*ep.spa(i7,i4) +
          ep.spb(i1,i7)*ep.spb(i5,i3)*ep.spa(i7,i5) + ep.spb(i1,i7)*ep.spb(i6,i3)*ep.spa(i7,i6))*
        ep.s(i3,i4,i5,i6)))/
   (ep.spb(i1,i2)*ep.spb(i2,i3)*ep.spb(i3,i4)*ep.spb(i4,i5)*ep.spb(i5,i6)*ep.spb(i6,i7)*
     ep.spb(i7,i1)));
}



template <class T> complex<T> R7g126(const eval_param<T>& ep,const mass_param_coll& mpc){return R7g1m<0,1,2,3,4,5,6>(ep,mpc);}  // 126
template <class T> complex<T> R7g1(const eval_param<T>& ep,const mass_param_coll& mpc){return R7g1p<0,1,2,3,4,5,6>(ep,mpc);}   //1
template <class T> complex<T> R7g0(const eval_param<T>& ep,const mass_param_coll& mpc){return Rallm(ep,mpc);}   //0
template <class T> complex<T> R7g127(const eval_param<T>& ep,const mass_param_coll& mpc){return Rallp(ep,mpc);}   //127


//#endif

template <class T> complex<T> (*R7g_Ptr_eval(int hc))(const eval_param<T>& ,const mass_param_coll&) {
	switch (hc) {
//#if _FAST_COMPILE == 0
	case 0  : return &R7g0;
	case 1  : return &R7g1;
	case 126: return &R7g126;
	case 127: return &R7g127;
//#endif
	default : return 0;
	}
}


template complex<R> (*R7g_Ptr_eval(int hc))(const eval_param<R>&,const mass_param_coll&) ;
 template complex<RHP> (*R7g_Ptr_eval(int hc))(const eval_param<RHP>&,const mass_param_coll&) ;
 template complex<RVHP> (*R7g_Ptr_eval(int hc))(const eval_param<RVHP>&,const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*R7g_Ptr_eval(int hc))(const eval_param<RGMP>&,const mass_param_coll&) ;
#endif

}

