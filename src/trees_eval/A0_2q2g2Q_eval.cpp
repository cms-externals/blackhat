/* The 2q 2Q amplitudes use a base 6 code to enumerate them
   m,qm,qp,p,Qm,Qp corresponds to 0,1,2,3,5,6.
   ep and em are left out of the code.
   For example, (qp,qm,m,p,Qm,Qp) = 2+ 1*6 + 0*36 + 3*216 + 4*1296 + 5*7776
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"


using namespace std;

#define _FAST_COMPILE_eval 0

namespace BH {

#if _FAST_COMPILE_eval==1

template <class T> complex<T> A2q2g2Qzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q2g2Q317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q1052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q1162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q1472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q1502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(-(ep.s(3,4,5)*ep.spa(1,
   5))-ep.spa(0,1)*ep.spab(5,
  ep.Sum(1,2),0))*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,1),2)*ep.spa(3,5)*
 ep.spb(2,0))/(pow(ep.spab(5,ep.Sum(3,4),
 2),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(0,1),
2)*(-(ep.s(3,4,5)*ep.spa(1,
   5))+ep.spa(1,2)*ep.spab(5,
  ep.Sum(0,1),2))*
 ((-ep.s(0,1)+ep.s(3,4,5))*
 ep.spab(5,ep.Sum(0,1),2)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),2)))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(5,
 ep.Sum(3,4),2),3)*pow(ep.spb(3,2),
2)*ep.spa(1,5)*ep.spb(4,2))/
(pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(3,4),2)*
 (ep.s(0,1,5)*ep.spab(5,ep.Sum(0,
   1),2)+(-ep.s(3,4)+
  ep.s(0,1,5))*ep.spab(5,
  ep.Sum(3,4),2))*
 (-(ep.s(0,1,5)*ep.spa(3,5))-
ep.spa(2,3)*ep.spab(5,ep.Sum(3,4),
  2)))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spa(2,5)*ep.spab(1,ep.Sum(5,2),
4)*ep.spab(5,ep.Sum(1,2),0)*
 ep.spb(4,3)*(-(ep.spa(1,2)*
  ep.spb(2,0))-ep.spa(1,5)*
 ep.spb(5,0)))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spab(5,ep.Sum(3,2),1),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spa(2,5)*ep.spab(3,ep.Sum(5,2),
4)*ep.spab(5,ep.Sum(3,2),4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
  ep.spb(2,0))+ep.spa(3,5)*
 ep.spb(5,0)))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),2),2)*
 pow(ep.spb(2,0),2)*ep.spa(1,3))/
(ep.spa(3,4)*ep.spab(1,ep.Sum(5,0),
2)*ep.spb(5,0)*
 (-(ep.spa(1,2)*ep.spb(2,0))-
ep.spa(1,5)*ep.spb(5,0))*
 (-(ep.spa(2,3)*ep.spb(2,0))+
ep.spa(3,5)*ep.spb(5,0))*
 ep.spb(5,2))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 pow(ep.spb(4,2),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spab(1,ep.Sum(5,2),
4)*ep.spab(3,ep.Sum(2,5),4)*
 ep.spab(3,ep.Sum(4,5),2)*
 ep.spb(5,2)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(1,0),2)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0)))-(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),
 3),2)*ep.spab(5,ep.Sum(0,1),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 pow(ep.spb(1,0),2)*ep.spa(3,5)*
 ep.spb(2,0))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2)*
 ep.spa(1,5)*ep.spab(5,ep.Sum(0,1),
4))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spab(1,ep.Sum(4,5),0)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1912_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.spab(2,ep.Sum(0,1),4))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q1935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 pow(ep.spb(2,0),3)*ep.spa(3,5))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),3)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(0,1),2))-
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),2)*
 ep.spab(5,ep.Sum(2,1),0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spab(2,ep.Sum(4,3),0))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(1,0),
2)*ep.spab(5,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
2)*ep.spab(5,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),2),2)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 (-(ep.spa(0,2)*ep.spb(4,0))-
ep.spa(2,5)*ep.spb(5,4)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2))/
(pow(ep.spa(0,1),2)*ep.s(2,3,4)*
ep.spb(3,2)*ep.spb(4,0)+
 ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(1,5)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),3),2)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 (-(ep.spa(0,2)*ep.spb(4,0))-
ep.spa(2,5)*ep.spb(5,4)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2))/
(pow(ep.spa(0,1),2)*ep.s(2,3,4)*
ep.spb(3,2)*ep.spb(4,0)+
 ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(1,5)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2662_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2747_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2767_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q2797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),3))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),3)*ep.spa(2,3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (-(ep.s(0,1,5)*ep.spa(3,5))-
ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
  4))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),3))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2)*ep.spb(2,0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(2,3))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*(-(ep.s(0,1,5)*ep.spa(3,
   5))-ep.spa(3,4)*ep.spab(5,
  ep.Sum(0,1),4))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(4,3),2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),3))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 pow(ep.spb(1,0),2)*ep.spb(2,0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(0,1,
   2))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(0,1,2)*ep.spab(5,ep.Sum(3,
   4),0)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),3)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(1,
ep.Sum(4,5),0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q2992_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(4,3),2))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2)*
 ep.spab(5,ep.Sum(0,1),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 pow(ep.spb(2,0),3))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(0,1,
   2))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(0,1,2)*ep.spab(5,ep.Sum(3,
   4),0)))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(1,
ep.Sum(4,5),0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,3),2)*ep.spa(2,3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (-(ep.s(0,1,5)*ep.spa(3,5))-
ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
  4)))-(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),
 0),3))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),3))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2)*ep.spa(1,5)*
 ep.spa(2,3))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*(-(ep.s(0,1,5)*ep.spa(3,
   5))-ep.spa(3,4)*ep.spab(5,
  ep.Sum(0,1),4)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),1),2)*
 ep.spab(5,ep.Sum(4,3),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),2),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3467_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3697_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q3805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q3835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,2),
3))/(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(3,4),0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(3,
ep.Sum(5,0),4))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,2),
3)*ep.spa(1,5))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2)*ep.spab(5,
ep.Sum(3,4),0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),3)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),0),3)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),2)*pow(ep.spb(4,3),
2)*ep.spab(5,ep.Sum(3,2),4)*
 ep.spb(4,2))/(pow(ep.spab(5,ep.Sum(0,1),
 4),2)*ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(3,2),
4)*(ep.s(0,1,5)*ep.spab(5,
  ep.Sum(0,1),4)+
(-ep.s(2,3)+ep.s(0,1,5))*
 ep.spab(5,ep.Sum(2,3),4))*
 (-(ep.s(0,1,5)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),3)*ep.spab(3,
ep.Sum(5,0),4))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5)*ep.spab(5,ep.Sum(3,4),
0))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),2)*pow(ep.spb(4,3),
2)*ep.spa(1,5)*ep.spab(5,ep.Sum(3,2),
4)*ep.spb(4,2))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(3,2),4)*
 (ep.s(0,1,5)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(0,1,5))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(0,1,5)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),
 4),2)*ep.spab(3,ep.Sum(5,0),4)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q4940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,2),
3))/(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(3,4),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(3,
ep.Sum(5,0),4))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q4945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,2),
3)*ep.spa(1,5))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2)*ep.spab(5,
ep.Sum(3,4),0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),3)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),0),3)*
 ep.spa(3,5))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,3),
2)*ep.spb(4,2))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),3)*ep.spab(3,
ep.Sum(5,0),4))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5)*ep.spab(5,ep.Sum(3,4),
0))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(5,
 ep.Sum(2,3),4),3)*pow(ep.spb(4,3),
2)*ep.spa(1,5)*ep.spb(4,2))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),
 4),2)*ep.spab(3,ep.Sum(5,0),4)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.spab(2,ep.Sum(0,5),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5647_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spab(2,ep.Sum(0,5),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(1,0),
2)*ep.spab(5,ep.Sum(1,2),0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q5915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
2)*ep.spab(5,ep.Sum(2,1),0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q5935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(0,1),
2)*(-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q6025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spab(2,ep.Sum(4,3),0))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6277_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6302_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(5,
 ep.Sum(1,2),0),3)*pow(ep.spb(1,0),
2)*ep.spa(3,5)*ep.spb(2,0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (-(ep.s(3,4,5)*ep.spa(1,5))-
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),3),2)*
 ep.spa(1,5)*ep.spab(5,ep.Sum(1,0),
4))/(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spab(1,ep.Sum(4,5),0)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(1,0),2)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0)))+(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),
 3),2)*ep.spab(5,ep.Sum(0,1),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q6337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(3,5)*
 ep.spab(5,ep.Sum(2,1),0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 (ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),0),3)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),3),2)*
 ep.spab(2,ep.Sum(0,1),4))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(4,3)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q6385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 pow(ep.spb(2,1),2)*ep.spa(3,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(2,0))/(pow(ep.spab(5,ep.Sum(3,4),
 2),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,0),
2)*(-(ep.s(3,4,5)*ep.spa(1,
   5))+ep.spa(1,2)*ep.spab(5,
  ep.Sum(0,1),2))*
 ((-ep.s(0,1)+ep.s(3,4,5))*
 ep.spab(5,ep.Sum(0,1),2)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),2)))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(5,
 ep.Sum(3,4),2),3)*pow(ep.spb(3,2),
2)*ep.spa(1,5)*ep.spb(4,2))/
(pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(3,4),2)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),2)+(-ep.s(3,4)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(3,4),2))*
 (ep.s(2,3,4)*ep.spa(3,5)+
ep.spa(2,3)*ep.spab(5,ep.Sum(3,4),
  2)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 2),2)*pow(ep.spb(2,0),2)*
 ep.spa(1,3))/(ep.spa(3,4)*
 ep.spab(1,ep.Sum(2,5),0)*
 ep.spab(1,ep.Sum(5,0),2)*
 ep.spab(3,ep.Sum(5,2),0)*
 ep.spb(5,0)*ep.spb(5,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spa(1,5)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,5)*
 ep.spab(1,ep.Sum(5,2),0)*
 ep.spab(5,ep.Sum(1,2),0)*
 ep.spb(4,3)*((ep.spa(1,2)*
  ep.spb(4,2))/ep.spa(1,5)-
ep.spb(5,4)))+(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spa(2,5)*ep.spab(3,ep.Sum(5,2),
0)*ep.spab(5,ep.Sum(2,3),4)*
 ep.spb(1,0)*(-(ep.spa(2,3)*
  ep.spb(4,2))-ep.spa(3,5)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(1,3)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spab(3,ep.Sum(4,5),
2)*ep.spb(5,2)*ep.spb(5,4)*
 (-(ep.spa(2,3)*ep.spb(4,2))-
ep.spa(3,5)*ep.spb(5,4))*
 (ep.spa(1,2)-(ep.spa(1,5)*
  ep.spb(5,4))/ep.spb(4,2)))
); }

template <class T> complex<T> A2q2g2Q6532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6562_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6772_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),
2)*ep.spab(5,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(0,1),4)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(3,1),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),
2))/(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q6952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,3),
2)*ep.spab(5,ep.Sum(3,2),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(5,ep.Sum(0,1),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q6975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.spab(5,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,3),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(0,1),
2)*(-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q7285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spab(5,ep.Sum(2,3),
4))/(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.spa(3,4)-(ep.s(2,3,4)*
  ep.spa(3,5))/ep.spab(5,ep.Sum(2,3),
  4)))-(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),
 1),2)*ep.spab(5,ep.Sum(3,4),0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(3,5),3)*pow(ep.spab(5,
 ep.Sum(1,2),0),3)*pow(ep.spb(1,0),
2)*ep.spb(2,0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (-(ep.s(3,4,5)*ep.spa(1,5))-
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),4),3)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(1,
ep.Sum(4,5),0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(1,0),2)*ep.spa(1,2))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0)))+(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),
 4),3))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(4,3)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q7417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spab(5,
ep.Sum(2,1),0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 (ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3)*ep.spab(1,
ep.Sum(4,5),0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2)*
 ep.spab(5,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(1,2))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 (ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))+(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(4,0),3))/(ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q7465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2)*ep.spb(2,0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(5,ep.Sum(2,3),4),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(2,3))/(pow(ep.spab(5,ep.Sum(0,1),
 4),2)*ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*(ep.s(2,3,4)*ep.spa(3,5)-
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7882_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q7972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q7982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(3,2),0),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(3,2),1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,5),3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(4,5),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,1),3)*ep.spa(0,4))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(1,0),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),3)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,2),2)*ep.spa(0,4)*
 ep.spb(3,1))/(pow(ep.spab(4,ep.Sum(5,0),
 3),2)*ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.s(1,2,3)*ep.spa(2,
   4))+ep.spa(2,3)*ep.spab(4,
  ep.Sum(1,2),3))*
 ((-ep.s(1,2)+ep.s(1,2,3))*
 ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4)*ep.spab(4,ep.Sum(2,3),
5))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(2,ep.Sum(4,5),3)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spab(1,ep.Sum(3,2),5))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),2),2)*
 ep.spab(1,ep.Sum(5,0),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(4,ep.Sum(2,3),1),3)*
 pow(ep.spb(2,1),2)*ep.spa(0,4)*
 ep.spb(3,1))/(pow(ep.spab(4,ep.Sum(5,0),
 1),2)*ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(2,3),
1)*(ep.s(1,2,3)*ep.spa(2,4)+
ep.spa(1,2)*ep.spab(4,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(2,3),1)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),1)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(0,2)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spab(2,ep.Sum(3,4),
1)*ep.spb(4,1)*ep.spb(4,3)*
 (-(ep.spa(1,2)*ep.spb(3,1))-
ep.spa(2,4)*ep.spb(4,3))*
 (ep.spa(0,1)-(ep.spa(0,4)*
  ep.spb(4,3))/ep.spb(3,1)))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(5,0),1),2)*
 pow(ep.spb(1,0),2)*ep.spa(2,4)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(5,1))/(pow(ep.spab(4,ep.Sum(2,3),
 1),2)*ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,5),
1)*(ep.s(2,3,4)*ep.spab(4,
  ep.Sum(2,3),1)+
(-ep.s(0,5)+ep.s(2,3,4))*
 ep.spab(4,ep.Sum(5,0),1))*
 (-(ep.s(2,3,4)*ep.spa(0,4))+
ep.spa(0,1)*ep.spab(4,ep.Sum(5,0),
  1)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),
 2),2)*ep.spa(0,4)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,4)*
 ep.spab(0,ep.Sum(4,1),5)*
 ep.spab(4,ep.Sum(0,1),5)*
 ep.spb(3,2)*((ep.spa(0,1)*
  ep.spb(3,1))/ep.spa(0,4)-
ep.spb(4,3)))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spa(1,4)*ep.spab(2,ep.Sum(4,1),
5)*ep.spab(4,ep.Sum(1,2),3)*
 (-(ep.spa(1,2)*ep.spb(3,1))-
ep.spa(2,4)*ep.spb(4,3))*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 pow(ep.spb(5,1),2)*ep.spa(0,2))/
(ep.spa(2,3)*ep.spab(0,ep.Sum(1,4),
5)*ep.spab(0,ep.Sum(4,5),1)*
 ep.spab(2,ep.Sum(4,1),5)*
 ep.spb(4,1)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8872_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8882_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q8932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q8990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spa(1,2)*
 ep.spab(4,ep.Sum(1,2),3))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (ep.spa(2,3)-(ep.s(1,2,3)*
  ep.spa(2,4))/ep.spab(4,ep.Sum(1,2),
  3)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),
 0),2)*ep.spab(4,ep.Sum(2,3),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q9002_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q9012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q9015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 ep.spab(4,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q9020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(4,ep.Sum(1,2),3),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(1,2))/(pow(ep.spab(4,ep.Sum(5,0),
 3),2)*ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*(ep.s(1,2,3)*ep.spa(2,4)-
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,0),3))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,1),2)*ep.spb(2,0))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),4),3)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10517_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*(ep.spa(3,5)*ep.spb(3,2)+
ep.spa(4,5)*ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2)*ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(4,0),3)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(-(ep.s(0,4,5)*ep.spa(3,
   5))+ep.spa(4,5)*ep.spab(3,
  ep.Sum(5,0),4))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10552_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*(ep.spa(3,5)*ep.spb(3,2)+
ep.spa(4,5)*ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q10575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(4,0),3)*ep.spa(0,5))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (-(ep.s(0,4,5)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,3),0),2)*
 ep.spa(2,4)*ep.spab(2,ep.Sum(4,3),
1))/(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(4,ep.Sum(1,2),3)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(4,3),2)*ep.spa(0,2)*
 ep.spb(5,3))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(ep.s(0,1,2)*ep.spab(2,
  ep.Sum(0,1),3)+
(-ep.s(4,5)+ep.s(0,1,2))*
 ep.spab(2,ep.Sum(4,5),3))*
 (ep.s(0,1,2)*ep.spa(2,4)-
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A2q2g2Q11050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,3),1),3)*
 ep.spa(2,4))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),3)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(4,5),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(4,3),2)*ep.spb(5,3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (ep.s(0,1,2)*ep.spa(2,4)-
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A2q2g2Q11153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4)*
 ep.spa(1,5))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(1,3)*
 ep.spa(2,5)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(1,4)*ep.spa(3,5))/
(ep.spa(0,1)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spa(2,5)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(5,3),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),0)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(0,4,5)*
 ep.spa(4,5)*(ep.s(0,4,5)*
 ep.spa(0,2)-ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),1),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.spab(4,ep.Sum(2,1),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11218_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,2),2)*pow(ep.spab(0,
 ep.Sum(1,2),3),3)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (ep.s(0,4,5)*ep.spa(0,2)+
ep.spa(2,3)*ep.spab(0,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(1,2),3)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),3))*ep.spab(4,ep.Sum(1,2),
3))-(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),
 3),2)*pow(ep.spb(3,1),2)*
 ep.spa(2,4))/(ep.spa(4,5)*
 ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(1,0)*ep.spb(3,0)*
 (-(ep.spa(0,2)*ep.spb(1,0))-
ep.spa(2,3)*ep.spb(3,1))*
 (ep.spa(0,4)*ep.spb(1,0)-
ep.spa(3,4)*ep.spb(3,1)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 pow(ep.spb(5,3),2)*ep.spa(2,4))/
(ep.spa(1,2)*ep.spab(2,ep.Sum(0,3),
5)*ep.spab(4,ep.Sum(3,0),5)*
 ep.spab(4,ep.Sum(5,0),3)*
 ep.spb(3,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spab(0,ep.Sum(4,3),2),2)*
 ep.spb(5,1))/(ep.spa(0,3)*
 ep.spa(3,4)*ep.spab(0,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(0,3),5)*
 ep.spb(2,1)*(ep.spa(0,4)*
 ep.spb(1,0)-ep.spa(3,4)*
 ep.spb(3,1)))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(0,
 ep.Sum(4,5),3),3)*pow(ep.spb(4,3),
2)*ep.spa(0,2)*ep.spb(5,3))/
(pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.s(0,1,2)*ep.spa(1,2)*
 (ep.s(0,1,2)*ep.spab(0,ep.Sum(1,
   2),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(0,
  ep.Sum(4,5),3))*
 (ep.s(0,1,2)*ep.spa(0,4)-
ep.spa(3,4)*ep.spab(0,ep.Sum(4,5),
  3))*ep.spab(2,ep.Sum(4,5),3))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(5,1))/(ep.spa(0,3)*
 ep.spa(2,3)*ep.spab(0,ep.Sum(2,3),
1)*ep.spab(2,ep.Sum(0,3),5)*
 (-(ep.spa(0,2)*ep.spb(1,0))-
ep.spa(2,3)*ep.spb(3,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11223_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11363_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(2,1),2)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(0,4,5)*
 ep.spa(4,5)*(-(ep.s(0,4,5)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(4,0),
3))/(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3)*
 ep.spab(5,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q11430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),3)*pow(ep.spb(2,1),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (-(ep.s(0,4,5)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),
 1),2)*ep.spab(2,ep.Sum(5,0),1)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spa(0,2)*ep.spab(0,ep.Sum(1,2),
5))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,3),2)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11578_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.spab(3,ep.Sum(1,2),5))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(4,0),
3)*ep.spa(1,3))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q11610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),3)*pow(ep.spb(3,1),
3)*ep.spa(0,4))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (-(ep.s(0,4,5)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),
 1),3)*ep.spb(5,1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q11613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q11765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q11770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q11825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q11840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(3,
 ep.Sum(5,0),4),2)*pow(ep.spb(4,0),
3)*ep.spa(1,3)*ep.spab(3,ep.Sum(0,5),
4))/(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(1,2,3)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),3)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q11860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q11870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,3),3)*pow(ep.spab(3,
 ep.Sum(5,0),4),2)*pow(ep.spb(4,0),
3)*ep.spab(3,ep.Sum(0,5),4))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(1,2,3)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q12125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(4,3),2)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(-(ep.s(3,4,5)*ep.spa(2,
   4))+ep.spa(3,4)*ep.spab(2,
  ep.Sum(4,5),3)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(1,0)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q12130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(4,3),2)*ep.spa(4,5))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (-(ep.s(3,4,5)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))+(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),
 1),3))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(1,0)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1)))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),3))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q12233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q12263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),
2))/(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(0,1),5))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),2)*
 ep.spab(0,ep.Sum(3,2),1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q12298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spab(3,ep.Sum(5,4),1))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12418_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12473_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12488_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12508_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(1,0),2)*
 ep.spab(2,ep.Sum(0,5),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),1)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12653_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12668_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),3))/(ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(4,0),3)*ep.spa(4,5))/
(ep.s(1,2,3)*ep.spa(2,3)*
 (ep.s(1,2,3)*ep.spa(1,5)+
ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
  0))*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(2,1),
2)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),1)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(5,ep.Sum(3,2),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12878_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),3)*ep.spa(1,3)*
 ep.spa(4,5))/(ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spa(1,5)+ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),0))*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
2)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),1)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q12903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12908_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q12977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q12992_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13037_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13152_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q13160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13572_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13703_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13733_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q13748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q13790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),
2))/(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2))-
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),3))/
(ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),3))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q13803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q13808_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q14042_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(1,0),2)*ep.spa(0,5))/
(ep.s(2,3,4)*ep.spa(3,4)*
 (-(ep.s(2,3,4)*ep.spa(0,2))-
ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),3))/
(ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q14060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),
3)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),3),3)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(2,
 ep.Sum(5,0),1),2)*pow(ep.spb(1,0),
2)*ep.spab(2,ep.Sum(0,5),1)*
 ep.spb(5,1))/(pow(ep.spab(2,ep.Sum(3,4),
 1),2)*ep.s(2,3,4)*
 ep.spa(3,4)*(ep.s(2,3,4)*
 ep.spab(2,ep.Sum(3,4),1)+
(-ep.s(0,5)+ep.s(2,3,4))*
 ep.spab(2,ep.Sum(5,0),1))*
 (-(ep.s(2,3,4)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(0,5),1))
); }

template <class T> complex<T> A2q2g2Q14112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q14130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),3),2)*
 ep.spa(0,2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q14150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),3)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spa(3,4),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(3,2),2)*ep.spb(4,2))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.s(2,3,4)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(2,3,
   4))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(2,3,4)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))+(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),
 0),3)*ep.spa(1,3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 (-(ep.spa(0,3)*ep.spb(5,0))-
ep.spa(1,3)*ep.spb(5,1)))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0)*(-(ep.spa(0,2)*
  ep.spb(5,0))-ep.spa(1,2)*
 ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
(ep.s(3,4,5)*ep.spa(0,2)*
ep.spa(1,2)*ep.spb(4,3)*
ep.spb(5,0)+pow(ep.spa(1,2),2)*
ep.s(3,4,5)*ep.spb(4,3)*
ep.spb(5,1))
); }

template <class T> complex<T> A2q2g2Q14163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q14168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q14272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14332_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q14420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(0,5),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q14980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q14998_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15028_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15038_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15068_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(1,0),2)*ep.spa(0,5)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(3,4)*(-(ep.s(2,3,4)*
  ep.spa(0,2))-ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),1))*
 ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),4),2)*
 ep.spab(2,ep.Sum(1,0),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q15140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(0,ep.Sum(2,3),1)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2)*ep.spab(2,ep.Sum(0,1),
3))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(2,
 ep.Sum(5,0),1),2)*pow(ep.spb(1,0),
2)*ep.spa(2,4)*ep.spab(2,ep.Sum(0,5),
1)*ep.spb(5,1))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(2,3,4)*ep.spa(3,4)*
 (ep.s(2,3,4)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(2,3,4))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(2,3,4)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(0,5),1))
); }

template <class T> complex<T> A2q2g2Q15192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q15210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),4),2)*
 ep.spa(0,2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spa(3,4),
2)*pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(4,2),3))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.s(2,3,4)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(2,3,
   4))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(2,3,4)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),3)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q15240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 (-(ep.spa(0,3)*ep.spb(5,0))-
ep.spa(1,3)*ep.spb(5,1)))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0)*(-(ep.spa(0,2)*
  ep.spb(5,0))-ep.spa(1,2)*
 ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2))/
(ep.s(3,4,5)*ep.spa(0,2)*
ep.spa(1,2)*ep.spb(4,3)*
ep.spb(5,0)+pow(ep.spa(1,2),2)*
ep.s(3,4,5)*ep.spb(4,3)*
ep.spb(5,1))
); }

template <class T> complex<T> A2q2g2Q15243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q15617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15727_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q15757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,0),
2)*ep.spab(4,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),
2)*ep.spab(4,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
(-(ep.s(1,2,3)*ep.spa(0,4)*
 ep.spa(0,5)*ep.spb(2,1)*
 ep.spb(4,3))-pow(ep.spa(0,5),2)*
ep.s(1,2,3)*ep.spb(2,1)*
ep.spb(5,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 (-(ep.spa(1,4)*ep.spb(4,3))-
ep.spa(1,5)*ep.spb(5,3)))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(1,2),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q15975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),2),2))/
(-(ep.s(1,2,3)*ep.spa(0,4)*
 ep.spa(0,5)*ep.spb(2,1)*
 ep.spb(4,3))-pow(ep.spa(0,5),2)*
ep.s(1,2,3)*ep.spb(2,1)*
ep.spb(5,3))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 (-(ep.spa(1,4)*ep.spb(4,3))-
ep.spa(1,5)*ep.spb(5,3)))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),
2))/(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(1,2),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,1),3))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(1,2),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,2),2)*ep.spb(3,1))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),5),3)*
 ep.spa(2,4))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),3)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16427_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),
2)*ep.spab(4,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*(-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,0),2)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5)))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 2),2)*ep.spab(4,ep.Sum(5,0),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(3,2)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16672_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q16765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.spab(4,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,2),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*(-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16777_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,0),2)*ep.spa(0,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5)))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 3),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(3,2)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3)))+(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16882_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,0),3)*ep.spa(3,5))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),
 2),3)*ep.spb(4,2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,1),2)*ep.spa(3,5)*
 ep.spb(2,0))/(pow(ep.spab(3,ep.Sum(4,5),
 2),2)*ep.s(0,1,2)*
 ep.spa(4,5)*(-(ep.s(0,1,2)*
  ep.spa(1,3))+ep.spa(1,2)*
 ep.spab(3,ep.Sum(0,1),2))*
 ((-ep.s(0,1)+ep.s(0,1,2))*
 ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),
 2),2)*ep.spab(1,ep.Sum(3,4),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3)*ep.spab(3,ep.Sum(1,2),
4))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q16997_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spab(0,ep.Sum(2,1),4))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(3,ep.Sum(1,2),0),3)*
 pow(ep.spb(1,0),2)*ep.spa(3,5)*
 ep.spb(2,0))/(pow(ep.spab(3,ep.Sum(4,5),
 0),2)*ep.s(0,1,2)*
 ep.spa(4,5)*(ep.s(0,1,2)*
 ep.spa(1,3)+ep.spa(0,1)*
 ep.spab(3,ep.Sum(1,2),0))*
 ((-ep.s(1,2)+ep.s(0,1,2))*
 ep.spab(3,ep.Sum(1,2),0)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),0))*ep.spab(5,ep.Sum(1,2),
0))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),
 0),2)*ep.spa(1,5)*ep.spb(2,0))/
(ep.spa(4,5)*ep.spab(1,ep.Sum(2,3),
0)*ep.spb(3,0)*ep.spb(3,2)*
 (-(ep.spa(0,1)*ep.spb(2,0))-
ep.spa(1,3)*ep.spb(3,2))*
 (-ep.spa(0,5)+(ep.spa(3,5)*
  ep.spb(3,2))/ep.spb(2,0)))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spa(4,5),2)*
 pow(ep.spab(3,ep.Sum(4,5),0),2)*
 pow(ep.spb(5,0),2)*ep.spa(1,3)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(4,0))/(pow(ep.spab(3,ep.Sum(1,2),
 0),2)*ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,4),
0)*(ep.s(1,2,3)*ep.spab(3,
  ep.Sum(1,2),0)+
(-ep.s(4,5)+ep.s(1,2,3))*
 ep.spab(3,ep.Sum(4,5),0))*
 (ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(0,5)*ep.spab(3,ep.Sum(4,5),
  0)))+(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),
 1),2)*ep.spa(3,5)*ep.spb(4,2))/
(ep.spa(0,3)*ep.spa(0,5)*
 ep.spab(3,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(3,0),4)*
 ep.spb(2,1)*((ep.spa(0,5)*
  ep.spb(2,0))/ep.spa(3,5)-
ep.spb(3,2)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 pow(ep.spb(4,0),2)*ep.spa(1,5))/
(ep.spa(1,2)*ep.spab(1,ep.Sum(3,0),
4)*ep.spab(5,ep.Sum(0,3),4)*
 ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(3,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spab(1,ep.Sum(3,0),
4)*ep.spab(3,ep.Sum(0,1),2)*
 (-(ep.spa(0,1)*ep.spb(2,0))-
ep.spa(1,3)*ep.spb(3,2))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17032_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,1),
2)*ep.spa(0,1)*ep.spa(3,5)*
 ep.spab(3,ep.Sum(0,1),2))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (ep.spa(1,2)-(ep.s(0,1,2)*
  ep.spa(1,3))/ep.spab(3,ep.Sum(0,1),
  2))*ep.spab(5,ep.Sum(0,1),2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(3,ep.Sum(0,1),2),2)*
 pow(ep.spb(2,0),3)*ep.spa(0,1)*
 ep.spa(3,5))/(pow(ep.spab(3,ep.Sum(4,5),
 2),2)*ep.s(0,1,2)*
 ep.spa(4,5)*(ep.s(0,1,2)*
 ep.spa(1,3)-ep.spa(1,2)*
 ep.spab(3,ep.Sum(0,1),2))*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spa(0,2)*
 ep.spab(2,ep.Sum(5,4),3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,0),2)*ep.spab(2,
ep.Sum(3,4),1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q17530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(4,5),
2)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spab(2,
ep.Sum(5,4),3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(2,
ep.Sum(3,4),1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q17633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2)*
 ep.spab(1,ep.Sum(5,4),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),0)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17663_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
3))/(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),3))/(ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),3)*
 ep.spa(3,4))/(ep.s(0,1,2)*
 ep.spa(1,2)*(ep.s(0,1,2)*
 ep.spa(0,4)-ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),5))*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spab(4,ep.Sum(2,1),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17698_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2)*ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3)*
 ep.spa(0,2)*ep.spa(3,4))/
(ep.s(0,1,2)*ep.spa(1,2)*
 (ep.s(0,1,2)*ep.spa(0,4)-
ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
  5))*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17703_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17742_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(4,ep.Sum(0,5),1))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),3))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17892_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),3)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),2),3)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),2)*pow(ep.spb(5,0),
2)*ep.spab(1,ep.Sum(5,4),0)*
 ep.spb(4,0))/(pow(ep.spab(1,ep.Sum(2,3),
 0),2)*ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spab(1,ep.Sum(2,3),0)+
(-ep.s(4,5)+ep.s(1,2,3))*
 ep.spab(1,ep.Sum(4,5),0))*
 (ep.s(1,2,3)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(5,4),0))
); }

template <class T> complex<T> A2q2g2Q17910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),3)*pow(ep.spb(2,1),
2)*ep.spb(3,1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(1,2,3)*ep.spa(4,5)*
 (-(ep.s(1,2,3)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(1,2,
   3))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(1,2,3)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(2,
ep.Sum(5,0),1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),3)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17913_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(2,5))/
(ep.spa(0,4)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 ep.spa(1,4)*ep.spa(2,5))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q17932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
3))/(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q17955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(5,4),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q18040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q18058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q18060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,2),2))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(4,ep.Sum(0,5),1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q18063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q18072_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q18075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(5,ep.Sum(1,2),0)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5)*ep.spab(1,ep.Sum(5,0),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),2)*pow(ep.spb(5,0),
2)*ep.spa(1,3)*ep.spab(1,ep.Sum(5,4),
0)*ep.spb(4,0))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(1,2,3)*ep.spa(2,3)*
 (ep.s(1,2,3)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(1,2,3))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(1,2,3)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(5,4),0))
); }

template <class T> complex<T> A2q2g2Q18090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),3)*pow(ep.spb(3,1),
3))/(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(1,2,3)*ep.spa(4,5)*
 (-(ep.s(1,2,3)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(1,2,
   3))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(1,2,3)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(2,
ep.Sum(5,0),1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q18093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,4),0),2)*
 ep.spa(3,5))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q19510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,4),1),2)*
 ep.spa(3,5))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q19595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q19615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(5,4),1),2)*
 ep.spa(3,5)*ep.spab(3,ep.Sum(5,4),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(5,4),
2)*ep.spa(1,3)*ep.spb(4,0))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),
 4),2)*ep.spab(5,ep.Sum(2,3),4)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q19630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q19645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(5,4),2),3)*
 ep.spa(3,5))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,3),3)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(5,4),
2)*ep.spb(4,0))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q19685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(5,4),0),2)*
 ep.spab(5,ep.Sum(3,4),1))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q19690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(0,1),5))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q19793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q19823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2)*ep.spa(3,4))/(ep.s(0,1,2)*
 ep.spa(1,2)*(ep.s(0,1,2)*
 ep.spa(0,4)-ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),5))*
 ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),1),3))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(5,1),3))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q19833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(0,1),2),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q19858_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q19860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,4),
2)*ep.spa(0,2)*ep.spa(3,4))/
(ep.s(0,1,2)*ep.spa(1,2)*
 (ep.s(0,1,2)*ep.spa(0,4)-
ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),2),2)*
 ep.spab(0,ep.Sum(5,4),1))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q19863_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20293_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(0,2)*
 ep.spa(2,5))/(ep.spa(0,3)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(2,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(0,4)*ep.spa(2,5))/
(ep.spa(0,3)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,4)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(0,4),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),1)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20363_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20423_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20443_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(ep.s(0,1,5)*
 ep.spa(1,3)-ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,4),2)*ep.spb(4,0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),2),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20493_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.spab(5,ep.Sum(3,2),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20578_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20653_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(1,2),
2)*pow(ep.spab(1,ep.Sum(5,0),4),3)*
 pow(ep.spb(5,4),2)*ep.spa(1,3)*
 ep.spb(4,0))/(pow(ep.spab(1,ep.Sum(2,3),
 4),2)*ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spab(1,ep.Sum(2,3),4)+
(-ep.s(0,5)+ep.s(1,2,3))*
 ep.spab(1,ep.Sum(5,0),4))*
 (ep.s(1,2,3)*ep.spa(1,5)-
ep.spa(4,5)*ep.spab(1,ep.Sum(5,0),
  4))*ep.spab(3,ep.Sum(5,0),4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 pow(ep.spb(4,0),2)*ep.spa(3,5))/
(ep.spa(2,3)*ep.spab(3,ep.Sum(1,4),
0)*ep.spab(5,ep.Sum(0,1),4)*
 ep.spab(5,ep.Sum(4,1),0)*
 ep.spb(1,0)*ep.spb(4,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(2,3),2)*
 pow(ep.spab(1,ep.Sum(2,3),4),3)*
 pow(ep.spb(4,3),2)*ep.spa(1,5)*
 ep.spb(4,2))/(pow(ep.spab(1,ep.Sum(5,0),
 4),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(ep.s(0,1,5)*
 ep.spa(1,3)+ep.spa(3,4)*
 ep.spab(1,ep.Sum(2,3),4))*
 ((-ep.s(2,3)+ep.s(0,1,5))*
 ep.spab(1,ep.Sum(2,3),4)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),4))*ep.spab(5,ep.Sum(2,3),
4))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spab(1,ep.Sum(5,4),3),2)*
 ep.spb(2,0))/(ep.spa(1,4)*
 ep.spa(4,5)*ep.spab(1,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(1,4),0)*
 ep.spb(3,2)*(ep.spa(1,5)*
 ep.spb(2,1)-ep.spa(4,5)*
 ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 pow(ep.spb(4,2),2)*ep.spa(3,5))/
(ep.spa(0,5)*ep.spab(3,ep.Sum(1,2),
4)*ep.spb(2,1)*ep.spb(4,1)*
 (-(ep.spa(1,3)*ep.spb(2,1))-
ep.spa(3,4)*ep.spb(4,2))*
 (ep.spa(1,5)*ep.spb(2,1)-
ep.spa(4,5)*ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(2,0))/(ep.spa(1,4)*
 ep.spa(3,4)*ep.spab(1,ep.Sum(3,4),
2)*ep.spab(3,ep.Sum(1,4),0)*
 (-(ep.spa(1,3)*ep.spb(2,1))-
ep.spa(3,4)*ep.spb(4,2))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),3),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20683_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q20747_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20767_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20777_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q20935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q20957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q20985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21083_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21093_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21443_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21503_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q21523_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q21565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(3,2),2)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(-(ep.s(0,1,5)*
  ep.spa(1,3))+ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q21583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q21817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(5,4),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q21835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(5,1),
3))/(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3)*
 ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q21852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q21870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21873_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q21925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spab(3,ep.Sum(0,1),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(3,2),2)*ep.spa(1,5)*
 ep.spb(4,2))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(-(ep.s(0,1,5)*
  ep.spa(1,3))+ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ((-ep.s(3,4)+ep.s(0,1,5))*
 ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))+(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),
 5),2)*ep.spa(1,3)*ep.spab(1,
ep.Sum(2,3),0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q21933_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q21943_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q22042_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22072_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22182_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q22195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,4),2)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22383_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.spa(2,4))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22738_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q22813_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q22825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(1,2),3),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.spab(4,ep.Sum(2,3),0))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q22843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q22897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spab(5,ep.Sum(1,0),3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q22915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(5,1),
3)*ep.spa(2,4))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q22932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q22950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q22953_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q23005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(4,2),3)*ep.spa(1,5))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(0,1,5)*ep.spa(0,5)*
 (-(ep.s(0,1,5)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(0,1,
   5))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))-(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),
 2),3)*ep.spb(2,0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q23010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q23013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q23023_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q23645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),2),2)*
 ep.spa(0,4)*ep.spab(4,ep.Sum(5,0),
3))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,0),
2)*ep.spa(2,4)*ep.spb(5,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(0,ep.Sum(3,4),5)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,1),
3)*ep.spa(2,4))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q23720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),1),2)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q23750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),3),3)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,0),
2)*ep.spb(5,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(0,1,
   5))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(0,1,5)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(0,
ep.Sum(3,4),5))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,1),
3))/(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(0,1,
   5))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(0,1,5)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(0,
ep.Sum(3,4),5))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q23935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q23950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q23965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),2),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
 ep.spa(0,4)*ep.spab(4,ep.Sum(0,5),
3))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,0),
2)*ep.spa(2,4)*ep.spb(5,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (-(ep.s(2,3,4)*ep.spa(0,4))+
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(0,ep.Sum(3,4),5)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(4,
 ep.Sum(0,1),5),2)*pow(ep.spb(5,1),
3)*ep.spa(2,4)*ep.spab(4,ep.Sum(1,0),
5))/(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 (ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,5),3),3)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(4,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,0),
2)*ep.spb(5,1))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 (-(ep.s(2,3,4)*ep.spa(0,4))+
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(0,
ep.Sum(3,4),5))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(4,
 ep.Sum(0,1),5),2)*pow(ep.spb(5,1),
3)*ep.spab(4,ep.Sum(1,0),5))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 (ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(0,
ep.Sum(3,4),5))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q24725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q24730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q24785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),1),2)*
 ep.spab(0,ep.Sum(4,5),2))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q24820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q24830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),
2))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
 ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q25085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(-(ep.s(3,4,5)*ep.spa(2,
   4))+ep.spa(3,4)*ep.spab(2,
  ep.Sum(4,5),3))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(1,0)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1)))
); }

template <class T> complex<T> A2q2g2Q25090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),3))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(1,0)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1)))+(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),3))/(ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spa(4,5))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (-(ep.s(3,4,5)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3))*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q25193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q25223_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(3,4),
2)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,3),3))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q25233_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(3,2),2)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),2)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q25258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,3),3)*ep.spa(0,2))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),3)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q25263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25433_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25448_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25478_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q25565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25583_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),
2)*ep.spa(4,5))/(ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spa(1,5)+ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),0))*
 ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),3))/
(ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),2),3))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q25680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q25683_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25688_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q25780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q25838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),
2)*ep.spa(1,3)*ep.spa(4,5))/
(ep.s(1,2,3)*ep.spa(2,3)*
 (ep.s(1,2,3)*ep.spa(1,5)+
ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),3),2)*
 ep.spab(1,ep.Sum(0,5),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q25860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q25863_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q25985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
(-(ep.s(0,1,2)*ep.spa(3,5)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(3,2))-pow(ep.spa(4,5),2)*
ep.s(0,1,2)*ep.spb(1,0)*
ep.spb(4,2))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 (-(ep.spa(0,3)*ep.spb(3,2))-
ep.spa(0,4)*ep.spb(4,2)))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(3,5)*ep.spb(3,2)+
ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q25990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2))/
(-(ep.s(0,1,2)*ep.spa(3,5)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(3,2))-pow(ep.spa(4,5),2)*
ep.s(0,1,2)*ep.spb(1,0)*
ep.spb(4,2))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 (-(ep.spa(0,3)*ep.spb(3,2))-
ep.spa(0,4)*ep.spb(4,2)))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 (ep.spa(3,5)*ep.spb(3,2)+
ep.spa(4,5)*ep.spb(4,2))*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(-(ep.s(0,4,5)*ep.spa(3,
   5))+ep.spa(4,5)*ep.spab(3,
  ep.Sum(5,0),4)))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spab(3,ep.Sum(4,5),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,5))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (-(ep.s(0,4,5)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),
 2),3))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(2,
 ep.Sum(0,1),5),3)*pow(ep.spb(5,0),
2)*ep.spa(2,4)*ep.spb(5,1))/
(pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spa(0,2)-
ep.spa(0,5)*ep.spab(2,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(0,1,
   5))*ep.spab(2,ep.Sum(0,1),5)+
ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),5))*ep.spab(4,ep.Sum(0,1),
5))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spa(2,5)*ep.spab(0,ep.Sum(2,5),
3)*ep.spab(2,ep.Sum(5,0),1)*
 ep.spb(4,3)*(-(ep.spa(0,2)*
  ep.spb(2,1))-ep.spa(0,5)*
 ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spa(2,4)*ep.spb(3,1))/
(ep.spa(2,5)*ep.spa(4,5)*
 ep.spab(2,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(2,5),3)*
 ep.spb(1,0)*(-ep.spb(2,1)+
(ep.spa(4,5)*ep.spb(5,1))/
 ep.spa(2,4)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 pow(ep.spb(5,3),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spab(0,ep.Sum(2,5),
3)*ep.spab(4,ep.Sum(2,3),5)*
 ep.spab(4,ep.Sum(5,2),3)*
 ep.spb(3,2)*ep.spb(5,2))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(0,4)*ep.spb(5,1))/
(ep.spa(3,4)*ep.spab(0,ep.Sum(1,2),
5)*ep.spb(2,1)*(ep.spa(4,5)-
(ep.spa(2,4)*ep.spb(2,1))/
 ep.spb(5,1))*(-(ep.spa(0,2)*
  ep.spb(2,1))-ep.spa(0,5)*
 ep.spb(5,1))*ep.spb(5,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(2,ep.Sum(3,4),5),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,2)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(5,3))/(pow(ep.spab(2,ep.Sum(0,1),
 5),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,3),
5)*(ep.s(0,1,2)*ep.spab(2,
  ep.Sum(0,1),5)+
(-ep.s(3,4)+ep.s(0,1,2))*
 ep.spab(2,ep.Sum(3,4),5))*
 (ep.s(0,1,2)*ep.spa(2,4)+
ep.spa(4,5)*ep.spab(2,ep.Sum(3,4),
  5)))
); }

template <class T> complex<T> A2q2g2Q26170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(2,ep.Sum(5,0),1),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,5)*
 ep.spa(2,4))/(pow(ep.spab(2,ep.Sum(3,4),
 1),2)*ep.s(0,1,5)*
 ep.spa(3,4)*(ep.s(0,1,5)*
 ep.spa(0,2)-ep.spa(0,1)*
 ep.spab(2,ep.Sum(5,0),1))*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26278_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),1),3)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(5,1),3)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(3,4),2)*
 pow(ep.spab(0,ep.Sum(3,4),5),2)*
 pow(ep.spb(5,4),2)*ep.spab(0,
ep.Sum(4,3),5)*ep.spb(5,3))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(0,1,2)*ep.spa(1,2)*
 (ep.s(0,1,2)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(0,1,2))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(0,1,2)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(4,3),5))
); }

template <class T> complex<T> A2q2g2Q26313_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,3)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(3,5))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(3,5)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(4,3),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),2)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q26338_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4)*ep.spab(0,ep.Sum(4,5),
1))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(4,ep.Sum(0,1),5)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(0,ep.Sum(3,4),5),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,2)*
 ep.spab(0,ep.Sum(4,3),5)*
 ep.spb(5,3))/(pow(ep.spab(0,ep.Sum(1,2),
 5),2)*ep.s(0,1,2)*
 ep.spa(1,2)*(ep.s(0,1,2)*
 ep.spab(0,ep.Sum(1,2),5)+
(-ep.s(3,4)+ep.s(0,1,2))*
 ep.spab(0,ep.Sum(3,4),5))*
 (ep.s(0,1,2)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(4,3),5))
); }

template <class T> complex<T> A2q2g2Q26343_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26633_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26758_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26773_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q26825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q26923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),
2))/(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(4,3),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q26970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q26973_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q26983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q27040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27058_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spab(4,ep.Sum(0,5),2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q27150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q27153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27593_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27608_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27628_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27638_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,4)*ep.spa(3,4)*
 ep.spa(3,5))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,3)*ep.spa(2,4)*
 ep.spa(3,5)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27718_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q27838_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q27853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28268_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28348_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28378_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28448_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3)*
 ep.spa(1,3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(3,5))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 ep.spa(0,3)*ep.spa(1,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,4)*ep.spa(2,3)*
 ep.spa(3,5)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.spab(0,ep.Sum(4,5),2))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q28580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),
2))/(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
 ep.spab(3,ep.Sum(4,5),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q28615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
2)*ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),2)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(1,5),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q28705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,1),
2)*ep.spab(3,ep.Sum(1,0),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),2)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q28745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q28760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,0),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*(ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0)))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),3))/
(ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q28823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q28833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q28868_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28928_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q28975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q28985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),
2)*ep.spa(1,3)*ep.spa(4,5)*
 ep.spab(1,ep.Sum(4,5),0))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (-ep.spa(0,5)+(ep.s(0,4,5)*
  ep.spa(1,5))/ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q29003_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q29013_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29063_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29083_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29213_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29348_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(ep.s(0,1,2)*ep.spa(2,4)-
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3))*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))
); }

template <class T> complex<T> A2q2g2Q29605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),3)*ep.spa(0,5))/
(ep.s(2,3,4)*ep.spa(3,4)*
 (-(ep.s(2,3,4)*ep.spa(0,2))-
ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
  1))*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),3))/(ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),3),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(3,4),
2)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,3),3))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q29643_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29648_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(3,2),2)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),2)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,3),3)*ep.spa(0,2))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),3)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q29673_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29683_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q29860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(-(ep.s(1,2,3)*ep.spa(3,
   5))+ep.spa(4,5)*ep.spab(3,
  ep.Sum(5,0),4)))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spab(3,ep.Sum(4,5),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),
2))/(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),3))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q29960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(4,3),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q29965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(4,3),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q30040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q30050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,4),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*(ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0)))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q30118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q30123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30158_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q30265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spab(4,ep.Sum(0,5),2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q30298_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q30303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30358_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30508_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30518_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30538_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q30608_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(2,3),
2)*pow(ep.spab(2,ep.Sum(0,1),5),3)*
 pow(ep.spb(5,0),2)*ep.spa(2,4)*
 ep.spb(5,1))/(pow(ep.spab(2,ep.Sum(3,4),
 5),2)*ep.s(2,3,4)*
 ep.spa(3,4)*(-(ep.s(2,3,4)*
  ep.spa(0,2))+ep.spa(0,5)*
 ep.spab(2,ep.Sum(0,1),5))*
 ((-ep.s(0,1)+ep.s(2,3,4))*
 ep.spab(2,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(2,ep.Sum(3,
   4),5))*ep.spab(4,ep.Sum(0,1),
5))+(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),
 5),2)*pow(ep.spb(5,1),2)*
 ep.spa(0,4))/(ep.spa(3,4)*
 ep.spab(0,ep.Sum(1,2),5)*
 ep.spab(0,ep.Sum(5,2),1)*
 ep.spab(4,ep.Sum(2,5),1)*
 ep.spb(2,1)*ep.spb(5,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(2,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,4),2)*ep.spa(0,2)*
 ep.spb(5,3))/(pow(ep.spab(2,ep.Sum(0,1),
 5),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(3,4),
5)*(ep.s(0,1,2)*ep.spab(2,
  ep.Sum(0,1),5)+
(-ep.s(3,4)+ep.s(0,1,2))*
 ep.spab(2,ep.Sum(3,4),5))*
 (ep.s(0,1,2)*ep.spa(2,4)+
ep.spa(4,5)*ep.spab(2,ep.Sum(3,4),
  5)))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spab(2,ep.Sum(0,5),4),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spa(2,5)*ep.spab(0,ep.Sum(2,5),
1)*ep.spab(2,ep.Sum(0,5),1)*
 ep.spb(4,3)*(-(ep.spa(0,2)*
  ep.spb(3,2))+ep.spa(0,5)*
 ep.spb(5,3)))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(3,1))/(ep.spa(2,5)*
 ep.spa(4,5)*ep.spab(2,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(2,5),1)*
 ep.spb(1,0)*(-(ep.spa(2,4)*
  ep.spb(3,2))-ep.spa(4,5)*
 ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 pow(ep.spb(5,3),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spab(4,ep.Sum(2,3),
5)*ep.spb(3,2)*ep.spb(5,2)*
 (-(ep.spa(0,2)*ep.spb(3,2))+
ep.spa(0,5)*ep.spb(5,3))*
 (-(ep.spa(2,4)*ep.spb(3,2))-
ep.spa(4,5)*ep.spb(5,3)))
); }

template <class T> complex<T> A2q2g2Q30685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,5)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(3,4)*(-(ep.s(2,3,4)*
  ep.spa(0,2))-ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),1))*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q30710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),4),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q30720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),1),3)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(5,1),3)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(3,4),2)*
 pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,4),2)*ep.spb(5,3))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))
); }

template <class T> complex<T> A2q2g2Q30723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,5)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 ep.spa(1,4)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,5)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30728_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),2)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q30750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4)*ep.spab(0,ep.Sum(4,5),
1))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(4,ep.Sum(0,1),5)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(0,ep.Sum(3,4),5),3)*
 pow(ep.spb(5,4),2)*ep.spa(0,2)*
 ep.spb(5,3))/(pow(ep.spab(0,ep.Sum(1,2),
 5),2)*ep.s(3,4,5)*
 ep.spa(1,2)*(ep.s(3,4,5)*
 ep.spab(0,ep.Sum(1,2),5)+
(-ep.s(3,4)+ep.s(3,4,5))*
 ep.spab(0,ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))
); }

template <class T> complex<T> A2q2g2Q30753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30788_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q30793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q31157_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31172_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31187_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31207_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31297_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(4,
 ep.Sum(2,3),1),3)*pow(ep.spb(2,1),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(2,3),1)*
 (-(ep.s(0,4,5)*ep.spa(2,4))-
ep.spa(1,2)*ep.spab(4,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(4,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(4,ep.Sum(5,
   0),1)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 pow(ep.spb(3,1),2)*ep.spa(0,2))/
(ep.spa(0,5)*ep.spab(0,ep.Sum(4,1),
3)*ep.spab(2,ep.Sum(1,4),3)*
 ep.spab(2,ep.Sum(3,4),1)*
 ep.spb(4,1)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(5,0),1),3)*
 pow(ep.spb(1,0),2)*ep.spa(2,4)*
 ep.spb(5,1))/(pow(ep.spab(4,ep.Sum(2,3),
 1),2)*ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(5,0),
1)*(ep.s(2,3,4)*ep.spab(4,
  ep.Sum(2,3),1)+
(-ep.s(0,5)+ep.s(2,3,4))*
 ep.spab(4,ep.Sum(5,0),1))*
 (-(ep.s(2,3,4)*ep.spa(0,4))+
ep.spa(0,1)*ep.spab(4,ep.Sum(5,0),
  1)))+(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spa(1,4)*ep.spab(0,ep.Sum(4,1),
3)*ep.spab(4,ep.Sum(0,1),5)*
 ep.spb(3,2)*(ep.spa(0,1)*
 ep.spb(5,1)+ep.spa(0,4)*
 ep.spb(5,4)))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spab(4,ep.Sum(2,1),0),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spa(1,4)*ep.spab(2,ep.Sum(4,1),
3)*ep.spab(4,ep.Sum(2,1),3)*
 ep.spb(5,0)*(ep.spa(1,2)*
 ep.spb(5,1)-ep.spa(2,4)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 pow(ep.spb(5,1),2)*ep.spa(0,2))/
(ep.spa(2,3)*ep.spab(0,ep.Sum(4,5),
1)*ep.spb(4,1)*ep.spb(5,4)*
 (ep.spa(0,1)*ep.spb(5,1)+
ep.spa(0,4)*ep.spb(5,4))*
 (ep.spa(1,2)*ep.spb(5,1)-
ep.spa(2,4)*ep.spb(5,4)))
); }

template <class T> complex<T> A2q2g2Q31422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.spab(1,ep.Sum(5,0),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),
2))/(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spab(1,ep.Sum(3,2),5))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31512_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31547_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(1,2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.s(0,4,5)*ep.spa(2,
   4))-ep.spa(2,3)*ep.spab(4,
  ep.Sum(5,0),3))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31602_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(3,2),1))-
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 ep.spab(4,ep.Sum(5,0),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spa(1,2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(0,4,5)*ep.spa(2,4))-
ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
  3)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(3,2),
 0),2)*ep.spab(4,ep.Sum(3,2),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,1),3)*ep.spa(0,4))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(1,0),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),3)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 pow(ep.spb(3,2),2)*ep.spa(0,4)*
 ep.spab(4,ep.Sum(2,1),3)*
 ep.spb(3,1))/(pow(ep.spab(4,ep.Sum(5,0),
 3),2)*ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(2,1),
3)*(-(ep.s(0,4,5)*ep.spa(2,
   4))+ep.spa(2,3)*ep.spab(4,
  ep.Sum(1,2),3))*
 ((-ep.s(1,2)+ep.s(0,4,5))*
 ep.spab(4,ep.Sum(1,2),3)+
ep.s(0,4,5)*ep.spab(4,ep.Sum(5,
   0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4)*ep.spab(4,ep.Sum(2,3),
5))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(2,ep.Sum(4,5),3)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q31940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q31945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32197_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(5,0),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(2,0),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,5),3),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32257_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,2),2)*ep.spab(4,
ep.Sum(2,1),3))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(4,ep.Sum(5,0),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q32417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),0),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q32525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q32837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),
3)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),3),3)*
 ep.spa(0,2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(1,0),
2)*ep.spb(5,1))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))
); }

template <class T> complex<T> A2q2g2Q32862_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(1,0),2))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*(ep.spa(2,4)*ep.spb(2,1)+
ep.spa(3,4)*ep.spb(3,1)))-
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),3))/
(ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q32945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),3)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spa(3,4),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(3,2),2)*ep.spb(4,2))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(0,1,5)*ep.spa(0,5)*
 (ep.s(0,1,5)*ep.spa(1,3)-
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(0,1,
   5))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))+(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),
 0),3)*ep.spa(1,3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32963_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,4)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,4))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,4)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q32970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.spab(3,ep.Sum(1,0),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q32973_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q33012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q33030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.spa(0,2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33033_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q33065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(3,2),2)*ep.spa(3,4))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.s(2,3,4)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*ep.spab(5,ep.Sum(3,4),2))+
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),3))/
(ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q33150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q33173_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q33188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33243_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33248_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33552_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33603_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33608_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33707_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33727_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33737_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q33815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33882_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q33917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(0,ep.Sum(2,3),1)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2)*ep.spab(2,ep.Sum(0,1),
3))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(1,0),
2)*ep.spa(2,4)*ep.spb(5,1))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))
); }

template <class T> complex<T> A2q2g2Q33942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q33945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(1,0),2)*ep.spa(0,5)*
 ep.spa(2,4)*ep.spab(2,ep.Sum(5,0),
1))/(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.spa(0,1)-(ep.s(0,1,5)*
  ep.spa(0,2))/ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q34025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spa(3,4),
2)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spab(1,
ep.Sum(4,3),2))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(0,1,5)*ep.spa(0,5)*
 (-(ep.s(0,1,5)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(0,1,
   5))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(4,3),
2))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),3)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34043_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spab(3,ep.Sum(1,0),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34053_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34092_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q34110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34113_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,3),0),2)*
 ep.spa(2,4))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(1,0),2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),1),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q34413_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q34483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34533_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34777_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34812_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34893_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q34903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(4,0),
3)*ep.spa(1,3))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(1,2,3)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),3)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q35060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,3),3)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(4,0),
3))/(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(0,4,5)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(0,4,5))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(0,4,5)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q35095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q35185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q35225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
3))/(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spa(3,4))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.s(2,3,4)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(4,0),3))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*(ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(2,1),
2)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),1)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35313_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35333_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35348_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35408_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 (ep.spa(2,5)*ep.spb(2,1)+
ep.spa(3,5)*ep.spb(3,1)))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(2,4)*ep.spb(2,1)+
ep.spa(3,4)*ep.spb(3,1))*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2))/
(ep.s(0,1,5)*ep.spa(2,4)*
ep.spa(3,4)*ep.spb(2,1)*
ep.spb(5,0)+pow(ep.spa(3,4),2)*
ep.s(0,1,5)*ep.spb(3,1)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),2)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(1,ep.Sum(4,5),0),2)*
 pow(ep.spb(4,0),3)*ep.spa(1,3)*
 ep.spa(4,5))/(pow(ep.spab(1,ep.Sum(2,3),
 0),2)*ep.s(0,4,5)*
 ep.spa(2,3)*(-(ep.s(0,4,5)*
  ep.spa(1,5))+ep.spa(0,5)*
 ep.spab(1,ep.Sum(4,5),0))*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
2)*ep.spab(0,ep.Sum(3,2),1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),1)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),1),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q35493_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35563_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35693_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q35828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q35833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(4,3),2)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(-(ep.s(0,1,2)*ep.spa(2,
   4))+ep.spa(3,4)*ep.spab(2,
  ep.Sum(4,5),3)))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q36085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),
2))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),3))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),3))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q36110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q36120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(0,
ep.Sum(2,3),1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(4,3),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))
); }

template <class T> complex<T> A2q2g2Q36123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36128_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q36150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),2),2)*
 ep.spab(3,ep.Sum(5,4),1))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q36153_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36188_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q36193_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q37592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37597_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2)*ep.spb(2,0))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(4,0),3)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4))*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(0,5),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37662_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),3)*ep.spa(0,1))/
(ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))-
ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
  2))*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2)*ep.spa(0,1))/(ep.s(3,4,5)*
 ep.spa(4,5)*(-(ep.s(3,4,5)*
  ep.spa(1,3))-ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),2))*
 ep.spab(5,ep.Sum(0,1),2))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,0),3))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,5),
3)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 pow(ep.spb(2,1),2)*ep.spab(3,
ep.Sum(1,0),2)*ep.spb(2,0))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(3,4,
   5))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(3,4,5)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(1,0),
2))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),3)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),4),3)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37812_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),4),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37872_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(4,0),
3))/(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3)*
 ep.spab(5,ep.Sum(1,2),0))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q37890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(2,3),
2)*pow(ep.spab(0,ep.Sum(2,3),1),3)*
 pow(ep.spb(2,1),2)*ep.spa(0,4)*
 ep.spb(3,1))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(0,4,5)*
 ep.spa(4,5)*(ep.s(0,4,5)*
 ep.spa(0,2)-ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ((-ep.s(2,3)+ep.s(0,4,5))*
 ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),
 1),2)*ep.spab(2,ep.Sum(5,0),1)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),4),2)*
 ep.spa(0,2)*ep.spab(0,ep.Sum(2,1),
5))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37893_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q37910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(2,1),2)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.s(1,2,3)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q37923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q37928_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38017_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38022_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,3),2)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(4,0),
3)*ep.spa(1,3))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q38070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),2)*pow(ep.spb(3,1),
3)*ep.spa(0,4)*ep.spab(0,ep.Sum(3,2),
1))/(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (-(ep.s(0,4,5)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(3,2),
1))-(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),
 1),3)*ep.spb(5,1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38073_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),4),2)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),2),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),4),2)*
 ep.spab(3,ep.Sum(1,2),5))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4)*ep.spab(2,ep.Sum(3,4),
1))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(4,ep.Sum(1,2),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(4,3),2)*ep.spa(0,2)*
 ep.spb(5,3))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(ep.s(0,1,2)*ep.spab(2,
  ep.Sum(0,1),3)+
(-ep.s(4,5)+ep.s(0,1,2))*
 ep.spab(2,ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A2q2g2Q38245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),1),3)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),3)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(4,5),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(4,3),2)*ep.spb(5,3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(3,4,5)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(3,4,5))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(3,4,5)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A2q2g2Q38270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.s(1,2,3)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38283_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),2)*
 (ep.spa(1,4)*ep.spb(1,0)+
ep.spa(2,4)*ep.spb(2,0)))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
(-(ep.s(0,4,5)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spb(1,0)*
 ep.spb(5,4))-pow(ep.spa(2,3),2)*
ep.s(0,4,5)*ep.spb(2,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2)*
 pow(ep.spb(3,1),2)*ep.spa(2,4))/
(ep.spa(4,5)*ep.spab(2,ep.Sum(0,1),
3)*ep.spab(2,ep.Sum(3,0),1)*
 ep.spab(4,ep.Sum(0,3),1)*
 ep.spb(1,0)*ep.spb(3,0))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(1,2),2)*
 pow(ep.spab(0,ep.Sum(1,2),3),2)*
 pow(ep.spb(3,2),2)*ep.spa(0,4)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(3,1))/(pow(ep.spab(0,ep.Sum(4,5),
 3),2)*ep.s(0,4,5)*
 ep.spa(4,5)*(ep.s(0,4,5)*
 ep.spa(0,2)+ep.spa(2,3)*
 ep.spab(0,ep.Sum(1,2),3))*
 ((-ep.s(1,2)+ep.s(0,4,5))*
 ep.spab(0,ep.Sum(1,2),3)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),3))*ep.spab(4,ep.Sum(2,1),
3))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(0,
 ep.Sum(4,5),3),3)*pow(ep.spb(4,3),
2)*ep.spa(0,2)*ep.spb(5,3))/
(pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),3)+(-ep.s(4,5)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(4,5),3))*
 (-(ep.s(3,4,5)*ep.spa(0,4))+
ep.spa(3,4)*ep.spab(0,ep.Sum(4,5),
  3))*ep.spab(2,ep.Sum(4,5),3))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(5,1))/(ep.spa(0,3)*
 ep.spa(3,4)*ep.spab(0,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(0,3),1)*
 ep.spb(2,1)*(-(ep.spa(0,4)*
  ep.spb(5,0))-ep.spa(3,4)*
 ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(2,4)*ep.spb(5,3))/
(ep.spa(1,2)*ep.spab(4,ep.Sum(5,0),
3)*ep.spb(3,0)*ep.spb(5,0)*
 (ep.spa(2,3)-(ep.spa(0,2)*
  ep.spb(5,0))/ep.spb(5,3))*
 (-(ep.spa(0,4)*ep.spb(5,0))-
ep.spa(3,4)*ep.spb(5,3)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spa(0,2)*ep.spb(5,1))/
(ep.spa(0,3)*ep.spa(2,3)*
 ep.spab(0,ep.Sum(2,3),1)*
 ep.spab(2,ep.Sum(0,3),1)*
 (ep.spb(5,0)-(ep.spa(2,3)*
  ep.spb(5,3))/ep.spa(0,2))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38313_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38323_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38348_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(0,3)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q38932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38942_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q38962_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q38977_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39032_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39037_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39122_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39172_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2)*ep.spb(3,1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(-(ep.s(2,3,4)*ep.spa(0,
   4))+ep.spa(0,5)*ep.spab(4,
  ep.Sum(0,1),5))*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,0),2)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 2),2)*ep.spab(4,ep.Sum(5,0),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),2)*
 ep.spab(4,ep.Sum(1,0),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39242_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),3)*ep.spa(1,2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(0,4,5)*ep.spa(2,4))-
ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
  3))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),3))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(3,2),1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2)*ep.spa(1,2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.s(0,4,5)*ep.spa(2,
   4))-ep.spa(2,3)*ep.spab(4,
  ep.Sum(5,0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(3,2),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39462_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(1,2),
2)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 pow(ep.spb(3,1),3))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(1,2),
2)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 pow(ep.spb(3,2),2)*ep.spab(4,
ep.Sum(2,1),3)*ep.spb(3,1))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(2,1),3)*
 (-(ep.s(0,4,5)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(0,4,
   5))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(0,4,5)*ep.spab(4,ep.Sum(5,
   0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),5),3)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,3),3)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q39680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39752_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(5,4),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(5,0),2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(5,4),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(3,4),5),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,0),
2)*ep.spab(4,ep.Sum(0,1),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q39865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),
2)*ep.spab(4,ep.Sum(1,0),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q40192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40202_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),5),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q40300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40332_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40612_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,4),3)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(5,1),
3))/(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3)*
 ep.spab(0,ep.Sum(2,3),1))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q40632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,0),2))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*(ep.spa(2,4)*ep.spb(2,1)+
ep.spa(3,4)*ep.spb(3,1)))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*(ep.spa(2,4)*
 ep.spb(2,1)+ep.spa(3,4)*
 ep.spb(3,1))*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q40720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spab(3,ep.Sum(0,1),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),3)*
 pow(ep.spb(3,2),2)*ep.spa(1,5)*
 ep.spb(4,2))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(ep.s(0,1,5)*
 ep.spa(1,3)-ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ((-ep.s(3,4)+ep.s(0,1,5))*
 ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))+(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),
 5),2)*ep.spa(1,3)*ep.spab(1,
ep.Sum(3,2),0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40738_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),
2)*ep.spab(0,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),5)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(4,2),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40743_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40752_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q40770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),
2))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40773_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(3,2),2)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.s(2,3,4)*
  ep.spa(1,3))+ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40918_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q40923_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40948_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q40958_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40983_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q40988_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41042_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41133_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41163_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41168_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41482_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41512_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,4),2)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q41590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41617_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q41692_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(2,
 ep.Sum(5,0),1),3)*pow(ep.spb(5,1),
3)*ep.spa(2,4))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),3)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q41712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spab(5,ep.Sum(1,0),3))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q41800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spab(1,ep.Sum(4,3),2))/
(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(0,1,5)*ep.spa(0,5)*
 (-(ep.s(0,1,5)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(0,1,
   5))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(4,3),
2))-(complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),
 2),3)*ep.spb(2,0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41818_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q41820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,4),
2)*ep.spab(0,ep.Sum(4,3),5))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),5)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41823_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q41832_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q41850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q41853_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q42130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,3),5),2)*
 ep.spa(2,4))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q42145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(1,2),3),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,3),5),2)*
 ep.spab(4,ep.Sum(2,3),0))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q42178_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q42180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q42183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q42238_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2g2Q42253_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42273_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42283_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42393_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42453_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42463_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spa(3,5)*ep.spab(3,ep.Sum(4,5),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(5,4),
2)*ep.spa(1,3)*ep.spb(4,0))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(1,2,3)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),
 4),2)*ep.spab(5,ep.Sum(2,3),4)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q42830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),3)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,3),3)*pow(ep.spab(3,
 ep.Sum(5,0),4),3)*pow(ep.spb(5,4),
2)*ep.spb(4,0))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 (ep.s(0,4,5)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(0,4,5))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(0,4,5)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q42865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q42920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.spa(3,5))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q42925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2g2Q43000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*(ep.spa(1,3)*
 ep.spb(1,0)+ep.spa(2,3)*
 ep.spb(2,0))*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(2,3,4)*
 ep.spa(0,5)*(-(ep.s(2,3,4)*
  ep.spa(1,3))+ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,4),2)*ep.spb(4,0))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43078_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43083_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43108_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43118_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43143_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43148_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2)*
 (ep.spa(2,5)*ep.spb(2,1)+
ep.spa(3,5)*ep.spb(3,1)))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 (ep.spa(2,4)*ep.spb(2,1)+
ep.spa(3,4)*ep.spb(3,1))*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),5),2))/
(ep.s(0,1,5)*ep.spa(2,4)*
ep.spa(3,4)*ep.spb(2,1)*
ep.spb(5,0)+pow(ep.spa(3,4),2)*
ep.s(0,1,5)*ep.spb(3,1)*
ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(1,2),2)*pow(ep.spab(1,
 ep.Sum(5,0),4),3)*pow(ep.spb(5,4),
2)*ep.spa(1,3)*ep.spb(4,0))/
(pow(ep.spab(1,ep.Sum(2,3),4),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),4)+(-ep.s(0,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(5,0),4))*
 (-(ep.s(0,4,5)*ep.spa(1,5))+
ep.spa(4,5)*ep.spab(1,ep.Sum(5,0),
  4))*ep.spab(3,ep.Sum(5,0),4))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(2,0))/(ep.spa(1,4)*
 ep.spa(4,5)*ep.spab(1,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(1,4),2)*
 ep.spb(3,2)*(ep.spa(1,5)*
 ep.spb(1,0)+ep.spa(4,5)*
 ep.spb(4,0)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 pow(ep.spb(4,2),2)*ep.spa(3,5))/
(ep.spa(0,5)*ep.spab(3,ep.Sum(1,2),
4)*ep.spab(3,ep.Sum(4,1),2)*
 ep.spab(5,ep.Sum(1,4),2)*
 ep.spb(2,1)*ep.spb(4,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(3,5)*ep.spb(4,0))/
(ep.spa(2,3)*ep.spab(5,ep.Sum(0,1),
4)*ep.spb(1,0)*(ep.spa(3,4)-
(ep.spa(1,3)*ep.spb(1,0))/
 ep.spb(4,0))*(ep.spa(1,5)*
 ep.spb(1,0)+ep.spa(4,5)*
 ep.spb(4,0))*ep.spb(4,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(2,3),2)*
 pow(ep.spab(1,ep.Sum(2,3),4),2)*
 pow(ep.spb(4,3),2)*ep.spa(1,5)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(4,2))/(pow(ep.spab(1,ep.Sum(5,0),
 4),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(ep.s(0,1,5)*
 ep.spa(1,3)+ep.spa(3,4)*
 ep.spab(1,ep.Sum(2,3),4))*
 ((-ep.s(2,3)+ep.s(0,1,5))*
 ep.spab(1,ep.Sum(2,3),4)+
ep.s(0,1,5)*ep.spab(1,ep.Sum(5,
   0),4))*ep.spab(5,ep.Sum(3,2),
4))-(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),
 5),2)*ep.spa(1,3)*ep.spb(2,0))/
(ep.spa(1,4)*ep.spa(3,4)*
 ep.spab(1,ep.Sum(3,4),2)*
 ep.spab(3,ep.Sum(1,4),2)*
 (-ep.spb(1,0)+(ep.spa(3,4)*
  ep.spb(4,0))/ep.spa(1,3))*
 ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43263_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43318_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43333_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43353_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43363_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43468_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43478_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43498_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,5))/
(ep.spa(0,5)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spa(2,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 ep.spa(1,4)*ep.spa(2,5))/
(ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(2,4)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43513_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g2Q43568_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.spab(5,ep.Sum(3,4),1))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q43645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q43670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q43680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,4),2))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))
); }

template <class T> complex<T> A2q2g2Q43683_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43688_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A2q2g2Q43710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,4),
2)*ep.spa(0,2)*ep.spa(3,4)*
 ep.spab(0,ep.Sum(3,4),5))/
(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.spa(4,5)+(ep.s(3,4,5)*
  ep.spa(0,4))/ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A2q2g2Q43713_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43748_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q43753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44072_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44077_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44102_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44112_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(3,
 ep.Sum(1,2),0),3)*pow(ep.spb(1,0),
2)*ep.spa(3,5)*ep.spb(2,0))/
(pow(ep.spab(3,ep.Sum(4,5),0),2)*
 ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))-
ep.spa(0,1)*ep.spab(3,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(3,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(3,ep.Sum(4,
   5),0))*ep.spab(5,ep.Sum(1,2),
0))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),
 0),2)*pow(ep.spb(2,0),2)*
 ep.spa(1,5))/(ep.spa(4,5)*
 ep.spab(1,ep.Sum(0,3),2)*
 ep.spab(1,ep.Sum(2,3),0)*
 ep.spab(5,ep.Sum(3,0),2)*
 ep.spb(3,0)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spa(4,5),2)*
 pow(ep.spab(3,ep.Sum(4,5),0),3)*
 pow(ep.spb(5,0),2)*ep.spa(1,3)*
 ep.spb(4,0))/(pow(ep.spab(3,ep.Sum(1,2),
 0),2)*ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(4,5),
0)*(ep.s(1,2,3)*ep.spab(3,
  ep.Sum(1,2),0)+
(-ep.s(4,5)+ep.s(1,2,3))*
 ep.spab(3,ep.Sum(4,5),0))*
 (ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(0,5)*ep.spab(3,ep.Sum(4,5),
  0)))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(4,2))/(ep.spa(0,3)*
 ep.spa(0,5)*ep.spab(3,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(3,0),2)*
 ep.spb(2,1)*(-(ep.spa(0,5)*
  ep.spb(4,0))-ep.spa(3,5)*
 ep.spb(4,3)))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 pow(ep.spb(4,0),2)*ep.spa(1,5))/
(ep.spa(1,2)*ep.spab(5,ep.Sum(3,4),
0)*ep.spb(3,0)*ep.spb(4,3)*
 (ep.spa(0,1)*ep.spb(4,0)-
ep.spa(1,3)*ep.spb(4,3))*
 (-(ep.spa(0,5)*ep.spb(4,0))-
ep.spa(3,5)*ep.spb(4,3)))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spab(3,ep.Sum(1,0),5),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spab(1,ep.Sum(3,0),
2)*ep.spab(3,ep.Sum(1,0),2)*
 (ep.spa(0,1)*ep.spb(4,0)-
ep.spa(1,3)*ep.spb(4,3))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spab(0,ep.Sum(2,1),4))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44137_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),3)*ep.spa(0,1)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(4,5)*(-(ep.s(3,4,5)*
  ep.spa(1,3))-ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),2))*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,1),
2)*ep.spa(0,1)*ep.spa(3,5))/
(ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))-
ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
  2))*ep.spab(5,ep.Sum(0,1),2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),5),2)*
 ep.spab(3,ep.Sum(2,1),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(3,ep.Sum(0,1),2),3)*
 pow(ep.spb(2,0),3)*ep.spa(3,5))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),
 2),3)*ep.spb(4,2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(3,4),
2)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 pow(ep.spb(2,1),2)*ep.spa(3,5)*
 ep.spab(3,ep.Sum(1,0),2)*
 ep.spb(2,0))/(pow(ep.spab(3,ep.Sum(4,5),
 2),2)*ep.s(3,4,5)*
 ep.spa(4,5)*(-(ep.s(3,4,5)*
  ep.spa(1,3))+ep.spa(1,2)*
 ep.spab(3,ep.Sum(0,1),2))*
 ((-ep.s(0,1)+ep.s(3,4,5))*
 ep.spab(3,ep.Sum(0,1),2)+
ep.s(3,4,5)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(1,0),
2))-(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),
 2),2)*ep.spab(1,ep.Sum(3,4),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3)*ep.spab(3,ep.Sum(1,2),
4))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44282_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44292_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),5),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),3)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),2),3)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(5,0),
2)*ep.spb(4,0))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))
); }

template <class T> complex<T> A2q2g2Q44370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(2,3),
2)*pow(ep.spab(0,ep.Sum(2,3),1),3)*
 pow(ep.spb(2,1),2)*ep.spb(3,1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (ep.s(0,4,5)*ep.spa(0,2)-
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(2,
ep.Sum(5,0),1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),5),3)*
 ep.spa(0,2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44373_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,3)*
 ep.spa(2,5))/(ep.spa(0,1)*
 ep.spa(0,3)*ep.spa(1,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 ep.spa(2,5)*ep.spa(3,5))/
(ep.spa(0,3)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(2,1),2)*ep.spa(2,3))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(1,2,3)*ep.spa(4,5)*
 (-(ep.s(1,2,3)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44403_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44408_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(5,ep.Sum(1,2),0)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5)*ep.spab(1,ep.Sum(5,0),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(1,
 ep.Sum(4,5),0),3)*pow(ep.spb(5,0),
2)*ep.spa(1,3)*ep.spb(4,0))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))
); }

template <class T> complex<T> A2q2g2Q44550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),2)*pow(ep.spb(3,1),
3)*ep.spab(0,ep.Sum(3,2),1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (-(ep.s(0,4,5)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(3,2),
1))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3)*ep.spab(2,
ep.Sum(5,0),1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44553_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),5),2)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44613_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44623_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(5,3),3)*ep.spa(0,2))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,0),2)*ep.spab(2,
ep.Sum(3,4),1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q44725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spa(4,5),
2)*pow(ep.spab(2,ep.Sum(4,5),3),3)*
 pow(ep.spb(5,3),3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(3,4,5)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(3,4,5))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(3,4,5)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(2,
ep.Sum(3,4),1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2g2Q44750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spa(2,3))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(1,2,3)*ep.spa(4,5)*
 (-(ep.s(1,2,3)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),3))/(ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44763_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44768_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 (ep.spa(1,4)*ep.spb(1,0)+
ep.spa(2,4)*ep.spb(2,0)))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2))/
(-(ep.s(0,4,5)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spb(1,0)*
 ep.spb(5,4))-pow(ep.spa(2,3),2)*
ep.s(0,4,5)*ep.spb(2,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2)*ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(0,ep.Sum(3,4),5),2)*
 pow(ep.spb(5,3),3)*ep.spa(0,2)*
 ep.spa(3,4))/(pow(ep.spab(0,ep.Sum(1,2),
 5),2)*ep.s(3,4,5)*
 ep.spa(1,2)*(-(ep.s(3,4,5)*
  ep.spa(0,4))-ep.spa(4,5)*
 ep.spab(0,ep.Sum(3,4),5))*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A2q2g2Q44793_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44803_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44828_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A2q2g2Q44833_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))

 ); }
#endif

template <class T> complex<T>  (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
#if _FAST_COMPILE_eval==1
case 317 : return &A2q2g2Q317_eval;
case 322 : return &A2q2g2Q322_eval;
case 377 : return &A2q2g2Q377_eval;
case 392 : return &A2q2g2Q392_eval;
case 412 : return &A2q2g2Q412_eval;
case 422 : return &A2q2g2Q422_eval;
case 497 : return &A2q2g2Q497_eval;
case 502 : return &A2q2g2Q502_eval;
case 587 : return &A2q2g2Q587_eval;
case 607 : return &A2q2g2Q607_eval;
case 622 : return &A2q2g2Q622_eval;
case 637 : return &A2q2g2Q637_eval;
case 917 : return &A2q2g2Q917_eval;
case 932 : return &A2q2g2Q932_eval;
case 947 : return &A2q2g2Q947_eval;
case 967 : return &A2q2g2Q967_eval;
case 1052 : return &A2q2g2Q1052_eval;
case 1057 : return &A2q2g2Q1057_eval;
case 1132 : return &A2q2g2Q1132_eval;
case 1142 : return &A2q2g2Q1142_eval;
case 1162 : return &A2q2g2Q1162_eval;
case 1177 : return &A2q2g2Q1177_eval;
case 1232 : return &A2q2g2Q1232_eval;
case 1237 : return &A2q2g2Q1237_eval;
case 1397 : return &A2q2g2Q1397_eval;
case 1402 : return &A2q2g2Q1402_eval;
case 1457 : return &A2q2g2Q1457_eval;
case 1472 : return &A2q2g2Q1472_eval;
case 1492 : return &A2q2g2Q1492_eval;
case 1502 : return &A2q2g2Q1502_eval;
case 1757 : return &A2q2g2Q1757_eval;
case 1762 : return &A2q2g2Q1762_eval;
case 1865 : return &A2q2g2Q1865_eval;
case 1870 : return &A2q2g2Q1870_eval;
case 1877 : return &A2q2g2Q1877_eval;
case 1895 : return &A2q2g2Q1895_eval;
case 1902 : return &A2q2g2Q1902_eval;
case 1905 : return &A2q2g2Q1905_eval;
case 1912 : return &A2q2g2Q1912_eval;
case 1930 : return &A2q2g2Q1930_eval;
case 1932 : return &A2q2g2Q1932_eval;
case 1935 : return &A2q2g2Q1935_eval;
case 2045 : return &A2q2g2Q2045_eval;
case 2050 : return &A2q2g2Q2050_eval;
case 2105 : return &A2q2g2Q2105_eval;
case 2120 : return &A2q2g2Q2120_eval;
case 2140 : return &A2q2g2Q2140_eval;
case 2150 : return &A2q2g2Q2150_eval;
case 2177 : return &A2q2g2Q2177_eval;
case 2192 : return &A2q2g2Q2192_eval;
case 2237 : return &A2q2g2Q2237_eval;
case 2255 : return &A2q2g2Q2255_eval;
case 2262 : return &A2q2g2Q2262_eval;
case 2265 : return &A2q2g2Q2265_eval;
case 2285 : return &A2q2g2Q2285_eval;
case 2300 : return &A2q2g2Q2300_eval;
case 2342 : return &A2q2g2Q2342_eval;
case 2352 : return &A2q2g2Q2352_eval;
case 2355 : return &A2q2g2Q2355_eval;
case 2360 : return &A2q2g2Q2360_eval;
case 2392 : return &A2q2g2Q2392_eval;
case 2402 : return &A2q2g2Q2402_eval;
case 2452 : return &A2q2g2Q2452_eval;
case 2470 : return &A2q2g2Q2470_eval;
case 2472 : return &A2q2g2Q2472_eval;
case 2475 : return &A2q2g2Q2475_eval;
case 2500 : return &A2q2g2Q2500_eval;
case 2510 : return &A2q2g2Q2510_eval;
case 2522 : return &A2q2g2Q2522_eval;
case 2532 : return &A2q2g2Q2532_eval;
case 2535 : return &A2q2g2Q2535_eval;
case 2540 : return &A2q2g2Q2540_eval;
case 2657 : return &A2q2g2Q2657_eval;
case 2662 : return &A2q2g2Q2662_eval;
case 2747 : return &A2q2g2Q2747_eval;
case 2767 : return &A2q2g2Q2767_eval;
case 2782 : return &A2q2g2Q2782_eval;
case 2797 : return &A2q2g2Q2797_eval;
case 2837 : return &A2q2g2Q2837_eval;
case 2842 : return &A2q2g2Q2842_eval;
case 2945 : return &A2q2g2Q2945_eval;
case 2950 : return &A2q2g2Q2950_eval;
case 2957 : return &A2q2g2Q2957_eval;
case 2975 : return &A2q2g2Q2975_eval;
case 2982 : return &A2q2g2Q2982_eval;
case 2985 : return &A2q2g2Q2985_eval;
case 2992 : return &A2q2g2Q2992_eval;
case 3010 : return &A2q2g2Q3010_eval;
case 3012 : return &A2q2g2Q3012_eval;
case 3015 : return &A2q2g2Q3015_eval;
case 3305 : return &A2q2g2Q3305_eval;
case 3310 : return &A2q2g2Q3310_eval;
case 3395 : return &A2q2g2Q3395_eval;
case 3415 : return &A2q2g2Q3415_eval;
case 3430 : return &A2q2g2Q3430_eval;
case 3445 : return &A2q2g2Q3445_eval;
case 3467 : return &A2q2g2Q3467_eval;
case 3487 : return &A2q2g2Q3487_eval;
case 3497 : return &A2q2g2Q3497_eval;
case 3515 : return &A2q2g2Q3515_eval;
case 3522 : return &A2q2g2Q3522_eval;
case 3525 : return &A2q2g2Q3525_eval;
case 3575 : return &A2q2g2Q3575_eval;
case 3595 : return &A2q2g2Q3595_eval;
case 3637 : return &A2q2g2Q3637_eval;
case 3642 : return &A2q2g2Q3642_eval;
case 3645 : return &A2q2g2Q3645_eval;
case 3655 : return &A2q2g2Q3655_eval;
case 3682 : return &A2q2g2Q3682_eval;
case 3697 : return &A2q2g2Q3697_eval;
case 3712 : return &A2q2g2Q3712_eval;
case 3730 : return &A2q2g2Q3730_eval;
case 3732 : return &A2q2g2Q3732_eval;
case 3735 : return &A2q2g2Q3735_eval;
case 3790 : return &A2q2g2Q3790_eval;
case 3805 : return &A2q2g2Q3805_eval;
case 3817 : return &A2q2g2Q3817_eval;
case 3822 : return &A2q2g2Q3822_eval;
case 3825 : return &A2q2g2Q3825_eval;
case 3835 : return &A2q2g2Q3835_eval;
case 4205 : return &A2q2g2Q4205_eval;
case 4210 : return &A2q2g2Q4210_eval;
case 4265 : return &A2q2g2Q4265_eval;
case 4280 : return &A2q2g2Q4280_eval;
case 4300 : return &A2q2g2Q4300_eval;
case 4310 : return &A2q2g2Q4310_eval;
case 4385 : return &A2q2g2Q4385_eval;
case 4390 : return &A2q2g2Q4390_eval;
case 4475 : return &A2q2g2Q4475_eval;
case 4495 : return &A2q2g2Q4495_eval;
case 4510 : return &A2q2g2Q4510_eval;
case 4525 : return &A2q2g2Q4525_eval;
case 4805 : return &A2q2g2Q4805_eval;
case 4820 : return &A2q2g2Q4820_eval;
case 4835 : return &A2q2g2Q4835_eval;
case 4855 : return &A2q2g2Q4855_eval;
case 4940 : return &A2q2g2Q4940_eval;
case 4945 : return &A2q2g2Q4945_eval;
case 5020 : return &A2q2g2Q5020_eval;
case 5030 : return &A2q2g2Q5030_eval;
case 5050 : return &A2q2g2Q5050_eval;
case 5065 : return &A2q2g2Q5065_eval;
case 5120 : return &A2q2g2Q5120_eval;
case 5125 : return &A2q2g2Q5125_eval;
case 5237 : return &A2q2g2Q5237_eval;
case 5252 : return &A2q2g2Q5252_eval;
case 5267 : return &A2q2g2Q5267_eval;
case 5287 : return &A2q2g2Q5287_eval;
case 5372 : return &A2q2g2Q5372_eval;
case 5377 : return &A2q2g2Q5377_eval;
case 5417 : return &A2q2g2Q5417_eval;
case 5432 : return &A2q2g2Q5432_eval;
case 5477 : return &A2q2g2Q5477_eval;
case 5495 : return &A2q2g2Q5495_eval;
case 5502 : return &A2q2g2Q5502_eval;
case 5505 : return &A2q2g2Q5505_eval;
case 5525 : return &A2q2g2Q5525_eval;
case 5540 : return &A2q2g2Q5540_eval;
case 5582 : return &A2q2g2Q5582_eval;
case 5592 : return &A2q2g2Q5592_eval;
case 5595 : return &A2q2g2Q5595_eval;
case 5600 : return &A2q2g2Q5600_eval;
case 5627 : return &A2q2g2Q5627_eval;
case 5647 : return &A2q2g2Q5647_eval;
case 5657 : return &A2q2g2Q5657_eval;
case 5675 : return &A2q2g2Q5675_eval;
case 5682 : return &A2q2g2Q5682_eval;
case 5685 : return &A2q2g2Q5685_eval;
case 5735 : return &A2q2g2Q5735_eval;
case 5755 : return &A2q2g2Q5755_eval;
case 5797 : return &A2q2g2Q5797_eval;
case 5802 : return &A2q2g2Q5802_eval;
case 5805 : return &A2q2g2Q5805_eval;
case 5815 : return &A2q2g2Q5815_eval;
case 5885 : return &A2q2g2Q5885_eval;
case 5900 : return &A2q2g2Q5900_eval;
case 5915 : return &A2q2g2Q5915_eval;
case 5935 : return &A2q2g2Q5935_eval;
case 6020 : return &A2q2g2Q6020_eval;
case 6025 : return &A2q2g2Q6025_eval;
case 6272 : return &A2q2g2Q6272_eval;
case 6277 : return &A2q2g2Q6277_eval;
case 6302 : return &A2q2g2Q6302_eval;
case 6312 : return &A2q2g2Q6312_eval;
case 6315 : return &A2q2g2Q6315_eval;
case 6320 : return &A2q2g2Q6320_eval;
case 6337 : return &A2q2g2Q6337_eval;
case 6342 : return &A2q2g2Q6342_eval;
case 6345 : return &A2q2g2Q6345_eval;
case 6355 : return &A2q2g2Q6355_eval;
case 6380 : return &A2q2g2Q6380_eval;
case 6385 : return &A2q2g2Q6385_eval;
case 6532 : return &A2q2g2Q6532_eval;
case 6542 : return &A2q2g2Q6542_eval;
case 6562 : return &A2q2g2Q6562_eval;
case 6577 : return &A2q2g2Q6577_eval;
case 6632 : return &A2q2g2Q6632_eval;
case 6637 : return &A2q2g2Q6637_eval;
case 6712 : return &A2q2g2Q6712_eval;
case 6722 : return &A2q2g2Q6722_eval;
case 6772 : return &A2q2g2Q6772_eval;
case 6790 : return &A2q2g2Q6790_eval;
case 6792 : return &A2q2g2Q6792_eval;
case 6795 : return &A2q2g2Q6795_eval;
case 6820 : return &A2q2g2Q6820_eval;
case 6830 : return &A2q2g2Q6830_eval;
case 6842 : return &A2q2g2Q6842_eval;
case 6852 : return &A2q2g2Q6852_eval;
case 6855 : return &A2q2g2Q6855_eval;
case 6860 : return &A2q2g2Q6860_eval;
case 6922 : return &A2q2g2Q6922_eval;
case 6937 : return &A2q2g2Q6937_eval;
case 6952 : return &A2q2g2Q6952_eval;
case 6970 : return &A2q2g2Q6970_eval;
case 6972 : return &A2q2g2Q6972_eval;
case 6975 : return &A2q2g2Q6975_eval;
case 7030 : return &A2q2g2Q7030_eval;
case 7045 : return &A2q2g2Q7045_eval;
case 7057 : return &A2q2g2Q7057_eval;
case 7062 : return &A2q2g2Q7062_eval;
case 7065 : return &A2q2g2Q7065_eval;
case 7075 : return &A2q2g2Q7075_eval;
case 7180 : return &A2q2g2Q7180_eval;
case 7190 : return &A2q2g2Q7190_eval;
case 7210 : return &A2q2g2Q7210_eval;
case 7225 : return &A2q2g2Q7225_eval;
case 7280 : return &A2q2g2Q7280_eval;
case 7285 : return &A2q2g2Q7285_eval;
case 7352 : return &A2q2g2Q7352_eval;
case 7357 : return &A2q2g2Q7357_eval;
case 7382 : return &A2q2g2Q7382_eval;
case 7392 : return &A2q2g2Q7392_eval;
case 7395 : return &A2q2g2Q7395_eval;
case 7400 : return &A2q2g2Q7400_eval;
case 7417 : return &A2q2g2Q7417_eval;
case 7422 : return &A2q2g2Q7422_eval;
case 7425 : return &A2q2g2Q7425_eval;
case 7435 : return &A2q2g2Q7435_eval;
case 7460 : return &A2q2g2Q7460_eval;
case 7465 : return &A2q2g2Q7465_eval;
case 7877 : return &A2q2g2Q7877_eval;
case 7882 : return &A2q2g2Q7882_eval;
case 7937 : return &A2q2g2Q7937_eval;
case 7952 : return &A2q2g2Q7952_eval;
case 7972 : return &A2q2g2Q7972_eval;
case 7982 : return &A2q2g2Q7982_eval;
case 8237 : return &A2q2g2Q8237_eval;
case 8242 : return &A2q2g2Q8242_eval;
case 8345 : return &A2q2g2Q8345_eval;
case 8350 : return &A2q2g2Q8350_eval;
case 8357 : return &A2q2g2Q8357_eval;
case 8375 : return &A2q2g2Q8375_eval;
case 8382 : return &A2q2g2Q8382_eval;
case 8385 : return &A2q2g2Q8385_eval;
case 8392 : return &A2q2g2Q8392_eval;
case 8410 : return &A2q2g2Q8410_eval;
case 8412 : return &A2q2g2Q8412_eval;
case 8415 : return &A2q2g2Q8415_eval;
case 8525 : return &A2q2g2Q8525_eval;
case 8530 : return &A2q2g2Q8530_eval;
case 8585 : return &A2q2g2Q8585_eval;
case 8600 : return &A2q2g2Q8600_eval;
case 8620 : return &A2q2g2Q8620_eval;
case 8630 : return &A2q2g2Q8630_eval;
case 8657 : return &A2q2g2Q8657_eval;
case 8672 : return &A2q2g2Q8672_eval;
case 8717 : return &A2q2g2Q8717_eval;
case 8735 : return &A2q2g2Q8735_eval;
case 8742 : return &A2q2g2Q8742_eval;
case 8745 : return &A2q2g2Q8745_eval;
case 8765 : return &A2q2g2Q8765_eval;
case 8780 : return &A2q2g2Q8780_eval;
case 8822 : return &A2q2g2Q8822_eval;
case 8832 : return &A2q2g2Q8832_eval;
case 8835 : return &A2q2g2Q8835_eval;
case 8840 : return &A2q2g2Q8840_eval;
case 8872 : return &A2q2g2Q8872_eval;
case 8882 : return &A2q2g2Q8882_eval;
case 8932 : return &A2q2g2Q8932_eval;
case 8950 : return &A2q2g2Q8950_eval;
case 8952 : return &A2q2g2Q8952_eval;
case 8955 : return &A2q2g2Q8955_eval;
case 8980 : return &A2q2g2Q8980_eval;
case 8990 : return &A2q2g2Q8990_eval;
case 9002 : return &A2q2g2Q9002_eval;
case 9012 : return &A2q2g2Q9012_eval;
case 9015 : return &A2q2g2Q9015_eval;
case 9020 : return &A2q2g2Q9020_eval;
case 10397 : return &A2q2g2Q10397_eval;
case 10402 : return &A2q2g2Q10402_eval;
case 10505 : return &A2q2g2Q10505_eval;
case 10510 : return &A2q2g2Q10510_eval;
case 10517 : return &A2q2g2Q10517_eval;
case 10535 : return &A2q2g2Q10535_eval;
case 10542 : return &A2q2g2Q10542_eval;
case 10545 : return &A2q2g2Q10545_eval;
case 10552 : return &A2q2g2Q10552_eval;
case 10570 : return &A2q2g2Q10570_eval;
case 10572 : return &A2q2g2Q10572_eval;
case 10575 : return &A2q2g2Q10575_eval;
case 11045 : return &A2q2g2Q11045_eval;
case 11050 : return &A2q2g2Q11050_eval;
case 11153 : return &A2q2g2Q11153_eval;
case 11158 : return &A2q2g2Q11158_eval;
case 11165 : return &A2q2g2Q11165_eval;
case 11183 : return &A2q2g2Q11183_eval;
case 11190 : return &A2q2g2Q11190_eval;
case 11193 : return &A2q2g2Q11193_eval;
case 11200 : return &A2q2g2Q11200_eval;
case 11218 : return &A2q2g2Q11218_eval;
case 11220 : return &A2q2g2Q11220_eval;
case 11223 : return &A2q2g2Q11223_eval;
case 11237 : return &A2q2g2Q11237_eval;
case 11255 : return &A2q2g2Q11255_eval;
case 11262 : return &A2q2g2Q11262_eval;
case 11265 : return &A2q2g2Q11265_eval;
case 11345 : return &A2q2g2Q11345_eval;
case 11363 : return &A2q2g2Q11363_eval;
case 11370 : return &A2q2g2Q11370_eval;
case 11373 : return &A2q2g2Q11373_eval;
case 11412 : return &A2q2g2Q11412_eval;
case 11415 : return &A2q2g2Q11415_eval;
case 11430 : return &A2q2g2Q11430_eval;
case 11433 : return &A2q2g2Q11433_eval;
case 11452 : return &A2q2g2Q11452_eval;
case 11470 : return &A2q2g2Q11470_eval;
case 11472 : return &A2q2g2Q11472_eval;
case 11475 : return &A2q2g2Q11475_eval;
case 11560 : return &A2q2g2Q11560_eval;
case 11578 : return &A2q2g2Q11578_eval;
case 11580 : return &A2q2g2Q11580_eval;
case 11583 : return &A2q2g2Q11583_eval;
case 11592 : return &A2q2g2Q11592_eval;
case 11595 : return &A2q2g2Q11595_eval;
case 11610 : return &A2q2g2Q11610_eval;
case 11613 : return &A2q2g2Q11613_eval;
case 11765 : return &A2q2g2Q11765_eval;
case 11770 : return &A2q2g2Q11770_eval;
case 11825 : return &A2q2g2Q11825_eval;
case 11840 : return &A2q2g2Q11840_eval;
case 11860 : return &A2q2g2Q11860_eval;
case 11870 : return &A2q2g2Q11870_eval;
case 12125 : return &A2q2g2Q12125_eval;
case 12130 : return &A2q2g2Q12130_eval;
case 12233 : return &A2q2g2Q12233_eval;
case 12238 : return &A2q2g2Q12238_eval;
case 12245 : return &A2q2g2Q12245_eval;
case 12263 : return &A2q2g2Q12263_eval;
case 12270 : return &A2q2g2Q12270_eval;
case 12273 : return &A2q2g2Q12273_eval;
case 12280 : return &A2q2g2Q12280_eval;
case 12298 : return &A2q2g2Q12298_eval;
case 12300 : return &A2q2g2Q12300_eval;
case 12303 : return &A2q2g2Q12303_eval;
case 12413 : return &A2q2g2Q12413_eval;
case 12418 : return &A2q2g2Q12418_eval;
case 12473 : return &A2q2g2Q12473_eval;
case 12488 : return &A2q2g2Q12488_eval;
case 12508 : return &A2q2g2Q12508_eval;
case 12518 : return &A2q2g2Q12518_eval;
case 12545 : return &A2q2g2Q12545_eval;
case 12560 : return &A2q2g2Q12560_eval;
case 12605 : return &A2q2g2Q12605_eval;
case 12623 : return &A2q2g2Q12623_eval;
case 12630 : return &A2q2g2Q12630_eval;
case 12633 : return &A2q2g2Q12633_eval;
case 12653 : return &A2q2g2Q12653_eval;
case 12668 : return &A2q2g2Q12668_eval;
case 12710 : return &A2q2g2Q12710_eval;
case 12720 : return &A2q2g2Q12720_eval;
case 12723 : return &A2q2g2Q12723_eval;
case 12728 : return &A2q2g2Q12728_eval;
case 12760 : return &A2q2g2Q12760_eval;
case 12770 : return &A2q2g2Q12770_eval;
case 12820 : return &A2q2g2Q12820_eval;
case 12838 : return &A2q2g2Q12838_eval;
case 12840 : return &A2q2g2Q12840_eval;
case 12843 : return &A2q2g2Q12843_eval;
case 12868 : return &A2q2g2Q12868_eval;
case 12878 : return &A2q2g2Q12878_eval;
case 12890 : return &A2q2g2Q12890_eval;
case 12900 : return &A2q2g2Q12900_eval;
case 12903 : return &A2q2g2Q12903_eval;
case 12908 : return &A2q2g2Q12908_eval;
case 12977 : return &A2q2g2Q12977_eval;
case 12992 : return &A2q2g2Q12992_eval;
case 13037 : return &A2q2g2Q13037_eval;
case 13055 : return &A2q2g2Q13055_eval;
case 13062 : return &A2q2g2Q13062_eval;
case 13065 : return &A2q2g2Q13065_eval;
case 13085 : return &A2q2g2Q13085_eval;
case 13100 : return &A2q2g2Q13100_eval;
case 13142 : return &A2q2g2Q13142_eval;
case 13152 : return &A2q2g2Q13152_eval;
case 13155 : return &A2q2g2Q13155_eval;
case 13160 : return &A2q2g2Q13160_eval;
case 13397 : return &A2q2g2Q13397_eval;
case 13415 : return &A2q2g2Q13415_eval;
case 13422 : return &A2q2g2Q13422_eval;
case 13425 : return &A2q2g2Q13425_eval;
case 13505 : return &A2q2g2Q13505_eval;
case 13523 : return &A2q2g2Q13523_eval;
case 13530 : return &A2q2g2Q13530_eval;
case 13533 : return &A2q2g2Q13533_eval;
case 13572 : return &A2q2g2Q13572_eval;
case 13575 : return &A2q2g2Q13575_eval;
case 13590 : return &A2q2g2Q13590_eval;
case 13593 : return &A2q2g2Q13593_eval;
case 13625 : return &A2q2g2Q13625_eval;
case 13640 : return &A2q2g2Q13640_eval;
case 13685 : return &A2q2g2Q13685_eval;
case 13703 : return &A2q2g2Q13703_eval;
case 13710 : return &A2q2g2Q13710_eval;
case 13713 : return &A2q2g2Q13713_eval;
case 13733 : return &A2q2g2Q13733_eval;
case 13748 : return &A2q2g2Q13748_eval;
case 13790 : return &A2q2g2Q13790_eval;
case 13800 : return &A2q2g2Q13800_eval;
case 13803 : return &A2q2g2Q13803_eval;
case 13808 : return &A2q2g2Q13808_eval;
case 14042 : return &A2q2g2Q14042_eval;
case 14052 : return &A2q2g2Q14052_eval;
case 14055 : return &A2q2g2Q14055_eval;
case 14060 : return &A2q2g2Q14060_eval;
case 14112 : return &A2q2g2Q14112_eval;
case 14115 : return &A2q2g2Q14115_eval;
case 14130 : return &A2q2g2Q14130_eval;
case 14133 : return &A2q2g2Q14133_eval;
case 14150 : return &A2q2g2Q14150_eval;
case 14160 : return &A2q2g2Q14160_eval;
case 14163 : return &A2q2g2Q14163_eval;
case 14168 : return &A2q2g2Q14168_eval;
case 14272 : return &A2q2g2Q14272_eval;
case 14282 : return &A2q2g2Q14282_eval;
case 14332 : return &A2q2g2Q14332_eval;
case 14350 : return &A2q2g2Q14350_eval;
case 14352 : return &A2q2g2Q14352_eval;
case 14355 : return &A2q2g2Q14355_eval;
case 14380 : return &A2q2g2Q14380_eval;
case 14390 : return &A2q2g2Q14390_eval;
case 14402 : return &A2q2g2Q14402_eval;
case 14412 : return &A2q2g2Q14412_eval;
case 14415 : return &A2q2g2Q14415_eval;
case 14420 : return &A2q2g2Q14420_eval;
case 14692 : return &A2q2g2Q14692_eval;
case 14710 : return &A2q2g2Q14710_eval;
case 14712 : return &A2q2g2Q14712_eval;
case 14715 : return &A2q2g2Q14715_eval;
case 14800 : return &A2q2g2Q14800_eval;
case 14818 : return &A2q2g2Q14818_eval;
case 14820 : return &A2q2g2Q14820_eval;
case 14823 : return &A2q2g2Q14823_eval;
case 14832 : return &A2q2g2Q14832_eval;
case 14835 : return &A2q2g2Q14835_eval;
case 14850 : return &A2q2g2Q14850_eval;
case 14853 : return &A2q2g2Q14853_eval;
case 14920 : return &A2q2g2Q14920_eval;
case 14930 : return &A2q2g2Q14930_eval;
case 14980 : return &A2q2g2Q14980_eval;
case 14998 : return &A2q2g2Q14998_eval;
case 15000 : return &A2q2g2Q15000_eval;
case 15003 : return &A2q2g2Q15003_eval;
case 15028 : return &A2q2g2Q15028_eval;
case 15038 : return &A2q2g2Q15038_eval;
case 15050 : return &A2q2g2Q15050_eval;
case 15060 : return &A2q2g2Q15060_eval;
case 15063 : return &A2q2g2Q15063_eval;
case 15068 : return &A2q2g2Q15068_eval;
case 15122 : return &A2q2g2Q15122_eval;
case 15132 : return &A2q2g2Q15132_eval;
case 15135 : return &A2q2g2Q15135_eval;
case 15140 : return &A2q2g2Q15140_eval;
case 15192 : return &A2q2g2Q15192_eval;
case 15195 : return &A2q2g2Q15195_eval;
case 15210 : return &A2q2g2Q15210_eval;
case 15213 : return &A2q2g2Q15213_eval;
case 15230 : return &A2q2g2Q15230_eval;
case 15240 : return &A2q2g2Q15240_eval;
case 15243 : return &A2q2g2Q15243_eval;
case 15248 : return &A2q2g2Q15248_eval;
case 15617 : return &A2q2g2Q15617_eval;
case 15622 : return &A2q2g2Q15622_eval;
case 15707 : return &A2q2g2Q15707_eval;
case 15727 : return &A2q2g2Q15727_eval;
case 15742 : return &A2q2g2Q15742_eval;
case 15757 : return &A2q2g2Q15757_eval;
case 15797 : return &A2q2g2Q15797_eval;
case 15802 : return &A2q2g2Q15802_eval;
case 15905 : return &A2q2g2Q15905_eval;
case 15910 : return &A2q2g2Q15910_eval;
case 15917 : return &A2q2g2Q15917_eval;
case 15935 : return &A2q2g2Q15935_eval;
case 15942 : return &A2q2g2Q15942_eval;
case 15945 : return &A2q2g2Q15945_eval;
case 15952 : return &A2q2g2Q15952_eval;
case 15970 : return &A2q2g2Q15970_eval;
case 15972 : return &A2q2g2Q15972_eval;
case 15975 : return &A2q2g2Q15975_eval;
case 16265 : return &A2q2g2Q16265_eval;
case 16270 : return &A2q2g2Q16270_eval;
case 16355 : return &A2q2g2Q16355_eval;
case 16375 : return &A2q2g2Q16375_eval;
case 16390 : return &A2q2g2Q16390_eval;
case 16405 : return &A2q2g2Q16405_eval;
case 16427 : return &A2q2g2Q16427_eval;
case 16447 : return &A2q2g2Q16447_eval;
case 16457 : return &A2q2g2Q16457_eval;
case 16475 : return &A2q2g2Q16475_eval;
case 16482 : return &A2q2g2Q16482_eval;
case 16485 : return &A2q2g2Q16485_eval;
case 16535 : return &A2q2g2Q16535_eval;
case 16555 : return &A2q2g2Q16555_eval;
case 16597 : return &A2q2g2Q16597_eval;
case 16602 : return &A2q2g2Q16602_eval;
case 16605 : return &A2q2g2Q16605_eval;
case 16615 : return &A2q2g2Q16615_eval;
case 16642 : return &A2q2g2Q16642_eval;
case 16657 : return &A2q2g2Q16657_eval;
case 16672 : return &A2q2g2Q16672_eval;
case 16690 : return &A2q2g2Q16690_eval;
case 16692 : return &A2q2g2Q16692_eval;
case 16695 : return &A2q2g2Q16695_eval;
case 16750 : return &A2q2g2Q16750_eval;
case 16765 : return &A2q2g2Q16765_eval;
case 16777 : return &A2q2g2Q16777_eval;
case 16782 : return &A2q2g2Q16782_eval;
case 16785 : return &A2q2g2Q16785_eval;
case 16795 : return &A2q2g2Q16795_eval;
case 16877 : return &A2q2g2Q16877_eval;
case 16882 : return &A2q2g2Q16882_eval;
case 16985 : return &A2q2g2Q16985_eval;
case 16990 : return &A2q2g2Q16990_eval;
case 16997 : return &A2q2g2Q16997_eval;
case 17015 : return &A2q2g2Q17015_eval;
case 17022 : return &A2q2g2Q17022_eval;
case 17025 : return &A2q2g2Q17025_eval;
case 17032 : return &A2q2g2Q17032_eval;
case 17050 : return &A2q2g2Q17050_eval;
case 17052 : return &A2q2g2Q17052_eval;
case 17055 : return &A2q2g2Q17055_eval;
case 17525 : return &A2q2g2Q17525_eval;
case 17530 : return &A2q2g2Q17530_eval;
case 17633 : return &A2q2g2Q17633_eval;
case 17638 : return &A2q2g2Q17638_eval;
case 17645 : return &A2q2g2Q17645_eval;
case 17663 : return &A2q2g2Q17663_eval;
case 17670 : return &A2q2g2Q17670_eval;
case 17673 : return &A2q2g2Q17673_eval;
case 17680 : return &A2q2g2Q17680_eval;
case 17698 : return &A2q2g2Q17698_eval;
case 17700 : return &A2q2g2Q17700_eval;
case 17703 : return &A2q2g2Q17703_eval;
case 17717 : return &A2q2g2Q17717_eval;
case 17735 : return &A2q2g2Q17735_eval;
case 17742 : return &A2q2g2Q17742_eval;
case 17745 : return &A2q2g2Q17745_eval;
case 17825 : return &A2q2g2Q17825_eval;
case 17843 : return &A2q2g2Q17843_eval;
case 17850 : return &A2q2g2Q17850_eval;
case 17853 : return &A2q2g2Q17853_eval;
case 17892 : return &A2q2g2Q17892_eval;
case 17895 : return &A2q2g2Q17895_eval;
case 17910 : return &A2q2g2Q17910_eval;
case 17913 : return &A2q2g2Q17913_eval;
case 17932 : return &A2q2g2Q17932_eval;
case 17950 : return &A2q2g2Q17950_eval;
case 17952 : return &A2q2g2Q17952_eval;
case 17955 : return &A2q2g2Q17955_eval;
case 18040 : return &A2q2g2Q18040_eval;
case 18058 : return &A2q2g2Q18058_eval;
case 18060 : return &A2q2g2Q18060_eval;
case 18063 : return &A2q2g2Q18063_eval;
case 18072 : return &A2q2g2Q18072_eval;
case 18075 : return &A2q2g2Q18075_eval;
case 18090 : return &A2q2g2Q18090_eval;
case 18093 : return &A2q2g2Q18093_eval;
case 19505 : return &A2q2g2Q19505_eval;
case 19510 : return &A2q2g2Q19510_eval;
case 19595 : return &A2q2g2Q19595_eval;
case 19615 : return &A2q2g2Q19615_eval;
case 19630 : return &A2q2g2Q19630_eval;
case 19645 : return &A2q2g2Q19645_eval;
case 19685 : return &A2q2g2Q19685_eval;
case 19690 : return &A2q2g2Q19690_eval;
case 19793 : return &A2q2g2Q19793_eval;
case 19798 : return &A2q2g2Q19798_eval;
case 19805 : return &A2q2g2Q19805_eval;
case 19823 : return &A2q2g2Q19823_eval;
case 19830 : return &A2q2g2Q19830_eval;
case 19833 : return &A2q2g2Q19833_eval;
case 19840 : return &A2q2g2Q19840_eval;
case 19858 : return &A2q2g2Q19858_eval;
case 19860 : return &A2q2g2Q19860_eval;
case 19863 : return &A2q2g2Q19863_eval;
case 20153 : return &A2q2g2Q20153_eval;
case 20158 : return &A2q2g2Q20158_eval;
case 20243 : return &A2q2g2Q20243_eval;
case 20263 : return &A2q2g2Q20263_eval;
case 20278 : return &A2q2g2Q20278_eval;
case 20293 : return &A2q2g2Q20293_eval;
case 20315 : return &A2q2g2Q20315_eval;
case 20335 : return &A2q2g2Q20335_eval;
case 20345 : return &A2q2g2Q20345_eval;
case 20363 : return &A2q2g2Q20363_eval;
case 20370 : return &A2q2g2Q20370_eval;
case 20373 : return &A2q2g2Q20373_eval;
case 20423 : return &A2q2g2Q20423_eval;
case 20443 : return &A2q2g2Q20443_eval;
case 20485 : return &A2q2g2Q20485_eval;
case 20490 : return &A2q2g2Q20490_eval;
case 20493 : return &A2q2g2Q20493_eval;
case 20503 : return &A2q2g2Q20503_eval;
case 20530 : return &A2q2g2Q20530_eval;
case 20545 : return &A2q2g2Q20545_eval;
case 20560 : return &A2q2g2Q20560_eval;
case 20578 : return &A2q2g2Q20578_eval;
case 20580 : return &A2q2g2Q20580_eval;
case 20583 : return &A2q2g2Q20583_eval;
case 20638 : return &A2q2g2Q20638_eval;
case 20653 : return &A2q2g2Q20653_eval;
case 20665 : return &A2q2g2Q20665_eval;
case 20670 : return &A2q2g2Q20670_eval;
case 20673 : return &A2q2g2Q20673_eval;
case 20683 : return &A2q2g2Q20683_eval;
case 20747 : return &A2q2g2Q20747_eval;
case 20767 : return &A2q2g2Q20767_eval;
case 20777 : return &A2q2g2Q20777_eval;
case 20795 : return &A2q2g2Q20795_eval;
case 20802 : return &A2q2g2Q20802_eval;
case 20805 : return &A2q2g2Q20805_eval;
case 20855 : return &A2q2g2Q20855_eval;
case 20875 : return &A2q2g2Q20875_eval;
case 20917 : return &A2q2g2Q20917_eval;
case 20922 : return &A2q2g2Q20922_eval;
case 20925 : return &A2q2g2Q20925_eval;
case 20935 : return &A2q2g2Q20935_eval;
case 20957 : return &A2q2g2Q20957_eval;
case 20975 : return &A2q2g2Q20975_eval;
case 20982 : return &A2q2g2Q20982_eval;
case 20985 : return &A2q2g2Q20985_eval;
case 21065 : return &A2q2g2Q21065_eval;
case 21083 : return &A2q2g2Q21083_eval;
case 21090 : return &A2q2g2Q21090_eval;
case 21093 : return &A2q2g2Q21093_eval;
case 21132 : return &A2q2g2Q21132_eval;
case 21135 : return &A2q2g2Q21135_eval;
case 21150 : return &A2q2g2Q21150_eval;
case 21153 : return &A2q2g2Q21153_eval;
case 21395 : return &A2q2g2Q21395_eval;
case 21415 : return &A2q2g2Q21415_eval;
case 21425 : return &A2q2g2Q21425_eval;
case 21443 : return &A2q2g2Q21443_eval;
case 21450 : return &A2q2g2Q21450_eval;
case 21453 : return &A2q2g2Q21453_eval;
case 21503 : return &A2q2g2Q21503_eval;
case 21523 : return &A2q2g2Q21523_eval;
case 21565 : return &A2q2g2Q21565_eval;
case 21570 : return &A2q2g2Q21570_eval;
case 21573 : return &A2q2g2Q21573_eval;
case 21583 : return &A2q2g2Q21583_eval;
case 21817 : return &A2q2g2Q21817_eval;
case 21822 : return &A2q2g2Q21822_eval;
case 21825 : return &A2q2g2Q21825_eval;
case 21835 : return &A2q2g2Q21835_eval;
case 21852 : return &A2q2g2Q21852_eval;
case 21855 : return &A2q2g2Q21855_eval;
case 21870 : return &A2q2g2Q21870_eval;
case 21873 : return &A2q2g2Q21873_eval;
case 21925 : return &A2q2g2Q21925_eval;
case 21930 : return &A2q2g2Q21930_eval;
case 21933 : return &A2q2g2Q21933_eval;
case 21943 : return &A2q2g2Q21943_eval;
case 22042 : return &A2q2g2Q22042_eval;
case 22057 : return &A2q2g2Q22057_eval;
case 22072 : return &A2q2g2Q22072_eval;
case 22090 : return &A2q2g2Q22090_eval;
case 22092 : return &A2q2g2Q22092_eval;
case 22095 : return &A2q2g2Q22095_eval;
case 22150 : return &A2q2g2Q22150_eval;
case 22165 : return &A2q2g2Q22165_eval;
case 22177 : return &A2q2g2Q22177_eval;
case 22182 : return &A2q2g2Q22182_eval;
case 22185 : return &A2q2g2Q22185_eval;
case 22195 : return &A2q2g2Q22195_eval;
case 22252 : return &A2q2g2Q22252_eval;
case 22270 : return &A2q2g2Q22270_eval;
case 22272 : return &A2q2g2Q22272_eval;
case 22275 : return &A2q2g2Q22275_eval;
case 22360 : return &A2q2g2Q22360_eval;
case 22378 : return &A2q2g2Q22378_eval;
case 22380 : return &A2q2g2Q22380_eval;
case 22383 : return &A2q2g2Q22383_eval;
case 22392 : return &A2q2g2Q22392_eval;
case 22395 : return &A2q2g2Q22395_eval;
case 22410 : return &A2q2g2Q22410_eval;
case 22413 : return &A2q2g2Q22413_eval;
case 22690 : return &A2q2g2Q22690_eval;
case 22705 : return &A2q2g2Q22705_eval;
case 22720 : return &A2q2g2Q22720_eval;
case 22738 : return &A2q2g2Q22738_eval;
case 22740 : return &A2q2g2Q22740_eval;
case 22743 : return &A2q2g2Q22743_eval;
case 22798 : return &A2q2g2Q22798_eval;
case 22813 : return &A2q2g2Q22813_eval;
case 22825 : return &A2q2g2Q22825_eval;
case 22830 : return &A2q2g2Q22830_eval;
case 22833 : return &A2q2g2Q22833_eval;
case 22843 : return &A2q2g2Q22843_eval;
case 22897 : return &A2q2g2Q22897_eval;
case 22902 : return &A2q2g2Q22902_eval;
case 22905 : return &A2q2g2Q22905_eval;
case 22915 : return &A2q2g2Q22915_eval;
case 22932 : return &A2q2g2Q22932_eval;
case 22935 : return &A2q2g2Q22935_eval;
case 22950 : return &A2q2g2Q22950_eval;
case 22953 : return &A2q2g2Q22953_eval;
case 23005 : return &A2q2g2Q23005_eval;
case 23010 : return &A2q2g2Q23010_eval;
case 23013 : return &A2q2g2Q23013_eval;
case 23023 : return &A2q2g2Q23023_eval;
case 23645 : return &A2q2g2Q23645_eval;
case 23650 : return &A2q2g2Q23650_eval;
case 23705 : return &A2q2g2Q23705_eval;
case 23720 : return &A2q2g2Q23720_eval;
case 23740 : return &A2q2g2Q23740_eval;
case 23750 : return &A2q2g2Q23750_eval;
case 23825 : return &A2q2g2Q23825_eval;
case 23830 : return &A2q2g2Q23830_eval;
case 23915 : return &A2q2g2Q23915_eval;
case 23935 : return &A2q2g2Q23935_eval;
case 23950 : return &A2q2g2Q23950_eval;
case 23965 : return &A2q2g2Q23965_eval;
case 24245 : return &A2q2g2Q24245_eval;
case 24260 : return &A2q2g2Q24260_eval;
case 24275 : return &A2q2g2Q24275_eval;
case 24295 : return &A2q2g2Q24295_eval;
case 24380 : return &A2q2g2Q24380_eval;
case 24385 : return &A2q2g2Q24385_eval;
case 24460 : return &A2q2g2Q24460_eval;
case 24470 : return &A2q2g2Q24470_eval;
case 24490 : return &A2q2g2Q24490_eval;
case 24505 : return &A2q2g2Q24505_eval;
case 24560 : return &A2q2g2Q24560_eval;
case 24565 : return &A2q2g2Q24565_eval;
case 24725 : return &A2q2g2Q24725_eval;
case 24730 : return &A2q2g2Q24730_eval;
case 24785 : return &A2q2g2Q24785_eval;
case 24800 : return &A2q2g2Q24800_eval;
case 24820 : return &A2q2g2Q24820_eval;
case 24830 : return &A2q2g2Q24830_eval;
case 25085 : return &A2q2g2Q25085_eval;
case 25090 : return &A2q2g2Q25090_eval;
case 25193 : return &A2q2g2Q25193_eval;
case 25198 : return &A2q2g2Q25198_eval;
case 25205 : return &A2q2g2Q25205_eval;
case 25223 : return &A2q2g2Q25223_eval;
case 25230 : return &A2q2g2Q25230_eval;
case 25233 : return &A2q2g2Q25233_eval;
case 25240 : return &A2q2g2Q25240_eval;
case 25258 : return &A2q2g2Q25258_eval;
case 25260 : return &A2q2g2Q25260_eval;
case 25263 : return &A2q2g2Q25263_eval;
case 25373 : return &A2q2g2Q25373_eval;
case 25378 : return &A2q2g2Q25378_eval;
case 25433 : return &A2q2g2Q25433_eval;
case 25448 : return &A2q2g2Q25448_eval;
case 25468 : return &A2q2g2Q25468_eval;
case 25478 : return &A2q2g2Q25478_eval;
case 25505 : return &A2q2g2Q25505_eval;
case 25520 : return &A2q2g2Q25520_eval;
case 25565 : return &A2q2g2Q25565_eval;
case 25583 : return &A2q2g2Q25583_eval;
case 25590 : return &A2q2g2Q25590_eval;
case 25593 : return &A2q2g2Q25593_eval;
case 25613 : return &A2q2g2Q25613_eval;
case 25628 : return &A2q2g2Q25628_eval;
case 25670 : return &A2q2g2Q25670_eval;
case 25680 : return &A2q2g2Q25680_eval;
case 25683 : return &A2q2g2Q25683_eval;
case 25688 : return &A2q2g2Q25688_eval;
case 25720 : return &A2q2g2Q25720_eval;
case 25730 : return &A2q2g2Q25730_eval;
case 25780 : return &A2q2g2Q25780_eval;
case 25798 : return &A2q2g2Q25798_eval;
case 25800 : return &A2q2g2Q25800_eval;
case 25803 : return &A2q2g2Q25803_eval;
case 25828 : return &A2q2g2Q25828_eval;
case 25838 : return &A2q2g2Q25838_eval;
case 25850 : return &A2q2g2Q25850_eval;
case 25860 : return &A2q2g2Q25860_eval;
case 25863 : return &A2q2g2Q25863_eval;
case 25868 : return &A2q2g2Q25868_eval;
case 25985 : return &A2q2g2Q25985_eval;
case 25990 : return &A2q2g2Q25990_eval;
case 26075 : return &A2q2g2Q26075_eval;
case 26095 : return &A2q2g2Q26095_eval;
case 26110 : return &A2q2g2Q26110_eval;
case 26125 : return &A2q2g2Q26125_eval;
case 26165 : return &A2q2g2Q26165_eval;
case 26170 : return &A2q2g2Q26170_eval;
case 26273 : return &A2q2g2Q26273_eval;
case 26278 : return &A2q2g2Q26278_eval;
case 26285 : return &A2q2g2Q26285_eval;
case 26303 : return &A2q2g2Q26303_eval;
case 26310 : return &A2q2g2Q26310_eval;
case 26313 : return &A2q2g2Q26313_eval;
case 26320 : return &A2q2g2Q26320_eval;
case 26338 : return &A2q2g2Q26338_eval;
case 26340 : return &A2q2g2Q26340_eval;
case 26343 : return &A2q2g2Q26343_eval;
case 26633 : return &A2q2g2Q26633_eval;
case 26638 : return &A2q2g2Q26638_eval;
case 26723 : return &A2q2g2Q26723_eval;
case 26743 : return &A2q2g2Q26743_eval;
case 26758 : return &A2q2g2Q26758_eval;
case 26773 : return &A2q2g2Q26773_eval;
case 26795 : return &A2q2g2Q26795_eval;
case 26815 : return &A2q2g2Q26815_eval;
case 26825 : return &A2q2g2Q26825_eval;
case 26843 : return &A2q2g2Q26843_eval;
case 26850 : return &A2q2g2Q26850_eval;
case 26853 : return &A2q2g2Q26853_eval;
case 26903 : return &A2q2g2Q26903_eval;
case 26923 : return &A2q2g2Q26923_eval;
case 26965 : return &A2q2g2Q26965_eval;
case 26970 : return &A2q2g2Q26970_eval;
case 26973 : return &A2q2g2Q26973_eval;
case 26983 : return &A2q2g2Q26983_eval;
case 27010 : return &A2q2g2Q27010_eval;
case 27025 : return &A2q2g2Q27025_eval;
case 27040 : return &A2q2g2Q27040_eval;
case 27058 : return &A2q2g2Q27058_eval;
case 27060 : return &A2q2g2Q27060_eval;
case 27063 : return &A2q2g2Q27063_eval;
case 27118 : return &A2q2g2Q27118_eval;
case 27133 : return &A2q2g2Q27133_eval;
case 27145 : return &A2q2g2Q27145_eval;
case 27150 : return &A2q2g2Q27150_eval;
case 27153 : return &A2q2g2Q27153_eval;
case 27163 : return &A2q2g2Q27163_eval;
case 27533 : return &A2q2g2Q27533_eval;
case 27538 : return &A2q2g2Q27538_eval;
case 27593 : return &A2q2g2Q27593_eval;
case 27608 : return &A2q2g2Q27608_eval;
case 27628 : return &A2q2g2Q27628_eval;
case 27638 : return &A2q2g2Q27638_eval;
case 27713 : return &A2q2g2Q27713_eval;
case 27718 : return &A2q2g2Q27718_eval;
case 27803 : return &A2q2g2Q27803_eval;
case 27823 : return &A2q2g2Q27823_eval;
case 27838 : return &A2q2g2Q27838_eval;
case 27853 : return &A2q2g2Q27853_eval;
case 28133 : return &A2q2g2Q28133_eval;
case 28148 : return &A2q2g2Q28148_eval;
case 28163 : return &A2q2g2Q28163_eval;
case 28183 : return &A2q2g2Q28183_eval;
case 28268 : return &A2q2g2Q28268_eval;
case 28273 : return &A2q2g2Q28273_eval;
case 28348 : return &A2q2g2Q28348_eval;
case 28358 : return &A2q2g2Q28358_eval;
case 28378 : return &A2q2g2Q28378_eval;
case 28393 : return &A2q2g2Q28393_eval;
case 28448 : return &A2q2g2Q28448_eval;
case 28453 : return &A2q2g2Q28453_eval;
case 28565 : return &A2q2g2Q28565_eval;
case 28580 : return &A2q2g2Q28580_eval;
case 28595 : return &A2q2g2Q28595_eval;
case 28615 : return &A2q2g2Q28615_eval;
case 28700 : return &A2q2g2Q28700_eval;
case 28705 : return &A2q2g2Q28705_eval;
case 28745 : return &A2q2g2Q28745_eval;
case 28760 : return &A2q2g2Q28760_eval;
case 28805 : return &A2q2g2Q28805_eval;
case 28823 : return &A2q2g2Q28823_eval;
case 28830 : return &A2q2g2Q28830_eval;
case 28833 : return &A2q2g2Q28833_eval;
case 28853 : return &A2q2g2Q28853_eval;
case 28868 : return &A2q2g2Q28868_eval;
case 28910 : return &A2q2g2Q28910_eval;
case 28920 : return &A2q2g2Q28920_eval;
case 28923 : return &A2q2g2Q28923_eval;
case 28928 : return &A2q2g2Q28928_eval;
case 28955 : return &A2q2g2Q28955_eval;
case 28975 : return &A2q2g2Q28975_eval;
case 28985 : return &A2q2g2Q28985_eval;
case 29003 : return &A2q2g2Q29003_eval;
case 29010 : return &A2q2g2Q29010_eval;
case 29013 : return &A2q2g2Q29013_eval;
case 29063 : return &A2q2g2Q29063_eval;
case 29083 : return &A2q2g2Q29083_eval;
case 29125 : return &A2q2g2Q29125_eval;
case 29130 : return &A2q2g2Q29130_eval;
case 29133 : return &A2q2g2Q29133_eval;
case 29143 : return &A2q2g2Q29143_eval;
case 29213 : return &A2q2g2Q29213_eval;
case 29228 : return &A2q2g2Q29228_eval;
case 29243 : return &A2q2g2Q29243_eval;
case 29263 : return &A2q2g2Q29263_eval;
case 29348 : return &A2q2g2Q29348_eval;
case 29353 : return &A2q2g2Q29353_eval;
case 29600 : return &A2q2g2Q29600_eval;
case 29605 : return &A2q2g2Q29605_eval;
case 29630 : return &A2q2g2Q29630_eval;
case 29640 : return &A2q2g2Q29640_eval;
case 29643 : return &A2q2g2Q29643_eval;
case 29648 : return &A2q2g2Q29648_eval;
case 29665 : return &A2q2g2Q29665_eval;
case 29670 : return &A2q2g2Q29670_eval;
case 29673 : return &A2q2g2Q29673_eval;
case 29683 : return &A2q2g2Q29683_eval;
case 29708 : return &A2q2g2Q29708_eval;
case 29713 : return &A2q2g2Q29713_eval;
case 29860 : return &A2q2g2Q29860_eval;
case 29870 : return &A2q2g2Q29870_eval;
case 29890 : return &A2q2g2Q29890_eval;
case 29905 : return &A2q2g2Q29905_eval;
case 29960 : return &A2q2g2Q29960_eval;
case 29965 : return &A2q2g2Q29965_eval;
case 30040 : return &A2q2g2Q30040_eval;
case 30050 : return &A2q2g2Q30050_eval;
case 30100 : return &A2q2g2Q30100_eval;
case 30118 : return &A2q2g2Q30118_eval;
case 30120 : return &A2q2g2Q30120_eval;
case 30123 : return &A2q2g2Q30123_eval;
case 30148 : return &A2q2g2Q30148_eval;
case 30158 : return &A2q2g2Q30158_eval;
case 30170 : return &A2q2g2Q30170_eval;
case 30180 : return &A2q2g2Q30180_eval;
case 30183 : return &A2q2g2Q30183_eval;
case 30188 : return &A2q2g2Q30188_eval;
case 30250 : return &A2q2g2Q30250_eval;
case 30265 : return &A2q2g2Q30265_eval;
case 30280 : return &A2q2g2Q30280_eval;
case 30298 : return &A2q2g2Q30298_eval;
case 30300 : return &A2q2g2Q30300_eval;
case 30303 : return &A2q2g2Q30303_eval;
case 30358 : return &A2q2g2Q30358_eval;
case 30373 : return &A2q2g2Q30373_eval;
case 30385 : return &A2q2g2Q30385_eval;
case 30390 : return &A2q2g2Q30390_eval;
case 30393 : return &A2q2g2Q30393_eval;
case 30403 : return &A2q2g2Q30403_eval;
case 30508 : return &A2q2g2Q30508_eval;
case 30518 : return &A2q2g2Q30518_eval;
case 30538 : return &A2q2g2Q30538_eval;
case 30553 : return &A2q2g2Q30553_eval;
case 30608 : return &A2q2g2Q30608_eval;
case 30613 : return &A2q2g2Q30613_eval;
case 30680 : return &A2q2g2Q30680_eval;
case 30685 : return &A2q2g2Q30685_eval;
case 30710 : return &A2q2g2Q30710_eval;
case 30720 : return &A2q2g2Q30720_eval;
case 30723 : return &A2q2g2Q30723_eval;
case 30728 : return &A2q2g2Q30728_eval;
case 30745 : return &A2q2g2Q30745_eval;
case 30750 : return &A2q2g2Q30750_eval;
case 30753 : return &A2q2g2Q30753_eval;
case 30763 : return &A2q2g2Q30763_eval;
case 30788 : return &A2q2g2Q30788_eval;
case 30793 : return &A2q2g2Q30793_eval;
case 31157 : return &A2q2g2Q31157_eval;
case 31172 : return &A2q2g2Q31172_eval;
case 31187 : return &A2q2g2Q31187_eval;
case 31207 : return &A2q2g2Q31207_eval;
case 31292 : return &A2q2g2Q31292_eval;
case 31297 : return &A2q2g2Q31297_eval;
case 31337 : return &A2q2g2Q31337_eval;
case 31352 : return &A2q2g2Q31352_eval;
case 31397 : return &A2q2g2Q31397_eval;
case 31415 : return &A2q2g2Q31415_eval;
case 31422 : return &A2q2g2Q31422_eval;
case 31425 : return &A2q2g2Q31425_eval;
case 31445 : return &A2q2g2Q31445_eval;
case 31460 : return &A2q2g2Q31460_eval;
case 31502 : return &A2q2g2Q31502_eval;
case 31512 : return &A2q2g2Q31512_eval;
case 31515 : return &A2q2g2Q31515_eval;
case 31520 : return &A2q2g2Q31520_eval;
case 31547 : return &A2q2g2Q31547_eval;
case 31567 : return &A2q2g2Q31567_eval;
case 31577 : return &A2q2g2Q31577_eval;
case 31595 : return &A2q2g2Q31595_eval;
case 31602 : return &A2q2g2Q31602_eval;
case 31605 : return &A2q2g2Q31605_eval;
case 31655 : return &A2q2g2Q31655_eval;
case 31675 : return &A2q2g2Q31675_eval;
case 31717 : return &A2q2g2Q31717_eval;
case 31722 : return &A2q2g2Q31722_eval;
case 31725 : return &A2q2g2Q31725_eval;
case 31735 : return &A2q2g2Q31735_eval;
case 31805 : return &A2q2g2Q31805_eval;
case 31820 : return &A2q2g2Q31820_eval;
case 31835 : return &A2q2g2Q31835_eval;
case 31855 : return &A2q2g2Q31855_eval;
case 31940 : return &A2q2g2Q31940_eval;
case 31945 : return &A2q2g2Q31945_eval;
case 32192 : return &A2q2g2Q32192_eval;
case 32197 : return &A2q2g2Q32197_eval;
case 32222 : return &A2q2g2Q32222_eval;
case 32232 : return &A2q2g2Q32232_eval;
case 32235 : return &A2q2g2Q32235_eval;
case 32240 : return &A2q2g2Q32240_eval;
case 32257 : return &A2q2g2Q32257_eval;
case 32262 : return &A2q2g2Q32262_eval;
case 32265 : return &A2q2g2Q32265_eval;
case 32275 : return &A2q2g2Q32275_eval;
case 32300 : return &A2q2g2Q32300_eval;
case 32305 : return &A2q2g2Q32305_eval;
case 32417 : return &A2q2g2Q32417_eval;
case 32432 : return &A2q2g2Q32432_eval;
case 32477 : return &A2q2g2Q32477_eval;
case 32495 : return &A2q2g2Q32495_eval;
case 32502 : return &A2q2g2Q32502_eval;
case 32505 : return &A2q2g2Q32505_eval;
case 32525 : return &A2q2g2Q32525_eval;
case 32540 : return &A2q2g2Q32540_eval;
case 32582 : return &A2q2g2Q32582_eval;
case 32592 : return &A2q2g2Q32592_eval;
case 32595 : return &A2q2g2Q32595_eval;
case 32600 : return &A2q2g2Q32600_eval;
case 32837 : return &A2q2g2Q32837_eval;
case 32855 : return &A2q2g2Q32855_eval;
case 32862 : return &A2q2g2Q32862_eval;
case 32865 : return &A2q2g2Q32865_eval;
case 32945 : return &A2q2g2Q32945_eval;
case 32963 : return &A2q2g2Q32963_eval;
case 32970 : return &A2q2g2Q32970_eval;
case 32973 : return &A2q2g2Q32973_eval;
case 33012 : return &A2q2g2Q33012_eval;
case 33015 : return &A2q2g2Q33015_eval;
case 33030 : return &A2q2g2Q33030_eval;
case 33033 : return &A2q2g2Q33033_eval;
case 33065 : return &A2q2g2Q33065_eval;
case 33080 : return &A2q2g2Q33080_eval;
case 33125 : return &A2q2g2Q33125_eval;
case 33143 : return &A2q2g2Q33143_eval;
case 33150 : return &A2q2g2Q33150_eval;
case 33153 : return &A2q2g2Q33153_eval;
case 33173 : return &A2q2g2Q33173_eval;
case 33188 : return &A2q2g2Q33188_eval;
case 33230 : return &A2q2g2Q33230_eval;
case 33240 : return &A2q2g2Q33240_eval;
case 33243 : return &A2q2g2Q33243_eval;
case 33248 : return &A2q2g2Q33248_eval;
case 33482 : return &A2q2g2Q33482_eval;
case 33492 : return &A2q2g2Q33492_eval;
case 33495 : return &A2q2g2Q33495_eval;
case 33500 : return &A2q2g2Q33500_eval;
case 33552 : return &A2q2g2Q33552_eval;
case 33555 : return &A2q2g2Q33555_eval;
case 33570 : return &A2q2g2Q33570_eval;
case 33573 : return &A2q2g2Q33573_eval;
case 33590 : return &A2q2g2Q33590_eval;
case 33600 : return &A2q2g2Q33600_eval;
case 33603 : return &A2q2g2Q33603_eval;
case 33608 : return &A2q2g2Q33608_eval;
case 33707 : return &A2q2g2Q33707_eval;
case 33727 : return &A2q2g2Q33727_eval;
case 33737 : return &A2q2g2Q33737_eval;
case 33755 : return &A2q2g2Q33755_eval;
case 33762 : return &A2q2g2Q33762_eval;
case 33765 : return &A2q2g2Q33765_eval;
case 33815 : return &A2q2g2Q33815_eval;
case 33835 : return &A2q2g2Q33835_eval;
case 33877 : return &A2q2g2Q33877_eval;
case 33882 : return &A2q2g2Q33882_eval;
case 33885 : return &A2q2g2Q33885_eval;
case 33895 : return &A2q2g2Q33895_eval;
case 33917 : return &A2q2g2Q33917_eval;
case 33935 : return &A2q2g2Q33935_eval;
case 33942 : return &A2q2g2Q33942_eval;
case 33945 : return &A2q2g2Q33945_eval;
case 34025 : return &A2q2g2Q34025_eval;
case 34043 : return &A2q2g2Q34043_eval;
case 34050 : return &A2q2g2Q34050_eval;
case 34053 : return &A2q2g2Q34053_eval;
case 34092 : return &A2q2g2Q34092_eval;
case 34095 : return &A2q2g2Q34095_eval;
case 34110 : return &A2q2g2Q34110_eval;
case 34113 : return &A2q2g2Q34113_eval;
case 34355 : return &A2q2g2Q34355_eval;
case 34375 : return &A2q2g2Q34375_eval;
case 34385 : return &A2q2g2Q34385_eval;
case 34403 : return &A2q2g2Q34403_eval;
case 34410 : return &A2q2g2Q34410_eval;
case 34413 : return &A2q2g2Q34413_eval;
case 34463 : return &A2q2g2Q34463_eval;
case 34483 : return &A2q2g2Q34483_eval;
case 34525 : return &A2q2g2Q34525_eval;
case 34530 : return &A2q2g2Q34530_eval;
case 34533 : return &A2q2g2Q34533_eval;
case 34543 : return &A2q2g2Q34543_eval;
case 34777 : return &A2q2g2Q34777_eval;
case 34782 : return &A2q2g2Q34782_eval;
case 34785 : return &A2q2g2Q34785_eval;
case 34795 : return &A2q2g2Q34795_eval;
case 34812 : return &A2q2g2Q34812_eval;
case 34815 : return &A2q2g2Q34815_eval;
case 34830 : return &A2q2g2Q34830_eval;
case 34833 : return &A2q2g2Q34833_eval;
case 34885 : return &A2q2g2Q34885_eval;
case 34890 : return &A2q2g2Q34890_eval;
case 34893 : return &A2q2g2Q34893_eval;
case 34903 : return &A2q2g2Q34903_eval;
case 35045 : return &A2q2g2Q35045_eval;
case 35060 : return &A2q2g2Q35060_eval;
case 35075 : return &A2q2g2Q35075_eval;
case 35095 : return &A2q2g2Q35095_eval;
case 35180 : return &A2q2g2Q35180_eval;
case 35185 : return &A2q2g2Q35185_eval;
case 35225 : return &A2q2g2Q35225_eval;
case 35240 : return &A2q2g2Q35240_eval;
case 35285 : return &A2q2g2Q35285_eval;
case 35303 : return &A2q2g2Q35303_eval;
case 35310 : return &A2q2g2Q35310_eval;
case 35313 : return &A2q2g2Q35313_eval;
case 35333 : return &A2q2g2Q35333_eval;
case 35348 : return &A2q2g2Q35348_eval;
case 35390 : return &A2q2g2Q35390_eval;
case 35400 : return &A2q2g2Q35400_eval;
case 35403 : return &A2q2g2Q35403_eval;
case 35408 : return &A2q2g2Q35408_eval;
case 35435 : return &A2q2g2Q35435_eval;
case 35455 : return &A2q2g2Q35455_eval;
case 35465 : return &A2q2g2Q35465_eval;
case 35483 : return &A2q2g2Q35483_eval;
case 35490 : return &A2q2g2Q35490_eval;
case 35493 : return &A2q2g2Q35493_eval;
case 35543 : return &A2q2g2Q35543_eval;
case 35563 : return &A2q2g2Q35563_eval;
case 35605 : return &A2q2g2Q35605_eval;
case 35610 : return &A2q2g2Q35610_eval;
case 35613 : return &A2q2g2Q35613_eval;
case 35623 : return &A2q2g2Q35623_eval;
case 35693 : return &A2q2g2Q35693_eval;
case 35708 : return &A2q2g2Q35708_eval;
case 35723 : return &A2q2g2Q35723_eval;
case 35743 : return &A2q2g2Q35743_eval;
case 35828 : return &A2q2g2Q35828_eval;
case 35833 : return &A2q2g2Q35833_eval;
case 36080 : return &A2q2g2Q36080_eval;
case 36085 : return &A2q2g2Q36085_eval;
case 36110 : return &A2q2g2Q36110_eval;
case 36120 : return &A2q2g2Q36120_eval;
case 36123 : return &A2q2g2Q36123_eval;
case 36128 : return &A2q2g2Q36128_eval;
case 36145 : return &A2q2g2Q36145_eval;
case 36150 : return &A2q2g2Q36150_eval;
case 36153 : return &A2q2g2Q36153_eval;
case 36163 : return &A2q2g2Q36163_eval;
case 36188 : return &A2q2g2Q36188_eval;
case 36193 : return &A2q2g2Q36193_eval;
case 37592 : return &A2q2g2Q37592_eval;
case 37597 : return &A2q2g2Q37597_eval;
case 37622 : return &A2q2g2Q37622_eval;
case 37632 : return &A2q2g2Q37632_eval;
case 37635 : return &A2q2g2Q37635_eval;
case 37640 : return &A2q2g2Q37640_eval;
case 37657 : return &A2q2g2Q37657_eval;
case 37662 : return &A2q2g2Q37662_eval;
case 37665 : return &A2q2g2Q37665_eval;
case 37675 : return &A2q2g2Q37675_eval;
case 37700 : return &A2q2g2Q37700_eval;
case 37705 : return &A2q2g2Q37705_eval;
case 37802 : return &A2q2g2Q37802_eval;
case 37812 : return &A2q2g2Q37812_eval;
case 37815 : return &A2q2g2Q37815_eval;
case 37820 : return &A2q2g2Q37820_eval;
case 37872 : return &A2q2g2Q37872_eval;
case 37875 : return &A2q2g2Q37875_eval;
case 37890 : return &A2q2g2Q37890_eval;
case 37893 : return &A2q2g2Q37893_eval;
case 37910 : return &A2q2g2Q37910_eval;
case 37920 : return &A2q2g2Q37920_eval;
case 37923 : return &A2q2g2Q37923_eval;
case 37928 : return &A2q2g2Q37928_eval;
case 38017 : return &A2q2g2Q38017_eval;
case 38022 : return &A2q2g2Q38022_eval;
case 38025 : return &A2q2g2Q38025_eval;
case 38035 : return &A2q2g2Q38035_eval;
case 38052 : return &A2q2g2Q38052_eval;
case 38055 : return &A2q2g2Q38055_eval;
case 38070 : return &A2q2g2Q38070_eval;
case 38073 : return &A2q2g2Q38073_eval;
case 38125 : return &A2q2g2Q38125_eval;
case 38130 : return &A2q2g2Q38130_eval;
case 38133 : return &A2q2g2Q38133_eval;
case 38143 : return &A2q2g2Q38143_eval;
case 38240 : return &A2q2g2Q38240_eval;
case 38245 : return &A2q2g2Q38245_eval;
case 38270 : return &A2q2g2Q38270_eval;
case 38280 : return &A2q2g2Q38280_eval;
case 38283 : return &A2q2g2Q38283_eval;
case 38288 : return &A2q2g2Q38288_eval;
case 38305 : return &A2q2g2Q38305_eval;
case 38310 : return &A2q2g2Q38310_eval;
case 38313 : return &A2q2g2Q38313_eval;
case 38323 : return &A2q2g2Q38323_eval;
case 38348 : return &A2q2g2Q38348_eval;
case 38353 : return &A2q2g2Q38353_eval;
case 38932 : return &A2q2g2Q38932_eval;
case 38942 : return &A2q2g2Q38942_eval;
case 38962 : return &A2q2g2Q38962_eval;
case 38977 : return &A2q2g2Q38977_eval;
case 39032 : return &A2q2g2Q39032_eval;
case 39037 : return &A2q2g2Q39037_eval;
case 39112 : return &A2q2g2Q39112_eval;
case 39122 : return &A2q2g2Q39122_eval;
case 39172 : return &A2q2g2Q39172_eval;
case 39190 : return &A2q2g2Q39190_eval;
case 39192 : return &A2q2g2Q39192_eval;
case 39195 : return &A2q2g2Q39195_eval;
case 39220 : return &A2q2g2Q39220_eval;
case 39230 : return &A2q2g2Q39230_eval;
case 39242 : return &A2q2g2Q39242_eval;
case 39252 : return &A2q2g2Q39252_eval;
case 39255 : return &A2q2g2Q39255_eval;
case 39260 : return &A2q2g2Q39260_eval;
case 39322 : return &A2q2g2Q39322_eval;
case 39337 : return &A2q2g2Q39337_eval;
case 39352 : return &A2q2g2Q39352_eval;
case 39370 : return &A2q2g2Q39370_eval;
case 39372 : return &A2q2g2Q39372_eval;
case 39375 : return &A2q2g2Q39375_eval;
case 39430 : return &A2q2g2Q39430_eval;
case 39445 : return &A2q2g2Q39445_eval;
case 39457 : return &A2q2g2Q39457_eval;
case 39462 : return &A2q2g2Q39462_eval;
case 39465 : return &A2q2g2Q39465_eval;
case 39475 : return &A2q2g2Q39475_eval;
case 39580 : return &A2q2g2Q39580_eval;
case 39590 : return &A2q2g2Q39590_eval;
case 39610 : return &A2q2g2Q39610_eval;
case 39625 : return &A2q2g2Q39625_eval;
case 39680 : return &A2q2g2Q39680_eval;
case 39685 : return &A2q2g2Q39685_eval;
case 39752 : return &A2q2g2Q39752_eval;
case 39757 : return &A2q2g2Q39757_eval;
case 39782 : return &A2q2g2Q39782_eval;
case 39792 : return &A2q2g2Q39792_eval;
case 39795 : return &A2q2g2Q39795_eval;
case 39800 : return &A2q2g2Q39800_eval;
case 39817 : return &A2q2g2Q39817_eval;
case 39822 : return &A2q2g2Q39822_eval;
case 39825 : return &A2q2g2Q39825_eval;
case 39835 : return &A2q2g2Q39835_eval;
case 39860 : return &A2q2g2Q39860_eval;
case 39865 : return &A2q2g2Q39865_eval;
case 40192 : return &A2q2g2Q40192_eval;
case 40202 : return &A2q2g2Q40202_eval;
case 40252 : return &A2q2g2Q40252_eval;
case 40270 : return &A2q2g2Q40270_eval;
case 40272 : return &A2q2g2Q40272_eval;
case 40275 : return &A2q2g2Q40275_eval;
case 40300 : return &A2q2g2Q40300_eval;
case 40310 : return &A2q2g2Q40310_eval;
case 40322 : return &A2q2g2Q40322_eval;
case 40332 : return &A2q2g2Q40332_eval;
case 40335 : return &A2q2g2Q40335_eval;
case 40340 : return &A2q2g2Q40340_eval;
case 40612 : return &A2q2g2Q40612_eval;
case 40630 : return &A2q2g2Q40630_eval;
case 40632 : return &A2q2g2Q40632_eval;
case 40635 : return &A2q2g2Q40635_eval;
case 40720 : return &A2q2g2Q40720_eval;
case 40738 : return &A2q2g2Q40738_eval;
case 40740 : return &A2q2g2Q40740_eval;
case 40743 : return &A2q2g2Q40743_eval;
case 40752 : return &A2q2g2Q40752_eval;
case 40755 : return &A2q2g2Q40755_eval;
case 40770 : return &A2q2g2Q40770_eval;
case 40773 : return &A2q2g2Q40773_eval;
case 40840 : return &A2q2g2Q40840_eval;
case 40850 : return &A2q2g2Q40850_eval;
case 40900 : return &A2q2g2Q40900_eval;
case 40918 : return &A2q2g2Q40918_eval;
case 40920 : return &A2q2g2Q40920_eval;
case 40923 : return &A2q2g2Q40923_eval;
case 40948 : return &A2q2g2Q40948_eval;
case 40958 : return &A2q2g2Q40958_eval;
case 40970 : return &A2q2g2Q40970_eval;
case 40980 : return &A2q2g2Q40980_eval;
case 40983 : return &A2q2g2Q40983_eval;
case 40988 : return &A2q2g2Q40988_eval;
case 41042 : return &A2q2g2Q41042_eval;
case 41052 : return &A2q2g2Q41052_eval;
case 41055 : return &A2q2g2Q41055_eval;
case 41060 : return &A2q2g2Q41060_eval;
case 41112 : return &A2q2g2Q41112_eval;
case 41115 : return &A2q2g2Q41115_eval;
case 41130 : return &A2q2g2Q41130_eval;
case 41133 : return &A2q2g2Q41133_eval;
case 41150 : return &A2q2g2Q41150_eval;
case 41160 : return &A2q2g2Q41160_eval;
case 41163 : return &A2q2g2Q41163_eval;
case 41168 : return &A2q2g2Q41168_eval;
case 41482 : return &A2q2g2Q41482_eval;
case 41497 : return &A2q2g2Q41497_eval;
case 41512 : return &A2q2g2Q41512_eval;
case 41530 : return &A2q2g2Q41530_eval;
case 41532 : return &A2q2g2Q41532_eval;
case 41535 : return &A2q2g2Q41535_eval;
case 41590 : return &A2q2g2Q41590_eval;
case 41605 : return &A2q2g2Q41605_eval;
case 41617 : return &A2q2g2Q41617_eval;
case 41622 : return &A2q2g2Q41622_eval;
case 41625 : return &A2q2g2Q41625_eval;
case 41635 : return &A2q2g2Q41635_eval;
case 41692 : return &A2q2g2Q41692_eval;
case 41710 : return &A2q2g2Q41710_eval;
case 41712 : return &A2q2g2Q41712_eval;
case 41715 : return &A2q2g2Q41715_eval;
case 41800 : return &A2q2g2Q41800_eval;
case 41818 : return &A2q2g2Q41818_eval;
case 41820 : return &A2q2g2Q41820_eval;
case 41823 : return &A2q2g2Q41823_eval;
case 41832 : return &A2q2g2Q41832_eval;
case 41835 : return &A2q2g2Q41835_eval;
case 41850 : return &A2q2g2Q41850_eval;
case 41853 : return &A2q2g2Q41853_eval;
case 42130 : return &A2q2g2Q42130_eval;
case 42145 : return &A2q2g2Q42145_eval;
case 42160 : return &A2q2g2Q42160_eval;
case 42178 : return &A2q2g2Q42178_eval;
case 42180 : return &A2q2g2Q42180_eval;
case 42183 : return &A2q2g2Q42183_eval;
case 42238 : return &A2q2g2Q42238_eval;
case 42253 : return &A2q2g2Q42253_eval;
case 42265 : return &A2q2g2Q42265_eval;
case 42270 : return &A2q2g2Q42270_eval;
case 42273 : return &A2q2g2Q42273_eval;
case 42283 : return &A2q2g2Q42283_eval;
case 42337 : return &A2q2g2Q42337_eval;
case 42342 : return &A2q2g2Q42342_eval;
case 42345 : return &A2q2g2Q42345_eval;
case 42355 : return &A2q2g2Q42355_eval;
case 42372 : return &A2q2g2Q42372_eval;
case 42375 : return &A2q2g2Q42375_eval;
case 42390 : return &A2q2g2Q42390_eval;
case 42393 : return &A2q2g2Q42393_eval;
case 42445 : return &A2q2g2Q42445_eval;
case 42450 : return &A2q2g2Q42450_eval;
case 42453 : return &A2q2g2Q42453_eval;
case 42463 : return &A2q2g2Q42463_eval;
case 42820 : return &A2q2g2Q42820_eval;
case 42830 : return &A2q2g2Q42830_eval;
case 42850 : return &A2q2g2Q42850_eval;
case 42865 : return &A2q2g2Q42865_eval;
case 42920 : return &A2q2g2Q42920_eval;
case 42925 : return &A2q2g2Q42925_eval;
case 43000 : return &A2q2g2Q43000_eval;
case 43010 : return &A2q2g2Q43010_eval;
case 43060 : return &A2q2g2Q43060_eval;
case 43078 : return &A2q2g2Q43078_eval;
case 43080 : return &A2q2g2Q43080_eval;
case 43083 : return &A2q2g2Q43083_eval;
case 43108 : return &A2q2g2Q43108_eval;
case 43118 : return &A2q2g2Q43118_eval;
case 43130 : return &A2q2g2Q43130_eval;
case 43140 : return &A2q2g2Q43140_eval;
case 43143 : return &A2q2g2Q43143_eval;
case 43148 : return &A2q2g2Q43148_eval;
case 43210 : return &A2q2g2Q43210_eval;
case 43225 : return &A2q2g2Q43225_eval;
case 43240 : return &A2q2g2Q43240_eval;
case 43258 : return &A2q2g2Q43258_eval;
case 43260 : return &A2q2g2Q43260_eval;
case 43263 : return &A2q2g2Q43263_eval;
case 43318 : return &A2q2g2Q43318_eval;
case 43333 : return &A2q2g2Q43333_eval;
case 43345 : return &A2q2g2Q43345_eval;
case 43350 : return &A2q2g2Q43350_eval;
case 43353 : return &A2q2g2Q43353_eval;
case 43363 : return &A2q2g2Q43363_eval;
case 43468 : return &A2q2g2Q43468_eval;
case 43478 : return &A2q2g2Q43478_eval;
case 43498 : return &A2q2g2Q43498_eval;
case 43513 : return &A2q2g2Q43513_eval;
case 43568 : return &A2q2g2Q43568_eval;
case 43573 : return &A2q2g2Q43573_eval;
case 43640 : return &A2q2g2Q43640_eval;
case 43645 : return &A2q2g2Q43645_eval;
case 43670 : return &A2q2g2Q43670_eval;
case 43680 : return &A2q2g2Q43680_eval;
case 43683 : return &A2q2g2Q43683_eval;
case 43688 : return &A2q2g2Q43688_eval;
case 43705 : return &A2q2g2Q43705_eval;
case 43710 : return &A2q2g2Q43710_eval;
case 43713 : return &A2q2g2Q43713_eval;
case 43723 : return &A2q2g2Q43723_eval;
case 43748 : return &A2q2g2Q43748_eval;
case 43753 : return &A2q2g2Q43753_eval;
case 44072 : return &A2q2g2Q44072_eval;
case 44077 : return &A2q2g2Q44077_eval;
case 44102 : return &A2q2g2Q44102_eval;
case 44112 : return &A2q2g2Q44112_eval;
case 44115 : return &A2q2g2Q44115_eval;
case 44120 : return &A2q2g2Q44120_eval;
case 44137 : return &A2q2g2Q44137_eval;
case 44142 : return &A2q2g2Q44142_eval;
case 44145 : return &A2q2g2Q44145_eval;
case 44155 : return &A2q2g2Q44155_eval;
case 44180 : return &A2q2g2Q44180_eval;
case 44185 : return &A2q2g2Q44185_eval;
case 44282 : return &A2q2g2Q44282_eval;
case 44292 : return &A2q2g2Q44292_eval;
case 44295 : return &A2q2g2Q44295_eval;
case 44300 : return &A2q2g2Q44300_eval;
case 44352 : return &A2q2g2Q44352_eval;
case 44355 : return &A2q2g2Q44355_eval;
case 44370 : return &A2q2g2Q44370_eval;
case 44373 : return &A2q2g2Q44373_eval;
case 44390 : return &A2q2g2Q44390_eval;
case 44400 : return &A2q2g2Q44400_eval;
case 44403 : return &A2q2g2Q44403_eval;
case 44408 : return &A2q2g2Q44408_eval;
case 44497 : return &A2q2g2Q44497_eval;
case 44502 : return &A2q2g2Q44502_eval;
case 44505 : return &A2q2g2Q44505_eval;
case 44515 : return &A2q2g2Q44515_eval;
case 44532 : return &A2q2g2Q44532_eval;
case 44535 : return &A2q2g2Q44535_eval;
case 44550 : return &A2q2g2Q44550_eval;
case 44553 : return &A2q2g2Q44553_eval;
case 44605 : return &A2q2g2Q44605_eval;
case 44610 : return &A2q2g2Q44610_eval;
case 44613 : return &A2q2g2Q44613_eval;
case 44623 : return &A2q2g2Q44623_eval;
case 44720 : return &A2q2g2Q44720_eval;
case 44725 : return &A2q2g2Q44725_eval;
case 44750 : return &A2q2g2Q44750_eval;
case 44760 : return &A2q2g2Q44760_eval;
case 44763 : return &A2q2g2Q44763_eval;
case 44768 : return &A2q2g2Q44768_eval;
case 44785 : return &A2q2g2Q44785_eval;
case 44790 : return &A2q2g2Q44790_eval;
case 44793 : return &A2q2g2Q44793_eval;
case 44803 : return &A2q2g2Q44803_eval;
case 44828 : return &A2q2g2Q44828_eval;
case 44833 : return &A2q2g2Q44833_eval;
#endif
default: return 0;

}
}



template complex<R>  (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2g2Q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
44532: m m qm qp Qm Qp
38052: m m qm qp Qp Qm
42372: m m qm Qm qp Qp
22932: m m qm Qm Qp qp
34812: m m qm Qp qp Qm
21852: m m qm Qp Qm qp
44352: m m qp qm Qm Qp
37872: m m qp qm Qp Qm
41112: m m qp Qm qm Qp
15192: m m qp Qm Qp qm
33552: m m qp Qp qm Qm
14112: m m qp Qp Qm qm
41832: m m Qm qm qp Qp
22392: m m Qm qm Qp qp
40752: m m Qm qp qm Qp
14832: m m Qm qp Qp qm
18072: m m Qm Qp qm qp
11592: m m Qm Qp qp qm
34092: m m Qp qm qp Qm
21132: m m Qp qm Qm qp
33012: m m Qp qp qm Qm
13572: m m Qp qp Qm qm
17892: m m Qp Qm qm qp
11412: m m Qp Qm qp qm
44502: m qm m qp Qm Qp
38022: m qm m qp Qp Qm
42342: m qm m Qm qp Qp
22902: m qm m Qm Qp qp
34782: m qm m Qp qp Qm
21822: m qm m Qp Qm qp
44142: m qm qp m Qm Qp
37662: m qm qp m Qp Qm
44790: m qm qp p Qm Qp
38310: m qm qp p Qp Qm
39822: m qm qp Qm m Qp
43710: m qm qp Qm p Qp
7422: m qm qp Qm Qp m
30750: m qm qp Qm Qp p
32262: m qm qp Qp m Qm
36150: m qm qp Qp p Qm
6342: m qm qp Qp Qm m
29670: m qm qp Qp Qm p
44610: m qm p qp Qm Qp
38130: m qm p qp Qp Qm
42450: m qm p Qm qp Qp
23010: m qm p Qm Qp qp
34890: m qm p Qp qp Qm
21930: m qm p Qp Qm qp
41622: m qm Qm m qp Qp
22182: m qm Qm m Qp qp
39462: m qm Qm qp m Qp
43350: m qm Qm qp p Qp
7062: m qm Qm qp Qp m
30390: m qm Qm qp Qp p
42270: m qm Qm p qp Qp
22830: m qm Qm p Qp qp
16782: m qm Qm Qp m qp
3822: m qm Qm Qp qp m
27150: m qm Qm Qp qp p
20670: m qm Qm Qp p qp
33882: m qm Qp m qp Qm
20922: m qm Qp m Qm qp
31722: m qm Qp qp m Qm
35610: m qm Qp qp p Qm
5802: m qm Qp qp Qm m
29130: m qm Qp qp Qm p
34530: m qm Qp p qp Qm
21570: m qm Qp p Qm qp
16602: m qm Qp Qm m qp
3642: m qm Qp Qm qp m
26970: m qm Qp Qm qp p
20490: m qm Qp Qm p qp
44292: m qp m qm Qm Qp
37812: m qp m qm Qp Qm
41052: m qp m Qm qm Qp
15132: m qp m Qm Qp qm
33492: m qp m Qp qm Qm
14052: m qp m Qp Qm qm
44112: m qp qm m Qm Qp
37632: m qp qm m Qp Qm
44760: m qp qm p Qm Qp
38280: m qp qm p Qp Qm
39792: m qp qm Qm m Qp
43680: m qp qm Qm p Qp
7392: m qp qm Qm Qp m
30720: m qp qm Qm Qp p
32232: m qp qm Qp m Qm
36120: m qp qm Qp p Qm
6312: m qp qm Qp Qm m
29640: m qp qm Qp Qm p
44400: m qp p qm Qm Qp
37920: m qp p qm Qp Qm
41160: m qp p Qm qm Qp
15240: m qp p Qm Qp qm
33600: m qp p Qp qm Qm
14160: m qp p Qp Qm qm
40332: m qp Qm m qm Qp
14412: m qp Qm m Qp qm
39252: m qp Qm qm m Qp
43140: m qp Qm qm p Qp
6852: m qp Qm qm Qp m
30180: m qp Qm qm Qp p
40980: m qp Qm p qm Qp
15060: m qp Qm p Qp qm
9012: m qp Qm Qp m qm
2532: m qp Qm Qp qm m
25860: m qp Qm Qp qm p
12900: m qp Qm Qp p qm
32592: m qp Qp m qm Qm
13152: m qp Qp m Qm qm
31512: m qp Qp qm m Qm
35400: m qp Qp qm p Qm
5592: m qp Qp qm Qm m
28920: m qp Qp qm Qm p
33240: m qp Qp p qm Qm
13800: m qp Qp p Qm qm
8832: m qp Qp Qm m qm
2352: m qp Qp Qm qm m
25680: m qp Qp Qm qm p
12720: m qp Qp Qm p qm
44550: m p qm qp Qm Qp
38070: m p qm qp Qp Qm
42390: m p qm Qm qp Qp
22950: m p qm Qm Qp qp
34830: m p qm Qp qp Qm
21870: m p qm Qp Qm qp
44370: m p qp qm Qm Qp
37890: m p qp qm Qp Qm
41130: m p qp Qm qm Qp
15210: m p qp Qm Qp qm
33570: m p qp Qp qm Qm
14130: m p qp Qp Qm qm
41850: m p Qm qm qp Qp
22410: m p Qm qm Qp qp
40770: m p Qm qp qm Qp
14850: m p Qm qp Qp qm
18090: m p Qm Qp qm qp
11610: m p Qm Qp qp qm
34110: m p Qp qm qp Qm
21150: m p Qp qm Qm qp
33030: m p Qp qp qm Qm
13590: m p Qp qp Qm qm
17910: m p Qp Qm qm qp
11430: m p Qp Qm qp qm
41712: m Qm m qm qp Qp
22272: m Qm m qm Qp qp
40632: m Qm m qp qm Qp
14712: m Qm m qp Qp qm
17952: m Qm m Qp qm qp
11472: m Qm m Qp qp qm
41532: m Qm qm m qp Qp
22092: m Qm qm m Qp qp
39372: m Qm qm qp m Qp
43260: m Qm qm qp p Qp
6972: m Qm qm qp Qp m
30300: m Qm qm qp Qp p
42180: m Qm qm p qp Qp
22740: m Qm qm p Qp qp
16692: m Qm qm Qp m qp
3732: m Qm qm Qp qp m
27060: m Qm qm Qp qp p
20580: m Qm qm Qp p qp
40272: m Qm qp m qm Qp
14352: m Qm qp m Qp qm
39192: m Qm qp qm m Qp
43080: m Qm qp qm p Qp
6792: m Qm qp qm Qp m
30120: m Qm qp qm Qp p
40920: m Qm qp p qm Qp
15000: m Qm qp p Qp qm
8952: m Qm qp Qp m qm
2472: m Qm qp Qp qm m
25800: m Qm qp Qp qm p
12840: m Qm qp Qp p qm
41820: m Qm p qm qp Qp
22380: m Qm p qm Qp qp
40740: m Qm p qp qm Qp
14820: m Qm p qp Qp qm
18060: m Qm p Qp qm qp
11580: m Qm p Qp qp qm
17052: m Qm Qp m qm qp
10572: m Qm Qp m qp qm
15972: m Qm Qp qm m qp
3012: m Qm Qp qm qp m
26340: m Qm Qp qm qp p
19860: m Qm Qp qm p qp
8412: m Qm Qp qp m qm
1932: m Qm Qp qp qm m
25260: m Qm Qp qp qm p
12300: m Qm Qp qp p qm
17700: m Qm Qp p qm qp
11220: m Qm Qp p qp qm
33942: m Qp m qm qp Qm
20982: m Qp m qm Qm qp
32862: m Qp m qp qm Qm
13422: m Qp m qp Qm qm
17742: m Qp m Qm qm qp
11262: m Qp m Qm qp qm
33762: m Qp qm m qp Qm
20802: m Qp qm m Qm qp
31602: m Qp qm qp m Qm
35490: m Qp qm qp p Qm
5682: m Qp qm qp Qm m
29010: m Qp qm qp Qm p
34410: m Qp qm p qp Qm
21450: m Qp qm p Qm qp
16482: m Qp qm Qm m qp
3522: m Qp qm Qm qp m
26850: m Qp qm Qm qp p
20370: m Qp qm Qm p qp
32502: m Qp qp m qm Qm
13062: m Qp qp m Qm qm
31422: m Qp qp qm m Qm
35310: m Qp qp qm p Qm
5502: m Qp qp qm Qm m
28830: m Qp qp qm Qm p
33150: m Qp qp p qm Qm
13710: m Qp qp p Qm qm
8742: m Qp qp Qm m qm
2262: m Qp qp Qm qm m
25590: m Qp qp Qm qm p
12630: m Qp qp Qm p qm
34050: m Qp p qm qp Qm
21090: m Qp p qm Qm qp
32970: m Qp p qp qm Qm
13530: m Qp p qp Qm qm
17850: m Qp p Qm qm qp
11370: m Qp p Qm qp qm
17022: m Qp Qm m qm qp
10542: m Qp Qm m qp qm
15942: m Qp Qm qm m qp
2982: m Qp Qm qm qp m
26310: m Qp Qm qm qp p
19830: m Qp Qm qm p qp
8382: m Qp Qm qp m qm
1902: m Qp Qm qp qm m
25230: m Qp Qm qp qm p
12270: m Qp Qm qp p qm
17670: m Qp Qm p qm qp
11190: m Qp Qm p qp qm
44497: qm m m qp Qm Qp
38017: qm m m qp Qp Qm
42337: qm m m Qm qp Qp
22897: qm m m Qm Qp qp
34777: qm m m Qp qp Qm
21817: qm m m Qp Qm qp
44137: qm m qp m Qm Qp
37657: qm m qp m Qp Qm
44785: qm m qp p Qm Qp
38305: qm m qp p Qp Qm
39817: qm m qp Qm m Qp
43705: qm m qp Qm p Qp
7417: qm m qp Qm Qp m
30745: qm m qp Qm Qp p
32257: qm m qp Qp m Qm
36145: qm m qp Qp p Qm
6337: qm m qp Qp Qm m
29665: qm m qp Qp Qm p
44605: qm m p qp Qm Qp
38125: qm m p qp Qp Qm
42445: qm m p Qm qp Qp
23005: qm m p Qm Qp qp
34885: qm m p Qp qp Qm
21925: qm m p Qp Qm qp
41617: qm m Qm m qp Qp
22177: qm m Qm m Qp qp
39457: qm m Qm qp m Qp
43345: qm m Qm qp p Qp
7057: qm m Qm qp Qp m
30385: qm m Qm qp Qp p
42265: qm m Qm p qp Qp
22825: qm m Qm p Qp qp
16777: qm m Qm Qp m qp
3817: qm m Qm Qp qp m
27145: qm m Qm Qp qp p
20665: qm m Qm Qp p qp
33877: qm m Qp m qp Qm
20917: qm m Qp m Qm qp
31717: qm m Qp qp m Qm
35605: qm m Qp qp p Qm
5797: qm m Qp qp Qm m
29125: qm m Qp qp Qm p
34525: qm m Qp p qp Qm
21565: qm m Qp p Qm qp
16597: qm m Qp Qm m qp
3637: qm m Qp Qm qp m
26965: qm m Qp Qm qp p
20485: qm m Qp Qm p qp
44077: qm qp m m Qm Qp
37597: qm qp m m Qp Qm
44725: qm qp m p Qm Qp
38245: qm qp m p Qp Qm
39757: qm qp m Qm m Qp
43645: qm qp m Qm p Qp
7357: qm qp m Qm Qp m
30685: qm qp m Qm Qp p
32197: qm qp m Qp m Qm
36085: qm qp m Qp p Qm
6277: qm qp m Qp Qm m
29605: qm qp m Qp Qm p
44185: qm qp p m Qm Qp
37705: qm qp p m Qp Qm
44833: qm qp p p Qm Qp
38353: qm qp p p Qp Qm
39865: qm qp p Qm m Qp
43753: qm qp p Qm p Qp
7465: qm qp p Qm Qp m
30793: qm qp p Qm Qp p
32305: qm qp p Qp m Qm
36193: qm qp p Qp p Qm
6385: qm qp p Qp Qm m
29713: qm qp p Qp Qm p
39037: qm qp Qm m m Qp
42925: qm qp Qm m p Qp
6637: qm qp Qm m Qp m
29965: qm qp Qm m Qp p
39685: qm qp Qm p m Qp
43573: qm qp Qm p p Qp
7285: qm qp Qm p Qp m
30613: qm qp Qm p Qp p
1237: qm qp Qm Qp m m
24565: qm qp Qm Qp m p
5125: qm qp Qm Qp p m
28453: qm qp Qm Qp p p
31297: qm qp Qp m m Qm
35185: qm qp Qp m p Qm
5377: qm qp Qp m Qm m
28705: qm qp Qp m Qm p
31945: qm qp Qp p m Qm
35833: qm qp Qp p p Qm
6025: qm qp Qp p Qm m
29353: qm qp Qp p Qm p
1057: qm qp Qp Qm m m
24385: qm qp Qp Qm m p
4945: qm qp Qp Qm p m
28273: qm qp Qp Qm p p
44515: qm p m qp Qm Qp
38035: qm p m qp Qp Qm
42355: qm p m Qm qp Qp
22915: qm p m Qm Qp qp
34795: qm p m Qp qp Qm
21835: qm p m Qp Qm qp
44155: qm p qp m Qm Qp
37675: qm p qp m Qp Qm
44803: qm p qp p Qm Qp
38323: qm p qp p Qp Qm
39835: qm p qp Qm m Qp
43723: qm p qp Qm p Qp
7435: qm p qp Qm Qp m
30763: qm p qp Qm Qp p
32275: qm p qp Qp m Qm
36163: qm p qp Qp p Qm
6355: qm p qp Qp Qm m
29683: qm p qp Qp Qm p
44623: qm p p qp Qm Qp
38143: qm p p qp Qp Qm
42463: qm p p Qm qp Qp
23023: qm p p Qm Qp qp
34903: qm p p Qp qp Qm
21943: qm p p Qp Qm qp
41635: qm p Qm m qp Qp
22195: qm p Qm m Qp qp
39475: qm p Qm qp m Qp
43363: qm p Qm qp p Qp
7075: qm p Qm qp Qp m
30403: qm p Qm qp Qp p
42283: qm p Qm p qp Qp
22843: qm p Qm p Qp qp
16795: qm p Qm Qp m qp
3835: qm p Qm Qp qp m
27163: qm p Qm Qp qp p
20683: qm p Qm Qp p qp
33895: qm p Qp m qp Qm
20935: qm p Qp m Qm qp
31735: qm p Qp qp m Qm
35623: qm p Qp qp p Qm
5815: qm p Qp qp Qm m
29143: qm p Qp qp Qm p
34543: qm p Qp p qp Qm
21583: qm p Qp p Qm qp
16615: qm p Qp Qm m qp
3655: qm p Qp Qm qp m
26983: qm p Qp Qm qp p
20503: qm p Qp Qm p qp
41497: qm Qm m m qp Qp
22057: qm Qm m m Qp qp
39337: qm Qm m qp m Qp
43225: qm Qm m qp p Qp
6937: qm Qm m qp Qp m
30265: qm Qm m qp Qp p
42145: qm Qm m p qp Qp
22705: qm Qm m p Qp qp
16657: qm Qm m Qp m qp
3697: qm Qm m Qp qp m
27025: qm Qm m Qp qp p
20545: qm Qm m Qp p qp
38977: qm Qm qp m m Qp
42865: qm Qm qp m p Qp
6577: qm Qm qp m Qp m
29905: qm Qm qp m Qp p
39625: qm Qm qp p m Qp
43513: qm Qm qp p p Qp
7225: qm Qm qp p Qp m
30553: qm Qm qp p Qp p
1177: qm Qm qp Qp m m
24505: qm Qm qp Qp m p
5065: qm Qm qp Qp p m
28393: qm Qm qp Qp p p
41605: qm Qm p m qp Qp
22165: qm Qm p m Qp qp
39445: qm Qm p qp m Qp
43333: qm Qm p qp p Qp
7045: qm Qm p qp Qp m
30373: qm Qm p qp Qp p
42253: qm Qm p p qp Qp
22813: qm Qm p p Qp qp
16765: qm Qm p Qp m qp
3805: qm Qm p Qp qp m
27133: qm Qm p Qp qp p
20653: qm Qm p Qp p qp
15757: qm Qm Qp m m qp
2797: qm Qm Qp m qp m
26125: qm Qm Qp m qp p
19645: qm Qm Qp m p qp
637: qm Qm Qp qp m m
23965: qm Qm Qp qp m p
4525: qm Qm Qp qp p m
27853: qm Qm Qp qp p p
16405: qm Qm Qp p m qp
3445: qm Qm Qp p qp m
26773: qm Qm Qp p qp p
20293: qm Qm Qp p p qp
33727: qm Qp m m qp Qm
20767: qm Qp m m Qm qp
31567: qm Qp m qp m Qm
35455: qm Qp m qp p Qm
5647: qm Qp m qp Qm m
28975: qm Qp m qp Qm p
34375: qm Qp m p qp Qm
21415: qm Qp m p Qm qp
16447: qm Qp m Qm m qp
3487: qm Qp m Qm qp m
26815: qm Qp m Qm qp p
20335: qm Qp m Qm p qp
31207: qm Qp qp m m Qm
35095: qm Qp qp m p Qm
5287: qm Qp qp m Qm m
28615: qm Qp qp m Qm p
31855: qm Qp qp p m Qm
35743: qm Qp qp p p Qm
5935: qm Qp qp p Qm m
29263: qm Qp qp p Qm p
967: qm Qp qp Qm m m
24295: qm Qp qp Qm m p
4855: qm Qp qp Qm p m
28183: qm Qp qp Qm p p
33835: qm Qp p m qp Qm
20875: qm Qp p m Qm qp
31675: qm Qp p qp m Qm
35563: qm Qp p qp p Qm
5755: qm Qp p qp Qm m
29083: qm Qp p qp Qm p
34483: qm Qp p p qp Qm
21523: qm Qp p p Qm qp
16555: qm Qp p Qm m qp
3595: qm Qp p Qm qp m
26923: qm Qp p Qm qp p
20443: qm Qp p Qm p qp
15727: qm Qp Qm m m qp
2767: qm Qp Qm m qp m
26095: qm Qp Qm m qp p
19615: qm Qp Qm m p qp
607: qm Qp Qm qp m m
23935: qm Qp Qm qp m p
4495: qm Qp Qm qp p m
27823: qm Qp Qm qp p p
16375: qm Qp Qm p m qp
3415: qm Qp Qm p qp m
26743: qm Qp Qm p qp p
20263: qm Qp Qm p p qp
44282: qp m m qm Qm Qp
37802: qp m m qm Qp Qm
41042: qp m m Qm qm Qp
15122: qp m m Qm Qp qm
33482: qp m m Qp qm Qm
14042: qp m m Qp Qm qm
44102: qp m qm m Qm Qp
37622: qp m qm m Qp Qm
44750: qp m qm p Qm Qp
38270: qp m qm p Qp Qm
39782: qp m qm Qm m Qp
43670: qp m qm Qm p Qp
7382: qp m qm Qm Qp m
30710: qp m qm Qm Qp p
32222: qp m qm Qp m Qm
36110: qp m qm Qp p Qm
6302: qp m qm Qp Qm m
29630: qp m qm Qp Qm p
44390: qp m p qm Qm Qp
37910: qp m p qm Qp Qm
41150: qp m p Qm qm Qp
15230: qp m p Qm Qp qm
33590: qp m p Qp qm Qm
14150: qp m p Qp Qm qm
40322: qp m Qm m qm Qp
14402: qp m Qm m Qp qm
39242: qp m Qm qm m Qp
43130: qp m Qm qm p Qp
6842: qp m Qm qm Qp m
30170: qp m Qm qm Qp p
40970: qp m Qm p qm Qp
15050: qp m Qm p Qp qm
9002: qp m Qm Qp m qm
2522: qp m Qm Qp qm m
25850: qp m Qm Qp qm p
12890: qp m Qm Qp p qm
32582: qp m Qp m qm Qm
13142: qp m Qp m Qm qm
31502: qp m Qp qm m Qm
35390: qp m Qp qm p Qm
5582: qp m Qp qm Qm m
28910: qp m Qp qm Qm p
33230: qp m Qp p qm Qm
13790: qp m Qp p Qm qm
8822: qp m Qp Qm m qm
2342: qp m Qp Qm qm m
25670: qp m Qp Qm qm p
12710: qp m Qp Qm p qm
44072: qp qm m m Qm Qp
37592: qp qm m m Qp Qm
44720: qp qm m p Qm Qp
38240: qp qm m p Qp Qm
39752: qp qm m Qm m Qp
43640: qp qm m Qm p Qp
7352: qp qm m Qm Qp m
30680: qp qm m Qm Qp p
32192: qp qm m Qp m Qm
36080: qp qm m Qp p Qm
6272: qp qm m Qp Qm m
29600: qp qm m Qp Qm p
44180: qp qm p m Qm Qp
37700: qp qm p m Qp Qm
44828: qp qm p p Qm Qp
38348: qp qm p p Qp Qm
39860: qp qm p Qm m Qp
43748: qp qm p Qm p Qp
7460: qp qm p Qm Qp m
30788: qp qm p Qm Qp p
32300: qp qm p Qp m Qm
36188: qp qm p Qp p Qm
6380: qp qm p Qp Qm m
29708: qp qm p Qp Qm p
39032: qp qm Qm m m Qp
42920: qp qm Qm m p Qp
6632: qp qm Qm m Qp m
29960: qp qm Qm m Qp p
39680: qp qm Qm p m Qp
43568: qp qm Qm p p Qp
7280: qp qm Qm p Qp m
30608: qp qm Qm p Qp p
1232: qp qm Qm Qp m m
24560: qp qm Qm Qp m p
5120: qp qm Qm Qp p m
28448: qp qm Qm Qp p p
31292: qp qm Qp m m Qm
35180: qp qm Qp m p Qm
5372: qp qm Qp m Qm m
28700: qp qm Qp m Qm p
31940: qp qm Qp p m Qm
35828: qp qm Qp p p Qm
6020: qp qm Qp p Qm m
29348: qp qm Qp p Qm p
1052: qp qm Qp Qm m m
24380: qp qm Qp Qm m p
4940: qp qm Qp Qm p m
28268: qp qm Qp Qm p p
44300: qp p m qm Qm Qp
37820: qp p m qm Qp Qm
41060: qp p m Qm qm Qp
15140: qp p m Qm Qp qm
33500: qp p m Qp qm Qm
14060: qp p m Qp Qm qm
44120: qp p qm m Qm Qp
37640: qp p qm m Qp Qm
44768: qp p qm p Qm Qp
38288: qp p qm p Qp Qm
39800: qp p qm Qm m Qp
43688: qp p qm Qm p Qp
7400: qp p qm Qm Qp m
30728: qp p qm Qm Qp p
32240: qp p qm Qp m Qm
36128: qp p qm Qp p Qm
6320: qp p qm Qp Qm m
29648: qp p qm Qp Qm p
44408: qp p p qm Qm Qp
37928: qp p p qm Qp Qm
41168: qp p p Qm qm Qp
15248: qp p p Qm Qp qm
33608: qp p p Qp qm Qm
14168: qp p p Qp Qm qm
40340: qp p Qm m qm Qp
14420: qp p Qm m Qp qm
39260: qp p Qm qm m Qp
43148: qp p Qm qm p Qp
6860: qp p Qm qm Qp m
30188: qp p Qm qm Qp p
40988: qp p Qm p qm Qp
15068: qp p Qm p Qp qm
9020: qp p Qm Qp m qm
2540: qp p Qm Qp qm m
25868: qp p Qm Qp qm p
12908: qp p Qm Qp p qm
32600: qp p Qp m qm Qm
13160: qp p Qp m Qm qm
31520: qp p Qp qm m Qm
35408: qp p Qp qm p Qm
5600: qp p Qp qm Qm m
28928: qp p Qp qm Qm p
33248: qp p Qp p qm Qm
13808: qp p Qp p Qm qm
8840: qp p Qp Qm m qm
2360: qp p Qp Qm qm m
25688: qp p Qp Qm qm p
12728: qp p Qp Qm p qm
40202: qp Qm m m qm Qp
14282: qp Qm m m Qp qm
39122: qp Qm m qm m Qp
43010: qp Qm m qm p Qp
6722: qp Qm m qm Qp m
30050: qp Qm m qm Qp p
40850: qp Qm m p qm Qp
14930: qp Qm m p Qp qm
8882: qp Qm m Qp m qm
2402: qp Qm m Qp qm m
25730: qp Qm m Qp qm p
12770: qp Qm m Qp p qm
38942: qp Qm qm m m Qp
42830: qp Qm qm m p Qp
6542: qp Qm qm m Qp m
29870: qp Qm qm m Qp p
39590: qp Qm qm p m Qp
43478: qp Qm qm p p Qp
7190: qp Qm qm p Qp m
30518: qp Qm qm p Qp p
1142: qp Qm qm Qp m m
24470: qp Qm qm Qp m p
5030: qp Qm qm Qp p m
28358: qp Qm qm Qp p p
40310: qp Qm p m qm Qp
14390: qp Qm p m Qp qm
39230: qp Qm p qm m Qp
43118: qp Qm p qm p Qp
6830: qp Qm p qm Qp m
30158: qp Qm p qm Qp p
40958: qp Qm p p qm Qp
15038: qp Qm p p Qp qm
8990: qp Qm p Qp m qm
2510: qp Qm p Qp qm m
25838: qp Qm p Qp qm p
12878: qp Qm p Qp p qm
7982: qp Qm Qp m m qm
1502: qp Qm Qp m qm m
24830: qp Qm Qp m qm p
11870: qp Qm Qp m p qm
422: qp Qm Qp qm m m
23750: qp Qm Qp qm m p
4310: qp Qm Qp qm p m
27638: qp Qm Qp qm p p
8630: qp Qm Qp p m qm
2150: qp Qm Qp p qm m
25478: qp Qm Qp p qm p
12518: qp Qm Qp p p qm
32432: qp Qp m m qm Qm
12992: qp Qp m m Qm qm
31352: qp Qp m qm m Qm
35240: qp Qp m qm p Qm
5432: qp Qp m qm Qm m
28760: qp Qp m qm Qm p
33080: qp Qp m p qm Qm
13640: qp Qp m p Qm qm
8672: qp Qp m Qm m qm
2192: qp Qp m Qm qm m
25520: qp Qp m Qm qm p
12560: qp Qp m Qm p qm
31172: qp Qp qm m m Qm
35060: qp Qp qm m p Qm
5252: qp Qp qm m Qm m
28580: qp Qp qm m Qm p
31820: qp Qp qm p m Qm
35708: qp Qp qm p p Qm
5900: qp Qp qm p Qm m
29228: qp Qp qm p Qm p
932: qp Qp qm Qm m m
24260: qp Qp qm Qm m p
4820: qp Qp qm Qm p m
28148: qp Qp qm Qm p p
32540: qp Qp p m qm Qm
13100: qp Qp p m Qm qm
31460: qp Qp p qm m Qm
35348: qp Qp p qm p Qm
5540: qp Qp p qm Qm m
28868: qp Qp p qm Qm p
33188: qp Qp p p qm Qm
13748: qp Qp p p Qm qm
8780: qp Qp p Qm m qm
2300: qp Qp p Qm qm m
25628: qp Qp p Qm qm p
12668: qp Qp p Qm p qm
7952: qp Qp Qm m m qm
1472: qp Qp Qm m qm m
24800: qp Qp Qm m qm p
11840: qp Qp Qm m p qm
392: qp Qp Qm qm m m
23720: qp Qp Qm qm m p
4280: qp Qp Qm qm p m
27608: qp Qp Qm qm p p
8600: qp Qp Qm p m qm
2120: qp Qp Qm p qm m
25448: qp Qp Qm p qm p
12488: qp Qp Qm p p qm
44535: p m qm qp Qm Qp
38055: p m qm qp Qp Qm
42375: p m qm Qm qp Qp
22935: p m qm Qm Qp qp
34815: p m qm Qp qp Qm
21855: p m qm Qp Qm qp
44355: p m qp qm Qm Qp
37875: p m qp qm Qp Qm
41115: p m qp Qm qm Qp
15195: p m qp Qm Qp qm
33555: p m qp Qp qm Qm
14115: p m qp Qp Qm qm
41835: p m Qm qm qp Qp
22395: p m Qm qm Qp qp
40755: p m Qm qp qm Qp
14835: p m Qm qp Qp qm
18075: p m Qm Qp qm qp
11595: p m Qm Qp qp qm
34095: p m Qp qm qp Qm
21135: p m Qp qm Qm qp
33015: p m Qp qp qm Qm
13575: p m Qp qp Qm qm
17895: p m Qp Qm qm qp
11415: p m Qp Qm qp qm
44505: p qm m qp Qm Qp
38025: p qm m qp Qp Qm
42345: p qm m Qm qp Qp
22905: p qm m Qm Qp qp
34785: p qm m Qp qp Qm
21825: p qm m Qp Qm qp
44145: p qm qp m Qm Qp
37665: p qm qp m Qp Qm
44793: p qm qp p Qm Qp
38313: p qm qp p Qp Qm
39825: p qm qp Qm m Qp
43713: p qm qp Qm p Qp
7425: p qm qp Qm Qp m
30753: p qm qp Qm Qp p
32265: p qm qp Qp m Qm
36153: p qm qp Qp p Qm
6345: p qm qp Qp Qm m
29673: p qm qp Qp Qm p
44613: p qm p qp Qm Qp
38133: p qm p qp Qp Qm
42453: p qm p Qm qp Qp
23013: p qm p Qm Qp qp
34893: p qm p Qp qp Qm
21933: p qm p Qp Qm qp
41625: p qm Qm m qp Qp
22185: p qm Qm m Qp qp
39465: p qm Qm qp m Qp
43353: p qm Qm qp p Qp
7065: p qm Qm qp Qp m
30393: p qm Qm qp Qp p
42273: p qm Qm p qp Qp
22833: p qm Qm p Qp qp
16785: p qm Qm Qp m qp
3825: p qm Qm Qp qp m
27153: p qm Qm Qp qp p
20673: p qm Qm Qp p qp
33885: p qm Qp m qp Qm
20925: p qm Qp m Qm qp
31725: p qm Qp qp m Qm
35613: p qm Qp qp p Qm
5805: p qm Qp qp Qm m
29133: p qm Qp qp Qm p
34533: p qm Qp p qp Qm
21573: p qm Qp p Qm qp
16605: p qm Qp Qm m qp
3645: p qm Qp Qm qp m
26973: p qm Qp Qm qp p
20493: p qm Qp Qm p qp
44295: p qp m qm Qm Qp
37815: p qp m qm Qp Qm
41055: p qp m Qm qm Qp
15135: p qp m Qm Qp qm
33495: p qp m Qp qm Qm
14055: p qp m Qp Qm qm
44115: p qp qm m Qm Qp
37635: p qp qm m Qp Qm
44763: p qp qm p Qm Qp
38283: p qp qm p Qp Qm
39795: p qp qm Qm m Qp
43683: p qp qm Qm p Qp
7395: p qp qm Qm Qp m
30723: p qp qm Qm Qp p
32235: p qp qm Qp m Qm
36123: p qp qm Qp p Qm
6315: p qp qm Qp Qm m
29643: p qp qm Qp Qm p
44403: p qp p qm Qm Qp
37923: p qp p qm Qp Qm
41163: p qp p Qm qm Qp
15243: p qp p Qm Qp qm
33603: p qp p Qp qm Qm
14163: p qp p Qp Qm qm
40335: p qp Qm m qm Qp
14415: p qp Qm m Qp qm
39255: p qp Qm qm m Qp
43143: p qp Qm qm p Qp
6855: p qp Qm qm Qp m
30183: p qp Qm qm Qp p
40983: p qp Qm p qm Qp
15063: p qp Qm p Qp qm
9015: p qp Qm Qp m qm
2535: p qp Qm Qp qm m
25863: p qp Qm Qp qm p
12903: p qp Qm Qp p qm
32595: p qp Qp m qm Qm
13155: p qp Qp m Qm qm
31515: p qp Qp qm m Qm
35403: p qp Qp qm p Qm
5595: p qp Qp qm Qm m
28923: p qp Qp qm Qm p
33243: p qp Qp p qm Qm
13803: p qp Qp p Qm qm
8835: p qp Qp Qm m qm
2355: p qp Qp Qm qm m
25683: p qp Qp Qm qm p
12723: p qp Qp Qm p qm
44553: p p qm qp Qm Qp
38073: p p qm qp Qp Qm
42393: p p qm Qm qp Qp
22953: p p qm Qm Qp qp
34833: p p qm Qp qp Qm
21873: p p qm Qp Qm qp
44373: p p qp qm Qm Qp
37893: p p qp qm Qp Qm
41133: p p qp Qm qm Qp
15213: p p qp Qm Qp qm
33573: p p qp Qp qm Qm
14133: p p qp Qp Qm qm
41853: p p Qm qm qp Qp
22413: p p Qm qm Qp qp
40773: p p Qm qp qm Qp
14853: p p Qm qp Qp qm
18093: p p Qm Qp qm qp
11613: p p Qm Qp qp qm
34113: p p Qp qm qp Qm
21153: p p Qp qm Qm qp
33033: p p Qp qp qm Qm
13593: p p Qp qp Qm qm
17913: p p Qp Qm qm qp
11433: p p Qp Qm qp qm
41715: p Qm m qm qp Qp
22275: p Qm m qm Qp qp
40635: p Qm m qp qm Qp
14715: p Qm m qp Qp qm
17955: p Qm m Qp qm qp
11475: p Qm m Qp qp qm
41535: p Qm qm m qp Qp
22095: p Qm qm m Qp qp
39375: p Qm qm qp m Qp
43263: p Qm qm qp p Qp
6975: p Qm qm qp Qp m
30303: p Qm qm qp Qp p
42183: p Qm qm p qp Qp
22743: p Qm qm p Qp qp
16695: p Qm qm Qp m qp
3735: p Qm qm Qp qp m
27063: p Qm qm Qp qp p
20583: p Qm qm Qp p qp
40275: p Qm qp m qm Qp
14355: p Qm qp m Qp qm
39195: p Qm qp qm m Qp
43083: p Qm qp qm p Qp
6795: p Qm qp qm Qp m
30123: p Qm qp qm Qp p
40923: p Qm qp p qm Qp
15003: p Qm qp p Qp qm
8955: p Qm qp Qp m qm
2475: p Qm qp Qp qm m
25803: p Qm qp Qp qm p
12843: p Qm qp Qp p qm
41823: p Qm p qm qp Qp
22383: p Qm p qm Qp qp
40743: p Qm p qp qm Qp
14823: p Qm p qp Qp qm
18063: p Qm p Qp qm qp
11583: p Qm p Qp qp qm
17055: p Qm Qp m qm qp
10575: p Qm Qp m qp qm
15975: p Qm Qp qm m qp
3015: p Qm Qp qm qp m
26343: p Qm Qp qm qp p
19863: p Qm Qp qm p qp
8415: p Qm Qp qp m qm
1935: p Qm Qp qp qm m
25263: p Qm Qp qp qm p
12303: p Qm Qp qp p qm
17703: p Qm Qp p qm qp
11223: p Qm Qp p qp qm
33945: p Qp m qm qp Qm
20985: p Qp m qm Qm qp
32865: p Qp m qp qm Qm
13425: p Qp m qp Qm qm
17745: p Qp m Qm qm qp
11265: p Qp m Qm qp qm
33765: p Qp qm m qp Qm
20805: p Qp qm m Qm qp
31605: p Qp qm qp m Qm
35493: p Qp qm qp p Qm
5685: p Qp qm qp Qm m
29013: p Qp qm qp Qm p
34413: p Qp qm p qp Qm
21453: p Qp qm p Qm qp
16485: p Qp qm Qm m qp
3525: p Qp qm Qm qp m
26853: p Qp qm Qm qp p
20373: p Qp qm Qm p qp
32505: p Qp qp m qm Qm
13065: p Qp qp m Qm qm
31425: p Qp qp qm m Qm
35313: p Qp qp qm p Qm
5505: p Qp qp qm Qm m
28833: p Qp qp qm Qm p
33153: p Qp qp p qm Qm
13713: p Qp qp p Qm qm
8745: p Qp qp Qm m qm
2265: p Qp qp Qm qm m
25593: p Qp qp Qm qm p
12633: p Qp qp Qm p qm
34053: p Qp p qm qp Qm
21093: p Qp p qm Qm qp
32973: p Qp p qp qm Qm
13533: p Qp p qp Qm qm
17853: p Qp p Qm qm qp
11373: p Qp p Qm qp qm
17025: p Qp Qm m qm qp
10545: p Qp Qm m qp qm
15945: p Qp Qm qm m qp
2985: p Qp Qm qm qp m
26313: p Qp Qm qm qp p
19833: p Qp Qm qm p qp
8385: p Qp Qm qp m qm
1905: p Qp Qm qp qm m
25233: p Qp Qm qp qm p
12273: p Qp Qm qp p qm
17673: p Qp Qm p qm qp
11193: p Qp Qm p qp qm
41692: Qm m m qm qp Qp
22252: Qm m m qm Qp qp
40612: Qm m m qp qm Qp
14692: Qm m m qp Qp qm
17932: Qm m m Qp qm qp
11452: Qm m m Qp qp qm
41512: Qm m qm m qp Qp
22072: Qm m qm m Qp qp
39352: Qm m qm qp m Qp
43240: Qm m qm qp p Qp
6952: Qm m qm qp Qp m
30280: Qm m qm qp Qp p
42160: Qm m qm p qp Qp
22720: Qm m qm p Qp qp
16672: Qm m qm Qp m qp
3712: Qm m qm Qp qp m
27040: Qm m qm Qp qp p
20560: Qm m qm Qp p qp
40252: Qm m qp m qm Qp
14332: Qm m qp m Qp qm
39172: Qm m qp qm m Qp
43060: Qm m qp qm p Qp
6772: Qm m qp qm Qp m
30100: Qm m qp qm Qp p
40900: Qm m qp p qm Qp
14980: Qm m qp p Qp qm
8932: Qm m qp Qp m qm
2452: Qm m qp Qp qm m
25780: Qm m qp Qp qm p
12820: Qm m qp Qp p qm
41800: Qm m p qm qp Qp
22360: Qm m p qm Qp qp
40720: Qm m p qp qm Qp
14800: Qm m p qp Qp qm
18040: Qm m p Qp qm qp
11560: Qm m p Qp qp qm
17032: Qm m Qp m qm qp
10552: Qm m Qp m qp qm
15952: Qm m Qp qm m qp
2992: Qm m Qp qm qp m
26320: Qm m Qp qm qp p
19840: Qm m Qp qm p qp
8392: Qm m Qp qp m qm
1912: Qm m Qp qp qm m
25240: Qm m Qp qp qm p
12280: Qm m Qp qp p qm
17680: Qm m Qp p qm qp
11200: Qm m Qp p qp qm
41482: Qm qm m m qp Qp
22042: Qm qm m m Qp qp
39322: Qm qm m qp m Qp
43210: Qm qm m qp p Qp
6922: Qm qm m qp Qp m
30250: Qm qm m qp Qp p
42130: Qm qm m p qp Qp
22690: Qm qm m p Qp qp
16642: Qm qm m Qp m qp
3682: Qm qm m Qp qp m
27010: Qm qm m Qp qp p
20530: Qm qm m Qp p qp
38962: Qm qm qp m m Qp
42850: Qm qm qp m p Qp
6562: Qm qm qp m Qp m
29890: Qm qm qp m Qp p
39610: Qm qm qp p m Qp
43498: Qm qm qp p p Qp
7210: Qm qm qp p Qp m
30538: Qm qm qp p Qp p
1162: Qm qm qp Qp m m
24490: Qm qm qp Qp m p
5050: Qm qm qp Qp p m
28378: Qm qm qp Qp p p
41590: Qm qm p m qp Qp
22150: Qm qm p m Qp qp
39430: Qm qm p qp m Qp
43318: Qm qm p qp p Qp
7030: Qm qm p qp Qp m
30358: Qm qm p qp Qp p
42238: Qm qm p p qp Qp
22798: Qm qm p p Qp qp
16750: Qm qm p Qp m qp
3790: Qm qm p Qp qp m
27118: Qm qm p Qp qp p
20638: Qm qm p Qp p qp
15742: Qm qm Qp m m qp
2782: Qm qm Qp m qp m
26110: Qm qm Qp m qp p
19630: Qm qm Qp m p qp
622: Qm qm Qp qp m m
23950: Qm qm Qp qp m p
4510: Qm qm Qp qp p m
27838: Qm qm Qp qp p p
16390: Qm qm Qp p m qp
3430: Qm qm Qp p qp m
26758: Qm qm Qp p qp p
20278: Qm qm Qp p p qp
40192: Qm qp m m qm Qp
14272: Qm qp m m Qp qm
39112: Qm qp m qm m Qp
43000: Qm qp m qm p Qp
6712: Qm qp m qm Qp m
30040: Qm qp m qm Qp p
40840: Qm qp m p qm Qp
14920: Qm qp m p Qp qm
8872: Qm qp m Qp m qm
2392: Qm qp m Qp qm m
25720: Qm qp m Qp qm p
12760: Qm qp m Qp p qm
38932: Qm qp qm m m Qp
42820: Qm qp qm m p Qp
6532: Qm qp qm m Qp m
29860: Qm qp qm m Qp p
39580: Qm qp qm p m Qp
43468: Qm qp qm p p Qp
7180: Qm qp qm p Qp m
30508: Qm qp qm p Qp p
1132: Qm qp qm Qp m m
24460: Qm qp qm Qp m p
5020: Qm qp qm Qp p m
28348: Qm qp qm Qp p p
40300: Qm qp p m qm Qp
14380: Qm qp p m Qp qm
39220: Qm qp p qm m Qp
43108: Qm qp p qm p Qp
6820: Qm qp p qm Qp m
30148: Qm qp p qm Qp p
40948: Qm qp p p qm Qp
15028: Qm qp p p Qp qm
8980: Qm qp p Qp m qm
2500: Qm qp p Qp qm m
25828: Qm qp p Qp qm p
12868: Qm qp p Qp p qm
7972: Qm qp Qp m m qm
1492: Qm qp Qp m qm m
24820: Qm qp Qp m qm p
11860: Qm qp Qp m p qm
412: Qm qp Qp qm m m
23740: Qm qp Qp qm m p
4300: Qm qp Qp qm p m
27628: Qm qp Qp qm p p
8620: Qm qp Qp p m qm
2140: Qm qp Qp p qm m
25468: Qm qp Qp p qm p
12508: Qm qp Qp p p qm
41710: Qm p m qm qp Qp
22270: Qm p m qm Qp qp
40630: Qm p m qp qm Qp
14710: Qm p m qp Qp qm
17950: Qm p m Qp qm qp
11470: Qm p m Qp qp qm
41530: Qm p qm m qp Qp
22090: Qm p qm m Qp qp
39370: Qm p qm qp m Qp
43258: Qm p qm qp p Qp
6970: Qm p qm qp Qp m
30298: Qm p qm qp Qp p
42178: Qm p qm p qp Qp
22738: Qm p qm p Qp qp
16690: Qm p qm Qp m qp
3730: Qm p qm Qp qp m
27058: Qm p qm Qp qp p
20578: Qm p qm Qp p qp
40270: Qm p qp m qm Qp
14350: Qm p qp m Qp qm
39190: Qm p qp qm m Qp
43078: Qm p qp qm p Qp
6790: Qm p qp qm Qp m
30118: Qm p qp qm Qp p
40918: Qm p qp p qm Qp
14998: Qm p qp p Qp qm
8950: Qm p qp Qp m qm
2470: Qm p qp Qp qm m
25798: Qm p qp Qp qm p
12838: Qm p qp Qp p qm
41818: Qm p p qm qp Qp
22378: Qm p p qm Qp qp
40738: Qm p p qp qm Qp
14818: Qm p p qp Qp qm
18058: Qm p p Qp qm qp
11578: Qm p p Qp qp qm
17050: Qm p Qp m qm qp
10570: Qm p Qp m qp qm
15970: Qm p Qp qm m qp
3010: Qm p Qp qm qp m
26338: Qm p Qp qm qp p
19858: Qm p Qp qm p qp
8410: Qm p Qp qp m qm
1930: Qm p Qp qp qm m
25258: Qm p Qp qp qm p
12298: Qm p Qp qp p qm
17698: Qm p Qp p qm qp
11218: Qm p Qp p qp qm
16882: Qm Qp m m qm qp
10402: Qm Qp m m qp qm
15802: Qm Qp m qm m qp
2842: Qm Qp m qm qp m
26170: Qm Qp m qm qp p
19690: Qm Qp m qm p qp
8242: Qm Qp m qp m qm
1762: Qm Qp m qp qm m
25090: Qm Qp m qp qm p
12130: Qm Qp m qp p qm
17530: Qm Qp m p qm qp
11050: Qm Qp m p qp qm
15622: Qm Qp qm m m qp
2662: Qm Qp qm m qp m
25990: Qm Qp qm m qp p
19510: Qm Qp qm m p qp
502: Qm Qp qm qp m m
23830: Qm Qp qm qp m p
4390: Qm Qp qm qp p m
27718: Qm Qp qm qp p p
16270: Qm Qp qm p m qp
3310: Qm Qp qm p qp m
26638: Qm Qp qm p qp p
20158: Qm Qp qm p p qp
7882: Qm Qp qp m m qm
1402: Qm Qp qp m qm m
24730: Qm Qp qp m qm p
11770: Qm Qp qp m p qm
322: Qm Qp qp qm m m
23650: Qm Qp qp qm m p
4210: Qm Qp qp qm p m
27538: Qm Qp qp qm p p
8530: Qm Qp qp p m qm
2050: Qm Qp qp p qm m
25378: Qm Qp qp p qm p
12418: Qm Qp qp p p qm
16990: Qm Qp p m qm qp
10510: Qm Qp p m qp qm
15910: Qm Qp p qm m qp
2950: Qm Qp p qm qp m
26278: Qm Qp p qm qp p
19798: Qm Qp p qm p qp
8350: Qm Qp p qp m qm
1870: Qm Qp p qp qm m
25198: Qm Qp p qp qm p
12238: Qm Qp p qp p qm
17638: Qm Qp p p qm qp
11158: Qm Qp p p qp qm
33917: Qp m m qm qp Qm
20957: Qp m m qm Qm qp
32837: Qp m m qp qm Qm
13397: Qp m m qp Qm qm
17717: Qp m m Qm qm qp
11237: Qp m m Qm qp qm
33737: Qp m qm m qp Qm
20777: Qp m qm m Qm qp
31577: Qp m qm qp m Qm
35465: Qp m qm qp p Qm
5657: Qp m qm qp Qm m
28985: Qp m qm qp Qm p
34385: Qp m qm p qp Qm
21425: Qp m qm p Qm qp
16457: Qp m qm Qm m qp
3497: Qp m qm Qm qp m
26825: Qp m qm Qm qp p
20345: Qp m qm Qm p qp
32477: Qp m qp m qm Qm
13037: Qp m qp m Qm qm
31397: Qp m qp qm m Qm
35285: Qp m qp qm p Qm
5477: Qp m qp qm Qm m
28805: Qp m qp qm Qm p
33125: Qp m qp p qm Qm
13685: Qp m qp p Qm qm
8717: Qp m qp Qm m qm
2237: Qp m qp Qm qm m
25565: Qp m qp Qm qm p
12605: Qp m qp Qm p qm
34025: Qp m p qm qp Qm
21065: Qp m p qm Qm qp
32945: Qp m p qp qm Qm
13505: Qp m p qp Qm qm
17825: Qp m p Qm qm qp
11345: Qp m p Qm qp qm
16997: Qp m Qm m qm qp
10517: Qp m Qm m qp qm
15917: Qp m Qm qm m qp
2957: Qp m Qm qm qp m
26285: Qp m Qm qm qp p
19805: Qp m Qm qm p qp
8357: Qp m Qm qp m qm
1877: Qp m Qm qp qm m
25205: Qp m Qm qp qm p
12245: Qp m Qm qp p qm
17645: Qp m Qm p qm qp
11165: Qp m Qm p qp qm
33707: Qp qm m m qp Qm
20747: Qp qm m m Qm qp
31547: Qp qm m qp m Qm
35435: Qp qm m qp p Qm
5627: Qp qm m qp Qm m
28955: Qp qm m qp Qm p
34355: Qp qm m p qp Qm
21395: Qp qm m p Qm qp
16427: Qp qm m Qm m qp
3467: Qp qm m Qm qp m
26795: Qp qm m Qm qp p
20315: Qp qm m Qm p qp
31187: Qp qm qp m m Qm
35075: Qp qm qp m p Qm
5267: Qp qm qp m Qm m
28595: Qp qm qp m Qm p
31835: Qp qm qp p m Qm
35723: Qp qm qp p p Qm
5915: Qp qm qp p Qm m
29243: Qp qm qp p Qm p
947: Qp qm qp Qm m m
24275: Qp qm qp Qm m p
4835: Qp qm qp Qm p m
28163: Qp qm qp Qm p p
33815: Qp qm p m qp Qm
20855: Qp qm p m Qm qp
31655: Qp qm p qp m Qm
35543: Qp qm p qp p Qm
5735: Qp qm p qp Qm m
29063: Qp qm p qp Qm p
34463: Qp qm p p qp Qm
21503: Qp qm p p Qm qp
16535: Qp qm p Qm m qp
3575: Qp qm p Qm qp m
26903: Qp qm p Qm qp p
20423: Qp qm p Qm p qp
15707: Qp qm Qm m m qp
2747: Qp qm Qm m qp m
26075: Qp qm Qm m qp p
19595: Qp qm Qm m p qp
587: Qp qm Qm qp m m
23915: Qp qm Qm qp m p
4475: Qp qm Qm qp p m
27803: Qp qm Qm qp p p
16355: Qp qm Qm p m qp
3395: Qp qm Qm p qp m
26723: Qp qm Qm p qp p
20243: Qp qm Qm p p qp
32417: Qp qp m m qm Qm
12977: Qp qp m m Qm qm
31337: Qp qp m qm m Qm
35225: Qp qp m qm p Qm
5417: Qp qp m qm Qm m
28745: Qp qp m qm Qm p
33065: Qp qp m p qm Qm
13625: Qp qp m p Qm qm
8657: Qp qp m Qm m qm
2177: Qp qp m Qm qm m
25505: Qp qp m Qm qm p
12545: Qp qp m Qm p qm
31157: Qp qp qm m m Qm
35045: Qp qp qm m p Qm
5237: Qp qp qm m Qm m
28565: Qp qp qm m Qm p
31805: Qp qp qm p m Qm
35693: Qp qp qm p p Qm
5885: Qp qp qm p Qm m
29213: Qp qp qm p Qm p
917: Qp qp qm Qm m m
24245: Qp qp qm Qm m p
4805: Qp qp qm Qm p m
28133: Qp qp qm Qm p p
32525: Qp qp p m qm Qm
13085: Qp qp p m Qm qm
31445: Qp qp p qm m Qm
35333: Qp qp p qm p Qm
5525: Qp qp p qm Qm m
28853: Qp qp p qm Qm p
33173: Qp qp p p qm Qm
13733: Qp qp p p Qm qm
8765: Qp qp p Qm m qm
2285: Qp qp p Qm qm m
25613: Qp qp p Qm qm p
12653: Qp qp p Qm p qm
7937: Qp qp Qm m m qm
1457: Qp qp Qm m qm m
24785: Qp qp Qm m qm p
11825: Qp qp Qm m p qm
377: Qp qp Qm qm m m
23705: Qp qp Qm qm m p
4265: Qp qp Qm qm p m
27593: Qp qp Qm qm p p
8585: Qp qp Qm p m qm
2105: Qp qp Qm p qm m
25433: Qp qp Qm p qm p
12473: Qp qp Qm p p qm
33935: Qp p m qm qp Qm
20975: Qp p m qm Qm qp
32855: Qp p m qp qm Qm
13415: Qp p m qp Qm qm
17735: Qp p m Qm qm qp
11255: Qp p m Qm qp qm
33755: Qp p qm m qp Qm
20795: Qp p qm m Qm qp
31595: Qp p qm qp m Qm
35483: Qp p qm qp p Qm
5675: Qp p qm qp Qm m
29003: Qp p qm qp Qm p
34403: Qp p qm p qp Qm
21443: Qp p qm p Qm qp
16475: Qp p qm Qm m qp
3515: Qp p qm Qm qp m
26843: Qp p qm Qm qp p
20363: Qp p qm Qm p qp
32495: Qp p qp m qm Qm
13055: Qp p qp m Qm qm
31415: Qp p qp qm m Qm
35303: Qp p qp qm p Qm
5495: Qp p qp qm Qm m
28823: Qp p qp qm Qm p
33143: Qp p qp p qm Qm
13703: Qp p qp p Qm qm
8735: Qp p qp Qm m qm
2255: Qp p qp Qm qm m
25583: Qp p qp Qm qm p
12623: Qp p qp Qm p qm
34043: Qp p p qm qp Qm
21083: Qp p p qm Qm qp
32963: Qp p p qp qm Qm
13523: Qp p p qp Qm qm
17843: Qp p p Qm qm qp
11363: Qp p p Qm qp qm
17015: Qp p Qm m qm qp
10535: Qp p Qm m qp qm
15935: Qp p Qm qm m qp
2975: Qp p Qm qm qp m
26303: Qp p Qm qm qp p
19823: Qp p Qm qm p qp
8375: Qp p Qm qp m qm
1895: Qp p Qm qp qm m
25223: Qp p Qm qp qm p
12263: Qp p Qm qp p qm
17663: Qp p Qm p qm qp
11183: Qp p Qm p qp qm
16877: Qp Qm m m qm qp
10397: Qp Qm m m qp qm
15797: Qp Qm m qm m qp
2837: Qp Qm m qm qp m
26165: Qp Qm m qm qp p
19685: Qp Qm m qm p qp
8237: Qp Qm m qp m qm
1757: Qp Qm m qp qm m
25085: Qp Qm m qp qm p
12125: Qp Qm m qp p qm
17525: Qp Qm m p qm qp
11045: Qp Qm m p qp qm
15617: Qp Qm qm m m qp
2657: Qp Qm qm m qp m
25985: Qp Qm qm m qp p
19505: Qp Qm qm m p qp
497: Qp Qm qm qp m m
23825: Qp Qm qm qp m p
4385: Qp Qm qm qp p m
27713: Qp Qm qm qp p p
16265: Qp Qm qm p m qp
3305: Qp Qm qm p qp m
26633: Qp Qm qm p qp p
20153: Qp Qm qm p p qp
7877: Qp Qm qp m m qm
1397: Qp Qm qp m qm m
24725: Qp Qm qp m qm p
11765: Qp Qm qp m p qm
317: Qp Qm qp qm m m
23645: Qp Qm qp qm m p
4205: Qp Qm qp qm p m
27533: Qp Qm qp qm p p
8525: Qp Qm qp p m qm
2045: Qp Qm qp p qm m
25373: Qp Qm qp p qm p
12413: Qp Qm qp p p qm
16985: Qp Qm p m qm qp
10505: Qp Qm p m qp qm
15905: Qp Qm p qm m qp
2945: Qp Qm p qm qp m
26273: Qp Qm p qm qp p
19793: Qp Qm p qm p qp
8345: Qp Qm p qp m qm
1865: Qp Qm p qp qm m
25193: Qp Qm p qp qm p
12233: Qp Qm p qp p qm
17633: Qp Qm p p qm qp
11153: Qp Qm p p qp qm
*/

