/* The 4q ng amplitudes use a base 6 code to enumerate them
   m,qm,qp,qbm,qbp,p corresponds to 0,1,2,3,4,5,6.
   ep and em are left out of the code.
   For example, (qp,qm,m,qbm,qbp,p) = 2 + 1*6 + 0*36 + 3*216 + 4*1296 + 5*7776
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

template <class T> complex<T> A4q2gzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A4q2g280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(3,0))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(3,0))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(3,0))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(3,0))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(4,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g1690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g1720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g1870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g1875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g1960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g1970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(3,1)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),3),2)*
 ep.spa(1,5)*ep.spab(5,ep.Sum(1,0),
4))/(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spab(1,ep.Sum(4,5),0)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
ep.spa(1,5)*ep.spb(5,4)))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 (-(ep.spa(0,2)*ep.spb(4,0))-
ep.spa(2,5)*ep.spb(5,4)))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),3),2))/
(pow(ep.spa(0,1),2)*ep.s(2,3,4)*
ep.spb(3,2)*ep.spb(4,0)+
 ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(1,5)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2)*ep.spb(4,2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A4q2g2355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A4q2g2535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g2650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(4,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g2955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g2985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(3,1)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3)*
 ep.spab(1,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g3835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g3950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(4,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g3985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(3,1)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spa(4,5),
2)*pow(ep.spab(5,ep.Sum(1,2),0),3)*
 pow(ep.spb(1,0),2)*ep.spa(3,5)*
 ep.spb(2,0))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(-(ep.s(3,4,5)*ep.spa(1,
   5))-ep.spa(0,1)*ep.spab(5,
  ep.Sum(1,2),0))*
 ((-ep.s(1,2)+ep.s(3,4,5))*
 ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),3),2)*
 ep.spa(1,5))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,0),3),2)*
 ep.spa(1,5)*ep.spab(5,ep.Sum(1,0),
4))/(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spab(1,ep.Sum(4,5),0)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(1,0),2)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0)))-(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(0,1),
 3),2)*ep.spab(5,ep.Sum(0,1),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),0),2)*
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

template <class T> complex<T> A4q2g4465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(5,ep.Sum(1,2),0),2)*
 pow(ep.spb(2,0),3)*ep.spa(1,2)*
 ep.spa(3,5))/(pow(ep.spab(5,ep.Sum(3,4),
 0),2)*ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*(ep.s(0,1,2)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(3,2),2)*ep.spb(4,2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3)*(-(ep.spa(0,1)*
  ep.spb(4,0))-ep.spa(1,5)*
 ep.spb(5,4)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A4q2g4525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g4825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g4940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2))/
 (ep.s(1,2,3)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g4945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
2)*ep.spab(5,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(1,2),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(3,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(0,1),2)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4)*
 (-(ep.spa(0,1)*ep.spb(4,0))-
ep.spa(1,5)*ep.spb(5,4)))
); }

template <class T> complex<T> A4q2g5065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(4,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1)*
ep.spb(4,0))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(3,1)*ep.spb(4,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,0))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(3,5),3)*pow(ep.spab(5,
 ep.Sum(1,2),0),2)*pow(ep.spb(2,0),
3)*ep.spab(5,ep.Sum(2,1),0))/
(pow(ep.spab(5,ep.Sum(3,4),0),2)*
 ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 (ep.s(3,4,5)*ep.spa(1,5)+
ep.spa(0,1)*ep.spab(5,ep.Sum(1,2),
  0))*((-ep.s(1,2)+ep.s(3,4,
   5))*ep.spab(5,ep.Sum(1,2),0)+
ep.s(3,4,5)*ep.spab(5,ep.Sum(3,
   4),0)))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(0,1),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(0,5),
4)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3)*
 ep.spab(1,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,1),2))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2))-
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,2),2)*ep.spab(5,
ep.Sum(2,3),4))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(5,0),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(1,2)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g5605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2)*ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(5,ep.Sum(2,3),4),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(2,3))/(pow(ep.spab(5,ep.Sum(0,1),
 4),2)*ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*(ep.s(2,3,4)*ep.spa(3,5)-
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4))*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spb(5,0)*
 ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 ep.spa(1,3)*ep.spb(4,0))/
(ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g5905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g5935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),3))/
(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(4,3),2))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),3))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),4),2))/
 (ep.s(0,4,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,3),
2)*ep.spa(1,5)*ep.spa(2,3)*
 ep.spab(5,ep.Sum(2,3),4))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.spa(3,4)-(ep.s(2,3,4)*
  ep.spa(3,5))/ep.spab(5,ep.Sum(2,3),
  4)))+(complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),
 1),2)*ep.spab(5,ep.Sum(3,4),0))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(4,3),1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(2,ep.Sum(4,5),0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spb(4,0))/(ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g6855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g6970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g6975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5)*ep.spab(5,ep.Sum(3,4),
0))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
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
 4),2)*ep.spb(4,0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(3,ep.Sum(5,0),4)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spa(2,3),
2)*pow(ep.spab(5,ep.Sum(2,3),4),3)*
 pow(ep.spb(4,2),3))/
(pow(ep.spab(5,ep.Sum(0,1),4),2)*
 ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 (ep.s(2,3,4)*ep.spab(5,ep.Sum(0,
   1),4)+(-ep.s(2,3)+
  ep.s(2,3,4))*ep.spab(5,
  ep.Sum(2,3),4))*
 (-(ep.s(2,3,4)*ep.spa(3,5))+
ep.spa(3,4)*ep.spab(5,ep.Sum(2,3),
  4)))+(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(5,
ep.Sum(3,4),0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 pow(ep.spb(4,0),3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3)*
 ep.spab(3,ep.Sum(5,0),4))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*pow(ep.spb(4,0),3))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5)*ep.spab(5,ep.Sum(3,4),
0))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(5,4),0)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
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
  4)))+(complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),
 4),2)*ep.spb(4,0))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(5,0),4),2)*
 ep.spab(3,ep.Sum(5,0),4)*
 ep.spb(4,0))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g7840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g7930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(5,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g7945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(5,3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(2,0)*ep.spb(5,3))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(5,3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(1,0),2))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1))-
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 ep.spab(4,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(5,0),3)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(2,0),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spab(4,ep.Sum(1,2),3),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(1,2))/(pow(ep.spab(4,ep.Sum(5,0),
 3),2)*ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*(ep.s(1,2,3)*ep.spa(2,4)-
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,5),3),2))/
 (ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),2))/
 (ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,2),2)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spab(4,ep.Sum(1,2),
3))/(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (ep.spa(2,3)-(ep.s(1,2,3)*
  ep.spa(2,4))/ep.spab(4,ep.Sum(1,2),
  3)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),
 0),2))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spab(4,ep.Sum(2,3),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g8835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g8980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g8990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g9010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(4,
 ep.Sum(1,2),3),3)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4)*ep.spab(4,ep.Sum(2,3),
5))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(2,ep.Sum(4,5),3)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g9015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g9220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g9240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g9245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g9250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g9960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g9965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g9970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g10110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g10390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g10480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g10510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g10515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g10530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g10535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g10545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g10570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g10575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g11190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g11655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g11765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g11770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g11860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g11870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g12120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g12125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g12245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g12270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g12770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g12840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g12845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g12950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g12970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g12975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g12980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g12985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(4,2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2)*ep.spb(2,0))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 ep.spb(4,3))-
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

template <class T> complex<T> A4q2g13030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(0,5),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),2)*ep.spab(3,
ep.Sum(0,5),4))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
2))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,0),2))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(5,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(1,0),
2)*ep.spa(0,5)*ep.spa(2,4))/
(ep.s(2,3,4)*ep.spa(3,4)*
 (-(ep.s(2,3,4)*ep.spa(0,2))-
ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),4),2)*
 ep.spab(2,ep.Sum(1,0),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),4),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(0,ep.Sum(2,3),1)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2)*ep.spab(2,ep.Sum(0,1),
3))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
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
  1))*ep.spab(4,ep.Sum(0,5),1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(5,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
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
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g13350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
1))-(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),
 4),2)*ep.spa(0,2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spab(2,ep.Sum(5,0),1)*
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

template <class T> complex<T> A4q2g13355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 ep.spa(2,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spa(3,4),2)*pow(ep.spab(1,
 ep.Sum(3,4),2),3)*pow(ep.spb(4,2),
3))/(pow(ep.spab(1,ep.Sum(5,0),2),2)*
 ep.s(2,3,4)*ep.spa(0,5)*
 (-(ep.s(2,3,4)*ep.spa(1,3))+
ep.spa(2,3)*ep.spab(1,ep.Sum(3,4),
  2))*((-ep.s(3,4)+ep.s(2,3,
   4))*ep.spab(1,ep.Sum(3,4),2)+
ep.s(2,3,4)*ep.spab(1,ep.Sum(5,
   0),2))*ep.spab(5,ep.Sum(3,4),
2))+(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 pow(ep.spb(2,0),3))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3)*
 ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(4,0),2))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
ep.spa(1,2)*ep.spb(5,1)))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 (-(ep.spa(0,3)*ep.spb(5,0))-
ep.spa(1,3)*ep.spb(5,1)))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0)*(-(ep.spa(0,2)*
  ep.spb(5,0))-ep.spa(1,2)*
 ep.spb(5,1)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2))/
(ep.s(3,4,5)*ep.spa(0,2)*
ep.spa(1,2)*ep.spb(4,3)*
ep.spb(5,0)+pow(ep.spa(1,2),2)*
ep.s(3,4,5)*ep.spb(4,3)*
ep.spb(5,1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 ep.spa(1,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(3,5),2))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g13825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g13865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g13895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4)*ep.spab(2,ep.Sum(3,4),
1))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(4,ep.Sum(1,2),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(2,
 ep.Sum(4,5),3),3)*pow(ep.spb(4,3),
2)*ep.spa(0,2)*ep.spb(5,3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A4q2g14055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,3),2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(5,ep.Sum(1,0),2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),2))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g14100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.s(1,2,3)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0))+
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

template <class T> complex<T> A4q2g14105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(1,5)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(0,3)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(0,3)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(1,3),
3)*pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(4,0),3)*ep.spab(3,
ep.Sum(0,5),4))/
(pow(ep.spab(3,ep.Sum(1,2),4),2)*
 ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 (ep.s(1,2,3)*ep.spab(3,ep.Sum(1,
   2),4)+(-ep.s(0,5)+
  ep.s(1,2,3))*ep.spab(3,
  ep.Sum(5,0),4))*
 (-(ep.s(1,2,3)*ep.spa(3,5))+
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3)*
 ep.spab(5,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g14415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g14420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g14425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g14500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g14895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g14930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g14980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spab(5,ep.Sum(3,2),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g15135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),3)*ep.spa(1,3)*
 ep.spa(4,5))/(ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spa(1,5)+ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),0))*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g15180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(0,
ep.Sum(2,3),1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(0,
ep.Sum(2,3),1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
1)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))
); }

template <class T> complex<T> A4q2g15185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(3,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(0,2)*ep.spa(3,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g15210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g15240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g15270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(3,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(2,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g15610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(5,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g15935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g15970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g15975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(5,3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(2,0)*ep.spb(5,3))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16290_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))-(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(5,3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,0),2)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5)))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 2),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),2),2)*
 ep.spab(4,ep.Sum(5,0),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3)))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(5,4),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(3,2)*
 (-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(3,4),5),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),
2)*ep.spab(4,ep.Sum(1,0),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*(-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spa(1,2),2)*pow(ep.spab(4,
 ep.Sum(1,2),3),3)*pow(ep.spb(3,1),
3))/(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g16785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g16795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g16960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g16995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g17025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g17855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g17895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g17915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g17925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g17950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g17955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g18060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g18065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g18070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g18075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g18105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g18130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g18135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g18165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g18285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g18975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g18990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g18995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g19005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g19245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g19355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g19365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g19425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g19450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(4,2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),3)*ep.spa(0,1)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(4,5)*(-(ep.s(3,4,5)*
  ep.spa(1,3))-ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),2))*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 ep.spb(4,3))-
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

template <class T> complex<T> A4q2g19545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,1),
2)*ep.spa(0,1)*ep.spa(3,5))/
(ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))-
ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
  2))*ep.spab(5,ep.Sum(0,1),2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),5),2)*
 ep.spab(3,ep.Sum(2,1),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
2))+(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),
 2),2)*ep.spb(4,2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(1,ep.Sum(3,4),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3)*ep.spab(3,ep.Sum(1,2),
4))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g19865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g19870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g19875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,3),2))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,3),2)*
 ep.spb(5,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,0),2))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(5,4),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),2)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spa(2,4),
3)*pow(ep.spab(2,ep.Sum(5,0),1),3)*
 pow(ep.spb(5,1),3))/
(pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.s(0,1,5)*ep.spa(3,4)*
 (ep.s(0,1,5)*ep.spab(2,ep.Sum(3,
   4),1)+(-ep.s(0,5)+
  ep.s(0,1,5))*ep.spab(2,
  ep.Sum(5,0),1))*
 (-(ep.s(0,1,5)*ep.spa(0,2))+
ep.spa(0,1)*ep.spab(2,ep.Sum(5,0),
  1))*ep.spab(4,ep.Sum(5,0),1))+
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3)*
 ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(2,ep.Sum(0,1),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spb(5,3),2)*
 ep.spb(5,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g19985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
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
  0))*ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g20010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
1))-(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(2,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spab(3,ep.Sum(0,1),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(3,ep.Sum(2,1),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
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
2))-(complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),
 5),2)*ep.spa(1,3)*ep.spab(1,
ep.Sum(2,3),0))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),5),2)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,2),2))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))-
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

template <class T> complex<T> A4q2g20075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(1,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,4),2))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3)*
 ep.spab(4,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(3,2),2)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(-(ep.s(0,1,5)*
  ep.spa(1,3))+ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 (ep.spa(1,4)*ep.spb(1,0)+
ep.spa(2,4)*ep.spb(2,0)))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2))/
(-(ep.s(0,4,5)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spb(1,0)*
 ep.spb(5,4))-pow(ep.spa(2,3),2)*
ep.s(0,4,5)*ep.spb(2,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2)*ep.spb(3,1))/(ep.s(0,4,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0))-
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

template <class T> complex<T> A4q2g20615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(1,5)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g20640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g20745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g20805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g20955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g20975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g20985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g21060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g21065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g21090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g21095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g21135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g21150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g21155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g21165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g21385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g21835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g21870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g21875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g21885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g21925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g21955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g22015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,4),1),2)*
 ep.spa(3,5)*ep.spab(3,ep.Sum(5,4),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(5,4),
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
 (ep.s(1,2,3)*ep.spa(3,5)-
ep.spa(4,5)*ep.spab(3,ep.Sum(5,0),
  4)))+(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),
 4),2)*ep.spb(4,2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spab(5,ep.Sum(2,3),4)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g22160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spab(2,ep.Sum(1,0),5))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(0,4),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),1)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g22700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g22740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 ep.spa(3,5)*ep.spb(2,0))/
(ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spab(1,ep.Sum(3,4),2),2)*
 pow(ep.spb(4,2),3)*ep.spa(1,5)*
 ep.spa(3,4))/(pow(ep.spab(1,ep.Sum(5,0),
 2),2)*ep.s(0,1,5)*
 ep.spa(0,5)*(ep.s(0,1,5)*
 ep.spa(1,3)-ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,4),2))*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,4),2)*ep.spb(4,0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g22770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),2),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g22775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(3,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,2)*ep.spa(3,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g22800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g22830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(3,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g22935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g22955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g22965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g23265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(2,0),2)*
 ep.spb(5,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g23750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(5,3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(4,5),3),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,0),2)*
 ep.spb(2,0)*ep.spb(5,3))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(3,1),2)*
 ep.spb(5,3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1))+
 (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,1),2)*
 ep.spab(4,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(5,0),3)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(1,2),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(2,0),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(4,5),3),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),3),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spa(1,2)*
 ep.spab(4,ep.Sum(1,2),3))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (ep.spa(2,3)-(ep.s(1,2,3)*
  ep.spa(2,4))/ep.spab(4,ep.Sum(1,2),
  3)))-(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),
 0),2))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spab(4,ep.Sum(2,3),5))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g23965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g23990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(4,
 ep.Sum(1,2),3),3)*pow(ep.spb(3,2),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.s(1,2,3)*ep.spa(2,4))+
ep.spa(2,3)*ep.spab(4,ep.Sum(1,2),
  3))*((-ep.s(1,2)+ep.s(1,2,
   3))*ep.spab(4,ep.Sum(1,2),3)+
ep.s(1,2,3)*ep.spab(4,ep.Sum(5,
   0),3)))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),0),2)*
 ep.spa(2,4)*ep.spab(4,ep.Sum(2,3),
5))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(4,5),3),2)*
 ep.spab(2,ep.Sum(4,5),3)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g24505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g24700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g24780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g24785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g24800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g24820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g24830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g25500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g25505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g25565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g25590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g25780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g25865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g25910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g25930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g25985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g25995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,0),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(4,2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(1,0),
2)*ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))-
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

template <class T> complex<T> A4q2g26060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(1,0),2))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(3,4),2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(0,5),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(3,ep.Sum(1,2),4)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(4,0),2)*ep.spab(3,
ep.Sum(0,5),4))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(0,5),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,5),3)*pow(ep.spab(3,
 ep.Sum(0,1),2),3)*pow(ep.spb(2,0),
3))/(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(0,1,2)*ep.spa(4,5)*
 (-(ep.s(0,1,2)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(0,1,
   2))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(0,1,2)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(0,1),
2))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3)*
 ep.spab(1,ep.Sum(3,4),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(3,ep.Sum(1,2),4))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,0),2))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,1),2)*
 ep.spb(5,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(1,0),2)*ep.spa(0,5)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(3,4)*(-(ep.s(2,3,4)*
  ep.spa(0,2))-ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),1))*
 ep.spab(4,ep.Sum(5,0),1))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
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
 ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(1,0),4),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),1),2)*
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
  1))*ep.spab(4,ep.Sum(0,5),1))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),4),2)*
 ep.spa(0,2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(4,2),2)*
 ep.spb(5,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spa(4,5),
2)*pow(ep.spab(1,ep.Sum(4,5),0),3)*
 pow(ep.spb(4,0),3))/
(pow(ep.spab(1,ep.Sum(2,3),0),2)*
 ep.s(0,4,5)*ep.spa(2,3)*
 (ep.s(0,4,5)*ep.spab(1,ep.Sum(2,
   3),0)+(-ep.s(4,5)+
  ep.s(0,4,5))*ep.spab(1,
  ep.Sum(4,5),0))*
 (ep.s(0,4,5)*ep.spa(1,5)-
ep.spa(0,5)*ep.spab(1,ep.Sum(4,5),
  0))*ep.spab(3,ep.Sum(4,5),0))+
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3)*
 ep.spab(5,ep.Sum(1,2),0))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(5,0),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,5),3)*
 pow(ep.spb(4,2),2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g26670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spa(2,3),2)*pow(ep.spab(0,
 ep.Sum(2,3),1),3)*pow(ep.spb(2,1),
2)*ep.spa(0,4)*ep.spb(3,1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (ep.s(0,4,5)*ep.spa(0,2)-
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(2,3),
1))+(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),
 4),2)*ep.spa(0,2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spab(2,ep.Sum(5,0),1)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,1),4),2)*
 ep.spa(0,2)*ep.spab(0,ep.Sum(2,1),
5))/(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 ep.spa(2,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 pow(ep.spb(2,0),3))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),3)*
 ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(4,0),2))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(2,1),2)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.s(1,2,3)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,0),1),2)*
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
ep.spb(5,1))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),4),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(3,5),2)*
 ep.spa(1,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(3,5),2))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g26825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g26975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g26995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4)*ep.spab(2,ep.Sum(3,4),
1))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(4,ep.Sum(1,2),3)*
 ep.spb(3,1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spa(4,5),2)*pow(ep.spab(2,
 ep.Sum(4,5),3),3)*pow(ep.spb(4,3),
2)*ep.spa(0,2)*ep.spb(5,3))/
(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(0,1,2)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(0,1,2))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(0,1,2)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))
); }

template <class T> complex<T> A4q2g27025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,3),
2))/(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(1,2),0))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
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
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(4,0),2))/(ep.s(1,2,3)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(4,5),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g27120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*
 pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spa(0,4)*
 ep.spa(2,3))/(pow(ep.spab(0,ep.Sum(4,5),
 1),2)*ep.s(1,2,3)*
 ep.spa(4,5)*(-(ep.s(1,2,3)*
  ep.spa(0,2))+ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,3),1))*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))-
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*(-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1))*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g27125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(1,5)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(0,3)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,4))/
(ep.spa(0,2)*ep.spa(1,2)*
 ep.spa(1,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 ep.spa(0,3)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,3)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g27230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27290_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(3,4),2),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g27320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g27725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(3,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g27750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g27830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g27845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g27860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,2)*
 ep.spa(1,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(4,5),2)*
 ep.spa(0,3)*ep.spa(2,5))/
(ep.spa(0,2)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(1,5)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(4,5),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(1,0),2))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(4,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(2,3),4),2))/
 (ep.s(2,3,4)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g28580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(4,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g28805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*pow(ep.spb(2,0),
3))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,5),3)*pow(ep.spb(4,2),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(1,0),
2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g28830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g28925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g28945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g28985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g29605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g29640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2))/
 (ep.s(2,3,4)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g29645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g29660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g29665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g29785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g29870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g29890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,0),
2))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g29895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g29900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 pow(ep.spb(2,0),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(3,5),3)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(4,5),2))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 pow(ep.spb(4,2),3))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(1,5),2)*pow(ep.spb(4,2),3)*
 ep.spab(5,ep.Sum(2,3),4))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g29905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*pow(ep.spb(2,1),
2))/(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(4,2),3))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g29930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g29960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g29965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(4,3),2)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(-(ep.s(0,1,2)*ep.spa(2,
   4))+ep.spa(3,4)*ep.spab(2,
  ep.Sum(4,5),3)))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spab(5,ep.Sum(3,2),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),3),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),0),2))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g30265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g30300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g30330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(5,ep.Sum(1,2),0),2)*
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

template <class T> complex<T> A4q2g30360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,5),2)*pow(ep.spb(3,1),
2)*ep.spab(0,ep.Sum(2,3),1))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,5),2)*
 pow(ep.spb(3,1),2)*ep.spab(0,
ep.Sum(2,3),1))/(ep.s(1,2,3)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
1)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(4,3),2))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))-
 (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(3,4),1),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0)*
 (-(ep.spa(0,2)*ep.spb(5,0))-
ep.spa(1,2)*ep.spb(5,1)))
); }

template <class T> complex<T> A4q2g30365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(1,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(3,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(2,5),2)*
 ep.spa(0,2)*ep.spa(3,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(5,ep.Sum(0,1),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g30390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(3,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g30725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g30950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,5),2)*
 ep.spa(2,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g30985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,5),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g31010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,5),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g31040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,5),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g31045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(5,2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(5,3))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(3,4),5),
 2))/(ep.s(3,4,5)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(2,ep.Sum(3,4),5),
 2))/(ep.s(0,1,2)*ep.spa(0,1)*
ep.spa(1,2)*ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(2,0)*ep.spb(5,3))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,3)*ep.spb(5,0)*
 ep.spb(5,4))+(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(5,3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,0),2)*ep.spa(0,1)*
 ep.spa(2,4))/(pow(ep.spab(4,ep.Sum(2,3),
 5),2)*ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*(ep.s(0,1,5)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5)))+(complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),
 2),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),2),2)*
 ep.spab(4,ep.Sum(5,0),3))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(5,4),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spb(4,3)*
 ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spa(0,2)*ep.spb(5,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(1,0),5))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(1,0),
5)*ep.spab(4,ep.Sum(2,3),5)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(0,1),5))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(2,1),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(5,0),
1)*(-(ep.spa(0,4)*ep.spb(4,3))-
ep.spa(0,5)*ep.spb(5,3)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),5),2)*
 ep.spab(1,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.spb(5,3))/(ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*(-(ep.spa(0,4)*
  ep.spb(4,3))-ep.spa(0,5)*
 ep.spb(5,3))*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2)*ep.spb(5,3))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0)*
ep.spb(5,3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(5,3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g31835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g31945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g31975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(4,3),
5)*ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,3),3))/
(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
   0),3)))-(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(4,
ep.Sum(2,3),5))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(4,3),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(5,0))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3)*
 ep.spab(2,ep.Sum(4,5),3))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g32275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g32305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,1),2)*
 ep.spb(4,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),2)*ep.spb(4,2)*
ep.spb(5,1))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(4,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,2),2)*
 ep.spb(4,2)*ep.spb(5,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3)*ep.spa(0,1)*ep.spa(3,5))/
(ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))-
ep.spa(1,2)*ep.spab(3,ep.Sum(4,5),
  2))*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spb(3,2)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spa(1,5)*ep.spb(4,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2)*
 ep.spb(4,0))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,1),2)*ep.spa(0,1)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(4,5)*(-(ep.s(3,4,5)*
  ep.spa(1,3))-ep.spa(1,2)*
 ep.spab(3,ep.Sum(4,5),2))*
 ep.spab(5,ep.Sum(0,1),2))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(2,3),4))/
(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,1)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(2,1),5),2)*
 ep.spab(3,ep.Sum(2,1),4))/
(ep.s(0,4,5)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spa(3,4),2)*pow(ep.spab(3,
 ep.Sum(0,1),2),2)*pow(ep.spb(2,1),
2)*ep.spa(3,5)*ep.spab(3,ep.Sum(1,0),
2)*ep.spb(2,0))/
(pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.s(3,4,5)*ep.spa(4,5)*
 (-(ep.s(3,4,5)*ep.spa(1,3))+
ep.spa(1,2)*ep.spab(3,ep.Sum(0,1),
  2))*((-ep.s(0,1)+ep.s(3,4,
   5))*ep.spab(3,ep.Sum(0,1),2)+
ep.s(3,4,5)*ep.spab(3,ep.Sum(4,
   5),2))*ep.spab(5,ep.Sum(1,0),
2))-(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),
 2),2)*ep.spb(4,2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(3,4),2),2)*
 ep.spab(1,ep.Sum(3,4),2)*
 ep.spb(4,2))/(ep.s(2,3,4)*
 ep.spa(0,1)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(1,2),5),2)*
 ep.spa(1,3)*ep.spab(3,ep.Sum(1,2),
4))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g32830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g32835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),
3))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(2,1),
3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,3),2)*ep.spb(5,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g32945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g32975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(2,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(5,2),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),
3))/(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,3),2))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,3),2)*
 ep.spb(5,1))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),
2))/(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(5,0),3),2)*
 ep.spab(5,ep.Sum(1,2),3))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(1,2),
3)*ep.spb(2,1)*ep.spb(3,2))+
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
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,3),2)*ep.spab(2,
ep.Sum(4,5),3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(2,ep.Sum(0,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
  1))*ep.spab(4,ep.Sum(5,0),1))-
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spb(2,1)*ep.spb(3,2))+
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
 ep.spb(4,3))-(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(2,1),3)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3)*
 ep.spb(5,0))+(complex<T>(0,1)*pow(ep.spb(5,3),2)*
 ep.spb(5,2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(5,0)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.s(0,1,2)*
 ep.spa(3,4)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(5,ep.Sum(1,2),0)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5)*ep.spab(1,ep.Sum(5,0),
2))/(ep.s(0,1,5)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
4)*ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(1,2),2)*
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
  0))*ep.spab(3,ep.Sum(4,5),0))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(5,0),3),2)*
 ep.spa(1,5))/(ep.s(0,1,5)*
 ep.spa(0,1)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g33330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spa(2,3),
2)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 pow(ep.spb(3,1),3)*ep.spab(0,
ep.Sum(3,2),1))/
(pow(ep.spab(0,ep.Sum(4,5),1),2)*
 ep.s(0,4,5)*ep.spa(4,5)*
 (-(ep.s(0,4,5)*ep.spa(0,2))+
ep.spa(1,2)*ep.spab(0,ep.Sum(2,3),
  1))*((-ep.s(2,3)+ep.s(0,4,
   5))*ep.spab(0,ep.Sum(2,3),1)+
ep.s(0,4,5)*ep.spab(0,ep.Sum(4,
   5),1))*ep.spab(4,ep.Sum(3,2),
1))+(complex<T>(0,1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,3),2))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3))/(ep.s(0,1,5)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(1,0),
5)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(1,2),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(1,2),
3)*ep.spab(2,ep.Sum(0,1),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(2,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.s(3,4,5)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(2,1),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
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
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,2),5),2)*
 ep.spa(1,3))/(ep.s(0,4,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2))/
(ep.s(0,1,2)*ep.spa(1,2)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(0,1),5))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(2,3)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),2)*
 ep.spab(0,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(1,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,4),2))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3)*
 pow(ep.spa(4,5),2)*pow(ep.spab(2,
 ep.Sum(4,5),3),3)*pow(ep.spb(5,3),
3))/(pow(ep.spab(2,ep.Sum(0,1),3),2)*
 ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(4,5),3)*
 (ep.s(3,4,5)*ep.spab(2,ep.Sum(0,
   1),3)+(-ep.s(4,5)+
  ep.s(3,4,5))*ep.spab(2,
  ep.Sum(4,5),3))*
 (-(ep.s(3,4,5)*ep.spa(2,4))+
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3)))+(complex<T>(0,1)*pow(ep.spa(2,4),3)*
 pow(ep.spb(5,1),2)*ep.spab(2,
ep.Sum(3,4),1))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,4),2)*pow(ep.spb(3,1),3)*
 ep.spab(4,ep.Sum(1,2),3))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g33520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,5)*
ep.spa(2,4))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),0),
3))/(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(1,2),
0)*ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),2))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),2),2)*
 (ep.spa(1,4)*ep.spb(1,0)+
ep.spa(2,4)*ep.spb(2,0)))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 (ep.spa(1,3)*ep.spb(1,0)+
ep.spa(2,3)*ep.spb(2,0))*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2)*
 ep.spab(1,ep.Sum(2,3),0))/
(ep.s(1,2,3)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),5),2))/
(-(ep.s(0,4,5)*ep.spa(1,3)*
 ep.spa(2,3)*ep.spb(1,0)*
 ep.spb(5,4))-pow(ep.spa(2,3),2)*
ep.s(0,4,5)*ep.spb(2,0)*
ep.spb(5,4))
); }

template <class T> complex<T> A4q2g33630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,2),2)*ep.spb(3,1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spb(1,0)*
 ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 ep.spa(2,4)*ep.spb(5,1))/
(ep.spa(2,3)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
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

template <class T> complex<T> A4q2g33635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(1,5)*ep.spa(2,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(2,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4))+(complex<T>(0,1)*pow(ep.spa(0,4),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g33705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g33735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g33810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g33815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g33835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g33885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g33915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g33945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g34525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g34555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g34785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g34890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g34895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g34915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g34965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g34975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g35000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g35095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g35425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,3),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g35495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g35605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g35615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g35635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*
 pow(ep.spb(3,1),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spab(4,ep.Sum(3,2),
1)*ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,1),2))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g36265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(2,3),4),2))/
 (ep.s(0,1,5)*ep.spa(0,1)*
ep.spa(0,5)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g36355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),3))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),2))/
(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(3,2),
4)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(5,4),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),
2))/(ep.s(0,1,2)*ep.spa(0,1)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),3))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g36690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(2,ep.Sum(3,4),
5)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),3))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g36695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g36715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g37375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(2,ep.Sum(0,1),5),2))/
 (ep.s(0,1,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g37415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g37435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q2g37555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q2g37630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spa(3,5))/(ep.s(3,4,5)*
 ep.spa(3,4)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(5,4),1),2)*
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
  4)))-(complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),
 4),2)*ep.spb(4,2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spab(5,ep.Sum(2,3),4)*
 ep.spb(4,2))/(ep.s(0,1,5)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g37665_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g37675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g37705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3)*
 pow(ep.spb(5,4),2))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(1,0),5))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(3,4),1))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spab(5,ep.Sum(0,4),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(2,3),
1)*ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spb(3,1))/(ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),2)*
 ep.spab(2,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(3,4),1)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(5,0))
); }

template <class T> complex<T> A4q2g37840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g37860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(4,5),2),
 2))/(ep.s(3,4,5)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g37890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(3,ep.Sum(1,2),0),
 2))/(ep.s(0,1,2)*ep.spa(3,4)*
ep.spa(4,5)*ep.spb(1,0)*ep.spb(2,1))
); }

template <class T> complex<T> A4q2g37920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g37925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spb(1,0)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
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

template <class T> complex<T> A4q2g37950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,2)*
 ep.spa(3,4)*ep.spab(0,ep.Sum(3,4),
5))/(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.spa(4,5)+(ep.s(3,4,5)*
  ep.spa(0,4))/ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(0,4,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(5,4),2),2))/
(ep.s(1,2,3)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(3,ep.Sum(5,0),1))/
(ep.s(0,1,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g37955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(1,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(3,5))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,2)*ep.spa(3,5))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g37995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(3,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2)*
ep.spa(3,5))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(3,5))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(3,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,5))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g38780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g38785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,5))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g38815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
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
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(0,ep.Sum(3,4),5)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))-
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3)*
 ep.spab(0,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39580_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
 ep.spa(0,4))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(0,5),2),2)*
 ep.spa(0,4)*ep.spab(4,ep.Sum(0,5),
3))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
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
   3),5)))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),5),2)*
 ep.spab(0,ep.Sum(3,4),5)*
 ep.spb(5,3))/(ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(5,4),
3)*ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g39805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spa(2,4),
3)*pow(ep.spab(4,ep.Sum(0,1),5),2)*
 pow(ep.spb(5,1),3)*ep.spab(4,
ep.Sum(1,0),5))/
(pow(ep.spab(4,ep.Sum(2,3),5),2)*
 ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(2,ep.Sum(1,0),5)*
 (ep.s(2,3,4)*ep.spa(0,4)-
ep.spa(0,5)*ep.spab(4,ep.Sum(0,1),
  5))*((-ep.s(0,1)+ep.s(2,3,
   4))*ep.spab(4,ep.Sum(0,1),5)+
ep.s(2,3,4)*ep.spab(4,ep.Sum(2,
   3),5)))+(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2)*ep.spab(4,
ep.Sum(5,0),3))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(4,5),3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 pow(ep.spb(5,3),3))/(ep.s(3,4,5)*
 ep.spa(1,2)*ep.spab(0,ep.Sum(5,4),
3)*ep.spb(4,3)*ep.spb(5,4))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3)*
 ep.spab(0,ep.Sum(3,4),5))/
(ep.s(3,4,5)*ep.spa(0,1)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spab(2,ep.Sum(3,4),5)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g39865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g39895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,2),
2))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(5,3),3))/
(ep.s(3,4,5)*ep.spa(1,2)*
 ep.spab(0,ep.Sum(5,4),3)*
 ep.spb(4,3)*ep.spb(5,4))
); }

template <class T> complex<T> A4q2g40240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,0),2))/(ep.s(1,2,3)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(3,2),
4)*ep.spab(3,ep.Sum(1,2),0))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
 ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),2)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0))+(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(1,5),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g40335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g40540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g40565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40630_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g40635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g40745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40775_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40810_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g40850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g40920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g40925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g40970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g40985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(4,ep.Sum(1,2),3),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
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

template <class T> complex<T> A4q2g41055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*
 pow(ep.spb(5,0),2)*ep.spa(1,3)*
 ep.spa(4,5))/(ep.s(1,2,3)*
 ep.spa(2,3)*(ep.s(1,2,3)*
 ep.spa(1,5)+ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),0))*
 ep.spab(3,ep.Sum(4,5),0))-
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(3,ep.Sum(1,2),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),3),2)*
 ep.spab(1,ep.Sum(0,5),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(0,5),4)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),3),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g41100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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
 ep.spb(2,1))-(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))-
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3)*
 ep.spab(4,ep.Sum(0,1),5))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g41105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(0,4)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(2,4))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(-(ep.s(0,4,5)*ep.spa(3,
   5))+ep.spa(4,5)*ep.spab(3,
  ep.Sum(5,0),4)))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spab(3,ep.Sum(4,5),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(4,3),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g41600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g41800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g41855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g41895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g41925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g42015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spb(5,1),3)*ep.spa(0,5)*
 ep.spa(2,4))/(ep.s(2,3,4)*
 ep.spa(3,4)*(-(ep.s(2,3,4)*
  ep.spa(0,2))-ep.spa(0,1)*
 ep.spab(2,ep.Sum(3,4),1))*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 ep.spb(3,2))-
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

template <class T> complex<T> A4q2g42140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*
 pow(ep.spb(5,4),2))/(ep.s(0,4,5)*
 ep.spa(2,3)*ep.spab(1,ep.Sum(2,3),
4)*ep.spab(3,ep.Sum(1,2),0))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),2),2)*
 ep.spb(2,0))/(ep.spa(4,5)*
 ep.spab(3,ep.Sum(2,1),0)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0)*ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spa(0,1),2)*pow(ep.spb(4,2),2)*
 ep.spab(1,ep.Sum(4,3),2))/
(ep.s(0,1,5)*ep.spa(0,5)*
 ep.spab(1,ep.Sum(5,0),4)*
 ep.spab(5,ep.Sum(4,3),2)*
 ep.spb(3,2))-(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),2)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
2)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g42210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4)*ep.spab(0,ep.Sum(4,5),
1))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))-
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

template <class T> complex<T> A4q2g42215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,4)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g42355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,4))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g42715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g42830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,0),
2))/(ep.s(1,2,3)*ep.spa(1,2)*
 ep.spab(1,ep.Sum(3,2),4)*
 ep.spab(3,ep.Sum(1,2),0))-
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(2,0),2)*
 ep.spab(3,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),2)*
 ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(1,0))-(complex<T>(0,1)*pow(ep.spa(3,4),2)*
 pow(ep.spb(2,0),2)*ep.spab(3,
ep.Sum(0,1),2))/(ep.s(0,1,2)*
 ep.spa(4,5)*ep.spab(3,ep.Sum(4,5),
0)*ep.spab(5,ep.Sum(0,1),2)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spab(0,ep.Sum(1,5),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(3,4),
2)*ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(3,4),2),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g42865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g42925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g43060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43080_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g43145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(3,2),1)*
 ep.spb(2,1)*ep.spb(3,2))-
 (complex<T>(0,1)*pow(ep.spa(1,2),2)*
 pow(ep.spab(2,ep.Sum(4,5),3),2)*
 pow(ep.spb(5,3),3)*ep.spa(0,2)*
 ep.spa(4,5))/(pow(ep.spab(2,ep.Sum(0,1),
 3),2)*ep.s(0,1,2)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(4,5),
3)*(ep.s(0,1,2)*ep.spa(2,4)-
ep.spa(3,4)*ep.spab(2,ep.Sum(4,5),
  3))*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spa(2,4),3)*pow(ep.spb(5,0),2)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(3,4),
5)*ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(1,0))
); }

template <class T> complex<T> A4q2g43225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43290_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*pow(ep.spb(5,0),
2)*ep.spa(1,3)*ep.spa(4,5))/
(ep.s(1,2,3)*ep.spa(2,3)*
 (ep.s(1,2,3)*ep.spa(1,5)+
ep.spa(0,5)*ep.spab(1,ep.Sum(2,3),
  0))*ep.spab(3,ep.Sum(4,5),0))+
 (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),0),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(0,1,2)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
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
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(1,ep.Sum(0,5),3),2))/
(ep.s(2,3,4)*ep.spa(0,5)*
 ep.spab(5,ep.Sum(1,0),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g43320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3)*
 pow(ep.spa(3,4),2)*pow(ep.spab(0,
 ep.Sum(3,4),5),3)*pow(ep.spb(5,3),
3))/(pow(ep.spab(0,ep.Sum(1,2),5),2)*
 ep.s(3,4,5)*ep.spa(1,2)*
 (ep.s(3,4,5)*ep.spab(0,ep.Sum(1,
   2),5)+(-ep.s(3,4)+
  ep.s(3,4,5))*ep.spab(0,
  ep.Sum(3,4),5))*
 (ep.s(3,4,5)*ep.spa(0,4)+
ep.spa(4,5)*ep.spab(0,ep.Sum(3,4),
  5))*ep.spab(2,ep.Sum(3,4),5))+
 (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),2)*
 ep.spab(0,ep.Sum(4,5),1))/
(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))+(complex<T>(0,1)*pow(ep.spa(0,4),3)*
 pow(ep.spb(3,1),2))/(ep.s(1,2,3)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 pow(ep.spb(5,1),3))/(ep.s(2,3,4)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3)*
 ep.spab(4,ep.Sum(0,1),5))/
(ep.s(2,3,4)*ep.spa(3,4)*
 ep.spab(2,ep.Sum(0,1),5)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g43325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(2,4),2)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(0,4)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(2,3),1),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g43490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(2,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(3,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g43505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(4,ep.Sum(1,2),3),2))/
 (ep.s(1,2,3)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g43645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g43680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,4),3)*pow(ep.spb(3,1),
2))/(ep.s(1,2,3)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(4,ep.Sum(0,5),
1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*pow(ep.spb(5,1),3))/
(ep.s(2,3,4)*ep.spa(2,3)*
 ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))
); }

template <class T> complex<T> A4q2g43685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g43910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,4),2)*
 ep.spa(2,4))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g43970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44000_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44110_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*
 pow(ep.spab(3,ep.Sum(5,0),4),2)*
 pow(ep.spb(5,4),2)*ep.spa(0,5)*
 ep.spa(1,3))/(pow(ep.spab(3,ep.Sum(1,2),
 4),2)*ep.s(0,4,5)*
 ep.spa(1,2)*ep.spab(1,ep.Sum(5,0),
4)*(-(ep.s(0,4,5)*ep.spa(3,
   5))+ep.spa(4,5)*ep.spab(3,
  ep.Sum(5,0),4)))-
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2))/
(ep.s(3,4,5)*ep.spa(4,5)*
 ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),1),2)*
 ep.spab(3,ep.Sum(4,5),2))/
(ep.s(0,1,2)*ep.spa(4,5)*
 ep.spab(3,ep.Sum(4,5),0)*
 ep.spb(2,1)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2)))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spab(0,ep.Sum(4,3),2))/
(ep.s(2,3,4)*ep.spa(0,1)*
 ep.spa(0,5)*ep.spab(5,ep.Sum(4,3),
2)*ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),4),2)*
 ep.spb(4,2))/(ep.spa(0,5)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spb(3,2)*(ep.spa(3,5)*
 ep.spb(3,2)+ep.spa(4,5)*
 ep.spb(4,2))*ep.spb(4,3))
); }

template <class T> complex<T> A4q2g44145_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44290_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*pow(ep.spb(5,1),
3)*ep.spa(0,5)*ep.spa(2,4))/
(ep.s(2,3,4)*ep.spa(3,4)*
 (-(ep.s(2,3,4)*ep.spa(0,2))-
ep.spa(0,1)*ep.spab(2,ep.Sum(3,4),
  1))*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spb(2,1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2)*
 ep.spa(0,4)*ep.spb(3,1))/
(ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(0,ep.Sum(1,2),3)*
 ep.spab(4,ep.Sum(2,3),1)*
 ep.spb(2,1)*ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spa(0,2),3)*pow(ep.spb(5,4),2)*
 ep.spb(5,3))/(ep.s(3,4,5)*
 ep.spa(0,1)*ep.spab(0,ep.Sum(2,1),
3)*ep.spab(2,ep.Sum(1,0),5)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g44320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3)*pow(ep.spb(5,4),
2))/(ep.s(0,4,5)*ep.spa(2,3)*
 ep.spab(1,ep.Sum(2,3),4)*
 ep.spab(3,ep.Sum(1,2),0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(4,5),2),2)*
 ep.spab(4,ep.Sum(0,1),2))/
(ep.s(3,4,5)*ep.spa(3,4)*
 ep.spa(4,5)*ep.spab(5,ep.Sum(0,1),
2)*ep.spb(1,0)*ep.spb(2,1))-
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
 ep.spb(3,2))+(complex<T>(0,1)*pow(ep.spa(0,1),2)*
 pow(ep.spb(4,2),2)*ep.spab(1,
ep.Sum(3,4),2))/(ep.s(2,3,4)*
 ep.spa(0,5)*ep.spab(1,ep.Sum(5,0),
2)*ep.spab(5,ep.Sum(3,4),2)*
 ep.spb(4,3))
); }

template <class T> complex<T> A4q2g44430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4)*ep.spab(0,ep.Sum(4,5),
1))/(ep.s(0,4,5)*ep.spa(0,5)*
 ep.spa(4,5)*ep.spab(0,ep.Sum(4,5),
3)*ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(2,1))+
 (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(4,5),2),2)*
 ep.spa(0,4))/(ep.s(0,4,5)*
 ep.spa(0,5)*ep.spa(4,5)*
 ep.spab(4,ep.Sum(0,5),1)*
 ep.spb(3,2))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spb(5,1))/(ep.s(0,1,5)*
 ep.spa(2,3)*ep.spab(4,ep.Sum(5,0),
1)*ep.spb(1,0)*ep.spb(5,0))+
 (complex<T>(0,1)*pow(ep.spab(3,ep.Sum(0,1),5),2)*
 ep.spab(4,ep.Sum(0,1),5)*
 ep.spb(5,1))/(ep.s(2,3,4)*
 ep.spa(3,4)*ep.spab(2,ep.Sum(0,1),
5)*ep.spab(4,ep.Sum(5,0),1)*
 ep.spb(1,0)*ep.spb(5,0))-
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

template <class T> complex<T> A4q2g44435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(1,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,4))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(1,2)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,4)*
ep.spa(1,3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44485_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,3),2)*
 ep.spa(0,4)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g44535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44615_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44695_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(2,3),1),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g44750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44780_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44785_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g44905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spab(0,ep.Sum(1,2),3),2))/
 (ep.s(0,4,5)*ep.spa(0,5)*
ep.spa(4,5)*ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A4q2g44965_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g44975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g44995_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45195_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,4))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(0,4))/
(ep.spa(0,1)*ep.spa(0,5)*
 ep.spa(2,3)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(2,4))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g45235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g45265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45295_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g45700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g45730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g45855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g45915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))-(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g45925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g45955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g46060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g46070_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(1,3),2)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46160_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g46270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g46275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*
ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
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

template <class T> complex<T> A4q2g46285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
 ep.spa(1,2)*ep.spa(3,4)*
 ep.spa(4,5))+(complex<T>(0,1)*pow(ep.spa(0,2),2)*
 ep.spa(0,3))/(ep.spa(0,1)*
 ep.spa(0,5)*ep.spa(2,3)*
 ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g46315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q2g46345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))
); }

template <class T> complex<T> A4q2g46375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,5)*
ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5))

 ); }
#endif


template <class T> complex<T>  (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
#if _FAST_COMPILE_eval==1
case 280 : return &A4q2g280_eval;
case 310 : return &A4q2g310_eval;
case 315 : return &A4q2g315_eval;
case 340 : return &A4q2g340_eval;
case 350 : return &A4q2g350_eval;
case 370 : return &A4q2g370_eval;
case 375 : return &A4q2g375_eval;
case 380 : return &A4q2g380_eval;
case 385 : return &A4q2g385_eval;
case 490 : return &A4q2g490_eval;
case 495 : return &A4q2g495_eval;
case 525 : return &A4q2g525_eval;
case 550 : return &A4q2g550_eval;
case 555 : return &A4q2g555_eval;
case 560 : return &A4q2g560_eval;
case 565 : return &A4q2g565_eval;
case 585 : return &A4q2g585_eval;
case 595 : return &A4q2g595_eval;
case 700 : return &A4q2g700_eval;
case 710 : return &A4q2g710_eval;
case 730 : return &A4q2g730_eval;
case 735 : return &A4q2g735_eval;
case 740 : return &A4q2g740_eval;
case 745 : return &A4q2g745_eval;
case 770 : return &A4q2g770_eval;
case 800 : return &A4q2g800_eval;
case 805 : return &A4q2g805_eval;
case 910 : return &A4q2g910_eval;
case 915 : return &A4q2g915_eval;
case 920 : return &A4q2g920_eval;
case 925 : return &A4q2g925_eval;
case 945 : return &A4q2g945_eval;
case 955 : return &A4q2g955_eval;
case 980 : return &A4q2g980_eval;
case 985 : return &A4q2g985_eval;
case 1015 : return &A4q2g1015_eval;
case 1360 : return &A4q2g1360_eval;
case 1390 : return &A4q2g1390_eval;
case 1395 : return &A4q2g1395_eval;
case 1420 : return &A4q2g1420_eval;
case 1430 : return &A4q2g1430_eval;
case 1450 : return &A4q2g1450_eval;
case 1455 : return &A4q2g1455_eval;
case 1460 : return &A4q2g1460_eval;
case 1465 : return &A4q2g1465_eval;
case 1540 : return &A4q2g1540_eval;
case 1660 : return &A4q2g1660_eval;
case 1680 : return &A4q2g1680_eval;
case 1685 : return &A4q2g1685_eval;
case 1690 : return &A4q2g1690_eval;
case 1720 : return &A4q2g1720_eval;
case 1750 : return &A4q2g1750_eval;
case 1755 : return &A4q2g1755_eval;
case 1840 : return &A4q2g1840_eval;
case 1860 : return &A4q2g1860_eval;
case 1865 : return &A4q2g1865_eval;
case 1870 : return &A4q2g1870_eval;
case 1875 : return &A4q2g1875_eval;
case 1890 : return &A4q2g1890_eval;
case 1895 : return &A4q2g1895_eval;
case 1905 : return &A4q2g1905_eval;
case 1930 : return &A4q2g1930_eval;
case 1935 : return &A4q2g1935_eval;
case 1960 : return &A4q2g1960_eval;
case 1970 : return &A4q2g1970_eval;
case 2020 : return &A4q2g2020_eval;
case 2040 : return &A4q2g2040_eval;
case 2045 : return &A4q2g2045_eval;
case 2050 : return &A4q2g2050_eval;
case 2090 : return &A4q2g2090_eval;
case 2100 : return &A4q2g2100_eval;
case 2105 : return &A4q2g2105_eval;
case 2120 : return &A4q2g2120_eval;
case 2140 : return &A4q2g2140_eval;
case 2150 : return &A4q2g2150_eval;
case 2170 : return &A4q2g2170_eval;
case 2175 : return &A4q2g2175_eval;
case 2180 : return &A4q2g2180_eval;
case 2185 : return &A4q2g2185_eval;
case 2200 : return &A4q2g2200_eval;
case 2220 : return &A4q2g2220_eval;
case 2225 : return &A4q2g2225_eval;
case 2230 : return &A4q2g2230_eval;
case 2235 : return &A4q2g2235_eval;
case 2250 : return &A4q2g2250_eval;
case 2255 : return &A4q2g2255_eval;
case 2265 : return &A4q2g2265_eval;
case 2270 : return &A4q2g2270_eval;
case 2280 : return &A4q2g2280_eval;
case 2285 : return &A4q2g2285_eval;
case 2300 : return &A4q2g2300_eval;
case 2305 : return &A4q2g2305_eval;
case 2310 : return &A4q2g2310_eval;
case 2315 : return &A4q2g2315_eval;
case 2335 : return &A4q2g2335_eval;
case 2350 : return &A4q2g2350_eval;
case 2355 : return &A4q2g2355_eval;
case 2360 : return &A4q2g2360_eval;
case 2365 : return &A4q2g2365_eval;
case 2440 : return &A4q2g2440_eval;
case 2470 : return &A4q2g2470_eval;
case 2475 : return &A4q2g2475_eval;
case 2500 : return &A4q2g2500_eval;
case 2510 : return &A4q2g2510_eval;
case 2530 : return &A4q2g2530_eval;
case 2535 : return &A4q2g2535_eval;
case 2540 : return &A4q2g2540_eval;
case 2545 : return &A4q2g2545_eval;
case 2650 : return &A4q2g2650_eval;
case 2655 : return &A4q2g2655_eval;
case 2685 : return &A4q2g2685_eval;
case 2710 : return &A4q2g2710_eval;
case 2715 : return &A4q2g2715_eval;
case 2720 : return &A4q2g2720_eval;
case 2725 : return &A4q2g2725_eval;
case 2745 : return &A4q2g2745_eval;
case 2755 : return &A4q2g2755_eval;
case 2830 : return &A4q2g2830_eval;
case 2835 : return &A4q2g2835_eval;
case 2920 : return &A4q2g2920_eval;
case 2940 : return &A4q2g2940_eval;
case 2945 : return &A4q2g2945_eval;
case 2950 : return &A4q2g2950_eval;
case 2955 : return &A4q2g2955_eval;
case 2970 : return &A4q2g2970_eval;
case 2975 : return &A4q2g2975_eval;
case 2985 : return &A4q2g2985_eval;
case 3010 : return &A4q2g3010_eval;
case 3015 : return &A4q2g3015_eval;
case 3045 : return &A4q2g3045_eval;
case 3135 : return &A4q2g3135_eval;
case 3150 : return &A4q2g3150_eval;
case 3155 : return &A4q2g3155_eval;
case 3165 : return &A4q2g3165_eval;
case 3225 : return &A4q2g3225_eval;
case 3250 : return &A4q2g3250_eval;
case 3255 : return &A4q2g3255_eval;
case 3260 : return &A4q2g3260_eval;
case 3265 : return &A4q2g3265_eval;
case 3280 : return &A4q2g3280_eval;
case 3300 : return &A4q2g3300_eval;
case 3305 : return &A4q2g3305_eval;
case 3310 : return &A4q2g3310_eval;
case 3315 : return &A4q2g3315_eval;
case 3330 : return &A4q2g3330_eval;
case 3335 : return &A4q2g3335_eval;
case 3345 : return &A4q2g3345_eval;
case 3350 : return &A4q2g3350_eval;
case 3360 : return &A4q2g3360_eval;
case 3365 : return &A4q2g3365_eval;
case 3380 : return &A4q2g3380_eval;
case 3385 : return &A4q2g3385_eval;
case 3390 : return &A4q2g3390_eval;
case 3395 : return &A4q2g3395_eval;
case 3415 : return &A4q2g3415_eval;
case 3430 : return &A4q2g3430_eval;
case 3435 : return &A4q2g3435_eval;
case 3440 : return &A4q2g3440_eval;
case 3445 : return &A4q2g3445_eval;
case 3465 : return &A4q2g3465_eval;
case 3475 : return &A4q2g3475_eval;
case 3495 : return &A4q2g3495_eval;
case 3510 : return &A4q2g3510_eval;
case 3515 : return &A4q2g3515_eval;
case 3525 : return &A4q2g3525_eval;
case 3565 : return &A4q2g3565_eval;
case 3570 : return &A4q2g3570_eval;
case 3575 : return &A4q2g3575_eval;
case 3595 : return &A4q2g3595_eval;
case 3645 : return &A4q2g3645_eval;
case 3655 : return &A4q2g3655_eval;
case 3730 : return &A4q2g3730_eval;
case 3735 : return &A4q2g3735_eval;
case 3765 : return &A4q2g3765_eval;
case 3790 : return &A4q2g3790_eval;
case 3795 : return &A4q2g3795_eval;
case 3800 : return &A4q2g3800_eval;
case 3805 : return &A4q2g3805_eval;
case 3825 : return &A4q2g3825_eval;
case 3835 : return &A4q2g3835_eval;
case 3940 : return &A4q2g3940_eval;
case 3950 : return &A4q2g3950_eval;
case 3970 : return &A4q2g3970_eval;
case 3975 : return &A4q2g3975_eval;
case 3980 : return &A4q2g3980_eval;
case 3985 : return &A4q2g3985_eval;
case 4010 : return &A4q2g4010_eval;
case 4040 : return &A4q2g4040_eval;
case 4045 : return &A4q2g4045_eval;
case 4120 : return &A4q2g4120_eval;
case 4130 : return &A4q2g4130_eval;
case 4180 : return &A4q2g4180_eval;
case 4200 : return &A4q2g4200_eval;
case 4205 : return &A4q2g4205_eval;
case 4210 : return &A4q2g4210_eval;
case 4250 : return &A4q2g4250_eval;
case 4260 : return &A4q2g4260_eval;
case 4265 : return &A4q2g4265_eval;
case 4280 : return &A4q2g4280_eval;
case 4300 : return &A4q2g4300_eval;
case 4310 : return &A4q2g4310_eval;
case 4330 : return &A4q2g4330_eval;
case 4335 : return &A4q2g4335_eval;
case 4340 : return &A4q2g4340_eval;
case 4345 : return &A4q2g4345_eval;
case 4360 : return &A4q2g4360_eval;
case 4380 : return &A4q2g4380_eval;
case 4385 : return &A4q2g4385_eval;
case 4390 : return &A4q2g4390_eval;
case 4395 : return &A4q2g4395_eval;
case 4410 : return &A4q2g4410_eval;
case 4415 : return &A4q2g4415_eval;
case 4425 : return &A4q2g4425_eval;
case 4430 : return &A4q2g4430_eval;
case 4440 : return &A4q2g4440_eval;
case 4445 : return &A4q2g4445_eval;
case 4460 : return &A4q2g4460_eval;
case 4465 : return &A4q2g4465_eval;
case 4470 : return &A4q2g4470_eval;
case 4475 : return &A4q2g4475_eval;
case 4495 : return &A4q2g4495_eval;
case 4510 : return &A4q2g4510_eval;
case 4515 : return &A4q2g4515_eval;
case 4520 : return &A4q2g4520_eval;
case 4525 : return &A4q2g4525_eval;
case 4550 : return &A4q2g4550_eval;
case 4610 : return &A4q2g4610_eval;
case 4620 : return &A4q2g4620_eval;
case 4625 : return &A4q2g4625_eval;
case 4640 : return &A4q2g4640_eval;
case 4730 : return &A4q2g4730_eval;
case 4760 : return &A4q2g4760_eval;
case 4765 : return &A4q2g4765_eval;
case 4790 : return &A4q2g4790_eval;
case 4800 : return &A4q2g4800_eval;
case 4805 : return &A4q2g4805_eval;
case 4820 : return &A4q2g4820_eval;
case 4825 : return &A4q2g4825_eval;
case 4830 : return &A4q2g4830_eval;
case 4835 : return &A4q2g4835_eval;
case 4855 : return &A4q2g4855_eval;
case 4940 : return &A4q2g4940_eval;
case 4945 : return &A4q2g4945_eval;
case 5020 : return &A4q2g5020_eval;
case 5030 : return &A4q2g5030_eval;
case 5050 : return &A4q2g5050_eval;
case 5055 : return &A4q2g5055_eval;
case 5060 : return &A4q2g5060_eval;
case 5065 : return &A4q2g5065_eval;
case 5090 : return &A4q2g5090_eval;
case 5120 : return &A4q2g5120_eval;
case 5125 : return &A4q2g5125_eval;
case 5230 : return &A4q2g5230_eval;
case 5235 : return &A4q2g5235_eval;
case 5240 : return &A4q2g5240_eval;
case 5245 : return &A4q2g5245_eval;
case 5265 : return &A4q2g5265_eval;
case 5275 : return &A4q2g5275_eval;
case 5300 : return &A4q2g5300_eval;
case 5305 : return &A4q2g5305_eval;
case 5335 : return &A4q2g5335_eval;
case 5410 : return &A4q2g5410_eval;
case 5415 : return &A4q2g5415_eval;
case 5420 : return &A4q2g5420_eval;
case 5425 : return &A4q2g5425_eval;
case 5440 : return &A4q2g5440_eval;
case 5460 : return &A4q2g5460_eval;
case 5465 : return &A4q2g5465_eval;
case 5470 : return &A4q2g5470_eval;
case 5475 : return &A4q2g5475_eval;
case 5490 : return &A4q2g5490_eval;
case 5495 : return &A4q2g5495_eval;
case 5505 : return &A4q2g5505_eval;
case 5510 : return &A4q2g5510_eval;
case 5520 : return &A4q2g5520_eval;
case 5525 : return &A4q2g5525_eval;
case 5540 : return &A4q2g5540_eval;
case 5545 : return &A4q2g5545_eval;
case 5550 : return &A4q2g5550_eval;
case 5555 : return &A4q2g5555_eval;
case 5575 : return &A4q2g5575_eval;
case 5590 : return &A4q2g5590_eval;
case 5595 : return &A4q2g5595_eval;
case 5600 : return &A4q2g5600_eval;
case 5605 : return &A4q2g5605_eval;
case 5625 : return &A4q2g5625_eval;
case 5635 : return &A4q2g5635_eval;
case 5655 : return &A4q2g5655_eval;
case 5670 : return &A4q2g5670_eval;
case 5675 : return &A4q2g5675_eval;
case 5685 : return &A4q2g5685_eval;
case 5725 : return &A4q2g5725_eval;
case 5730 : return &A4q2g5730_eval;
case 5735 : return &A4q2g5735_eval;
case 5755 : return &A4q2g5755_eval;
case 5805 : return &A4q2g5805_eval;
case 5815 : return &A4q2g5815_eval;
case 5840 : return &A4q2g5840_eval;
case 5845 : return &A4q2g5845_eval;
case 5870 : return &A4q2g5870_eval;
case 5880 : return &A4q2g5880_eval;
case 5885 : return &A4q2g5885_eval;
case 5900 : return &A4q2g5900_eval;
case 5905 : return &A4q2g5905_eval;
case 5910 : return &A4q2g5910_eval;
case 5915 : return &A4q2g5915_eval;
case 5935 : return &A4q2g5935_eval;
case 6020 : return &A4q2g6020_eval;
case 6025 : return &A4q2g6025_eval;
case 6055 : return &A4q2g6055_eval;
case 6085 : return &A4q2g6085_eval;
case 6090 : return &A4q2g6090_eval;
case 6095 : return &A4q2g6095_eval;
case 6115 : return &A4q2g6115_eval;
case 6235 : return &A4q2g6235_eval;
case 6310 : return &A4q2g6310_eval;
case 6315 : return &A4q2g6315_eval;
case 6320 : return &A4q2g6320_eval;
case 6325 : return &A4q2g6325_eval;
case 6345 : return &A4q2g6345_eval;
case 6355 : return &A4q2g6355_eval;
case 6380 : return &A4q2g6380_eval;
case 6385 : return &A4q2g6385_eval;
case 6415 : return &A4q2g6415_eval;
case 6760 : return &A4q2g6760_eval;
case 6790 : return &A4q2g6790_eval;
case 6795 : return &A4q2g6795_eval;
case 6820 : return &A4q2g6820_eval;
case 6830 : return &A4q2g6830_eval;
case 6850 : return &A4q2g6850_eval;
case 6855 : return &A4q2g6855_eval;
case 6860 : return &A4q2g6860_eval;
case 6865 : return &A4q2g6865_eval;
case 6970 : return &A4q2g6970_eval;
case 6975 : return &A4q2g6975_eval;
case 7005 : return &A4q2g7005_eval;
case 7030 : return &A4q2g7030_eval;
case 7035 : return &A4q2g7035_eval;
case 7040 : return &A4q2g7040_eval;
case 7045 : return &A4q2g7045_eval;
case 7065 : return &A4q2g7065_eval;
case 7075 : return &A4q2g7075_eval;
case 7180 : return &A4q2g7180_eval;
case 7190 : return &A4q2g7190_eval;
case 7210 : return &A4q2g7210_eval;
case 7215 : return &A4q2g7215_eval;
case 7220 : return &A4q2g7220_eval;
case 7225 : return &A4q2g7225_eval;
case 7250 : return &A4q2g7250_eval;
case 7280 : return &A4q2g7280_eval;
case 7285 : return &A4q2g7285_eval;
case 7390 : return &A4q2g7390_eval;
case 7395 : return &A4q2g7395_eval;
case 7400 : return &A4q2g7400_eval;
case 7405 : return &A4q2g7405_eval;
case 7425 : return &A4q2g7425_eval;
case 7435 : return &A4q2g7435_eval;
case 7460 : return &A4q2g7460_eval;
case 7465 : return &A4q2g7465_eval;
case 7495 : return &A4q2g7495_eval;
case 7840 : return &A4q2g7840_eval;
case 7870 : return &A4q2g7870_eval;
case 7875 : return &A4q2g7875_eval;
case 7900 : return &A4q2g7900_eval;
case 7910 : return &A4q2g7910_eval;
case 7930 : return &A4q2g7930_eval;
case 7935 : return &A4q2g7935_eval;
case 7940 : return &A4q2g7940_eval;
case 7945 : return &A4q2g7945_eval;
case 8020 : return &A4q2g8020_eval;
case 8140 : return &A4q2g8140_eval;
case 8160 : return &A4q2g8160_eval;
case 8165 : return &A4q2g8165_eval;
case 8170 : return &A4q2g8170_eval;
case 8200 : return &A4q2g8200_eval;
case 8230 : return &A4q2g8230_eval;
case 8235 : return &A4q2g8235_eval;
case 8320 : return &A4q2g8320_eval;
case 8340 : return &A4q2g8340_eval;
case 8345 : return &A4q2g8345_eval;
case 8350 : return &A4q2g8350_eval;
case 8355 : return &A4q2g8355_eval;
case 8370 : return &A4q2g8370_eval;
case 8375 : return &A4q2g8375_eval;
case 8385 : return &A4q2g8385_eval;
case 8410 : return &A4q2g8410_eval;
case 8415 : return &A4q2g8415_eval;
case 8440 : return &A4q2g8440_eval;
case 8450 : return &A4q2g8450_eval;
case 8500 : return &A4q2g8500_eval;
case 8520 : return &A4q2g8520_eval;
case 8525 : return &A4q2g8525_eval;
case 8530 : return &A4q2g8530_eval;
case 8570 : return &A4q2g8570_eval;
case 8580 : return &A4q2g8580_eval;
case 8585 : return &A4q2g8585_eval;
case 8600 : return &A4q2g8600_eval;
case 8620 : return &A4q2g8620_eval;
case 8630 : return &A4q2g8630_eval;
case 8650 : return &A4q2g8650_eval;
case 8655 : return &A4q2g8655_eval;
case 8660 : return &A4q2g8660_eval;
case 8665 : return &A4q2g8665_eval;
case 8680 : return &A4q2g8680_eval;
case 8700 : return &A4q2g8700_eval;
case 8705 : return &A4q2g8705_eval;
case 8710 : return &A4q2g8710_eval;
case 8715 : return &A4q2g8715_eval;
case 8730 : return &A4q2g8730_eval;
case 8735 : return &A4q2g8735_eval;
case 8745 : return &A4q2g8745_eval;
case 8750 : return &A4q2g8750_eval;
case 8760 : return &A4q2g8760_eval;
case 8765 : return &A4q2g8765_eval;
case 8780 : return &A4q2g8780_eval;
case 8785 : return &A4q2g8785_eval;
case 8790 : return &A4q2g8790_eval;
case 8795 : return &A4q2g8795_eval;
case 8815 : return &A4q2g8815_eval;
case 8830 : return &A4q2g8830_eval;
case 8835 : return &A4q2g8835_eval;
case 8840 : return &A4q2g8840_eval;
case 8845 : return &A4q2g8845_eval;
case 8920 : return &A4q2g8920_eval;
case 8950 : return &A4q2g8950_eval;
case 8955 : return &A4q2g8955_eval;
case 8980 : return &A4q2g8980_eval;
case 8990 : return &A4q2g8990_eval;
case 9010 : return &A4q2g9010_eval;
case 9015 : return &A4q2g9015_eval;
case 9020 : return &A4q2g9020_eval;
case 9025 : return &A4q2g9025_eval;
case 9100 : return &A4q2g9100_eval;
case 9220 : return &A4q2g9220_eval;
case 9240 : return &A4q2g9240_eval;
case 9245 : return &A4q2g9245_eval;
case 9250 : return &A4q2g9250_eval;
case 9280 : return &A4q2g9280_eval;
case 9940 : return &A4q2g9940_eval;
case 9960 : return &A4q2g9960_eval;
case 9965 : return &A4q2g9965_eval;
case 9970 : return &A4q2g9970_eval;
case 10080 : return &A4q2g10080_eval;
case 10085 : return &A4q2g10085_eval;
case 10110 : return &A4q2g10110_eval;
case 10115 : return &A4q2g10115_eval;
case 10120 : return &A4q2g10120_eval;
case 10140 : return &A4q2g10140_eval;
case 10145 : return &A4q2g10145_eval;
case 10150 : return &A4q2g10150_eval;
case 10180 : return &A4q2g10180_eval;
case 10300 : return &A4q2g10300_eval;
case 10320 : return &A4q2g10320_eval;
case 10325 : return &A4q2g10325_eval;
case 10330 : return &A4q2g10330_eval;
case 10360 : return &A4q2g10360_eval;
case 10390 : return &A4q2g10390_eval;
case 10395 : return &A4q2g10395_eval;
case 10480 : return &A4q2g10480_eval;
case 10500 : return &A4q2g10500_eval;
case 10505 : return &A4q2g10505_eval;
case 10510 : return &A4q2g10510_eval;
case 10515 : return &A4q2g10515_eval;
case 10530 : return &A4q2g10530_eval;
case 10535 : return &A4q2g10535_eval;
case 10545 : return &A4q2g10545_eval;
case 10570 : return &A4q2g10570_eval;
case 10575 : return &A4q2g10575_eval;
case 11020 : return &A4q2g11020_eval;
case 11040 : return &A4q2g11040_eval;
case 11045 : return &A4q2g11045_eval;
case 11050 : return &A4q2g11050_eval;
case 11160 : return &A4q2g11160_eval;
case 11165 : return &A4q2g11165_eval;
case 11190 : return &A4q2g11190_eval;
case 11195 : return &A4q2g11195_eval;
case 11200 : return &A4q2g11200_eval;
case 11220 : return &A4q2g11220_eval;
case 11225 : return &A4q2g11225_eval;
case 11230 : return &A4q2g11230_eval;
case 11235 : return &A4q2g11235_eval;
case 11250 : return &A4q2g11250_eval;
case 11255 : return &A4q2g11255_eval;
case 11265 : return &A4q2g11265_eval;
case 11340 : return &A4q2g11340_eval;
case 11345 : return &A4q2g11345_eval;
case 11370 : return &A4q2g11370_eval;
case 11375 : return &A4q2g11375_eval;
case 11415 : return &A4q2g11415_eval;
case 11430 : return &A4q2g11430_eval;
case 11435 : return &A4q2g11435_eval;
case 11445 : return &A4q2g11445_eval;
case 11470 : return &A4q2g11470_eval;
case 11475 : return &A4q2g11475_eval;
case 11560 : return &A4q2g11560_eval;
case 11580 : return &A4q2g11580_eval;
case 11585 : return &A4q2g11585_eval;
case 11590 : return &A4q2g11590_eval;
case 11595 : return &A4q2g11595_eval;
case 11610 : return &A4q2g11610_eval;
case 11615 : return &A4q2g11615_eval;
case 11625 : return &A4q2g11625_eval;
case 11650 : return &A4q2g11650_eval;
case 11655 : return &A4q2g11655_eval;
case 11680 : return &A4q2g11680_eval;
case 11690 : return &A4q2g11690_eval;
case 11740 : return &A4q2g11740_eval;
case 11760 : return &A4q2g11760_eval;
case 11765 : return &A4q2g11765_eval;
case 11770 : return &A4q2g11770_eval;
case 11810 : return &A4q2g11810_eval;
case 11820 : return &A4q2g11820_eval;
case 11825 : return &A4q2g11825_eval;
case 11840 : return &A4q2g11840_eval;
case 11860 : return &A4q2g11860_eval;
case 11870 : return &A4q2g11870_eval;
case 12100 : return &A4q2g12100_eval;
case 12120 : return &A4q2g12120_eval;
case 12125 : return &A4q2g12125_eval;
case 12130 : return &A4q2g12130_eval;
case 12240 : return &A4q2g12240_eval;
case 12245 : return &A4q2g12245_eval;
case 12270 : return &A4q2g12270_eval;
case 12275 : return &A4q2g12275_eval;
case 12280 : return &A4q2g12280_eval;
case 12300 : return &A4q2g12300_eval;
case 12305 : return &A4q2g12305_eval;
case 12310 : return &A4q2g12310_eval;
case 12530 : return &A4q2g12530_eval;
case 12540 : return &A4q2g12540_eval;
case 12545 : return &A4q2g12545_eval;
case 12560 : return &A4q2g12560_eval;
case 12600 : return &A4q2g12600_eval;
case 12605 : return &A4q2g12605_eval;
case 12630 : return &A4q2g12630_eval;
case 12635 : return &A4q2g12635_eval;
case 12710 : return &A4q2g12710_eval;
case 12720 : return &A4q2g12720_eval;
case 12725 : return &A4q2g12725_eval;
case 12740 : return &A4q2g12740_eval;
case 12760 : return &A4q2g12760_eval;
case 12770 : return &A4q2g12770_eval;
case 12820 : return &A4q2g12820_eval;
case 12840 : return &A4q2g12840_eval;
case 12845 : return &A4q2g12845_eval;
case 12850 : return &A4q2g12850_eval;
case 12890 : return &A4q2g12890_eval;
case 12900 : return &A4q2g12900_eval;
case 12905 : return &A4q2g12905_eval;
case 12920 : return &A4q2g12920_eval;
case 12940 : return &A4q2g12940_eval;
case 12950 : return &A4q2g12950_eval;
case 12970 : return &A4q2g12970_eval;
case 12975 : return &A4q2g12975_eval;
case 12980 : return &A4q2g12980_eval;
case 12985 : return &A4q2g12985_eval;
case 13000 : return &A4q2g13000_eval;
case 13020 : return &A4q2g13020_eval;
case 13025 : return &A4q2g13025_eval;
case 13030 : return &A4q2g13030_eval;
case 13035 : return &A4q2g13035_eval;
case 13050 : return &A4q2g13050_eval;
case 13055 : return &A4q2g13055_eval;
case 13065 : return &A4q2g13065_eval;
case 13070 : return &A4q2g13070_eval;
case 13080 : return &A4q2g13080_eval;
case 13085 : return &A4q2g13085_eval;
case 13100 : return &A4q2g13100_eval;
case 13105 : return &A4q2g13105_eval;
case 13110 : return &A4q2g13110_eval;
case 13115 : return &A4q2g13115_eval;
case 13135 : return &A4q2g13135_eval;
case 13150 : return &A4q2g13150_eval;
case 13155 : return &A4q2g13155_eval;
case 13160 : return &A4q2g13160_eval;
case 13165 : return &A4q2g13165_eval;
case 13180 : return &A4q2g13180_eval;
case 13200 : return &A4q2g13200_eval;
case 13205 : return &A4q2g13205_eval;
case 13210 : return &A4q2g13210_eval;
case 13320 : return &A4q2g13320_eval;
case 13325 : return &A4q2g13325_eval;
case 13350 : return &A4q2g13350_eval;
case 13355 : return &A4q2g13355_eval;
case 13360 : return &A4q2g13360_eval;
case 13380 : return &A4q2g13380_eval;
case 13385 : return &A4q2g13385_eval;
case 13390 : return &A4q2g13390_eval;
case 13395 : return &A4q2g13395_eval;
case 13410 : return &A4q2g13410_eval;
case 13415 : return &A4q2g13415_eval;
case 13425 : return &A4q2g13425_eval;
case 13500 : return &A4q2g13500_eval;
case 13505 : return &A4q2g13505_eval;
case 13530 : return &A4q2g13530_eval;
case 13535 : return &A4q2g13535_eval;
case 13575 : return &A4q2g13575_eval;
case 13590 : return &A4q2g13590_eval;
case 13595 : return &A4q2g13595_eval;
case 13605 : return &A4q2g13605_eval;
case 13610 : return &A4q2g13610_eval;
case 13620 : return &A4q2g13620_eval;
case 13625 : return &A4q2g13625_eval;
case 13640 : return &A4q2g13640_eval;
case 13680 : return &A4q2g13680_eval;
case 13685 : return &A4q2g13685_eval;
case 13710 : return &A4q2g13710_eval;
case 13715 : return &A4q2g13715_eval;
case 13790 : return &A4q2g13790_eval;
case 13800 : return &A4q2g13800_eval;
case 13805 : return &A4q2g13805_eval;
case 13820 : return &A4q2g13820_eval;
case 13825 : return &A4q2g13825_eval;
case 13830 : return &A4q2g13830_eval;
case 13835 : return &A4q2g13835_eval;
case 13855 : return &A4q2g13855_eval;
case 13860 : return &A4q2g13860_eval;
case 13865 : return &A4q2g13865_eval;
case 13890 : return &A4q2g13890_eval;
case 13895 : return &A4q2g13895_eval;
case 14005 : return &A4q2g14005_eval;
case 14010 : return &A4q2g14010_eval;
case 14015 : return &A4q2g14015_eval;
case 14035 : return &A4q2g14035_eval;
case 14050 : return &A4q2g14050_eval;
case 14055 : return &A4q2g14055_eval;
case 14060 : return &A4q2g14060_eval;
case 14065 : return &A4q2g14065_eval;
case 14080 : return &A4q2g14080_eval;
case 14100 : return &A4q2g14100_eval;
case 14105 : return &A4q2g14105_eval;
case 14110 : return &A4q2g14110_eval;
case 14115 : return &A4q2g14115_eval;
case 14130 : return &A4q2g14130_eval;
case 14135 : return &A4q2g14135_eval;
case 14145 : return &A4q2g14145_eval;
case 14150 : return &A4q2g14150_eval;
case 14160 : return &A4q2g14160_eval;
case 14165 : return &A4q2g14165_eval;
case 14180 : return &A4q2g14180_eval;
case 14185 : return &A4q2g14185_eval;
case 14190 : return &A4q2g14190_eval;
case 14195 : return &A4q2g14195_eval;
case 14215 : return &A4q2g14215_eval;
case 14230 : return &A4q2g14230_eval;
case 14235 : return &A4q2g14235_eval;
case 14240 : return &A4q2g14240_eval;
case 14245 : return &A4q2g14245_eval;
case 14320 : return &A4q2g14320_eval;
case 14350 : return &A4q2g14350_eval;
case 14355 : return &A4q2g14355_eval;
case 14380 : return &A4q2g14380_eval;
case 14390 : return &A4q2g14390_eval;
case 14410 : return &A4q2g14410_eval;
case 14415 : return &A4q2g14415_eval;
case 14420 : return &A4q2g14420_eval;
case 14425 : return &A4q2g14425_eval;
case 14500 : return &A4q2g14500_eval;
case 14620 : return &A4q2g14620_eval;
case 14640 : return &A4q2g14640_eval;
case 14645 : return &A4q2g14645_eval;
case 14650 : return &A4q2g14650_eval;
case 14680 : return &A4q2g14680_eval;
case 14710 : return &A4q2g14710_eval;
case 14715 : return &A4q2g14715_eval;
case 14800 : return &A4q2g14800_eval;
case 14820 : return &A4q2g14820_eval;
case 14825 : return &A4q2g14825_eval;
case 14830 : return &A4q2g14830_eval;
case 14835 : return &A4q2g14835_eval;
case 14850 : return &A4q2g14850_eval;
case 14855 : return &A4q2g14855_eval;
case 14865 : return &A4q2g14865_eval;
case 14890 : return &A4q2g14890_eval;
case 14895 : return &A4q2g14895_eval;
case 14920 : return &A4q2g14920_eval;
case 14930 : return &A4q2g14930_eval;
case 14980 : return &A4q2g14980_eval;
case 15000 : return &A4q2g15000_eval;
case 15005 : return &A4q2g15005_eval;
case 15010 : return &A4q2g15010_eval;
case 15050 : return &A4q2g15050_eval;
case 15060 : return &A4q2g15060_eval;
case 15065 : return &A4q2g15065_eval;
case 15080 : return &A4q2g15080_eval;
case 15100 : return &A4q2g15100_eval;
case 15110 : return &A4q2g15110_eval;
case 15130 : return &A4q2g15130_eval;
case 15135 : return &A4q2g15135_eval;
case 15140 : return &A4q2g15140_eval;
case 15145 : return &A4q2g15145_eval;
case 15160 : return &A4q2g15160_eval;
case 15180 : return &A4q2g15180_eval;
case 15185 : return &A4q2g15185_eval;
case 15190 : return &A4q2g15190_eval;
case 15195 : return &A4q2g15195_eval;
case 15210 : return &A4q2g15210_eval;
case 15215 : return &A4q2g15215_eval;
case 15225 : return &A4q2g15225_eval;
case 15230 : return &A4q2g15230_eval;
case 15240 : return &A4q2g15240_eval;
case 15245 : return &A4q2g15245_eval;
case 15260 : return &A4q2g15260_eval;
case 15265 : return &A4q2g15265_eval;
case 15270 : return &A4q2g15270_eval;
case 15275 : return &A4q2g15275_eval;
case 15295 : return &A4q2g15295_eval;
case 15310 : return &A4q2g15310_eval;
case 15315 : return &A4q2g15315_eval;
case 15320 : return &A4q2g15320_eval;
case 15325 : return &A4q2g15325_eval;
case 15400 : return &A4q2g15400_eval;
case 15430 : return &A4q2g15430_eval;
case 15435 : return &A4q2g15435_eval;
case 15460 : return &A4q2g15460_eval;
case 15470 : return &A4q2g15470_eval;
case 15490 : return &A4q2g15490_eval;
case 15495 : return &A4q2g15495_eval;
case 15500 : return &A4q2g15500_eval;
case 15505 : return &A4q2g15505_eval;
case 15610 : return &A4q2g15610_eval;
case 15615 : return &A4q2g15615_eval;
case 15645 : return &A4q2g15645_eval;
case 15670 : return &A4q2g15670_eval;
case 15675 : return &A4q2g15675_eval;
case 15680 : return &A4q2g15680_eval;
case 15685 : return &A4q2g15685_eval;
case 15705 : return &A4q2g15705_eval;
case 15715 : return &A4q2g15715_eval;
case 15790 : return &A4q2g15790_eval;
case 15795 : return &A4q2g15795_eval;
case 15880 : return &A4q2g15880_eval;
case 15900 : return &A4q2g15900_eval;
case 15905 : return &A4q2g15905_eval;
case 15910 : return &A4q2g15910_eval;
case 15915 : return &A4q2g15915_eval;
case 15930 : return &A4q2g15930_eval;
case 15935 : return &A4q2g15935_eval;
case 15945 : return &A4q2g15945_eval;
case 15970 : return &A4q2g15970_eval;
case 15975 : return &A4q2g15975_eval;
case 16005 : return &A4q2g16005_eval;
case 16095 : return &A4q2g16095_eval;
case 16110 : return &A4q2g16110_eval;
case 16115 : return &A4q2g16115_eval;
case 16125 : return &A4q2g16125_eval;
case 16185 : return &A4q2g16185_eval;
case 16210 : return &A4q2g16210_eval;
case 16215 : return &A4q2g16215_eval;
case 16220 : return &A4q2g16220_eval;
case 16225 : return &A4q2g16225_eval;
case 16240 : return &A4q2g16240_eval;
case 16260 : return &A4q2g16260_eval;
case 16265 : return &A4q2g16265_eval;
case 16270 : return &A4q2g16270_eval;
case 16275 : return &A4q2g16275_eval;
case 16290 : return &A4q2g16290_eval;
case 16295 : return &A4q2g16295_eval;
case 16305 : return &A4q2g16305_eval;
case 16310 : return &A4q2g16310_eval;
case 16320 : return &A4q2g16320_eval;
case 16325 : return &A4q2g16325_eval;
case 16340 : return &A4q2g16340_eval;
case 16345 : return &A4q2g16345_eval;
case 16350 : return &A4q2g16350_eval;
case 16355 : return &A4q2g16355_eval;
case 16375 : return &A4q2g16375_eval;
case 16390 : return &A4q2g16390_eval;
case 16395 : return &A4q2g16395_eval;
case 16400 : return &A4q2g16400_eval;
case 16405 : return &A4q2g16405_eval;
case 16425 : return &A4q2g16425_eval;
case 16435 : return &A4q2g16435_eval;
case 16455 : return &A4q2g16455_eval;
case 16470 : return &A4q2g16470_eval;
case 16475 : return &A4q2g16475_eval;
case 16485 : return &A4q2g16485_eval;
case 16525 : return &A4q2g16525_eval;
case 16530 : return &A4q2g16530_eval;
case 16535 : return &A4q2g16535_eval;
case 16555 : return &A4q2g16555_eval;
case 16605 : return &A4q2g16605_eval;
case 16615 : return &A4q2g16615_eval;
case 16690 : return &A4q2g16690_eval;
case 16695 : return &A4q2g16695_eval;
case 16725 : return &A4q2g16725_eval;
case 16750 : return &A4q2g16750_eval;
case 16755 : return &A4q2g16755_eval;
case 16760 : return &A4q2g16760_eval;
case 16765 : return &A4q2g16765_eval;
case 16785 : return &A4q2g16785_eval;
case 16795 : return &A4q2g16795_eval;
case 16870 : return &A4q2g16870_eval;
case 16875 : return &A4q2g16875_eval;
case 16960 : return &A4q2g16960_eval;
case 16980 : return &A4q2g16980_eval;
case 16985 : return &A4q2g16985_eval;
case 16990 : return &A4q2g16990_eval;
case 16995 : return &A4q2g16995_eval;
case 17010 : return &A4q2g17010_eval;
case 17015 : return &A4q2g17015_eval;
case 17025 : return &A4q2g17025_eval;
case 17050 : return &A4q2g17050_eval;
case 17055 : return &A4q2g17055_eval;
case 17500 : return &A4q2g17500_eval;
case 17520 : return &A4q2g17520_eval;
case 17525 : return &A4q2g17525_eval;
case 17530 : return &A4q2g17530_eval;
case 17640 : return &A4q2g17640_eval;
case 17645 : return &A4q2g17645_eval;
case 17670 : return &A4q2g17670_eval;
case 17675 : return &A4q2g17675_eval;
case 17680 : return &A4q2g17680_eval;
case 17700 : return &A4q2g17700_eval;
case 17705 : return &A4q2g17705_eval;
case 17710 : return &A4q2g17710_eval;
case 17715 : return &A4q2g17715_eval;
case 17730 : return &A4q2g17730_eval;
case 17735 : return &A4q2g17735_eval;
case 17745 : return &A4q2g17745_eval;
case 17820 : return &A4q2g17820_eval;
case 17825 : return &A4q2g17825_eval;
case 17850 : return &A4q2g17850_eval;
case 17855 : return &A4q2g17855_eval;
case 17895 : return &A4q2g17895_eval;
case 17910 : return &A4q2g17910_eval;
case 17915 : return &A4q2g17915_eval;
case 17925 : return &A4q2g17925_eval;
case 17950 : return &A4q2g17950_eval;
case 17955 : return &A4q2g17955_eval;
case 18040 : return &A4q2g18040_eval;
case 18060 : return &A4q2g18060_eval;
case 18065 : return &A4q2g18065_eval;
case 18070 : return &A4q2g18070_eval;
case 18075 : return &A4q2g18075_eval;
case 18090 : return &A4q2g18090_eval;
case 18095 : return &A4q2g18095_eval;
case 18105 : return &A4q2g18105_eval;
case 18130 : return &A4q2g18130_eval;
case 18135 : return &A4q2g18135_eval;
case 18165 : return &A4q2g18165_eval;
case 18255 : return &A4q2g18255_eval;
case 18270 : return &A4q2g18270_eval;
case 18275 : return &A4q2g18275_eval;
case 18285 : return &A4q2g18285_eval;
case 18345 : return &A4q2g18345_eval;
case 18795 : return &A4q2g18795_eval;
case 18810 : return &A4q2g18810_eval;
case 18815 : return &A4q2g18815_eval;
case 18825 : return &A4q2g18825_eval;
case 18900 : return &A4q2g18900_eval;
case 18905 : return &A4q2g18905_eval;
case 18930 : return &A4q2g18930_eval;
case 18935 : return &A4q2g18935_eval;
case 18975 : return &A4q2g18975_eval;
case 18990 : return &A4q2g18990_eval;
case 18995 : return &A4q2g18995_eval;
case 19005 : return &A4q2g19005_eval;
case 19245 : return &A4q2g19245_eval;
case 19335 : return &A4q2g19335_eval;
case 19350 : return &A4q2g19350_eval;
case 19355 : return &A4q2g19355_eval;
case 19365 : return &A4q2g19365_eval;
case 19425 : return &A4q2g19425_eval;
case 19450 : return &A4q2g19450_eval;
case 19455 : return &A4q2g19455_eval;
case 19460 : return &A4q2g19460_eval;
case 19465 : return &A4q2g19465_eval;
case 19480 : return &A4q2g19480_eval;
case 19500 : return &A4q2g19500_eval;
case 19505 : return &A4q2g19505_eval;
case 19510 : return &A4q2g19510_eval;
case 19515 : return &A4q2g19515_eval;
case 19530 : return &A4q2g19530_eval;
case 19535 : return &A4q2g19535_eval;
case 19545 : return &A4q2g19545_eval;
case 19550 : return &A4q2g19550_eval;
case 19560 : return &A4q2g19560_eval;
case 19565 : return &A4q2g19565_eval;
case 19580 : return &A4q2g19580_eval;
case 19585 : return &A4q2g19585_eval;
case 19590 : return &A4q2g19590_eval;
case 19595 : return &A4q2g19595_eval;
case 19615 : return &A4q2g19615_eval;
case 19630 : return &A4q2g19630_eval;
case 19635 : return &A4q2g19635_eval;
case 19640 : return &A4q2g19640_eval;
case 19645 : return &A4q2g19645_eval;
case 19660 : return &A4q2g19660_eval;
case 19680 : return &A4q2g19680_eval;
case 19685 : return &A4q2g19685_eval;
case 19690 : return &A4q2g19690_eval;
case 19800 : return &A4q2g19800_eval;
case 19805 : return &A4q2g19805_eval;
case 19830 : return &A4q2g19830_eval;
case 19835 : return &A4q2g19835_eval;
case 19840 : return &A4q2g19840_eval;
case 19860 : return &A4q2g19860_eval;
case 19865 : return &A4q2g19865_eval;
case 19870 : return &A4q2g19870_eval;
case 19875 : return &A4q2g19875_eval;
case 19890 : return &A4q2g19890_eval;
case 19895 : return &A4q2g19895_eval;
case 19905 : return &A4q2g19905_eval;
case 19980 : return &A4q2g19980_eval;
case 19985 : return &A4q2g19985_eval;
case 20010 : return &A4q2g20010_eval;
case 20015 : return &A4q2g20015_eval;
case 20055 : return &A4q2g20055_eval;
case 20070 : return &A4q2g20070_eval;
case 20075 : return &A4q2g20075_eval;
case 20085 : return &A4q2g20085_eval;
case 20090 : return &A4q2g20090_eval;
case 20100 : return &A4q2g20100_eval;
case 20105 : return &A4q2g20105_eval;
case 20120 : return &A4q2g20120_eval;
case 20160 : return &A4q2g20160_eval;
case 20165 : return &A4q2g20165_eval;
case 20190 : return &A4q2g20190_eval;
case 20195 : return &A4q2g20195_eval;
case 20270 : return &A4q2g20270_eval;
case 20280 : return &A4q2g20280_eval;
case 20285 : return &A4q2g20285_eval;
case 20300 : return &A4q2g20300_eval;
case 20305 : return &A4q2g20305_eval;
case 20310 : return &A4q2g20310_eval;
case 20315 : return &A4q2g20315_eval;
case 20335 : return &A4q2g20335_eval;
case 20340 : return &A4q2g20340_eval;
case 20345 : return &A4q2g20345_eval;
case 20370 : return &A4q2g20370_eval;
case 20375 : return &A4q2g20375_eval;
case 20485 : return &A4q2g20485_eval;
case 20490 : return &A4q2g20490_eval;
case 20495 : return &A4q2g20495_eval;
case 20515 : return &A4q2g20515_eval;
case 20530 : return &A4q2g20530_eval;
case 20535 : return &A4q2g20535_eval;
case 20540 : return &A4q2g20540_eval;
case 20545 : return &A4q2g20545_eval;
case 20560 : return &A4q2g20560_eval;
case 20580 : return &A4q2g20580_eval;
case 20585 : return &A4q2g20585_eval;
case 20590 : return &A4q2g20590_eval;
case 20595 : return &A4q2g20595_eval;
case 20610 : return &A4q2g20610_eval;
case 20615 : return &A4q2g20615_eval;
case 20625 : return &A4q2g20625_eval;
case 20630 : return &A4q2g20630_eval;
case 20640 : return &A4q2g20640_eval;
case 20645 : return &A4q2g20645_eval;
case 20660 : return &A4q2g20660_eval;
case 20665 : return &A4q2g20665_eval;
case 20670 : return &A4q2g20670_eval;
case 20675 : return &A4q2g20675_eval;
case 20695 : return &A4q2g20695_eval;
case 20710 : return &A4q2g20710_eval;
case 20715 : return &A4q2g20715_eval;
case 20720 : return &A4q2g20720_eval;
case 20725 : return &A4q2g20725_eval;
case 20745 : return &A4q2g20745_eval;
case 20755 : return &A4q2g20755_eval;
case 20775 : return &A4q2g20775_eval;
case 20790 : return &A4q2g20790_eval;
case 20795 : return &A4q2g20795_eval;
case 20805 : return &A4q2g20805_eval;
case 20845 : return &A4q2g20845_eval;
case 20850 : return &A4q2g20850_eval;
case 20855 : return &A4q2g20855_eval;
case 20875 : return &A4q2g20875_eval;
case 20925 : return &A4q2g20925_eval;
case 20935 : return &A4q2g20935_eval;
case 20955 : return &A4q2g20955_eval;
case 20970 : return &A4q2g20970_eval;
case 20975 : return &A4q2g20975_eval;
case 20985 : return &A4q2g20985_eval;
case 21060 : return &A4q2g21060_eval;
case 21065 : return &A4q2g21065_eval;
case 21090 : return &A4q2g21090_eval;
case 21095 : return &A4q2g21095_eval;
case 21135 : return &A4q2g21135_eval;
case 21150 : return &A4q2g21150_eval;
case 21155 : return &A4q2g21155_eval;
case 21165 : return &A4q2g21165_eval;
case 21385 : return &A4q2g21385_eval;
case 21390 : return &A4q2g21390_eval;
case 21395 : return &A4q2g21395_eval;
case 21415 : return &A4q2g21415_eval;
case 21420 : return &A4q2g21420_eval;
case 21425 : return &A4q2g21425_eval;
case 21450 : return &A4q2g21450_eval;
case 21455 : return &A4q2g21455_eval;
case 21565 : return &A4q2g21565_eval;
case 21570 : return &A4q2g21570_eval;
case 21575 : return &A4q2g21575_eval;
case 21595 : return &A4q2g21595_eval;
case 21825 : return &A4q2g21825_eval;
case 21835 : return &A4q2g21835_eval;
case 21855 : return &A4q2g21855_eval;
case 21870 : return &A4q2g21870_eval;
case 21875 : return &A4q2g21875_eval;
case 21885 : return &A4q2g21885_eval;
case 21925 : return &A4q2g21925_eval;
case 21930 : return &A4q2g21930_eval;
case 21935 : return &A4q2g21935_eval;
case 21955 : return &A4q2g21955_eval;
case 22005 : return &A4q2g22005_eval;
case 22015 : return &A4q2g22015_eval;
case 22090 : return &A4q2g22090_eval;
case 22095 : return &A4q2g22095_eval;
case 22125 : return &A4q2g22125_eval;
case 22150 : return &A4q2g22150_eval;
case 22155 : return &A4q2g22155_eval;
case 22160 : return &A4q2g22160_eval;
case 22165 : return &A4q2g22165_eval;
case 22185 : return &A4q2g22185_eval;
case 22195 : return &A4q2g22195_eval;
case 22270 : return &A4q2g22270_eval;
case 22275 : return &A4q2g22275_eval;
case 22360 : return &A4q2g22360_eval;
case 22380 : return &A4q2g22380_eval;
case 22385 : return &A4q2g22385_eval;
case 22390 : return &A4q2g22390_eval;
case 22395 : return &A4q2g22395_eval;
case 22410 : return &A4q2g22410_eval;
case 22415 : return &A4q2g22415_eval;
case 22425 : return &A4q2g22425_eval;
case 22450 : return &A4q2g22450_eval;
case 22455 : return &A4q2g22455_eval;
case 22485 : return &A4q2g22485_eval;
case 22575 : return &A4q2g22575_eval;
case 22590 : return &A4q2g22590_eval;
case 22595 : return &A4q2g22595_eval;
case 22605 : return &A4q2g22605_eval;
case 22665 : return &A4q2g22665_eval;
case 22690 : return &A4q2g22690_eval;
case 22695 : return &A4q2g22695_eval;
case 22700 : return &A4q2g22700_eval;
case 22705 : return &A4q2g22705_eval;
case 22720 : return &A4q2g22720_eval;
case 22740 : return &A4q2g22740_eval;
case 22745 : return &A4q2g22745_eval;
case 22750 : return &A4q2g22750_eval;
case 22755 : return &A4q2g22755_eval;
case 22770 : return &A4q2g22770_eval;
case 22775 : return &A4q2g22775_eval;
case 22785 : return &A4q2g22785_eval;
case 22790 : return &A4q2g22790_eval;
case 22800 : return &A4q2g22800_eval;
case 22805 : return &A4q2g22805_eval;
case 22820 : return &A4q2g22820_eval;
case 22825 : return &A4q2g22825_eval;
case 22830 : return &A4q2g22830_eval;
case 22835 : return &A4q2g22835_eval;
case 22855 : return &A4q2g22855_eval;
case 22870 : return &A4q2g22870_eval;
case 22875 : return &A4q2g22875_eval;
case 22880 : return &A4q2g22880_eval;
case 22885 : return &A4q2g22885_eval;
case 22905 : return &A4q2g22905_eval;
case 22915 : return &A4q2g22915_eval;
case 22935 : return &A4q2g22935_eval;
case 22950 : return &A4q2g22950_eval;
case 22955 : return &A4q2g22955_eval;
case 22965 : return &A4q2g22965_eval;
case 23005 : return &A4q2g23005_eval;
case 23010 : return &A4q2g23010_eval;
case 23015 : return &A4q2g23015_eval;
case 23035 : return &A4q2g23035_eval;
case 23085 : return &A4q2g23085_eval;
case 23095 : return &A4q2g23095_eval;
case 23170 : return &A4q2g23170_eval;
case 23175 : return &A4q2g23175_eval;
case 23205 : return &A4q2g23205_eval;
case 23230 : return &A4q2g23230_eval;
case 23235 : return &A4q2g23235_eval;
case 23240 : return &A4q2g23240_eval;
case 23245 : return &A4q2g23245_eval;
case 23265 : return &A4q2g23265_eval;
case 23275 : return &A4q2g23275_eval;
case 23380 : return &A4q2g23380_eval;
case 23390 : return &A4q2g23390_eval;
case 23410 : return &A4q2g23410_eval;
case 23415 : return &A4q2g23415_eval;
case 23420 : return &A4q2g23420_eval;
case 23425 : return &A4q2g23425_eval;
case 23450 : return &A4q2g23450_eval;
case 23480 : return &A4q2g23480_eval;
case 23485 : return &A4q2g23485_eval;
case 23560 : return &A4q2g23560_eval;
case 23570 : return &A4q2g23570_eval;
case 23620 : return &A4q2g23620_eval;
case 23640 : return &A4q2g23640_eval;
case 23645 : return &A4q2g23645_eval;
case 23650 : return &A4q2g23650_eval;
case 23690 : return &A4q2g23690_eval;
case 23700 : return &A4q2g23700_eval;
case 23705 : return &A4q2g23705_eval;
case 23720 : return &A4q2g23720_eval;
case 23740 : return &A4q2g23740_eval;
case 23750 : return &A4q2g23750_eval;
case 23770 : return &A4q2g23770_eval;
case 23775 : return &A4q2g23775_eval;
case 23780 : return &A4q2g23780_eval;
case 23785 : return &A4q2g23785_eval;
case 23800 : return &A4q2g23800_eval;
case 23820 : return &A4q2g23820_eval;
case 23825 : return &A4q2g23825_eval;
case 23830 : return &A4q2g23830_eval;
case 23835 : return &A4q2g23835_eval;
case 23850 : return &A4q2g23850_eval;
case 23855 : return &A4q2g23855_eval;
case 23865 : return &A4q2g23865_eval;
case 23870 : return &A4q2g23870_eval;
case 23880 : return &A4q2g23880_eval;
case 23885 : return &A4q2g23885_eval;
case 23900 : return &A4q2g23900_eval;
case 23905 : return &A4q2g23905_eval;
case 23910 : return &A4q2g23910_eval;
case 23915 : return &A4q2g23915_eval;
case 23935 : return &A4q2g23935_eval;
case 23950 : return &A4q2g23950_eval;
case 23955 : return &A4q2g23955_eval;
case 23960 : return &A4q2g23960_eval;
case 23965 : return &A4q2g23965_eval;
case 23990 : return &A4q2g23990_eval;
case 24050 : return &A4q2g24050_eval;
case 24060 : return &A4q2g24060_eval;
case 24065 : return &A4q2g24065_eval;
case 24080 : return &A4q2g24080_eval;
case 24170 : return &A4q2g24170_eval;
case 24200 : return &A4q2g24200_eval;
case 24205 : return &A4q2g24205_eval;
case 24230 : return &A4q2g24230_eval;
case 24240 : return &A4q2g24240_eval;
case 24245 : return &A4q2g24245_eval;
case 24260 : return &A4q2g24260_eval;
case 24265 : return &A4q2g24265_eval;
case 24270 : return &A4q2g24270_eval;
case 24275 : return &A4q2g24275_eval;
case 24295 : return &A4q2g24295_eval;
case 24380 : return &A4q2g24380_eval;
case 24385 : return &A4q2g24385_eval;
case 24460 : return &A4q2g24460_eval;
case 24470 : return &A4q2g24470_eval;
case 24490 : return &A4q2g24490_eval;
case 24495 : return &A4q2g24495_eval;
case 24500 : return &A4q2g24500_eval;
case 24505 : return &A4q2g24505_eval;
case 24530 : return &A4q2g24530_eval;
case 24560 : return &A4q2g24560_eval;
case 24565 : return &A4q2g24565_eval;
case 24640 : return &A4q2g24640_eval;
case 24650 : return &A4q2g24650_eval;
case 24700 : return &A4q2g24700_eval;
case 24720 : return &A4q2g24720_eval;
case 24725 : return &A4q2g24725_eval;
case 24730 : return &A4q2g24730_eval;
case 24770 : return &A4q2g24770_eval;
case 24780 : return &A4q2g24780_eval;
case 24785 : return &A4q2g24785_eval;
case 24800 : return &A4q2g24800_eval;
case 24820 : return &A4q2g24820_eval;
case 24830 : return &A4q2g24830_eval;
case 25060 : return &A4q2g25060_eval;
case 25080 : return &A4q2g25080_eval;
case 25085 : return &A4q2g25085_eval;
case 25090 : return &A4q2g25090_eval;
case 25200 : return &A4q2g25200_eval;
case 25205 : return &A4q2g25205_eval;
case 25230 : return &A4q2g25230_eval;
case 25235 : return &A4q2g25235_eval;
case 25240 : return &A4q2g25240_eval;
case 25260 : return &A4q2g25260_eval;
case 25265 : return &A4q2g25265_eval;
case 25270 : return &A4q2g25270_eval;
case 25490 : return &A4q2g25490_eval;
case 25500 : return &A4q2g25500_eval;
case 25505 : return &A4q2g25505_eval;
case 25520 : return &A4q2g25520_eval;
case 25560 : return &A4q2g25560_eval;
case 25565 : return &A4q2g25565_eval;
case 25590 : return &A4q2g25590_eval;
case 25595 : return &A4q2g25595_eval;
case 25670 : return &A4q2g25670_eval;
case 25680 : return &A4q2g25680_eval;
case 25685 : return &A4q2g25685_eval;
case 25700 : return &A4q2g25700_eval;
case 25720 : return &A4q2g25720_eval;
case 25730 : return &A4q2g25730_eval;
case 25780 : return &A4q2g25780_eval;
case 25800 : return &A4q2g25800_eval;
case 25805 : return &A4q2g25805_eval;
case 25810 : return &A4q2g25810_eval;
case 25850 : return &A4q2g25850_eval;
case 25860 : return &A4q2g25860_eval;
case 25865 : return &A4q2g25865_eval;
case 25880 : return &A4q2g25880_eval;
case 25900 : return &A4q2g25900_eval;
case 25910 : return &A4q2g25910_eval;
case 25930 : return &A4q2g25930_eval;
case 25935 : return &A4q2g25935_eval;
case 25940 : return &A4q2g25940_eval;
case 25945 : return &A4q2g25945_eval;
case 25960 : return &A4q2g25960_eval;
case 25980 : return &A4q2g25980_eval;
case 25985 : return &A4q2g25985_eval;
case 25990 : return &A4q2g25990_eval;
case 25995 : return &A4q2g25995_eval;
case 26010 : return &A4q2g26010_eval;
case 26015 : return &A4q2g26015_eval;
case 26025 : return &A4q2g26025_eval;
case 26030 : return &A4q2g26030_eval;
case 26040 : return &A4q2g26040_eval;
case 26045 : return &A4q2g26045_eval;
case 26060 : return &A4q2g26060_eval;
case 26065 : return &A4q2g26065_eval;
case 26070 : return &A4q2g26070_eval;
case 26075 : return &A4q2g26075_eval;
case 26095 : return &A4q2g26095_eval;
case 26110 : return &A4q2g26110_eval;
case 26115 : return &A4q2g26115_eval;
case 26120 : return &A4q2g26120_eval;
case 26125 : return &A4q2g26125_eval;
case 26140 : return &A4q2g26140_eval;
case 26160 : return &A4q2g26160_eval;
case 26165 : return &A4q2g26165_eval;
case 26170 : return &A4q2g26170_eval;
case 26280 : return &A4q2g26280_eval;
case 26285 : return &A4q2g26285_eval;
case 26310 : return &A4q2g26310_eval;
case 26315 : return &A4q2g26315_eval;
case 26320 : return &A4q2g26320_eval;
case 26340 : return &A4q2g26340_eval;
case 26345 : return &A4q2g26345_eval;
case 26350 : return &A4q2g26350_eval;
case 26355 : return &A4q2g26355_eval;
case 26370 : return &A4q2g26370_eval;
case 26375 : return &A4q2g26375_eval;
case 26385 : return &A4q2g26385_eval;
case 26460 : return &A4q2g26460_eval;
case 26465 : return &A4q2g26465_eval;
case 26490 : return &A4q2g26490_eval;
case 26495 : return &A4q2g26495_eval;
case 26535 : return &A4q2g26535_eval;
case 26550 : return &A4q2g26550_eval;
case 26555 : return &A4q2g26555_eval;
case 26565 : return &A4q2g26565_eval;
case 26570 : return &A4q2g26570_eval;
case 26580 : return &A4q2g26580_eval;
case 26585 : return &A4q2g26585_eval;
case 26600 : return &A4q2g26600_eval;
case 26640 : return &A4q2g26640_eval;
case 26645 : return &A4q2g26645_eval;
case 26670 : return &A4q2g26670_eval;
case 26675 : return &A4q2g26675_eval;
case 26750 : return &A4q2g26750_eval;
case 26760 : return &A4q2g26760_eval;
case 26765 : return &A4q2g26765_eval;
case 26780 : return &A4q2g26780_eval;
case 26785 : return &A4q2g26785_eval;
case 26790 : return &A4q2g26790_eval;
case 26795 : return &A4q2g26795_eval;
case 26815 : return &A4q2g26815_eval;
case 26820 : return &A4q2g26820_eval;
case 26825 : return &A4q2g26825_eval;
case 26850 : return &A4q2g26850_eval;
case 26855 : return &A4q2g26855_eval;
case 26965 : return &A4q2g26965_eval;
case 26970 : return &A4q2g26970_eval;
case 26975 : return &A4q2g26975_eval;
case 26995 : return &A4q2g26995_eval;
case 27010 : return &A4q2g27010_eval;
case 27015 : return &A4q2g27015_eval;
case 27020 : return &A4q2g27020_eval;
case 27025 : return &A4q2g27025_eval;
case 27040 : return &A4q2g27040_eval;
case 27060 : return &A4q2g27060_eval;
case 27065 : return &A4q2g27065_eval;
case 27070 : return &A4q2g27070_eval;
case 27075 : return &A4q2g27075_eval;
case 27090 : return &A4q2g27090_eval;
case 27095 : return &A4q2g27095_eval;
case 27105 : return &A4q2g27105_eval;
case 27110 : return &A4q2g27110_eval;
case 27120 : return &A4q2g27120_eval;
case 27125 : return &A4q2g27125_eval;
case 27140 : return &A4q2g27140_eval;
case 27145 : return &A4q2g27145_eval;
case 27150 : return &A4q2g27150_eval;
case 27155 : return &A4q2g27155_eval;
case 27175 : return &A4q2g27175_eval;
case 27190 : return &A4q2g27190_eval;
case 27195 : return &A4q2g27195_eval;
case 27200 : return &A4q2g27200_eval;
case 27205 : return &A4q2g27205_eval;
case 27230 : return &A4q2g27230_eval;
case 27290 : return &A4q2g27290_eval;
case 27300 : return &A4q2g27300_eval;
case 27305 : return &A4q2g27305_eval;
case 27320 : return &A4q2g27320_eval;
case 27410 : return &A4q2g27410_eval;
case 27650 : return &A4q2g27650_eval;
case 27660 : return &A4q2g27660_eval;
case 27665 : return &A4q2g27665_eval;
case 27680 : return &A4q2g27680_eval;
case 27720 : return &A4q2g27720_eval;
case 27725 : return &A4q2g27725_eval;
case 27750 : return &A4q2g27750_eval;
case 27755 : return &A4q2g27755_eval;
case 27830 : return &A4q2g27830_eval;
case 27840 : return &A4q2g27840_eval;
case 27845 : return &A4q2g27845_eval;
case 27860 : return &A4q2g27860_eval;
case 28310 : return &A4q2g28310_eval;
case 28370 : return &A4q2g28370_eval;
case 28380 : return &A4q2g28380_eval;
case 28385 : return &A4q2g28385_eval;
case 28400 : return &A4q2g28400_eval;
case 28490 : return &A4q2g28490_eval;
case 28520 : return &A4q2g28520_eval;
case 28525 : return &A4q2g28525_eval;
case 28550 : return &A4q2g28550_eval;
case 28560 : return &A4q2g28560_eval;
case 28565 : return &A4q2g28565_eval;
case 28580 : return &A4q2g28580_eval;
case 28585 : return &A4q2g28585_eval;
case 28590 : return &A4q2g28590_eval;
case 28595 : return &A4q2g28595_eval;
case 28615 : return &A4q2g28615_eval;
case 28700 : return &A4q2g28700_eval;
case 28705 : return &A4q2g28705_eval;
case 28730 : return &A4q2g28730_eval;
case 28740 : return &A4q2g28740_eval;
case 28745 : return &A4q2g28745_eval;
case 28760 : return &A4q2g28760_eval;
case 28800 : return &A4q2g28800_eval;
case 28805 : return &A4q2g28805_eval;
case 28830 : return &A4q2g28830_eval;
case 28835 : return &A4q2g28835_eval;
case 28910 : return &A4q2g28910_eval;
case 28920 : return &A4q2g28920_eval;
case 28925 : return &A4q2g28925_eval;
case 28940 : return &A4q2g28940_eval;
case 28945 : return &A4q2g28945_eval;
case 28950 : return &A4q2g28950_eval;
case 28955 : return &A4q2g28955_eval;
case 28975 : return &A4q2g28975_eval;
case 28980 : return &A4q2g28980_eval;
case 28985 : return &A4q2g28985_eval;
case 29010 : return &A4q2g29010_eval;
case 29015 : return &A4q2g29015_eval;
case 29125 : return &A4q2g29125_eval;
case 29130 : return &A4q2g29130_eval;
case 29135 : return &A4q2g29135_eval;
case 29155 : return &A4q2g29155_eval;
case 29600 : return &A4q2g29600_eval;
case 29605 : return &A4q2g29605_eval;
case 29630 : return &A4q2g29630_eval;
case 29640 : return &A4q2g29640_eval;
case 29645 : return &A4q2g29645_eval;
case 29660 : return &A4q2g29660_eval;
case 29665 : return &A4q2g29665_eval;
case 29670 : return &A4q2g29670_eval;
case 29675 : return &A4q2g29675_eval;
case 29695 : return &A4q2g29695_eval;
case 29780 : return &A4q2g29780_eval;
case 29785 : return &A4q2g29785_eval;
case 29860 : return &A4q2g29860_eval;
case 29870 : return &A4q2g29870_eval;
case 29890 : return &A4q2g29890_eval;
case 29895 : return &A4q2g29895_eval;
case 29900 : return &A4q2g29900_eval;
case 29905 : return &A4q2g29905_eval;
case 29930 : return &A4q2g29930_eval;
case 29960 : return &A4q2g29960_eval;
case 29965 : return &A4q2g29965_eval;
case 30040 : return &A4q2g30040_eval;
case 30050 : return &A4q2g30050_eval;
case 30100 : return &A4q2g30100_eval;
case 30120 : return &A4q2g30120_eval;
case 30125 : return &A4q2g30125_eval;
case 30130 : return &A4q2g30130_eval;
case 30170 : return &A4q2g30170_eval;
case 30180 : return &A4q2g30180_eval;
case 30185 : return &A4q2g30185_eval;
case 30200 : return &A4q2g30200_eval;
case 30220 : return &A4q2g30220_eval;
case 30230 : return &A4q2g30230_eval;
case 30250 : return &A4q2g30250_eval;
case 30255 : return &A4q2g30255_eval;
case 30260 : return &A4q2g30260_eval;
case 30265 : return &A4q2g30265_eval;
case 30280 : return &A4q2g30280_eval;
case 30300 : return &A4q2g30300_eval;
case 30305 : return &A4q2g30305_eval;
case 30310 : return &A4q2g30310_eval;
case 30315 : return &A4q2g30315_eval;
case 30330 : return &A4q2g30330_eval;
case 30335 : return &A4q2g30335_eval;
case 30345 : return &A4q2g30345_eval;
case 30350 : return &A4q2g30350_eval;
case 30360 : return &A4q2g30360_eval;
case 30365 : return &A4q2g30365_eval;
case 30380 : return &A4q2g30380_eval;
case 30385 : return &A4q2g30385_eval;
case 30390 : return &A4q2g30390_eval;
case 30395 : return &A4q2g30395_eval;
case 30415 : return &A4q2g30415_eval;
case 30430 : return &A4q2g30430_eval;
case 30435 : return &A4q2g30435_eval;
case 30440 : return &A4q2g30440_eval;
case 30445 : return &A4q2g30445_eval;
case 30470 : return &A4q2g30470_eval;
case 30530 : return &A4q2g30530_eval;
case 30540 : return &A4q2g30540_eval;
case 30545 : return &A4q2g30545_eval;
case 30560 : return &A4q2g30560_eval;
case 30650 : return &A4q2g30650_eval;
case 30680 : return &A4q2g30680_eval;
case 30685 : return &A4q2g30685_eval;
case 30710 : return &A4q2g30710_eval;
case 30720 : return &A4q2g30720_eval;
case 30725 : return &A4q2g30725_eval;
case 30740 : return &A4q2g30740_eval;
case 30745 : return &A4q2g30745_eval;
case 30750 : return &A4q2g30750_eval;
case 30755 : return &A4q2g30755_eval;
case 30775 : return &A4q2g30775_eval;
case 30860 : return &A4q2g30860_eval;
case 30865 : return &A4q2g30865_eval;
case 30940 : return &A4q2g30940_eval;
case 30950 : return &A4q2g30950_eval;
case 30970 : return &A4q2g30970_eval;
case 30975 : return &A4q2g30975_eval;
case 30980 : return &A4q2g30980_eval;
case 30985 : return &A4q2g30985_eval;
case 31010 : return &A4q2g31010_eval;
case 31040 : return &A4q2g31040_eval;
case 31045 : return &A4q2g31045_eval;
case 31150 : return &A4q2g31150_eval;
case 31155 : return &A4q2g31155_eval;
case 31160 : return &A4q2g31160_eval;
case 31165 : return &A4q2g31165_eval;
case 31185 : return &A4q2g31185_eval;
case 31195 : return &A4q2g31195_eval;
case 31220 : return &A4q2g31220_eval;
case 31225 : return &A4q2g31225_eval;
case 31255 : return &A4q2g31255_eval;
case 31330 : return &A4q2g31330_eval;
case 31335 : return &A4q2g31335_eval;
case 31340 : return &A4q2g31340_eval;
case 31345 : return &A4q2g31345_eval;
case 31360 : return &A4q2g31360_eval;
case 31380 : return &A4q2g31380_eval;
case 31385 : return &A4q2g31385_eval;
case 31390 : return &A4q2g31390_eval;
case 31395 : return &A4q2g31395_eval;
case 31410 : return &A4q2g31410_eval;
case 31415 : return &A4q2g31415_eval;
case 31425 : return &A4q2g31425_eval;
case 31430 : return &A4q2g31430_eval;
case 31440 : return &A4q2g31440_eval;
case 31445 : return &A4q2g31445_eval;
case 31460 : return &A4q2g31460_eval;
case 31465 : return &A4q2g31465_eval;
case 31470 : return &A4q2g31470_eval;
case 31475 : return &A4q2g31475_eval;
case 31495 : return &A4q2g31495_eval;
case 31510 : return &A4q2g31510_eval;
case 31515 : return &A4q2g31515_eval;
case 31520 : return &A4q2g31520_eval;
case 31525 : return &A4q2g31525_eval;
case 31545 : return &A4q2g31545_eval;
case 31555 : return &A4q2g31555_eval;
case 31575 : return &A4q2g31575_eval;
case 31590 : return &A4q2g31590_eval;
case 31595 : return &A4q2g31595_eval;
case 31605 : return &A4q2g31605_eval;
case 31645 : return &A4q2g31645_eval;
case 31650 : return &A4q2g31650_eval;
case 31655 : return &A4q2g31655_eval;
case 31675 : return &A4q2g31675_eval;
case 31725 : return &A4q2g31725_eval;
case 31735 : return &A4q2g31735_eval;
case 31760 : return &A4q2g31760_eval;
case 31765 : return &A4q2g31765_eval;
case 31790 : return &A4q2g31790_eval;
case 31800 : return &A4q2g31800_eval;
case 31805 : return &A4q2g31805_eval;
case 31820 : return &A4q2g31820_eval;
case 31825 : return &A4q2g31825_eval;
case 31830 : return &A4q2g31830_eval;
case 31835 : return &A4q2g31835_eval;
case 31855 : return &A4q2g31855_eval;
case 31940 : return &A4q2g31940_eval;
case 31945 : return &A4q2g31945_eval;
case 31975 : return &A4q2g31975_eval;
case 32005 : return &A4q2g32005_eval;
case 32010 : return &A4q2g32010_eval;
case 32015 : return &A4q2g32015_eval;
case 32035 : return &A4q2g32035_eval;
case 32155 : return &A4q2g32155_eval;
case 32230 : return &A4q2g32230_eval;
case 32235 : return &A4q2g32235_eval;
case 32240 : return &A4q2g32240_eval;
case 32245 : return &A4q2g32245_eval;
case 32265 : return &A4q2g32265_eval;
case 32275 : return &A4q2g32275_eval;
case 32300 : return &A4q2g32300_eval;
case 32305 : return &A4q2g32305_eval;
case 32335 : return &A4q2g32335_eval;
case 32410 : return &A4q2g32410_eval;
case 32415 : return &A4q2g32415_eval;
case 32420 : return &A4q2g32420_eval;
case 32425 : return &A4q2g32425_eval;
case 32440 : return &A4q2g32440_eval;
case 32460 : return &A4q2g32460_eval;
case 32465 : return &A4q2g32465_eval;
case 32470 : return &A4q2g32470_eval;
case 32475 : return &A4q2g32475_eval;
case 32490 : return &A4q2g32490_eval;
case 32495 : return &A4q2g32495_eval;
case 32505 : return &A4q2g32505_eval;
case 32510 : return &A4q2g32510_eval;
case 32520 : return &A4q2g32520_eval;
case 32525 : return &A4q2g32525_eval;
case 32540 : return &A4q2g32540_eval;
case 32545 : return &A4q2g32545_eval;
case 32550 : return &A4q2g32550_eval;
case 32555 : return &A4q2g32555_eval;
case 32575 : return &A4q2g32575_eval;
case 32590 : return &A4q2g32590_eval;
case 32595 : return &A4q2g32595_eval;
case 32600 : return &A4q2g32600_eval;
case 32605 : return &A4q2g32605_eval;
case 32620 : return &A4q2g32620_eval;
case 32640 : return &A4q2g32640_eval;
case 32645 : return &A4q2g32645_eval;
case 32650 : return &A4q2g32650_eval;
case 32760 : return &A4q2g32760_eval;
case 32765 : return &A4q2g32765_eval;
case 32790 : return &A4q2g32790_eval;
case 32795 : return &A4q2g32795_eval;
case 32800 : return &A4q2g32800_eval;
case 32820 : return &A4q2g32820_eval;
case 32825 : return &A4q2g32825_eval;
case 32830 : return &A4q2g32830_eval;
case 32835 : return &A4q2g32835_eval;
case 32850 : return &A4q2g32850_eval;
case 32855 : return &A4q2g32855_eval;
case 32865 : return &A4q2g32865_eval;
case 32940 : return &A4q2g32940_eval;
case 32945 : return &A4q2g32945_eval;
case 32970 : return &A4q2g32970_eval;
case 32975 : return &A4q2g32975_eval;
case 33015 : return &A4q2g33015_eval;
case 33030 : return &A4q2g33030_eval;
case 33035 : return &A4q2g33035_eval;
case 33045 : return &A4q2g33045_eval;
case 33050 : return &A4q2g33050_eval;
case 33060 : return &A4q2g33060_eval;
case 33065 : return &A4q2g33065_eval;
case 33080 : return &A4q2g33080_eval;
case 33120 : return &A4q2g33120_eval;
case 33125 : return &A4q2g33125_eval;
case 33150 : return &A4q2g33150_eval;
case 33155 : return &A4q2g33155_eval;
case 33230 : return &A4q2g33230_eval;
case 33240 : return &A4q2g33240_eval;
case 33245 : return &A4q2g33245_eval;
case 33260 : return &A4q2g33260_eval;
case 33265 : return &A4q2g33265_eval;
case 33270 : return &A4q2g33270_eval;
case 33275 : return &A4q2g33275_eval;
case 33295 : return &A4q2g33295_eval;
case 33300 : return &A4q2g33300_eval;
case 33305 : return &A4q2g33305_eval;
case 33330 : return &A4q2g33330_eval;
case 33335 : return &A4q2g33335_eval;
case 33445 : return &A4q2g33445_eval;
case 33450 : return &A4q2g33450_eval;
case 33455 : return &A4q2g33455_eval;
case 33475 : return &A4q2g33475_eval;
case 33490 : return &A4q2g33490_eval;
case 33495 : return &A4q2g33495_eval;
case 33500 : return &A4q2g33500_eval;
case 33505 : return &A4q2g33505_eval;
case 33520 : return &A4q2g33520_eval;
case 33540 : return &A4q2g33540_eval;
case 33545 : return &A4q2g33545_eval;
case 33550 : return &A4q2g33550_eval;
case 33555 : return &A4q2g33555_eval;
case 33570 : return &A4q2g33570_eval;
case 33575 : return &A4q2g33575_eval;
case 33585 : return &A4q2g33585_eval;
case 33590 : return &A4q2g33590_eval;
case 33600 : return &A4q2g33600_eval;
case 33605 : return &A4q2g33605_eval;
case 33620 : return &A4q2g33620_eval;
case 33625 : return &A4q2g33625_eval;
case 33630 : return &A4q2g33630_eval;
case 33635 : return &A4q2g33635_eval;
case 33655 : return &A4q2g33655_eval;
case 33670 : return &A4q2g33670_eval;
case 33675 : return &A4q2g33675_eval;
case 33680 : return &A4q2g33680_eval;
case 33685 : return &A4q2g33685_eval;
case 33705 : return &A4q2g33705_eval;
case 33715 : return &A4q2g33715_eval;
case 33735 : return &A4q2g33735_eval;
case 33750 : return &A4q2g33750_eval;
case 33755 : return &A4q2g33755_eval;
case 33765 : return &A4q2g33765_eval;
case 33805 : return &A4q2g33805_eval;
case 33810 : return &A4q2g33810_eval;
case 33815 : return &A4q2g33815_eval;
case 33835 : return &A4q2g33835_eval;
case 33885 : return &A4q2g33885_eval;
case 33895 : return &A4q2g33895_eval;
case 33915 : return &A4q2g33915_eval;
case 33930 : return &A4q2g33930_eval;
case 33935 : return &A4q2g33935_eval;
case 33945 : return &A4q2g33945_eval;
case 34020 : return &A4q2g34020_eval;
case 34025 : return &A4q2g34025_eval;
case 34050 : return &A4q2g34050_eval;
case 34055 : return &A4q2g34055_eval;
case 34095 : return &A4q2g34095_eval;
case 34110 : return &A4q2g34110_eval;
case 34115 : return &A4q2g34115_eval;
case 34125 : return &A4q2g34125_eval;
case 34345 : return &A4q2g34345_eval;
case 34350 : return &A4q2g34350_eval;
case 34355 : return &A4q2g34355_eval;
case 34375 : return &A4q2g34375_eval;
case 34380 : return &A4q2g34380_eval;
case 34385 : return &A4q2g34385_eval;
case 34410 : return &A4q2g34410_eval;
case 34415 : return &A4q2g34415_eval;
case 34525 : return &A4q2g34525_eval;
case 34530 : return &A4q2g34530_eval;
case 34535 : return &A4q2g34535_eval;
case 34555 : return &A4q2g34555_eval;
case 34785 : return &A4q2g34785_eval;
case 34795 : return &A4q2g34795_eval;
case 34815 : return &A4q2g34815_eval;
case 34830 : return &A4q2g34830_eval;
case 34835 : return &A4q2g34835_eval;
case 34845 : return &A4q2g34845_eval;
case 34885 : return &A4q2g34885_eval;
case 34890 : return &A4q2g34890_eval;
case 34895 : return &A4q2g34895_eval;
case 34915 : return &A4q2g34915_eval;
case 34965 : return &A4q2g34965_eval;
case 34975 : return &A4q2g34975_eval;
case 35000 : return &A4q2g35000_eval;
case 35005 : return &A4q2g35005_eval;
case 35030 : return &A4q2g35030_eval;
case 35040 : return &A4q2g35040_eval;
case 35045 : return &A4q2g35045_eval;
case 35060 : return &A4q2g35060_eval;
case 35065 : return &A4q2g35065_eval;
case 35070 : return &A4q2g35070_eval;
case 35075 : return &A4q2g35075_eval;
case 35095 : return &A4q2g35095_eval;
case 35180 : return &A4q2g35180_eval;
case 35185 : return &A4q2g35185_eval;
case 35210 : return &A4q2g35210_eval;
case 35220 : return &A4q2g35220_eval;
case 35225 : return &A4q2g35225_eval;
case 35240 : return &A4q2g35240_eval;
case 35280 : return &A4q2g35280_eval;
case 35285 : return &A4q2g35285_eval;
case 35310 : return &A4q2g35310_eval;
case 35315 : return &A4q2g35315_eval;
case 35390 : return &A4q2g35390_eval;
case 35400 : return &A4q2g35400_eval;
case 35405 : return &A4q2g35405_eval;
case 35420 : return &A4q2g35420_eval;
case 35425 : return &A4q2g35425_eval;
case 35430 : return &A4q2g35430_eval;
case 35435 : return &A4q2g35435_eval;
case 35455 : return &A4q2g35455_eval;
case 35460 : return &A4q2g35460_eval;
case 35465 : return &A4q2g35465_eval;
case 35490 : return &A4q2g35490_eval;
case 35495 : return &A4q2g35495_eval;
case 35605 : return &A4q2g35605_eval;
case 35610 : return &A4q2g35610_eval;
case 35615 : return &A4q2g35615_eval;
case 35635 : return &A4q2g35635_eval;
case 36080 : return &A4q2g36080_eval;
case 36085 : return &A4q2g36085_eval;
case 36110 : return &A4q2g36110_eval;
case 36120 : return &A4q2g36120_eval;
case 36125 : return &A4q2g36125_eval;
case 36140 : return &A4q2g36140_eval;
case 36145 : return &A4q2g36145_eval;
case 36150 : return &A4q2g36150_eval;
case 36155 : return &A4q2g36155_eval;
case 36175 : return &A4q2g36175_eval;
case 36260 : return &A4q2g36260_eval;
case 36265 : return &A4q2g36265_eval;
case 36295 : return &A4q2g36295_eval;
case 36325 : return &A4q2g36325_eval;
case 36330 : return &A4q2g36330_eval;
case 36335 : return &A4q2g36335_eval;
case 36355 : return &A4q2g36355_eval;
case 36475 : return &A4q2g36475_eval;
case 36505 : return &A4q2g36505_eval;
case 36510 : return &A4q2g36510_eval;
case 36515 : return &A4q2g36515_eval;
case 36535 : return &A4q2g36535_eval;
case 36540 : return &A4q2g36540_eval;
case 36545 : return &A4q2g36545_eval;
case 36570 : return &A4q2g36570_eval;
case 36575 : return &A4q2g36575_eval;
case 36685 : return &A4q2g36685_eval;
case 36690 : return &A4q2g36690_eval;
case 36695 : return &A4q2g36695_eval;
case 36715 : return &A4q2g36715_eval;
case 37375 : return &A4q2g37375_eval;
case 37405 : return &A4q2g37405_eval;
case 37410 : return &A4q2g37410_eval;
case 37415 : return &A4q2g37415_eval;
case 37435 : return &A4q2g37435_eval;
case 37555 : return &A4q2g37555_eval;
case 37630 : return &A4q2g37630_eval;
case 37635 : return &A4q2g37635_eval;
case 37640 : return &A4q2g37640_eval;
case 37645 : return &A4q2g37645_eval;
case 37665 : return &A4q2g37665_eval;
case 37675 : return &A4q2g37675_eval;
case 37700 : return &A4q2g37700_eval;
case 37705 : return &A4q2g37705_eval;
case 37735 : return &A4q2g37735_eval;
case 37810 : return &A4q2g37810_eval;
case 37815 : return &A4q2g37815_eval;
case 37820 : return &A4q2g37820_eval;
case 37825 : return &A4q2g37825_eval;
case 37840 : return &A4q2g37840_eval;
case 37860 : return &A4q2g37860_eval;
case 37865 : return &A4q2g37865_eval;
case 37870 : return &A4q2g37870_eval;
case 37875 : return &A4q2g37875_eval;
case 37890 : return &A4q2g37890_eval;
case 37895 : return &A4q2g37895_eval;
case 37905 : return &A4q2g37905_eval;
case 37910 : return &A4q2g37910_eval;
case 37920 : return &A4q2g37920_eval;
case 37925 : return &A4q2g37925_eval;
case 37940 : return &A4q2g37940_eval;
case 37945 : return &A4q2g37945_eval;
case 37950 : return &A4q2g37950_eval;
case 37955 : return &A4q2g37955_eval;
case 37975 : return &A4q2g37975_eval;
case 37990 : return &A4q2g37990_eval;
case 37995 : return &A4q2g37995_eval;
case 38000 : return &A4q2g38000_eval;
case 38005 : return &A4q2g38005_eval;
case 38025 : return &A4q2g38025_eval;
case 38035 : return &A4q2g38035_eval;
case 38055 : return &A4q2g38055_eval;
case 38070 : return &A4q2g38070_eval;
case 38075 : return &A4q2g38075_eval;
case 38085 : return &A4q2g38085_eval;
case 38125 : return &A4q2g38125_eval;
case 38130 : return &A4q2g38130_eval;
case 38135 : return &A4q2g38135_eval;
case 38155 : return &A4q2g38155_eval;
case 38205 : return &A4q2g38205_eval;
case 38215 : return &A4q2g38215_eval;
case 38240 : return &A4q2g38240_eval;
case 38245 : return &A4q2g38245_eval;
case 38270 : return &A4q2g38270_eval;
case 38280 : return &A4q2g38280_eval;
case 38285 : return &A4q2g38285_eval;
case 38300 : return &A4q2g38300_eval;
case 38305 : return &A4q2g38305_eval;
case 38310 : return &A4q2g38310_eval;
case 38315 : return &A4q2g38315_eval;
case 38335 : return &A4q2g38335_eval;
case 38420 : return &A4q2g38420_eval;
case 38425 : return &A4q2g38425_eval;
case 38455 : return &A4q2g38455_eval;
case 38485 : return &A4q2g38485_eval;
case 38490 : return &A4q2g38490_eval;
case 38495 : return &A4q2g38495_eval;
case 38515 : return &A4q2g38515_eval;
case 38635 : return &A4q2g38635_eval;
case 38710 : return &A4q2g38710_eval;
case 38715 : return &A4q2g38715_eval;
case 38720 : return &A4q2g38720_eval;
case 38725 : return &A4q2g38725_eval;
case 38745 : return &A4q2g38745_eval;
case 38755 : return &A4q2g38755_eval;
case 38780 : return &A4q2g38780_eval;
case 38785 : return &A4q2g38785_eval;
case 38815 : return &A4q2g38815_eval;
case 39160 : return &A4q2g39160_eval;
case 39190 : return &A4q2g39190_eval;
case 39195 : return &A4q2g39195_eval;
case 39220 : return &A4q2g39220_eval;
case 39230 : return &A4q2g39230_eval;
case 39250 : return &A4q2g39250_eval;
case 39255 : return &A4q2g39255_eval;
case 39260 : return &A4q2g39260_eval;
case 39265 : return &A4q2g39265_eval;
case 39370 : return &A4q2g39370_eval;
case 39375 : return &A4q2g39375_eval;
case 39405 : return &A4q2g39405_eval;
case 39430 : return &A4q2g39430_eval;
case 39435 : return &A4q2g39435_eval;
case 39440 : return &A4q2g39440_eval;
case 39445 : return &A4q2g39445_eval;
case 39465 : return &A4q2g39465_eval;
case 39475 : return &A4q2g39475_eval;
case 39580 : return &A4q2g39580_eval;
case 39590 : return &A4q2g39590_eval;
case 39610 : return &A4q2g39610_eval;
case 39615 : return &A4q2g39615_eval;
case 39620 : return &A4q2g39620_eval;
case 39625 : return &A4q2g39625_eval;
case 39650 : return &A4q2g39650_eval;
case 39680 : return &A4q2g39680_eval;
case 39685 : return &A4q2g39685_eval;
case 39790 : return &A4q2g39790_eval;
case 39795 : return &A4q2g39795_eval;
case 39800 : return &A4q2g39800_eval;
case 39805 : return &A4q2g39805_eval;
case 39825 : return &A4q2g39825_eval;
case 39835 : return &A4q2g39835_eval;
case 39860 : return &A4q2g39860_eval;
case 39865 : return &A4q2g39865_eval;
case 39895 : return &A4q2g39895_eval;
case 40240 : return &A4q2g40240_eval;
case 40270 : return &A4q2g40270_eval;
case 40275 : return &A4q2g40275_eval;
case 40300 : return &A4q2g40300_eval;
case 40310 : return &A4q2g40310_eval;
case 40330 : return &A4q2g40330_eval;
case 40335 : return &A4q2g40335_eval;
case 40340 : return &A4q2g40340_eval;
case 40345 : return &A4q2g40345_eval;
case 40420 : return &A4q2g40420_eval;
case 40540 : return &A4q2g40540_eval;
case 40560 : return &A4q2g40560_eval;
case 40565 : return &A4q2g40565_eval;
case 40570 : return &A4q2g40570_eval;
case 40600 : return &A4q2g40600_eval;
case 40630 : return &A4q2g40630_eval;
case 40635 : return &A4q2g40635_eval;
case 40720 : return &A4q2g40720_eval;
case 40740 : return &A4q2g40740_eval;
case 40745 : return &A4q2g40745_eval;
case 40750 : return &A4q2g40750_eval;
case 40755 : return &A4q2g40755_eval;
case 40770 : return &A4q2g40770_eval;
case 40775 : return &A4q2g40775_eval;
case 40785 : return &A4q2g40785_eval;
case 40810 : return &A4q2g40810_eval;
case 40815 : return &A4q2g40815_eval;
case 40840 : return &A4q2g40840_eval;
case 40850 : return &A4q2g40850_eval;
case 40900 : return &A4q2g40900_eval;
case 40920 : return &A4q2g40920_eval;
case 40925 : return &A4q2g40925_eval;
case 40930 : return &A4q2g40930_eval;
case 40970 : return &A4q2g40970_eval;
case 40980 : return &A4q2g40980_eval;
case 40985 : return &A4q2g40985_eval;
case 41000 : return &A4q2g41000_eval;
case 41020 : return &A4q2g41020_eval;
case 41030 : return &A4q2g41030_eval;
case 41050 : return &A4q2g41050_eval;
case 41055 : return &A4q2g41055_eval;
case 41060 : return &A4q2g41060_eval;
case 41065 : return &A4q2g41065_eval;
case 41080 : return &A4q2g41080_eval;
case 41100 : return &A4q2g41100_eval;
case 41105 : return &A4q2g41105_eval;
case 41110 : return &A4q2g41110_eval;
case 41115 : return &A4q2g41115_eval;
case 41130 : return &A4q2g41130_eval;
case 41135 : return &A4q2g41135_eval;
case 41145 : return &A4q2g41145_eval;
case 41150 : return &A4q2g41150_eval;
case 41160 : return &A4q2g41160_eval;
case 41165 : return &A4q2g41165_eval;
case 41180 : return &A4q2g41180_eval;
case 41185 : return &A4q2g41185_eval;
case 41190 : return &A4q2g41190_eval;
case 41195 : return &A4q2g41195_eval;
case 41215 : return &A4q2g41215_eval;
case 41230 : return &A4q2g41230_eval;
case 41235 : return &A4q2g41235_eval;
case 41240 : return &A4q2g41240_eval;
case 41245 : return &A4q2g41245_eval;
case 41320 : return &A4q2g41320_eval;
case 41350 : return &A4q2g41350_eval;
case 41355 : return &A4q2g41355_eval;
case 41380 : return &A4q2g41380_eval;
case 41390 : return &A4q2g41390_eval;
case 41410 : return &A4q2g41410_eval;
case 41415 : return &A4q2g41415_eval;
case 41420 : return &A4q2g41420_eval;
case 41425 : return &A4q2g41425_eval;
case 41530 : return &A4q2g41530_eval;
case 41535 : return &A4q2g41535_eval;
case 41565 : return &A4q2g41565_eval;
case 41590 : return &A4q2g41590_eval;
case 41595 : return &A4q2g41595_eval;
case 41600 : return &A4q2g41600_eval;
case 41605 : return &A4q2g41605_eval;
case 41625 : return &A4q2g41625_eval;
case 41635 : return &A4q2g41635_eval;
case 41710 : return &A4q2g41710_eval;
case 41715 : return &A4q2g41715_eval;
case 41800 : return &A4q2g41800_eval;
case 41820 : return &A4q2g41820_eval;
case 41825 : return &A4q2g41825_eval;
case 41830 : return &A4q2g41830_eval;
case 41835 : return &A4q2g41835_eval;
case 41850 : return &A4q2g41850_eval;
case 41855 : return &A4q2g41855_eval;
case 41865 : return &A4q2g41865_eval;
case 41890 : return &A4q2g41890_eval;
case 41895 : return &A4q2g41895_eval;
case 41925 : return &A4q2g41925_eval;
case 42015 : return &A4q2g42015_eval;
case 42030 : return &A4q2g42030_eval;
case 42035 : return &A4q2g42035_eval;
case 42045 : return &A4q2g42045_eval;
case 42105 : return &A4q2g42105_eval;
case 42130 : return &A4q2g42130_eval;
case 42135 : return &A4q2g42135_eval;
case 42140 : return &A4q2g42140_eval;
case 42145 : return &A4q2g42145_eval;
case 42160 : return &A4q2g42160_eval;
case 42180 : return &A4q2g42180_eval;
case 42185 : return &A4q2g42185_eval;
case 42190 : return &A4q2g42190_eval;
case 42195 : return &A4q2g42195_eval;
case 42210 : return &A4q2g42210_eval;
case 42215 : return &A4q2g42215_eval;
case 42225 : return &A4q2g42225_eval;
case 42230 : return &A4q2g42230_eval;
case 42240 : return &A4q2g42240_eval;
case 42245 : return &A4q2g42245_eval;
case 42260 : return &A4q2g42260_eval;
case 42265 : return &A4q2g42265_eval;
case 42270 : return &A4q2g42270_eval;
case 42275 : return &A4q2g42275_eval;
case 42295 : return &A4q2g42295_eval;
case 42310 : return &A4q2g42310_eval;
case 42315 : return &A4q2g42315_eval;
case 42320 : return &A4q2g42320_eval;
case 42325 : return &A4q2g42325_eval;
case 42345 : return &A4q2g42345_eval;
case 42355 : return &A4q2g42355_eval;
case 42375 : return &A4q2g42375_eval;
case 42390 : return &A4q2g42390_eval;
case 42395 : return &A4q2g42395_eval;
case 42405 : return &A4q2g42405_eval;
case 42445 : return &A4q2g42445_eval;
case 42450 : return &A4q2g42450_eval;
case 42455 : return &A4q2g42455_eval;
case 42475 : return &A4q2g42475_eval;
case 42525 : return &A4q2g42525_eval;
case 42535 : return &A4q2g42535_eval;
case 42610 : return &A4q2g42610_eval;
case 42615 : return &A4q2g42615_eval;
case 42645 : return &A4q2g42645_eval;
case 42670 : return &A4q2g42670_eval;
case 42675 : return &A4q2g42675_eval;
case 42680 : return &A4q2g42680_eval;
case 42685 : return &A4q2g42685_eval;
case 42705 : return &A4q2g42705_eval;
case 42715 : return &A4q2g42715_eval;
case 42820 : return &A4q2g42820_eval;
case 42830 : return &A4q2g42830_eval;
case 42850 : return &A4q2g42850_eval;
case 42855 : return &A4q2g42855_eval;
case 42860 : return &A4q2g42860_eval;
case 42865 : return &A4q2g42865_eval;
case 42890 : return &A4q2g42890_eval;
case 42920 : return &A4q2g42920_eval;
case 42925 : return &A4q2g42925_eval;
case 43000 : return &A4q2g43000_eval;
case 43010 : return &A4q2g43010_eval;
case 43060 : return &A4q2g43060_eval;
case 43080 : return &A4q2g43080_eval;
case 43085 : return &A4q2g43085_eval;
case 43090 : return &A4q2g43090_eval;
case 43130 : return &A4q2g43130_eval;
case 43140 : return &A4q2g43140_eval;
case 43145 : return &A4q2g43145_eval;
case 43160 : return &A4q2g43160_eval;
case 43180 : return &A4q2g43180_eval;
case 43190 : return &A4q2g43190_eval;
case 43210 : return &A4q2g43210_eval;
case 43215 : return &A4q2g43215_eval;
case 43220 : return &A4q2g43220_eval;
case 43225 : return &A4q2g43225_eval;
case 43240 : return &A4q2g43240_eval;
case 43260 : return &A4q2g43260_eval;
case 43265 : return &A4q2g43265_eval;
case 43270 : return &A4q2g43270_eval;
case 43275 : return &A4q2g43275_eval;
case 43290 : return &A4q2g43290_eval;
case 43295 : return &A4q2g43295_eval;
case 43305 : return &A4q2g43305_eval;
case 43310 : return &A4q2g43310_eval;
case 43320 : return &A4q2g43320_eval;
case 43325 : return &A4q2g43325_eval;
case 43340 : return &A4q2g43340_eval;
case 43345 : return &A4q2g43345_eval;
case 43350 : return &A4q2g43350_eval;
case 43355 : return &A4q2g43355_eval;
case 43375 : return &A4q2g43375_eval;
case 43390 : return &A4q2g43390_eval;
case 43395 : return &A4q2g43395_eval;
case 43400 : return &A4q2g43400_eval;
case 43405 : return &A4q2g43405_eval;
case 43430 : return &A4q2g43430_eval;
case 43490 : return &A4q2g43490_eval;
case 43500 : return &A4q2g43500_eval;
case 43505 : return &A4q2g43505_eval;
case 43520 : return &A4q2g43520_eval;
case 43610 : return &A4q2g43610_eval;
case 43640 : return &A4q2g43640_eval;
case 43645 : return &A4q2g43645_eval;
case 43670 : return &A4q2g43670_eval;
case 43680 : return &A4q2g43680_eval;
case 43685 : return &A4q2g43685_eval;
case 43700 : return &A4q2g43700_eval;
case 43705 : return &A4q2g43705_eval;
case 43710 : return &A4q2g43710_eval;
case 43715 : return &A4q2g43715_eval;
case 43735 : return &A4q2g43735_eval;
case 43820 : return &A4q2g43820_eval;
case 43825 : return &A4q2g43825_eval;
case 43900 : return &A4q2g43900_eval;
case 43910 : return &A4q2g43910_eval;
case 43930 : return &A4q2g43930_eval;
case 43935 : return &A4q2g43935_eval;
case 43940 : return &A4q2g43940_eval;
case 43945 : return &A4q2g43945_eval;
case 43970 : return &A4q2g43970_eval;
case 44000 : return &A4q2g44000_eval;
case 44005 : return &A4q2g44005_eval;
case 44110 : return &A4q2g44110_eval;
case 44115 : return &A4q2g44115_eval;
case 44120 : return &A4q2g44120_eval;
case 44125 : return &A4q2g44125_eval;
case 44145 : return &A4q2g44145_eval;
case 44155 : return &A4q2g44155_eval;
case 44180 : return &A4q2g44180_eval;
case 44185 : return &A4q2g44185_eval;
case 44215 : return &A4q2g44215_eval;
case 44290 : return &A4q2g44290_eval;
case 44295 : return &A4q2g44295_eval;
case 44300 : return &A4q2g44300_eval;
case 44305 : return &A4q2g44305_eval;
case 44320 : return &A4q2g44320_eval;
case 44340 : return &A4q2g44340_eval;
case 44345 : return &A4q2g44345_eval;
case 44350 : return &A4q2g44350_eval;
case 44355 : return &A4q2g44355_eval;
case 44370 : return &A4q2g44370_eval;
case 44375 : return &A4q2g44375_eval;
case 44385 : return &A4q2g44385_eval;
case 44390 : return &A4q2g44390_eval;
case 44400 : return &A4q2g44400_eval;
case 44405 : return &A4q2g44405_eval;
case 44420 : return &A4q2g44420_eval;
case 44425 : return &A4q2g44425_eval;
case 44430 : return &A4q2g44430_eval;
case 44435 : return &A4q2g44435_eval;
case 44455 : return &A4q2g44455_eval;
case 44470 : return &A4q2g44470_eval;
case 44475 : return &A4q2g44475_eval;
case 44480 : return &A4q2g44480_eval;
case 44485 : return &A4q2g44485_eval;
case 44505 : return &A4q2g44505_eval;
case 44515 : return &A4q2g44515_eval;
case 44535 : return &A4q2g44535_eval;
case 44550 : return &A4q2g44550_eval;
case 44555 : return &A4q2g44555_eval;
case 44565 : return &A4q2g44565_eval;
case 44605 : return &A4q2g44605_eval;
case 44610 : return &A4q2g44610_eval;
case 44615 : return &A4q2g44615_eval;
case 44635 : return &A4q2g44635_eval;
case 44685 : return &A4q2g44685_eval;
case 44695 : return &A4q2g44695_eval;
case 44720 : return &A4q2g44720_eval;
case 44725 : return &A4q2g44725_eval;
case 44750 : return &A4q2g44750_eval;
case 44760 : return &A4q2g44760_eval;
case 44765 : return &A4q2g44765_eval;
case 44780 : return &A4q2g44780_eval;
case 44785 : return &A4q2g44785_eval;
case 44790 : return &A4q2g44790_eval;
case 44795 : return &A4q2g44795_eval;
case 44815 : return &A4q2g44815_eval;
case 44900 : return &A4q2g44900_eval;
case 44905 : return &A4q2g44905_eval;
case 44935 : return &A4q2g44935_eval;
case 44965 : return &A4q2g44965_eval;
case 44970 : return &A4q2g44970_eval;
case 44975 : return &A4q2g44975_eval;
case 44995 : return &A4q2g44995_eval;
case 45115 : return &A4q2g45115_eval;
case 45190 : return &A4q2g45190_eval;
case 45195 : return &A4q2g45195_eval;
case 45200 : return &A4q2g45200_eval;
case 45205 : return &A4q2g45205_eval;
case 45225 : return &A4q2g45225_eval;
case 45235 : return &A4q2g45235_eval;
case 45260 : return &A4q2g45260_eval;
case 45265 : return &A4q2g45265_eval;
case 45295 : return &A4q2g45295_eval;
case 45640 : return &A4q2g45640_eval;
case 45670 : return &A4q2g45670_eval;
case 45675 : return &A4q2g45675_eval;
case 45700 : return &A4q2g45700_eval;
case 45710 : return &A4q2g45710_eval;
case 45730 : return &A4q2g45730_eval;
case 45735 : return &A4q2g45735_eval;
case 45740 : return &A4q2g45740_eval;
case 45745 : return &A4q2g45745_eval;
case 45850 : return &A4q2g45850_eval;
case 45855 : return &A4q2g45855_eval;
case 45885 : return &A4q2g45885_eval;
case 45910 : return &A4q2g45910_eval;
case 45915 : return &A4q2g45915_eval;
case 45920 : return &A4q2g45920_eval;
case 45925 : return &A4q2g45925_eval;
case 45945 : return &A4q2g45945_eval;
case 45955 : return &A4q2g45955_eval;
case 46060 : return &A4q2g46060_eval;
case 46070 : return &A4q2g46070_eval;
case 46090 : return &A4q2g46090_eval;
case 46095 : return &A4q2g46095_eval;
case 46100 : return &A4q2g46100_eval;
case 46105 : return &A4q2g46105_eval;
case 46130 : return &A4q2g46130_eval;
case 46160 : return &A4q2g46160_eval;
case 46165 : return &A4q2g46165_eval;
case 46270 : return &A4q2g46270_eval;
case 46275 : return &A4q2g46275_eval;
case 46280 : return &A4q2g46280_eval;
case 46285 : return &A4q2g46285_eval;
case 46305 : return &A4q2g46305_eval;
case 46315 : return &A4q2g46315_eval;
case 46340 : return &A4q2g46340_eval;
case 46345 : return &A4q2g46345_eval;
case 46375 : return &A4q2g46375_eval;
#endif
default: return 0;
		throw BHerror("case missing for tree amplitude!");

}
}



template complex<R>  (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A4q2g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
36540: m m qm qm qbp qbp
35460: m m qm qp qbm qbp
28980: m m qm qp qbp qbm
34380: m m qm qbm qp qbp
21420: m m qm qbm qbp qp
33300: m m qm qbp qm qbp
26820: m m qm qbp qp qbm
20340: m m qm qbp qbm qp
13860: m m qm qbp qbp qm
35280: m m qp qm qbm qbp
28800: m m qp qm qbp qbm
27720: m m qp qp qbm qbm
33120: m m qp qbm qm qbp
26640: m m qp qbm qp qbm
20160: m m qp qbm qbm qp
13680: m m qp qbm qbp qm
25560: m m qp qbp qm qbm
12600: m m qp qbp qbm qm
34020: m m qbm qm qp qbp
21060: m m qbm qm qbp qp
32940: m m qbm qp qm qbp
26460: m m qbm qp qp qbm
19980: m m qbm qp qbm qp
13500: m m qbm qp qbp qm
18900: m m qbm qbm qp qp
17820: m m qbm qbp qm qp
11340: m m qbm qbp qp qm
32760: m m qbp qm qm qbp
26280: m m qbp qm qp qbm
19800: m m qbp qm qbm qp
13320: m m qbp qm qbp qm
25200: m m qbp qp qm qbm
12240: m m qbp qp qbm qm
17640: m m qbp qbm qm qp
11160: m m qbp qbm qp qm
10080: m m qbp qbp qm qm
36510: m qm m qm qbp qbp
35430: m qm m qp qbm qbp
28950: m qm m qp qbp qbm
34350: m qm m qbm qp qbp
21390: m qm m qbm qbp qp
33270: m qm m qbp qm qbp
26790: m qm m qbp qp qbm
20310: m qm m qbp qbm qp
13830: m qm m qbp qbp qm
36330: m qm qm m qbp qbp
32010: m qm qm qbp m qbp
6090: m qm qm qbp qbp m
44970: m qm qm qbp qbp p
38490: m qm qm qbp p qbp
37410: m qm qm p qbp qbp
35070: m qm qp m qbm qbp
28590: m qm qp m qbp qbm
31830: m qm qp qbm m qbp
5910: m qm qp qbm qbp m
44790: m qm qp qbm qbp p
38310: m qm qp qbm p qbp
24270: m qm qp qbp m qbm
4830: m qm qp qbp qbm m
43710: m qm qp qbp qbm p
30750: m qm qp qbp p qbm
36150: m qm qp p qbm qbp
29670: m qm qp p qbp qbm
33810: m qm qbm m qp qbp
20850: m qm qbm m qbp qp
31650: m qm qbm qp m qbp
5730: m qm qbm qp qbp m
44610: m qm qbm qp qbp p
38130: m qm qbm qp p qbp
16530: m qm qbm qbp m qp
3570: m qm qbm qbp qp m
42450: m qm qbm qbp qp p
23010: m qm qbm qbp p qp
34890: m qm qbm p qp qbp
21930: m qm qbm p qbp qp
32550: m qm qbp m qm qbp
26070: m qm qbp m qp qbm
19590: m qm qbp m qbm qp
13110: m qm qbp m qbp qm
31470: m qm qbp qm m qbp
5550: m qm qbp qm qbp m
44430: m qm qbp qm qbp p
37950: m qm qbp qm p qbp
23910: m qm qbp qp m qbm
4470: m qm qbp qp qbm m
43350: m qm qbp qp qbm p
30390: m qm qbp qp p qbm
16350: m qm qbp qbm m qp
3390: m qm qbp qbm qp m
42270: m qm qbp qbm qp p
22830: m qm qbp qbm p qp
8790: m qm qbp qbp m qm
2310: m qm qbp qbp qm m
41190: m qm qbp qbp qm p
15270: m qm qbp qbp p qm
33630: m qm qbp p qm qbp
27150: m qm qbp p qp qbm
20670: m qm qbp p qbm qp
14190: m qm qbp p qbp qm
36690: m qm p qm qbp qbp
35610: m qm p qp qbm qbp
29130: m qm p qp qbp qbm
34530: m qm p qbm qp qbp
21570: m qm p qbm qbp qp
33450: m qm p qbp qm qbp
26970: m qm p qbp qp qbm
20490: m qm p qbp qbm qp
14010: m qm p qbp qbp qm
35220: m qp m qm qbm qbp
28740: m qp m qm qbp qbm
27660: m qp m qp qbm qbm
33060: m qp m qbm qm qbp
26580: m qp m qbm qp qbm
20100: m qp m qbm qbm qp
13620: m qp m qbm qbp qm
25500: m qp m qbp qm qbm
12540: m qp m qbp qbm qm
35040: m qp qm m qbm qbp
28560: m qp qm m qbp qbm
31800: m qp qm qbm m qbp
5880: m qp qm qbm qbp m
44760: m qp qm qbm qbp p
38280: m qp qm qbm p qbp
24240: m qp qm qbp m qbm
4800: m qp qm qbp qbm m
43680: m qp qm qbp qbm p
30720: m qp qm qbp p qbm
36120: m qp qm p qbm qbp
29640: m qp qm p qbp qbm
27300: m qp qp m qbm qbm
24060: m qp qp qbm m qbm
4620: m qp qp qbm qbm m
43500: m qp qp qbm qbm p
30540: m qp qp qbm p qbm
28380: m qp qp p qbm qbm
32520: m qp qbm m qm qbp
26040: m qp qbm m qp qbm
19560: m qp qbm m qbm qp
13080: m qp qbm m qbp qm
31440: m qp qbm qm m qbp
5520: m qp qbm qm qbp m
44400: m qp qbm qm qbp p
37920: m qp qbm qm p qbp
23880: m qp qbm qp m qbm
4440: m qp qbm qp qbm m
43320: m qp qbm qp qbm p
30360: m qp qbm qp p qbm
16320: m qp qbm qbm m qp
3360: m qp qbm qbm qp m
42240: m qp qbm qbm qp p
22800: m qp qbm qbm p qp
8760: m qp qbm qbp m qm
2280: m qp qbm qbp qm m
41160: m qp qbm qbp qm p
15240: m qp qbm qbp p qm
33600: m qp qbm p qm qbp
27120: m qp qbm p qp qbm
20640: m qp qbm p qbm qp
14160: m qp qbm p qbp qm
24780: m qp qbp m qm qbm
11820: m qp qbp m qbm qm
23700: m qp qbp qm m qbm
4260: m qp qbp qm qbm m
43140: m qp qbp qm qbm p
30180: m qp qbp qm p qbm
8580: m qp qbp qbm m qm
2100: m qp qbp qbm qm m
40980: m qp qbp qbm qm p
15060: m qp qbp qbm p qm
25860: m qp qbp p qm qbm
12900: m qp qbp p qbm qm
35400: m qp p qm qbm qbp
28920: m qp p qm qbp qbm
27840: m qp p qp qbm qbm
33240: m qp p qbm qm qbp
26760: m qp p qbm qp qbm
20280: m qp p qbm qbm qp
13800: m qp p qbm qbp qm
25680: m qp p qbp qm qbm
12720: m qp p qbp qbm qm
33930: m qbm m qm qp qbp
20970: m qbm m qm qbp qp
32850: m qbm m qp qm qbp
26370: m qbm m qp qp qbm
19890: m qbm m qp qbm qp
13410: m qbm m qp qbp qm
18810: m qbm m qbm qp qp
17730: m qbm m qbp qm qp
11250: m qbm m qbp qp qm
33750: m qbm qm m qp qbp
20790: m qbm qm m qbp qp
31590: m qbm qm qp m qbp
5670: m qbm qm qp qbp m
44550: m qbm qm qp qbp p
38070: m qbm qm qp p qbp
16470: m qbm qm qbp m qp
3510: m qbm qm qbp qp m
42390: m qbm qm qbp qp p
22950: m qbm qm qbp p qp
34830: m qbm qm p qp qbp
21870: m qbm qm p qbp qp
32490: m qbm qp m qm qbp
26010: m qbm qp m qp qbm
19530: m qbm qp m qbm qp
13050: m qbm qp m qbp qm
31410: m qbm qp qm m qbp
5490: m qbm qp qm qbp m
44370: m qbm qp qm qbp p
37890: m qbm qp qm p qbp
23850: m qbm qp qp m qbm
4410: m qbm qp qp qbm m
43290: m qbm qp qp qbm p
30330: m qbm qp qp p qbm
16290: m qbm qp qbm m qp
3330: m qbm qp qbm qp m
42210: m qbm qp qbm qp p
22770: m qbm qp qbm p qp
8730: m qbm qp qbp m qm
2250: m qbm qp qbp qm m
41130: m qbm qp qbp qm p
15210: m qbm qp qbp p qm
33570: m qbm qp p qm qbp
27090: m qbm qp p qp qbm
20610: m qbm qp p qbm qp
14130: m qbm qp p qbp qm
18270: m qbm qbm m qp qp
16110: m qbm qbm qp m qp
3150: m qbm qbm qp qp m
42030: m qbm qbm qp qp p
22590: m qbm qbm qp p qp
19350: m qbm qbm p qp qp
17010: m qbm qbp m qm qp
10530: m qbm qbp m qp qm
15930: m qbm qbp qm m qp
2970: m qbm qbp qm qp m
41850: m qbm qbp qm qp p
22410: m qbm qbp qm p qp
8370: m qbm qbp qp m qm
1890: m qbm qbp qp qm m
40770: m qbm qbp qp qm p
14850: m qbm qbp qp p qm
18090: m qbm qbp p qm qp
11610: m qbm qbp p qp qm
34110: m qbm p qm qp qbp
21150: m qbm p qm qbp qp
33030: m qbm p qp qm qbp
26550: m qbm p qp qp qbm
20070: m qbm p qp qbm qp
13590: m qbm p qp qbp qm
18990: m qbm p qbm qp qp
17910: m qbm p qbp qm qp
11430: m qbm p qbp qp qm
32640: m qbp m qm qm qbp
26160: m qbp m qm qp qbm
19680: m qbp m qm qbm qp
13200: m qbp m qm qbp qm
25080: m qbp m qp qm qbm
12120: m qbp m qp qbm qm
17520: m qbp m qbm qm qp
11040: m qbp m qbm qp qm
9960: m qbp m qbp qm qm
32460: m qbp qm m qm qbp
25980: m qbp qm m qp qbm
19500: m qbp qm m qbm qp
13020: m qbp qm m qbp qm
31380: m qbp qm qm m qbp
5460: m qbp qm qm qbp m
44340: m qbp qm qm qbp p
37860: m qbp qm qm p qbp
23820: m qbp qm qp m qbm
4380: m qbp qm qp qbm m
43260: m qbp qm qp qbm p
30300: m qbp qm qp p qbm
16260: m qbp qm qbm m qp
3300: m qbp qm qbm qp m
42180: m qbp qm qbm qp p
22740: m qbp qm qbm p qp
8700: m qbp qm qbp m qm
2220: m qbp qm qbp qm m
41100: m qbp qm qbp qm p
15180: m qbp qm qbp p qm
33540: m qbp qm p qm qbp
27060: m qbp qm p qp qbm
20580: m qbp qm p qbm qp
14100: m qbp qm p qbp qm
24720: m qbp qp m qm qbm
11760: m qbp qp m qbm qm
23640: m qbp qp qm m qbm
4200: m qbp qp qm qbm m
43080: m qbp qp qm qbm p
30120: m qbp qp qm p qbm
8520: m qbp qp qbm m qm
2040: m qbp qp qbm qm m
40920: m qbp qp qbm qm p
15000: m qbp qp qbm p qm
25800: m qbp qp p qm qbm
12840: m qbp qp p qbm qm
16980: m qbp qbm m qm qp
10500: m qbp qbm m qp qm
15900: m qbp qbm qm m qp
2940: m qbp qbm qm qp m
41820: m qbp qbm qm qp p
22380: m qbp qbm qm p qp
8340: m qbp qbm qp m qm
1860: m qbp qbm qp qm m
40740: m qbp qbm qp qm p
14820: m qbp qbm qp p qm
18060: m qbp qbm p qm qp
11580: m qbp qbm p qp qm
9240: m qbp qbp m qm qm
8160: m qbp qbp qm m qm
1680: m qbp qbp qm qm m
40560: m qbp qbp qm qm p
14640: m qbp qbp qm p qm
10320: m qbp qbp p qm qm
32820: m qbp p qm qm qbp
26340: m qbp p qm qp qbm
19860: m qbp p qm qbm qp
13380: m qbp p qm qbp qm
25260: m qbp p qp qm qbm
12300: m qbp p qp qbm qm
17700: m qbp p qbm qm qp
11220: m qbp p qbm qp qm
10140: m qbp p qbp qm qm
36570: m p qm qm qbp qbp
35490: m p qm qp qbm qbp
29010: m p qm qp qbp qbm
34410: m p qm qbm qp qbp
21450: m p qm qbm qbp qp
33330: m p qm qbp qm qbp
26850: m p qm qbp qp qbm
20370: m p qm qbp qbm qp
13890: m p qm qbp qbp qm
35310: m p qp qm qbm qbp
28830: m p qp qm qbp qbm
27750: m p qp qp qbm qbm
33150: m p qp qbm qm qbp
26670: m p qp qbm qp qbm
20190: m p qp qbm qbm qp
13710: m p qp qbm qbp qm
25590: m p qp qbp qm qbm
12630: m p qp qbp qbm qm
34050: m p qbm qm qp qbp
21090: m p qbm qm qbp qp
32970: m p qbm qp qm qbp
26490: m p qbm qp qp qbm
20010: m p qbm qp qbm qp
13530: m p qbm qp qbp qm
18930: m p qbm qbm qp qp
17850: m p qbm qbp qm qp
11370: m p qbm qbp qp qm
32790: m p qbp qm qm qbp
26310: m p qbp qm qp qbm
19830: m p qbp qm qbm qp
13350: m p qbp qm qbp qm
25230: m p qbp qp qm qbm
12270: m p qbp qp qbm qm
17670: m p qbp qbm qm qp
11190: m p qbp qbm qp qm
10110: m p qbp qbp qm qm
36505: qm m m qm qbp qbp
35425: qm m m qp qbm qbp
28945: qm m m qp qbp qbm
34345: qm m m qbm qp qbp
21385: qm m m qbm qbp qp
33265: qm m m qbp qm qbp
26785: qm m m qbp qp qbm
20305: qm m m qbp qbm qp
13825: qm m m qbp qbp qm
36325: qm m qm m qbp qbp
32005: qm m qm qbp m qbp
6085: qm m qm qbp qbp m
44965: qm m qm qbp qbp p
38485: qm m qm qbp p qbp
37405: qm m qm p qbp qbp
35065: qm m qp m qbm qbp
28585: qm m qp m qbp qbm
31825: qm m qp qbm m qbp
5905: qm m qp qbm qbp m
44785: qm m qp qbm qbp p
38305: qm m qp qbm p qbp
24265: qm m qp qbp m qbm
4825: qm m qp qbp qbm m
43705: qm m qp qbp qbm p
30745: qm m qp qbp p qbm
36145: qm m qp p qbm qbp
29665: qm m qp p qbp qbm
33805: qm m qbm m qp qbp
20845: qm m qbm m qbp qp
31645: qm m qbm qp m qbp
5725: qm m qbm qp qbp m
44605: qm m qbm qp qbp p
38125: qm m qbm qp p qbp
16525: qm m qbm qbp m qp
3565: qm m qbm qbp qp m
42445: qm m qbm qbp qp p
23005: qm m qbm qbp p qp
34885: qm m qbm p qp qbp
21925: qm m qbm p qbp qp
32545: qm m qbp m qm qbp
26065: qm m qbp m qp qbm
19585: qm m qbp m qbm qp
13105: qm m qbp m qbp qm
31465: qm m qbp qm m qbp
5545: qm m qbp qm qbp m
44425: qm m qbp qm qbp p
37945: qm m qbp qm p qbp
23905: qm m qbp qp m qbm
4465: qm m qbp qp qbm m
43345: qm m qbp qp qbm p
30385: qm m qbp qp p qbm
16345: qm m qbp qbm m qp
3385: qm m qbp qbm qp m
42265: qm m qbp qbm qp p
22825: qm m qbp qbm p qp
8785: qm m qbp qbp m qm
2305: qm m qbp qbp qm m
41185: qm m qbp qbp qm p
15265: qm m qbp qbp p qm
33625: qm m qbp p qm qbp
27145: qm m qbp p qp qbm
20665: qm m qbp p qbm qp
14185: qm m qbp p qbp qm
36685: qm m p qm qbp qbp
35605: qm m p qp qbm qbp
29125: qm m p qp qbp qbm
34525: qm m p qbm qp qbp
21565: qm m p qbm qbp qp
33445: qm m p qbp qm qbp
26965: qm m p qbp qp qbm
20485: qm m p qbp qbm qp
14005: qm m p qbp qbp qm
36295: qm qm m m qbp qbp
31975: qm qm m qbp m qbp
6055: qm qm m qbp qbp m
44935: qm qm m qbp qbp p
38455: qm qm m qbp p qbp
37375: qm qm m p qbp qbp
31255: qm qm qbp m m qbp
5335: qm qm qbp m qbp m
44215: qm qm qbp m qbp p
37735: qm qm qbp m p qbp
1015: qm qm qbp qbp m m
39895: qm qm qbp qbp m p
7495: qm qm qbp qbp p m
46375: qm qm qbp qbp p p
32335: qm qm qbp p m qbp
6415: qm qm qbp p qbp m
45295: qm qm qbp p qbp p
38815: qm qm qbp p p qbp
36475: qm qm p m qbp qbp
32155: qm qm p qbp m qbp
6235: qm qm p qbp qbp m
45115: qm qm p qbp qbp p
38635: qm qm p qbp p qbp
37555: qm qm p p qbp qbp
35005: qm qp m m qbm qbp
28525: qm qp m m qbp qbm
31765: qm qp m qbm m qbp
5845: qm qp m qbm qbp m
44725: qm qp m qbm qbp p
38245: qm qp m qbm p qbp
24205: qm qp m qbp m qbm
4765: qm qp m qbp qbm m
43645: qm qp m qbp qbm p
30685: qm qp m qbp p qbm
36085: qm qp m p qbm qbp
29605: qm qp m p qbp qbm
31225: qm qp qbm m m qbp
5305: qm qp qbm m qbp m
44185: qm qp qbm m qbp p
37705: qm qp qbm m p qbp
985: qm qp qbm qbp m m
39865: qm qp qbm qbp m p
7465: qm qp qbm qbp p m
46345: qm qp qbm qbp p p
32305: qm qp qbm p m qbp
6385: qm qp qbm p qbp m
45265: qm qp qbm p qbp p
38785: qm qp qbm p p qbp
23485: qm qp qbp m m qbm
4045: qm qp qbp m qbm m
42925: qm qp qbp m qbm p
29965: qm qp qbp m p qbm
805: qm qp qbp qbm m m
39685: qm qp qbp qbm m p
7285: qm qp qbp qbm p m
46165: qm qp qbp qbm p p
24565: qm qp qbp p m qbm
5125: qm qp qbp p qbm m
44005: qm qp qbp p qbm p
31045: qm qp qbp p p qbm
35185: qm qp p m qbm qbp
28705: qm qp p m qbp qbm
31945: qm qp p qbm m qbp
6025: qm qp p qbm qbp m
44905: qm qp p qbm qbp p
38425: qm qp p qbm p qbp
24385: qm qp p qbp m qbm
4945: qm qp p qbp qbm m
43825: qm qp p qbp qbm p
30865: qm qp p qbp p qbm
36265: qm qp p p qbm qbp
29785: qm qp p p qbp qbm
33715: qm qbm m m qp qbp
20755: qm qbm m m qbp qp
31555: qm qbm m qp m qbp
5635: qm qbm m qp qbp m
44515: qm qbm m qp qbp p
38035: qm qbm m qp p qbp
16435: qm qbm m qbp m qp
3475: qm qbm m qbp qp m
42355: qm qbm m qbp qp p
22915: qm qbm m qbp p qp
34795: qm qbm m p qp qbp
21835: qm qbm m p qbp qp
31195: qm qbm qp m m qbp
5275: qm qbm qp m qbp m
44155: qm qbm qp m qbp p
37675: qm qbm qp m p qbp
955: qm qbm qp qbp m m
39835: qm qbm qp qbp m p
7435: qm qbm qp qbp p m
46315: qm qbm qp qbp p p
32275: qm qbm qp p m qbp
6355: qm qbm qp p qbp m
45235: qm qbm qp p qbp p
38755: qm qbm qp p p qbp
15715: qm qbm qbp m m qp
2755: qm qbm qbp m qp m
41635: qm qbm qbp m qp p
22195: qm qbm qbp m p qp
595: qm qbm qbp qp m m
39475: qm qbm qbp qp m p
7075: qm qbm qbp qp p m
45955: qm qbm qbp qp p p
16795: qm qbm qbp p m qp
3835: qm qbm qbp p qp m
42715: qm qbm qbp p qp p
23275: qm qbm qbp p p qp
33895: qm qbm p m qp qbp
20935: qm qbm p m qbp qp
31735: qm qbm p qp m qbp
5815: qm qbm p qp qbp m
44695: qm qbm p qp qbp p
38215: qm qbm p qp p qbp
16615: qm qbm p qbp m qp
3655: qm qbm p qbp qp m
42535: qm qbm p qbp qp p
23095: qm qbm p qbp p qp
34975: qm qbm p p qp qbp
22015: qm qbm p p qbp qp
32425: qm qbp m m qm qbp
25945: qm qbp m m qp qbm
19465: qm qbp m m qbm qp
12985: qm qbp m m qbp qm
31345: qm qbp m qm m qbp
5425: qm qbp m qm qbp m
44305: qm qbp m qm qbp p
37825: qm qbp m qm p qbp
23785: qm qbp m qp m qbm
4345: qm qbp m qp qbm m
43225: qm qbp m qp qbm p
30265: qm qbp m qp p qbm
16225: qm qbp m qbm m qp
3265: qm qbp m qbm qp m
42145: qm qbp m qbm qp p
22705: qm qbp m qbm p qp
8665: qm qbp m qbp m qm
2185: qm qbp m qbp qm m
41065: qm qbp m qbp qm p
15145: qm qbp m qbp p qm
33505: qm qbp m p qm qbp
27025: qm qbp m p qp qbm
20545: qm qbp m p qbm qp
14065: qm qbp m p qbp qm
31165: qm qbp qm m m qbp
5245: qm qbp qm m qbp m
44125: qm qbp qm m qbp p
37645: qm qbp qm m p qbp
925: qm qbp qm qbp m m
39805: qm qbp qm qbp m p
7405: qm qbp qm qbp p m
46285: qm qbp qm qbp p p
32245: qm qbp qm p m qbp
6325: qm qbp qm p qbp m
45205: qm qbp qm p qbp p
38725: qm qbp qm p p qbp
23425: qm qbp qp m m qbm
3985: qm qbp qp m qbm m
42865: qm qbp qp m qbm p
29905: qm qbp qp m p qbm
745: qm qbp qp qbm m m
39625: qm qbp qp qbm m p
7225: qm qbp qp qbm p m
46105: qm qbp qp qbm p p
24505: qm qbp qp p m qbm
5065: qm qbp qp p qbm m
43945: qm qbp qp p qbm p
30985: qm qbp qp p p qbm
15685: qm qbp qbm m m qp
2725: qm qbp qbm m qp m
41605: qm qbp qbm m qp p
22165: qm qbp qbm m p qp
565: qm qbp qbm qp m m
39445: qm qbp qbm qp m p
7045: qm qbp qbm qp p m
45925: qm qbp qbm qp p p
16765: qm qbp qbm p m qp
3805: qm qbp qbm p qp m
42685: qm qbp qbm p qp p
23245: qm qbp qbm p p qp
7945: qm qbp qbp m m qm
1465: qm qbp qbp m qm m
40345: qm qbp qbp m qm p
14425: qm qbp qbp m p qm
385: qm qbp qbp qm m m
39265: qm qbp qbp qm m p
6865: qm qbp qbp qm p m
45745: qm qbp qbp qm p p
9025: qm qbp qbp p m qm
2545: qm qbp qbp p qm m
41425: qm qbp qbp p qm p
15505: qm qbp qbp p p qm
32605: qm qbp p m qm qbp
26125: qm qbp p m qp qbm
19645: qm qbp p m qbm qp
13165: qm qbp p m qbp qm
31525: qm qbp p qm m qbp
5605: qm qbp p qm qbp m
44485: qm qbp p qm qbp p
38005: qm qbp p qm p qbp
23965: qm qbp p qp m qbm
4525: qm qbp p qp qbm m
43405: qm qbp p qp qbm p
30445: qm qbp p qp p qbm
16405: qm qbp p qbm m qp
3445: qm qbp p qbm qp m
42325: qm qbp p qbm qp p
22885: qm qbp p qbm p qp
8845: qm qbp p qbp m qm
2365: qm qbp p qbp qm m
41245: qm qbp p qbp qm p
15325: qm qbp p qbp p qm
33685: qm qbp p p qm qbp
27205: qm qbp p p qp qbm
20725: qm qbp p p qbm qp
14245: qm qbp p p qbp qm
36535: qm p m qm qbp qbp
35455: qm p m qp qbm qbp
28975: qm p m qp qbp qbm
34375: qm p m qbm qp qbp
21415: qm p m qbm qbp qp
33295: qm p m qbp qm qbp
26815: qm p m qbp qp qbm
20335: qm p m qbp qbm qp
13855: qm p m qbp qbp qm
36355: qm p qm m qbp qbp
32035: qm p qm qbp m qbp
6115: qm p qm qbp qbp m
44995: qm p qm qbp qbp p
38515: qm p qm qbp p qbp
37435: qm p qm p qbp qbp
35095: qm p qp m qbm qbp
28615: qm p qp m qbp qbm
31855: qm p qp qbm m qbp
5935: qm p qp qbm qbp m
44815: qm p qp qbm qbp p
38335: qm p qp qbm p qbp
24295: qm p qp qbp m qbm
4855: qm p qp qbp qbm m
43735: qm p qp qbp qbm p
30775: qm p qp qbp p qbm
36175: qm p qp p qbm qbp
29695: qm p qp p qbp qbm
33835: qm p qbm m qp qbp
20875: qm p qbm m qbp qp
31675: qm p qbm qp m qbp
5755: qm p qbm qp qbp m
44635: qm p qbm qp qbp p
38155: qm p qbm qp p qbp
16555: qm p qbm qbp m qp
3595: qm p qbm qbp qp m
42475: qm p qbm qbp qp p
23035: qm p qbm qbp p qp
34915: qm p qbm p qp qbp
21955: qm p qbm p qbp qp
32575: qm p qbp m qm qbp
26095: qm p qbp m qp qbm
19615: qm p qbp m qbm qp
13135: qm p qbp m qbp qm
31495: qm p qbp qm m qbp
5575: qm p qbp qm qbp m
44455: qm p qbp qm qbp p
37975: qm p qbp qm p qbp
23935: qm p qbp qp m qbm
4495: qm p qbp qp qbm m
43375: qm p qbp qp qbm p
30415: qm p qbp qp p qbm
16375: qm p qbp qbm m qp
3415: qm p qbp qbm qp m
42295: qm p qbp qbm qp p
22855: qm p qbp qbm p qp
8815: qm p qbp qbp m qm
2335: qm p qbp qbp qm m
41215: qm p qbp qbp qm p
15295: qm p qbp qbp p qm
33655: qm p qbp p qm qbp
27175: qm p qbp p qp qbm
20695: qm p qbp p qbm qp
14215: qm p qbp p qbp qm
36715: qm p p qm qbp qbp
35635: qm p p qp qbm qbp
29155: qm p p qp qbp qbm
34555: qm p p qbm qp qbp
21595: qm p p qbm qbp qp
33475: qm p p qbp qm qbp
26995: qm p p qbp qp qbm
20515: qm p p qbp qbm qp
14035: qm p p qbp qbp qm
35210: qp m m qm qbm qbp
28730: qp m m qm qbp qbm
27650: qp m m qp qbm qbm
33050: qp m m qbm qm qbp
26570: qp m m qbm qp qbm
20090: qp m m qbm qbm qp
13610: qp m m qbm qbp qm
25490: qp m m qbp qm qbm
12530: qp m m qbp qbm qm
35030: qp m qm m qbm qbp
28550: qp m qm m qbp qbm
31790: qp m qm qbm m qbp
5870: qp m qm qbm qbp m
44750: qp m qm qbm qbp p
38270: qp m qm qbm p qbp
24230: qp m qm qbp m qbm
4790: qp m qm qbp qbm m
43670: qp m qm qbp qbm p
30710: qp m qm qbp p qbm
36110: qp m qm p qbm qbp
29630: qp m qm p qbp qbm
27290: qp m qp m qbm qbm
24050: qp m qp qbm m qbm
4610: qp m qp qbm qbm m
43490: qp m qp qbm qbm p
30530: qp m qp qbm p qbm
28370: qp m qp p qbm qbm
32510: qp m qbm m qm qbp
26030: qp m qbm m qp qbm
19550: qp m qbm m qbm qp
13070: qp m qbm m qbp qm
31430: qp m qbm qm m qbp
5510: qp m qbm qm qbp m
44390: qp m qbm qm qbp p
37910: qp m qbm qm p qbp
23870: qp m qbm qp m qbm
4430: qp m qbm qp qbm m
43310: qp m qbm qp qbm p
30350: qp m qbm qp p qbm
16310: qp m qbm qbm m qp
3350: qp m qbm qbm qp m
42230: qp m qbm qbm qp p
22790: qp m qbm qbm p qp
8750: qp m qbm qbp m qm
2270: qp m qbm qbp qm m
41150: qp m qbm qbp qm p
15230: qp m qbm qbp p qm
33590: qp m qbm p qm qbp
27110: qp m qbm p qp qbm
20630: qp m qbm p qbm qp
14150: qp m qbm p qbp qm
24770: qp m qbp m qm qbm
11810: qp m qbp m qbm qm
23690: qp m qbp qm m qbm
4250: qp m qbp qm qbm m
43130: qp m qbp qm qbm p
30170: qp m qbp qm p qbm
8570: qp m qbp qbm m qm
2090: qp m qbp qbm qm m
40970: qp m qbp qbm qm p
15050: qp m qbp qbm p qm
25850: qp m qbp p qm qbm
12890: qp m qbp p qbm qm
35390: qp m p qm qbm qbp
28910: qp m p qm qbp qbm
27830: qp m p qp qbm qbm
33230: qp m p qbm qm qbp
26750: qp m p qbm qp qbm
20270: qp m p qbm qbm qp
13790: qp m p qbm qbp qm
25670: qp m p qbp qm qbm
12710: qp m p qbp qbm qm
35000: qp qm m m qbm qbp
28520: qp qm m m qbp qbm
31760: qp qm m qbm m qbp
5840: qp qm m qbm qbp m
44720: qp qm m qbm qbp p
38240: qp qm m qbm p qbp
24200: qp qm m qbp m qbm
4760: qp qm m qbp qbm m
43640: qp qm m qbp qbm p
30680: qp qm m qbp p qbm
36080: qp qm m p qbm qbp
29600: qp qm m p qbp qbm
31220: qp qm qbm m m qbp
5300: qp qm qbm m qbp m
44180: qp qm qbm m qbp p
37700: qp qm qbm m p qbp
980: qp qm qbm qbp m m
39860: qp qm qbm qbp m p
7460: qp qm qbm qbp p m
46340: qp qm qbm qbp p p
32300: qp qm qbm p m qbp
6380: qp qm qbm p qbp m
45260: qp qm qbm p qbp p
38780: qp qm qbm p p qbp
23480: qp qm qbp m m qbm
4040: qp qm qbp m qbm m
42920: qp qm qbp m qbm p
29960: qp qm qbp m p qbm
800: qp qm qbp qbm m m
39680: qp qm qbp qbm m p
7280: qp qm qbp qbm p m
46160: qp qm qbp qbm p p
24560: qp qm qbp p m qbm
5120: qp qm qbp p qbm m
44000: qp qm qbp p qbm p
31040: qp qm qbp p p qbm
35180: qp qm p m qbm qbp
28700: qp qm p m qbp qbm
31940: qp qm p qbm m qbp
6020: qp qm p qbm qbp m
44900: qp qm p qbm qbp p
38420: qp qm p qbm p qbp
24380: qp qm p qbp m qbm
4940: qp qm p qbp qbm m
43820: qp qm p qbp qbm p
30860: qp qm p qbp p qbm
36260: qp qm p p qbm qbp
29780: qp qm p p qbp qbm
27230: qp qp m m qbm qbm
23990: qp qp m qbm m qbm
4550: qp qp m qbm qbm m
43430: qp qp m qbm qbm p
30470: qp qp m qbm p qbm
28310: qp qp m p qbm qbm
23450: qp qp qbm m m qbm
4010: qp qp qbm m qbm m
42890: qp qp qbm m qbm p
29930: qp qp qbm m p qbm
770: qp qp qbm qbm m m
39650: qp qp qbm qbm m p
7250: qp qp qbm qbm p m
46130: qp qp qbm qbm p p
24530: qp qp qbm p m qbm
5090: qp qp qbm p qbm m
43970: qp qp qbm p qbm p
31010: qp qp qbm p p qbm
27410: qp qp p m qbm qbm
24170: qp qp p qbm m qbm
4730: qp qp p qbm qbm m
43610: qp qp p qbm qbm p
30650: qp qp p qbm p qbm
28490: qp qp p p qbm qbm
32420: qp qbm m m qm qbp
25940: qp qbm m m qp qbm
19460: qp qbm m m qbm qp
12980: qp qbm m m qbp qm
31340: qp qbm m qm m qbp
5420: qp qbm m qm qbp m
44300: qp qbm m qm qbp p
37820: qp qbm m qm p qbp
23780: qp qbm m qp m qbm
4340: qp qbm m qp qbm m
43220: qp qbm m qp qbm p
30260: qp qbm m qp p qbm
16220: qp qbm m qbm m qp
3260: qp qbm m qbm qp m
42140: qp qbm m qbm qp p
22700: qp qbm m qbm p qp
8660: qp qbm m qbp m qm
2180: qp qbm m qbp qm m
41060: qp qbm m qbp qm p
15140: qp qbm m qbp p qm
33500: qp qbm m p qm qbp
27020: qp qbm m p qp qbm
20540: qp qbm m p qbm qp
14060: qp qbm m p qbp qm
31160: qp qbm qm m m qbp
5240: qp qbm qm m qbp m
44120: qp qbm qm m qbp p
37640: qp qbm qm m p qbp
920: qp qbm qm qbp m m
39800: qp qbm qm qbp m p
7400: qp qbm qm qbp p m
46280: qp qbm qm qbp p p
32240: qp qbm qm p m qbp
6320: qp qbm qm p qbp m
45200: qp qbm qm p qbp p
38720: qp qbm qm p p qbp
23420: qp qbm qp m m qbm
3980: qp qbm qp m qbm m
42860: qp qbm qp m qbm p
29900: qp qbm qp m p qbm
740: qp qbm qp qbm m m
39620: qp qbm qp qbm m p
7220: qp qbm qp qbm p m
46100: qp qbm qp qbm p p
24500: qp qbm qp p m qbm
5060: qp qbm qp p qbm m
43940: qp qbm qp p qbm p
30980: qp qbm qp p p qbm
15680: qp qbm qbm m m qp
2720: qp qbm qbm m qp m
41600: qp qbm qbm m qp p
22160: qp qbm qbm m p qp
560: qp qbm qbm qp m m
39440: qp qbm qbm qp m p
7040: qp qbm qbm qp p m
45920: qp qbm qbm qp p p
16760: qp qbm qbm p m qp
3800: qp qbm qbm p qp m
42680: qp qbm qbm p qp p
23240: qp qbm qbm p p qp
7940: qp qbm qbp m m qm
1460: qp qbm qbp m qm m
40340: qp qbm qbp m qm p
14420: qp qbm qbp m p qm
380: qp qbm qbp qm m m
39260: qp qbm qbp qm m p
6860: qp qbm qbp qm p m
45740: qp qbm qbp qm p p
9020: qp qbm qbp p m qm
2540: qp qbm qbp p qm m
41420: qp qbm qbp p qm p
15500: qp qbm qbp p p qm
32600: qp qbm p m qm qbp
26120: qp qbm p m qp qbm
19640: qp qbm p m qbm qp
13160: qp qbm p m qbp qm
31520: qp qbm p qm m qbp
5600: qp qbm p qm qbp m
44480: qp qbm p qm qbp p
38000: qp qbm p qm p qbp
23960: qp qbm p qp m qbm
4520: qp qbm p qp qbm m
43400: qp qbm p qp qbm p
30440: qp qbm p qp p qbm
16400: qp qbm p qbm m qp
3440: qp qbm p qbm qp m
42320: qp qbm p qbm qp p
22880: qp qbm p qbm p qp
8840: qp qbm p qbp m qm
2360: qp qbm p qbp qm m
41240: qp qbm p qbp qm p
15320: qp qbm p qbp p qm
33680: qp qbm p p qm qbp
27200: qp qbm p p qp qbm
20720: qp qbm p p qbm qp
14240: qp qbm p p qbp qm
24650: qp qbp m m qm qbm
11690: qp qbp m m qbm qm
23570: qp qbp m qm m qbm
4130: qp qbp m qm qbm m
43010: qp qbp m qm qbm p
30050: qp qbp m qm p qbm
8450: qp qbp m qbm m qm
1970: qp qbp m qbm qm m
40850: qp qbp m qbm qm p
14930: qp qbp m qbm p qm
25730: qp qbp m p qm qbm
12770: qp qbp m p qbm qm
23390: qp qbp qm m m qbm
3950: qp qbp qm m qbm m
42830: qp qbp qm m qbm p
29870: qp qbp qm m p qbm
710: qp qbp qm qbm m m
39590: qp qbp qm qbm m p
7190: qp qbp qm qbm p m
46070: qp qbp qm qbm p p
24470: qp qbp qm p m qbm
5030: qp qbp qm p qbm m
43910: qp qbp qm p qbm p
30950: qp qbp qm p p qbm
7910: qp qbp qbm m m qm
1430: qp qbp qbm m qm m
40310: qp qbp qbm m qm p
14390: qp qbp qbm m p qm
350: qp qbp qbm qm m m
39230: qp qbp qbm qm m p
6830: qp qbp qbm qm p m
45710: qp qbp qbm qm p p
8990: qp qbp qbm p m qm
2510: qp qbp qbm p qm m
41390: qp qbp qbm p qm p
15470: qp qbp qbm p p qm
24830: qp qbp p m qm qbm
11870: qp qbp p m qbm qm
23750: qp qbp p qm m qbm
4310: qp qbp p qm qbm m
43190: qp qbp p qm qbm p
30230: qp qbp p qm p qbm
8630: qp qbp p qbm m qm
2150: qp qbp p qbm qm m
41030: qp qbp p qbm qm p
15110: qp qbp p qbm p qm
25910: qp qbp p p qm qbm
12950: qp qbp p p qbm qm
35240: qp p m qm qbm qbp
28760: qp p m qm qbp qbm
27680: qp p m qp qbm qbm
33080: qp p m qbm qm qbp
26600: qp p m qbm qp qbm
20120: qp p m qbm qbm qp
13640: qp p m qbm qbp qm
25520: qp p m qbp qm qbm
12560: qp p m qbp qbm qm
35060: qp p qm m qbm qbp
28580: qp p qm m qbp qbm
31820: qp p qm qbm m qbp
5900: qp p qm qbm qbp m
44780: qp p qm qbm qbp p
38300: qp p qm qbm p qbp
24260: qp p qm qbp m qbm
4820: qp p qm qbp qbm m
43700: qp p qm qbp qbm p
30740: qp p qm qbp p qbm
36140: qp p qm p qbm qbp
29660: qp p qm p qbp qbm
27320: qp p qp m qbm qbm
24080: qp p qp qbm m qbm
4640: qp p qp qbm qbm m
43520: qp p qp qbm qbm p
30560: qp p qp qbm p qbm
28400: qp p qp p qbm qbm
32540: qp p qbm m qm qbp
26060: qp p qbm m qp qbm
19580: qp p qbm m qbm qp
13100: qp p qbm m qbp qm
31460: qp p qbm qm m qbp
5540: qp p qbm qm qbp m
44420: qp p qbm qm qbp p
37940: qp p qbm qm p qbp
23900: qp p qbm qp m qbm
4460: qp p qbm qp qbm m
43340: qp p qbm qp qbm p
30380: qp p qbm qp p qbm
16340: qp p qbm qbm m qp
3380: qp p qbm qbm qp m
42260: qp p qbm qbm qp p
22820: qp p qbm qbm p qp
8780: qp p qbm qbp m qm
2300: qp p qbm qbp qm m
41180: qp p qbm qbp qm p
15260: qp p qbm qbp p qm
33620: qp p qbm p qm qbp
27140: qp p qbm p qp qbm
20660: qp p qbm p qbm qp
14180: qp p qbm p qbp qm
24800: qp p qbp m qm qbm
11840: qp p qbp m qbm qm
23720: qp p qbp qm m qbm
4280: qp p qbp qm qbm m
43160: qp p qbp qm qbm p
30200: qp p qbp qm p qbm
8600: qp p qbp qbm m qm
2120: qp p qbp qbm qm m
41000: qp p qbp qbm qm p
15080: qp p qbp qbm p qm
25880: qp p qbp p qm qbm
12920: qp p qbp p qbm qm
35420: qp p p qm qbm qbp
28940: qp p p qm qbp qbm
27860: qp p p qp qbm qbm
33260: qp p p qbm qm qbp
26780: qp p p qbm qp qbm
20300: qp p p qbm qbm qp
13820: qp p p qbm qbp qm
25700: qp p p qbp qm qbm
12740: qp p p qbp qbm qm
33915: qbm m m qm qp qbp
20955: qbm m m qm qbp qp
32835: qbm m m qp qm qbp
26355: qbm m m qp qp qbm
19875: qbm m m qp qbm qp
13395: qbm m m qp qbp qm
18795: qbm m m qbm qp qp
17715: qbm m m qbp qm qp
11235: qbm m m qbp qp qm
33735: qbm m qm m qp qbp
20775: qbm m qm m qbp qp
31575: qbm m qm qp m qbp
5655: qbm m qm qp qbp m
44535: qbm m qm qp qbp p
38055: qbm m qm qp p qbp
16455: qbm m qm qbp m qp
3495: qbm m qm qbp qp m
42375: qbm m qm qbp qp p
22935: qbm m qm qbp p qp
34815: qbm m qm p qp qbp
21855: qbm m qm p qbp qp
32475: qbm m qp m qm qbp
25995: qbm m qp m qp qbm
19515: qbm m qp m qbm qp
13035: qbm m qp m qbp qm
31395: qbm m qp qm m qbp
5475: qbm m qp qm qbp m
44355: qbm m qp qm qbp p
37875: qbm m qp qm p qbp
23835: qbm m qp qp m qbm
4395: qbm m qp qp qbm m
43275: qbm m qp qp qbm p
30315: qbm m qp qp p qbm
16275: qbm m qp qbm m qp
3315: qbm m qp qbm qp m
42195: qbm m qp qbm qp p
22755: qbm m qp qbm p qp
8715: qbm m qp qbp m qm
2235: qbm m qp qbp qm m
41115: qbm m qp qbp qm p
15195: qbm m qp qbp p qm
33555: qbm m qp p qm qbp
27075: qbm m qp p qp qbm
20595: qbm m qp p qbm qp
14115: qbm m qp p qbp qm
18255: qbm m qbm m qp qp
16095: qbm m qbm qp m qp
3135: qbm m qbm qp qp m
42015: qbm m qbm qp qp p
22575: qbm m qbm qp p qp
19335: qbm m qbm p qp qp
16995: qbm m qbp m qm qp
10515: qbm m qbp m qp qm
15915: qbm m qbp qm m qp
2955: qbm m qbp qm qp m
41835: qbm m qbp qm qp p
22395: qbm m qbp qm p qp
8355: qbm m qbp qp m qm
1875: qbm m qbp qp qm m
40755: qbm m qbp qp qm p
14835: qbm m qbp qp p qm
18075: qbm m qbp p qm qp
11595: qbm m qbp p qp qm
34095: qbm m p qm qp qbp
21135: qbm m p qm qbp qp
33015: qbm m p qp qm qbp
26535: qbm m p qp qp qbm
20055: qbm m p qp qbm qp
13575: qbm m p qp qbp qm
18975: qbm m p qbm qp qp
17895: qbm m p qbp qm qp
11415: qbm m p qbp qp qm
33705: qbm qm m m qp qbp
20745: qbm qm m m qbp qp
31545: qbm qm m qp m qbp
5625: qbm qm m qp qbp m
44505: qbm qm m qp qbp p
38025: qbm qm m qp p qbp
16425: qbm qm m qbp m qp
3465: qbm qm m qbp qp m
42345: qbm qm m qbp qp p
22905: qbm qm m qbp p qp
34785: qbm qm m p qp qbp
21825: qbm qm m p qbp qp
31185: qbm qm qp m m qbp
5265: qbm qm qp m qbp m
44145: qbm qm qp m qbp p
37665: qbm qm qp m p qbp
945: qbm qm qp qbp m m
39825: qbm qm qp qbp m p
7425: qbm qm qp qbp p m
46305: qbm qm qp qbp p p
32265: qbm qm qp p m qbp
6345: qbm qm qp p qbp m
45225: qbm qm qp p qbp p
38745: qbm qm qp p p qbp
15705: qbm qm qbp m m qp
2745: qbm qm qbp m qp m
41625: qbm qm qbp m qp p
22185: qbm qm qbp m p qp
585: qbm qm qbp qp m m
39465: qbm qm qbp qp m p
7065: qbm qm qbp qp p m
45945: qbm qm qbp qp p p
16785: qbm qm qbp p m qp
3825: qbm qm qbp p qp m
42705: qbm qm qbp p qp p
23265: qbm qm qbp p p qp
33885: qbm qm p m qp qbp
20925: qbm qm p m qbp qp
31725: qbm qm p qp m qbp
5805: qbm qm p qp qbp m
44685: qbm qm p qp qbp p
38205: qbm qm p qp p qbp
16605: qbm qm p qbp m qp
3645: qbm qm p qbp qp m
42525: qbm qm p qbp qp p
23085: qbm qm p qbp p qp
34965: qbm qm p p qp qbp
22005: qbm qm p p qbp qp
32415: qbm qp m m qm qbp
25935: qbm qp m m qp qbm
19455: qbm qp m m qbm qp
12975: qbm qp m m qbp qm
31335: qbm qp m qm m qbp
5415: qbm qp m qm qbp m
44295: qbm qp m qm qbp p
37815: qbm qp m qm p qbp
23775: qbm qp m qp m qbm
4335: qbm qp m qp qbm m
43215: qbm qp m qp qbm p
30255: qbm qp m qp p qbm
16215: qbm qp m qbm m qp
3255: qbm qp m qbm qp m
42135: qbm qp m qbm qp p
22695: qbm qp m qbm p qp
8655: qbm qp m qbp m qm
2175: qbm qp m qbp qm m
41055: qbm qp m qbp qm p
15135: qbm qp m qbp p qm
33495: qbm qp m p qm qbp
27015: qbm qp m p qp qbm
20535: qbm qp m p qbm qp
14055: qbm qp m p qbp qm
31155: qbm qp qm m m qbp
5235: qbm qp qm m qbp m
44115: qbm qp qm m qbp p
37635: qbm qp qm m p qbp
915: qbm qp qm qbp m m
39795: qbm qp qm qbp m p
7395: qbm qp qm qbp p m
46275: qbm qp qm qbp p p
32235: qbm qp qm p m qbp
6315: qbm qp qm p qbp m
45195: qbm qp qm p qbp p
38715: qbm qp qm p p qbp
23415: qbm qp qp m m qbm
3975: qbm qp qp m qbm m
42855: qbm qp qp m qbm p
29895: qbm qp qp m p qbm
735: qbm qp qp qbm m m
39615: qbm qp qp qbm m p
7215: qbm qp qp qbm p m
46095: qbm qp qp qbm p p
24495: qbm qp qp p m qbm
5055: qbm qp qp p qbm m
43935: qbm qp qp p qbm p
30975: qbm qp qp p p qbm
15675: qbm qp qbm m m qp
2715: qbm qp qbm m qp m
41595: qbm qp qbm m qp p
22155: qbm qp qbm m p qp
555: qbm qp qbm qp m m
39435: qbm qp qbm qp m p
7035: qbm qp qbm qp p m
45915: qbm qp qbm qp p p
16755: qbm qp qbm p m qp
3795: qbm qp qbm p qp m
42675: qbm qp qbm p qp p
23235: qbm qp qbm p p qp
7935: qbm qp qbp m m qm
1455: qbm qp qbp m qm m
40335: qbm qp qbp m qm p
14415: qbm qp qbp m p qm
375: qbm qp qbp qm m m
39255: qbm qp qbp qm m p
6855: qbm qp qbp qm p m
45735: qbm qp qbp qm p p
9015: qbm qp qbp p m qm
2535: qbm qp qbp p qm m
41415: qbm qp qbp p qm p
15495: qbm qp qbp p p qm
32595: qbm qp p m qm qbp
26115: qbm qp p m qp qbm
19635: qbm qp p m qbm qp
13155: qbm qp p m qbp qm
31515: qbm qp p qm m qbp
5595: qbm qp p qm qbp m
44475: qbm qp p qm qbp p
37995: qbm qp p qm p qbp
23955: qbm qp p qp m qbm
4515: qbm qp p qp qbm m
43395: qbm qp p qp qbm p
30435: qbm qp p qp p qbm
16395: qbm qp p qbm m qp
3435: qbm qp p qbm qp m
42315: qbm qp p qbm qp p
22875: qbm qp p qbm p qp
8835: qbm qp p qbp m qm
2355: qbm qp p qbp qm m
41235: qbm qp p qbp qm p
15315: qbm qp p qbp p qm
33675: qbm qp p p qm qbp
27195: qbm qp p p qp qbm
20715: qbm qp p p qbm qp
14235: qbm qp p p qbp qm
18165: qbm qbm m m qp qp
16005: qbm qbm m qp m qp
3045: qbm qbm m qp qp m
41925: qbm qbm m qp qp p
22485: qbm qbm m qp p qp
19245: qbm qbm m p qp qp
15645: qbm qbm qp m m qp
2685: qbm qbm qp m qp m
41565: qbm qbm qp m qp p
22125: qbm qbm qp m p qp
525: qbm qbm qp qp m m
39405: qbm qbm qp qp m p
7005: qbm qbm qp qp p m
45885: qbm qbm qp qp p p
16725: qbm qbm qp p m qp
3765: qbm qbm qp p qp m
42645: qbm qbm qp p qp p
23205: qbm qbm qp p p qp
18345: qbm qbm p m qp qp
16185: qbm qbm p qp m qp
3225: qbm qbm p qp qp m
42105: qbm qbm p qp qp p
22665: qbm qbm p qp p qp
19425: qbm qbm p p qp qp
16875: qbm qbp m m qm qp
10395: qbm qbp m m qp qm
15795: qbm qbp m qm m qp
2835: qbm qbp m qm qp m
41715: qbm qbp m qm qp p
22275: qbm qbp m qm p qp
8235: qbm qbp m qp m qm
1755: qbm qbp m qp qm m
40635: qbm qbp m qp qm p
14715: qbm qbp m qp p qm
17955: qbm qbp m p qm qp
11475: qbm qbp m p qp qm
15615: qbm qbp qm m m qp
2655: qbm qbp qm m qp m
41535: qbm qbp qm m qp p
22095: qbm qbp qm m p qp
495: qbm qbp qm qp m m
39375: qbm qbp qm qp m p
6975: qbm qbp qm qp p m
45855: qbm qbp qm qp p p
16695: qbm qbp qm p m qp
3735: qbm qbp qm p qp m
42615: qbm qbp qm p qp p
23175: qbm qbp qm p p qp
7875: qbm qbp qp m m qm
1395: qbm qbp qp m qm m
40275: qbm qbp qp m qm p
14355: qbm qbp qp m p qm
315: qbm qbp qp qm m m
39195: qbm qbp qp qm m p
6795: qbm qbp qp qm p m
45675: qbm qbp qp qm p p
8955: qbm qbp qp p m qm
2475: qbm qbp qp p qm m
41355: qbm qbp qp p qm p
15435: qbm qbp qp p p qm
17055: qbm qbp p m qm qp
10575: qbm qbp p m qp qm
15975: qbm qbp p qm m qp
3015: qbm qbp p qm qp m
41895: qbm qbp p qm qp p
22455: qbm qbp p qm p qp
8415: qbm qbp p qp m qm
1935: qbm qbp p qp qm m
40815: qbm qbp p qp qm p
14895: qbm qbp p qp p qm
18135: qbm qbp p p qm qp
11655: qbm qbp p p qp qm
33945: qbm p m qm qp qbp
20985: qbm p m qm qbp qp
32865: qbm p m qp qm qbp
26385: qbm p m qp qp qbm
19905: qbm p m qp qbm qp
13425: qbm p m qp qbp qm
18825: qbm p m qbm qp qp
17745: qbm p m qbp qm qp
11265: qbm p m qbp qp qm
33765: qbm p qm m qp qbp
20805: qbm p qm m qbp qp
31605: qbm p qm qp m qbp
5685: qbm p qm qp qbp m
44565: qbm p qm qp qbp p
38085: qbm p qm qp p qbp
16485: qbm p qm qbp m qp
3525: qbm p qm qbp qp m
42405: qbm p qm qbp qp p
22965: qbm p qm qbp p qp
34845: qbm p qm p qp qbp
21885: qbm p qm p qbp qp
32505: qbm p qp m qm qbp
26025: qbm p qp m qp qbm
19545: qbm p qp m qbm qp
13065: qbm p qp m qbp qm
31425: qbm p qp qm m qbp
5505: qbm p qp qm qbp m
44385: qbm p qp qm qbp p
37905: qbm p qp qm p qbp
23865: qbm p qp qp m qbm
4425: qbm p qp qp qbm m
43305: qbm p qp qp qbm p
30345: qbm p qp qp p qbm
16305: qbm p qp qbm m qp
3345: qbm p qp qbm qp m
42225: qbm p qp qbm qp p
22785: qbm p qp qbm p qp
8745: qbm p qp qbp m qm
2265: qbm p qp qbp qm m
41145: qbm p qp qbp qm p
15225: qbm p qp qbp p qm
33585: qbm p qp p qm qbp
27105: qbm p qp p qp qbm
20625: qbm p qp p qbm qp
14145: qbm p qp p qbp qm
18285: qbm p qbm m qp qp
16125: qbm p qbm qp m qp
3165: qbm p qbm qp qp m
42045: qbm p qbm qp qp p
22605: qbm p qbm qp p qp
19365: qbm p qbm p qp qp
17025: qbm p qbp m qm qp
10545: qbm p qbp m qp qm
15945: qbm p qbp qm m qp
2985: qbm p qbp qm qp m
41865: qbm p qbp qm qp p
22425: qbm p qbp qm p qp
8385: qbm p qbp qp m qm
1905: qbm p qbp qp qm m
40785: qbm p qbp qp qm p
14865: qbm p qbp qp p qm
18105: qbm p qbp p qm qp
11625: qbm p qbp p qp qm
34125: qbm p p qm qp qbp
21165: qbm p p qm qbp qp
33045: qbm p p qp qm qbp
26565: qbm p p qp qp qbm
20085: qbm p p qp qbm qp
13605: qbm p p qp qbp qm
19005: qbm p p qbm qp qp
17925: qbm p p qbp qm qp
11445: qbm p p qbp qp qm
32620: qbp m m qm qm qbp
26140: qbp m m qm qp qbm
19660: qbp m m qm qbm qp
13180: qbp m m qm qbp qm
25060: qbp m m qp qm qbm
12100: qbp m m qp qbm qm
17500: qbp m m qbm qm qp
11020: qbp m m qbm qp qm
9940: qbp m m qbp qm qm
32440: qbp m qm m qm qbp
25960: qbp m qm m qp qbm
19480: qbp m qm m qbm qp
13000: qbp m qm m qbp qm
31360: qbp m qm qm m qbp
5440: qbp m qm qm qbp m
44320: qbp m qm qm qbp p
37840: qbp m qm qm p qbp
23800: qbp m qm qp m qbm
4360: qbp m qm qp qbm m
43240: qbp m qm qp qbm p
30280: qbp m qm qp p qbm
16240: qbp m qm qbm m qp
3280: qbp m qm qbm qp m
42160: qbp m qm qbm qp p
22720: qbp m qm qbm p qp
8680: qbp m qm qbp m qm
2200: qbp m qm qbp qm m
41080: qbp m qm qbp qm p
15160: qbp m qm qbp p qm
33520: qbp m qm p qm qbp
27040: qbp m qm p qp qbm
20560: qbp m qm p qbm qp
14080: qbp m qm p qbp qm
24700: qbp m qp m qm qbm
11740: qbp m qp m qbm qm
23620: qbp m qp qm m qbm
4180: qbp m qp qm qbm m
43060: qbp m qp qm qbm p
30100: qbp m qp qm p qbm
8500: qbp m qp qbm m qm
2020: qbp m qp qbm qm m
40900: qbp m qp qbm qm p
14980: qbp m qp qbm p qm
25780: qbp m qp p qm qbm
12820: qbp m qp p qbm qm
16960: qbp m qbm m qm qp
10480: qbp m qbm m qp qm
15880: qbp m qbm qm m qp
2920: qbp m qbm qm qp m
41800: qbp m qbm qm qp p
22360: qbp m qbm qm p qp
8320: qbp m qbm qp m qm
1840: qbp m qbm qp qm m
40720: qbp m qbm qp qm p
14800: qbp m qbm qp p qm
18040: qbp m qbm p qm qp
11560: qbp m qbm p qp qm
9220: qbp m qbp m qm qm
8140: qbp m qbp qm m qm
1660: qbp m qbp qm qm m
40540: qbp m qbp qm qm p
14620: qbp m qbp qm p qm
10300: qbp m qbp p qm qm
32800: qbp m p qm qm qbp
26320: qbp m p qm qp qbm
19840: qbp m p qm qbm qp
13360: qbp m p qm qbp qm
25240: qbp m p qp qm qbm
12280: qbp m p qp qbm qm
17680: qbp m p qbm qm qp
11200: qbp m p qbm qp qm
10120: qbp m p qbp qm qm
32410: qbp qm m m qm qbp
25930: qbp qm m m qp qbm
19450: qbp qm m m qbm qp
12970: qbp qm m m qbp qm
31330: qbp qm m qm m qbp
5410: qbp qm m qm qbp m
44290: qbp qm m qm qbp p
37810: qbp qm m qm p qbp
23770: qbp qm m qp m qbm
4330: qbp qm m qp qbm m
43210: qbp qm m qp qbm p
30250: qbp qm m qp p qbm
16210: qbp qm m qbm m qp
3250: qbp qm m qbm qp m
42130: qbp qm m qbm qp p
22690: qbp qm m qbm p qp
8650: qbp qm m qbp m qm
2170: qbp qm m qbp qm m
41050: qbp qm m qbp qm p
15130: qbp qm m qbp p qm
33490: qbp qm m p qm qbp
27010: qbp qm m p qp qbm
20530: qbp qm m p qbm qp
14050: qbp qm m p qbp qm
31150: qbp qm qm m m qbp
5230: qbp qm qm m qbp m
44110: qbp qm qm m qbp p
37630: qbp qm qm m p qbp
910: qbp qm qm qbp m m
39790: qbp qm qm qbp m p
7390: qbp qm qm qbp p m
46270: qbp qm qm qbp p p
32230: qbp qm qm p m qbp
6310: qbp qm qm p qbp m
45190: qbp qm qm p qbp p
38710: qbp qm qm p p qbp
23410: qbp qm qp m m qbm
3970: qbp qm qp m qbm m
42850: qbp qm qp m qbm p
29890: qbp qm qp m p qbm
730: qbp qm qp qbm m m
39610: qbp qm qp qbm m p
7210: qbp qm qp qbm p m
46090: qbp qm qp qbm p p
24490: qbp qm qp p m qbm
5050: qbp qm qp p qbm m
43930: qbp qm qp p qbm p
30970: qbp qm qp p p qbm
15670: qbp qm qbm m m qp
2710: qbp qm qbm m qp m
41590: qbp qm qbm m qp p
22150: qbp qm qbm m p qp
550: qbp qm qbm qp m m
39430: qbp qm qbm qp m p
7030: qbp qm qbm qp p m
45910: qbp qm qbm qp p p
16750: qbp qm qbm p m qp
3790: qbp qm qbm p qp m
42670: qbp qm qbm p qp p
23230: qbp qm qbm p p qp
7930: qbp qm qbp m m qm
1450: qbp qm qbp m qm m
40330: qbp qm qbp m qm p
14410: qbp qm qbp m p qm
370: qbp qm qbp qm m m
39250: qbp qm qbp qm m p
6850: qbp qm qbp qm p m
45730: qbp qm qbp qm p p
9010: qbp qm qbp p m qm
2530: qbp qm qbp p qm m
41410: qbp qm qbp p qm p
15490: qbp qm qbp p p qm
32590: qbp qm p m qm qbp
26110: qbp qm p m qp qbm
19630: qbp qm p m qbm qp
13150: qbp qm p m qbp qm
31510: qbp qm p qm m qbp
5590: qbp qm p qm qbp m
44470: qbp qm p qm qbp p
37990: qbp qm p qm p qbp
23950: qbp qm p qp m qbm
4510: qbp qm p qp qbm m
43390: qbp qm p qp qbm p
30430: qbp qm p qp p qbm
16390: qbp qm p qbm m qp
3430: qbp qm p qbm qp m
42310: qbp qm p qbm qp p
22870: qbp qm p qbm p qp
8830: qbp qm p qbp m qm
2350: qbp qm p qbp qm m
41230: qbp qm p qbp qm p
15310: qbp qm p qbp p qm
33670: qbp qm p p qm qbp
27190: qbp qm p p qp qbm
20710: qbp qm p p qbm qp
14230: qbp qm p p qbp qm
24640: qbp qp m m qm qbm
11680: qbp qp m m qbm qm
23560: qbp qp m qm m qbm
4120: qbp qp m qm qbm m
43000: qbp qp m qm qbm p
30040: qbp qp m qm p qbm
8440: qbp qp m qbm m qm
1960: qbp qp m qbm qm m
40840: qbp qp m qbm qm p
14920: qbp qp m qbm p qm
25720: qbp qp m p qm qbm
12760: qbp qp m p qbm qm
23380: qbp qp qm m m qbm
3940: qbp qp qm m qbm m
42820: qbp qp qm m qbm p
29860: qbp qp qm m p qbm
700: qbp qp qm qbm m m
39580: qbp qp qm qbm m p
7180: qbp qp qm qbm p m
46060: qbp qp qm qbm p p
24460: qbp qp qm p m qbm
5020: qbp qp qm p qbm m
43900: qbp qp qm p qbm p
30940: qbp qp qm p p qbm
7900: qbp qp qbm m m qm
1420: qbp qp qbm m qm m
40300: qbp qp qbm m qm p
14380: qbp qp qbm m p qm
340: qbp qp qbm qm m m
39220: qbp qp qbm qm m p
6820: qbp qp qbm qm p m
45700: qbp qp qbm qm p p
8980: qbp qp qbm p m qm
2500: qbp qp qbm p qm m
41380: qbp qp qbm p qm p
15460: qbp qp qbm p p qm
24820: qbp qp p m qm qbm
11860: qbp qp p m qbm qm
23740: qbp qp p qm m qbm
4300: qbp qp p qm qbm m
43180: qbp qp p qm qbm p
30220: qbp qp p qm p qbm
8620: qbp qp p qbm m qm
2140: qbp qp p qbm qm m
41020: qbp qp p qbm qm p
15100: qbp qp p qbm p qm
25900: qbp qp p p qm qbm
12940: qbp qp p p qbm qm
16870: qbp qbm m m qm qp
10390: qbp qbm m m qp qm
15790: qbp qbm m qm m qp
2830: qbp qbm m qm qp m
41710: qbp qbm m qm qp p
22270: qbp qbm m qm p qp
8230: qbp qbm m qp m qm
1750: qbp qbm m qp qm m
40630: qbp qbm m qp qm p
14710: qbp qbm m qp p qm
17950: qbp qbm m p qm qp
11470: qbp qbm m p qp qm
15610: qbp qbm qm m m qp
2650: qbp qbm qm m qp m
41530: qbp qbm qm m qp p
22090: qbp qbm qm m p qp
490: qbp qbm qm qp m m
39370: qbp qbm qm qp m p
6970: qbp qbm qm qp p m
45850: qbp qbm qm qp p p
16690: qbp qbm qm p m qp
3730: qbp qbm qm p qp m
42610: qbp qbm qm p qp p
23170: qbp qbm qm p p qp
7870: qbp qbm qp m m qm
1390: qbp qbm qp m qm m
40270: qbp qbm qp m qm p
14350: qbp qbm qp m p qm
310: qbp qbm qp qm m m
39190: qbp qbm qp qm m p
6790: qbp qbm qp qm p m
45670: qbp qbm qp qm p p
8950: qbp qbm qp p m qm
2470: qbp qbm qp p qm m
41350: qbp qbm qp p qm p
15430: qbp qbm qp p p qm
17050: qbp qbm p m qm qp
10570: qbp qbm p m qp qm
15970: qbp qbm p qm m qp
3010: qbp qbm p qm qp m
41890: qbp qbm p qm qp p
22450: qbp qbm p qm p qp
8410: qbp qbm p qp m qm
1930: qbp qbm p qp qm m
40810: qbp qbm p qp qm p
14890: qbp qbm p qp p qm
18130: qbp qbm p p qm qp
11650: qbp qbm p p qp qm
9100: qbp qbp m m qm qm
8020: qbp qbp m qm m qm
1540: qbp qbp m qm qm m
40420: qbp qbp m qm qm p
14500: qbp qbp m qm p qm
10180: qbp qbp m p qm qm
7840: qbp qbp qm m m qm
1360: qbp qbp qm m qm m
40240: qbp qbp qm m qm p
14320: qbp qbp qm m p qm
280: qbp qbp qm qm m m
39160: qbp qbp qm qm m p
6760: qbp qbp qm qm p m
45640: qbp qbp qm qm p p
8920: qbp qbp qm p m qm
2440: qbp qbp qm p qm m
41320: qbp qbp qm p qm p
15400: qbp qbp qm p p qm
9280: qbp qbp p m qm qm
8200: qbp qbp p qm m qm
1720: qbp qbp p qm qm m
40600: qbp qbp p qm qm p
14680: qbp qbp p qm p qm
10360: qbp qbp p p qm qm
32650: qbp p m qm qm qbp
26170: qbp p m qm qp qbm
19690: qbp p m qm qbm qp
13210: qbp p m qm qbp qm
25090: qbp p m qp qm qbm
12130: qbp p m qp qbm qm
17530: qbp p m qbm qm qp
11050: qbp p m qbm qp qm
9970: qbp p m qbp qm qm
32470: qbp p qm m qm qbp
25990: qbp p qm m qp qbm
19510: qbp p qm m qbm qp
13030: qbp p qm m qbp qm
31390: qbp p qm qm m qbp
5470: qbp p qm qm qbp m
44350: qbp p qm qm qbp p
37870: qbp p qm qm p qbp
23830: qbp p qm qp m qbm
4390: qbp p qm qp qbm m
43270: qbp p qm qp qbm p
30310: qbp p qm qp p qbm
16270: qbp p qm qbm m qp
3310: qbp p qm qbm qp m
42190: qbp p qm qbm qp p
22750: qbp p qm qbm p qp
8710: qbp p qm qbp m qm
2230: qbp p qm qbp qm m
41110: qbp p qm qbp qm p
15190: qbp p qm qbp p qm
33550: qbp p qm p qm qbp
27070: qbp p qm p qp qbm
20590: qbp p qm p qbm qp
14110: qbp p qm p qbp qm
24730: qbp p qp m qm qbm
11770: qbp p qp m qbm qm
23650: qbp p qp qm m qbm
4210: qbp p qp qm qbm m
43090: qbp p qp qm qbm p
30130: qbp p qp qm p qbm
8530: qbp p qp qbm m qm
2050: qbp p qp qbm qm m
40930: qbp p qp qbm qm p
15010: qbp p qp qbm p qm
25810: qbp p qp p qm qbm
12850: qbp p qp p qbm qm
16990: qbp p qbm m qm qp
10510: qbp p qbm m qp qm
15910: qbp p qbm qm m qp
2950: qbp p qbm qm qp m
41830: qbp p qbm qm qp p
22390: qbp p qbm qm p qp
8350: qbp p qbm qp m qm
1870: qbp p qbm qp qm m
40750: qbp p qbm qp qm p
14830: qbp p qbm qp p qm
18070: qbp p qbm p qm qp
11590: qbp p qbm p qp qm
9250: qbp p qbp m qm qm
8170: qbp p qbp qm m qm
1690: qbp p qbp qm qm m
40570: qbp p qbp qm qm p
14650: qbp p qbp qm p qm
10330: qbp p qbp p qm qm
32830: qbp p p qm qm qbp
26350: qbp p p qm qp qbm
19870: qbp p p qm qbm qp
13390: qbp p p qm qbp qm
25270: qbp p p qp qm qbm
12310: qbp p p qp qbm qm
17710: qbp p p qbm qm qp
11230: qbp p p qbm qp qm
10150: qbp p p qbp qm qm
36545: p m qm qm qbp qbp
35465: p m qm qp qbm qbp
28985: p m qm qp qbp qbm
34385: p m qm qbm qp qbp
21425: p m qm qbm qbp qp
33305: p m qm qbp qm qbp
26825: p m qm qbp qp qbm
20345: p m qm qbp qbm qp
13865: p m qm qbp qbp qm
35285: p m qp qm qbm qbp
28805: p m qp qm qbp qbm
27725: p m qp qp qbm qbm
33125: p m qp qbm qm qbp
26645: p m qp qbm qp qbm
20165: p m qp qbm qbm qp
13685: p m qp qbm qbp qm
25565: p m qp qbp qm qbm
12605: p m qp qbp qbm qm
34025: p m qbm qm qp qbp
21065: p m qbm qm qbp qp
32945: p m qbm qp qm qbp
26465: p m qbm qp qp qbm
19985: p m qbm qp qbm qp
13505: p m qbm qp qbp qm
18905: p m qbm qbm qp qp
17825: p m qbm qbp qm qp
11345: p m qbm qbp qp qm
32765: p m qbp qm qm qbp
26285: p m qbp qm qp qbm
19805: p m qbp qm qbm qp
13325: p m qbp qm qbp qm
25205: p m qbp qp qm qbm
12245: p m qbp qp qbm qm
17645: p m qbp qbm qm qp
11165: p m qbp qbm qp qm
10085: p m qbp qbp qm qm
36515: p qm m qm qbp qbp
35435: p qm m qp qbm qbp
28955: p qm m qp qbp qbm
34355: p qm m qbm qp qbp
21395: p qm m qbm qbp qp
33275: p qm m qbp qm qbp
26795: p qm m qbp qp qbm
20315: p qm m qbp qbm qp
13835: p qm m qbp qbp qm
36335: p qm qm m qbp qbp
32015: p qm qm qbp m qbp
6095: p qm qm qbp qbp m
44975: p qm qm qbp qbp p
38495: p qm qm qbp p qbp
37415: p qm qm p qbp qbp
35075: p qm qp m qbm qbp
28595: p qm qp m qbp qbm
31835: p qm qp qbm m qbp
5915: p qm qp qbm qbp m
44795: p qm qp qbm qbp p
38315: p qm qp qbm p qbp
24275: p qm qp qbp m qbm
4835: p qm qp qbp qbm m
43715: p qm qp qbp qbm p
30755: p qm qp qbp p qbm
36155: p qm qp p qbm qbp
29675: p qm qp p qbp qbm
33815: p qm qbm m qp qbp
20855: p qm qbm m qbp qp
31655: p qm qbm qp m qbp
5735: p qm qbm qp qbp m
44615: p qm qbm qp qbp p
38135: p qm qbm qp p qbp
16535: p qm qbm qbp m qp
3575: p qm qbm qbp qp m
42455: p qm qbm qbp qp p
23015: p qm qbm qbp p qp
34895: p qm qbm p qp qbp
21935: p qm qbm p qbp qp
32555: p qm qbp m qm qbp
26075: p qm qbp m qp qbm
19595: p qm qbp m qbm qp
13115: p qm qbp m qbp qm
31475: p qm qbp qm m qbp
5555: p qm qbp qm qbp m
44435: p qm qbp qm qbp p
37955: p qm qbp qm p qbp
23915: p qm qbp qp m qbm
4475: p qm qbp qp qbm m
43355: p qm qbp qp qbm p
30395: p qm qbp qp p qbm
16355: p qm qbp qbm m qp
3395: p qm qbp qbm qp m
42275: p qm qbp qbm qp p
22835: p qm qbp qbm p qp
8795: p qm qbp qbp m qm
2315: p qm qbp qbp qm m
41195: p qm qbp qbp qm p
15275: p qm qbp qbp p qm
33635: p qm qbp p qm qbp
27155: p qm qbp p qp qbm
20675: p qm qbp p qbm qp
14195: p qm qbp p qbp qm
36695: p qm p qm qbp qbp
35615: p qm p qp qbm qbp
29135: p qm p qp qbp qbm
34535: p qm p qbm qp qbp
21575: p qm p qbm qbp qp
33455: p qm p qbp qm qbp
26975: p qm p qbp qp qbm
20495: p qm p qbp qbm qp
14015: p qm p qbp qbp qm
35225: p qp m qm qbm qbp
28745: p qp m qm qbp qbm
27665: p qp m qp qbm qbm
33065: p qp m qbm qm qbp
26585: p qp m qbm qp qbm
20105: p qp m qbm qbm qp
13625: p qp m qbm qbp qm
25505: p qp m qbp qm qbm
12545: p qp m qbp qbm qm
35045: p qp qm m qbm qbp
28565: p qp qm m qbp qbm
31805: p qp qm qbm m qbp
5885: p qp qm qbm qbp m
44765: p qp qm qbm qbp p
38285: p qp qm qbm p qbp
24245: p qp qm qbp m qbm
4805: p qp qm qbp qbm m
43685: p qp qm qbp qbm p
30725: p qp qm qbp p qbm
36125: p qp qm p qbm qbp
29645: p qp qm p qbp qbm
27305: p qp qp m qbm qbm
24065: p qp qp qbm m qbm
4625: p qp qp qbm qbm m
43505: p qp qp qbm qbm p
30545: p qp qp qbm p qbm
28385: p qp qp p qbm qbm
32525: p qp qbm m qm qbp
26045: p qp qbm m qp qbm
19565: p qp qbm m qbm qp
13085: p qp qbm m qbp qm
31445: p qp qbm qm m qbp
5525: p qp qbm qm qbp m
44405: p qp qbm qm qbp p
37925: p qp qbm qm p qbp
23885: p qp qbm qp m qbm
4445: p qp qbm qp qbm m
43325: p qp qbm qp qbm p
30365: p qp qbm qp p qbm
16325: p qp qbm qbm m qp
3365: p qp qbm qbm qp m
42245: p qp qbm qbm qp p
22805: p qp qbm qbm p qp
8765: p qp qbm qbp m qm
2285: p qp qbm qbp qm m
41165: p qp qbm qbp qm p
15245: p qp qbm qbp p qm
33605: p qp qbm p qm qbp
27125: p qp qbm p qp qbm
20645: p qp qbm p qbm qp
14165: p qp qbm p qbp qm
24785: p qp qbp m qm qbm
11825: p qp qbp m qbm qm
23705: p qp qbp qm m qbm
4265: p qp qbp qm qbm m
43145: p qp qbp qm qbm p
30185: p qp qbp qm p qbm
8585: p qp qbp qbm m qm
2105: p qp qbp qbm qm m
40985: p qp qbp qbm qm p
15065: p qp qbp qbm p qm
25865: p qp qbp p qm qbm
12905: p qp qbp p qbm qm
35405: p qp p qm qbm qbp
28925: p qp p qm qbp qbm
27845: p qp p qp qbm qbm
33245: p qp p qbm qm qbp
26765: p qp p qbm qp qbm
20285: p qp p qbm qbm qp
13805: p qp p qbm qbp qm
25685: p qp p qbp qm qbm
12725: p qp p qbp qbm qm
33935: p qbm m qm qp qbp
20975: p qbm m qm qbp qp
32855: p qbm m qp qm qbp
26375: p qbm m qp qp qbm
19895: p qbm m qp qbm qp
13415: p qbm m qp qbp qm
18815: p qbm m qbm qp qp
17735: p qbm m qbp qm qp
11255: p qbm m qbp qp qm
33755: p qbm qm m qp qbp
20795: p qbm qm m qbp qp
31595: p qbm qm qp m qbp
5675: p qbm qm qp qbp m
44555: p qbm qm qp qbp p
38075: p qbm qm qp p qbp
16475: p qbm qm qbp m qp
3515: p qbm qm qbp qp m
42395: p qbm qm qbp qp p
22955: p qbm qm qbp p qp
34835: p qbm qm p qp qbp
21875: p qbm qm p qbp qp
32495: p qbm qp m qm qbp
26015: p qbm qp m qp qbm
19535: p qbm qp m qbm qp
13055: p qbm qp m qbp qm
31415: p qbm qp qm m qbp
5495: p qbm qp qm qbp m
44375: p qbm qp qm qbp p
37895: p qbm qp qm p qbp
23855: p qbm qp qp m qbm
4415: p qbm qp qp qbm m
43295: p qbm qp qp qbm p
30335: p qbm qp qp p qbm
16295: p qbm qp qbm m qp
3335: p qbm qp qbm qp m
42215: p qbm qp qbm qp p
22775: p qbm qp qbm p qp
8735: p qbm qp qbp m qm
2255: p qbm qp qbp qm m
41135: p qbm qp qbp qm p
15215: p qbm qp qbp p qm
33575: p qbm qp p qm qbp
27095: p qbm qp p qp qbm
20615: p qbm qp p qbm qp
14135: p qbm qp p qbp qm
18275: p qbm qbm m qp qp
16115: p qbm qbm qp m qp
3155: p qbm qbm qp qp m
42035: p qbm qbm qp qp p
22595: p qbm qbm qp p qp
19355: p qbm qbm p qp qp
17015: p qbm qbp m qm qp
10535: p qbm qbp m qp qm
15935: p qbm qbp qm m qp
2975: p qbm qbp qm qp m
41855: p qbm qbp qm qp p
22415: p qbm qbp qm p qp
8375: p qbm qbp qp m qm
1895: p qbm qbp qp qm m
40775: p qbm qbp qp qm p
14855: p qbm qbp qp p qm
18095: p qbm qbp p qm qp
11615: p qbm qbp p qp qm
34115: p qbm p qm qp qbp
21155: p qbm p qm qbp qp
33035: p qbm p qp qm qbp
26555: p qbm p qp qp qbm
20075: p qbm p qp qbm qp
13595: p qbm p qp qbp qm
18995: p qbm p qbm qp qp
17915: p qbm p qbp qm qp
11435: p qbm p qbp qp qm
32645: p qbp m qm qm qbp
26165: p qbp m qm qp qbm
19685: p qbp m qm qbm qp
13205: p qbp m qm qbp qm
25085: p qbp m qp qm qbm
12125: p qbp m qp qbm qm
17525: p qbp m qbm qm qp
11045: p qbp m qbm qp qm
9965: p qbp m qbp qm qm
32465: p qbp qm m qm qbp
25985: p qbp qm m qp qbm
19505: p qbp qm m qbm qp
13025: p qbp qm m qbp qm
31385: p qbp qm qm m qbp
5465: p qbp qm qm qbp m
44345: p qbp qm qm qbp p
37865: p qbp qm qm p qbp
23825: p qbp qm qp m qbm
4385: p qbp qm qp qbm m
43265: p qbp qm qp qbm p
30305: p qbp qm qp p qbm
16265: p qbp qm qbm m qp
3305: p qbp qm qbm qp m
42185: p qbp qm qbm qp p
22745: p qbp qm qbm p qp
8705: p qbp qm qbp m qm
2225: p qbp qm qbp qm m
41105: p qbp qm qbp qm p
15185: p qbp qm qbp p qm
33545: p qbp qm p qm qbp
27065: p qbp qm p qp qbm
20585: p qbp qm p qbm qp
14105: p qbp qm p qbp qm
24725: p qbp qp m qm qbm
11765: p qbp qp m qbm qm
23645: p qbp qp qm m qbm
4205: p qbp qp qm qbm m
43085: p qbp qp qm qbm p
30125: p qbp qp qm p qbm
8525: p qbp qp qbm m qm
2045: p qbp qp qbm qm m
40925: p qbp qp qbm qm p
15005: p qbp qp qbm p qm
25805: p qbp qp p qm qbm
12845: p qbp qp p qbm qm
16985: p qbp qbm m qm qp
10505: p qbp qbm m qp qm
15905: p qbp qbm qm m qp
2945: p qbp qbm qm qp m
41825: p qbp qbm qm qp p
22385: p qbp qbm qm p qp
8345: p qbp qbm qp m qm
1865: p qbp qbm qp qm m
40745: p qbp qbm qp qm p
14825: p qbp qbm qp p qm
18065: p qbp qbm p qm qp
11585: p qbp qbm p qp qm
9245: p qbp qbp m qm qm
8165: p qbp qbp qm m qm
1685: p qbp qbp qm qm m
40565: p qbp qbp qm qm p
14645: p qbp qbp qm p qm
10325: p qbp qbp p qm qm
32825: p qbp p qm qm qbp
26345: p qbp p qm qp qbm
19865: p qbp p qm qbm qp
13385: p qbp p qm qbp qm
25265: p qbp p qp qm qbm
12305: p qbp p qp qbm qm
17705: p qbp p qbm qm qp
11225: p qbp p qbm qp qm
10145: p qbp p qbp qm qm
36575: p p qm qm qbp qbp
35495: p p qm qp qbm qbp
29015: p p qm qp qbp qbm
34415: p p qm qbm qp qbp
21455: p p qm qbm qbp qp
33335: p p qm qbp qm qbp
26855: p p qm qbp qp qbm
20375: p p qm qbp qbm qp
13895: p p qm qbp qbp qm
35315: p p qp qm qbm qbp
28835: p p qp qm qbp qbm
27755: p p qp qp qbm qbm
33155: p p qp qbm qm qbp
26675: p p qp qbm qp qbm
20195: p p qp qbm qbm qp
13715: p p qp qbm qbp qm
25595: p p qp qbp qm qbm
12635: p p qp qbp qbm qm
34055: p p qbm qm qp qbp
21095: p p qbm qm qbp qp
32975: p p qbm qp qm qbp
26495: p p qbm qp qp qbm
20015: p p qbm qp qbm qp
13535: p p qbm qp qbp qm
18935: p p qbm qbm qp qp
17855: p p qbm qbp qm qp
11375: p p qbm qbp qp qm
32795: p p qbp qm qm qbp
26315: p p qbp qm qp qbm
19835: p p qbp qm qbm qp
13355: p p qbp qm qbp qm
25235: p p qbp qp qm qbm
12275: p p qbp qp qbm qm
17675: p p qbp qbm qm qp
11195: p p qbp qbm qp qm
10115: p p qbp qbp qm qm
*/

