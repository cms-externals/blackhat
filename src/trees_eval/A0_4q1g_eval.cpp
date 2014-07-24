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

namespace BH {


template <class T> complex<T> A4q1gzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A4q1g280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
(ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1455_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g1540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1660_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1680_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1690_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1750_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1875_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1890_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g1935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g1960_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g1970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g2040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g2045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g2050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g2090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2100_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g2150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2170_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2175_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2185_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2530_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g2650_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2715_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3))-
 (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
(ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g2755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g2920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g2955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g2970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g2975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g2985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3155_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3165_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0))-
 (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3365_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))-
 (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g3515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3570_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g3655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g3835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g3950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
(ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g3985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4130_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4200_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4260_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4330_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g4475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g4550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4610_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4620_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4640_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g4805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g4940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g4945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5060_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5230_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5240_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5245_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
(ep.spb(2,1)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5275_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5335_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5410_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5420_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5440_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5520_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5545_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0))+
 (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
(ep.spb(1,0)*ep.spb(2,1)*
 ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A4q1g5555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(0,4)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
 ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5590_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5605_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(2,3))+
 (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
(ep.spa(0,1)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g5625_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5635_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5670_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5725_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g5755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g5805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g5840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5845_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5880_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g5905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g5915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g5935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q1g6020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g6055_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g6085_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A4q1g6090_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A4q1g6095_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g6115_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g6235_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A4q1g6310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6325_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,1)*
 ep.spa(2,3)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
(ep.spa(0,4)*ep.spa(1,2)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6760_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6850_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g6970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g6975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7035_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))-
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7040_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g7180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g7190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7215_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7220_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7250_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g7390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7405_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
 ep.spa(1,2)*ep.spa(3,4))+
 (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
(ep.spa(0,1)*ep.spa(0,4)*
 ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g7435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q1g7465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A4q1g7495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
 ); }


template <class T> complex<T>  (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 280 : return &A4q1g280_eval;
case 310 : return &A4q1g310_eval;
case 315 : return &A4q1g315_eval;
case 340 : return &A4q1g340_eval;
case 350 : return &A4q1g350_eval;
case 370 : return &A4q1g370_eval;
case 375 : return &A4q1g375_eval;
case 380 : return &A4q1g380_eval;
case 385 : return &A4q1g385_eval;
case 490 : return &A4q1g490_eval;
case 495 : return &A4q1g495_eval;
case 525 : return &A4q1g525_eval;
case 550 : return &A4q1g550_eval;
case 555 : return &A4q1g555_eval;
case 560 : return &A4q1g560_eval;
case 565 : return &A4q1g565_eval;
case 585 : return &A4q1g585_eval;
case 595 : return &A4q1g595_eval;
case 700 : return &A4q1g700_eval;
case 710 : return &A4q1g710_eval;
case 730 : return &A4q1g730_eval;
case 735 : return &A4q1g735_eval;
case 740 : return &A4q1g740_eval;
case 745 : return &A4q1g745_eval;
case 770 : return &A4q1g770_eval;
case 800 : return &A4q1g800_eval;
case 805 : return &A4q1g805_eval;
case 910 : return &A4q1g910_eval;
case 915 : return &A4q1g915_eval;
case 920 : return &A4q1g920_eval;
case 925 : return &A4q1g925_eval;
case 945 : return &A4q1g945_eval;
case 955 : return &A4q1g955_eval;
case 980 : return &A4q1g980_eval;
case 985 : return &A4q1g985_eval;
case 1015 : return &A4q1g1015_eval;
case 1360 : return &A4q1g1360_eval;
case 1390 : return &A4q1g1390_eval;
case 1395 : return &A4q1g1395_eval;
case 1420 : return &A4q1g1420_eval;
case 1430 : return &A4q1g1430_eval;
case 1450 : return &A4q1g1450_eval;
case 1455 : return &A4q1g1455_eval;
case 1460 : return &A4q1g1460_eval;
case 1465 : return &A4q1g1465_eval;
case 1540 : return &A4q1g1540_eval;
case 1660 : return &A4q1g1660_eval;
case 1680 : return &A4q1g1680_eval;
case 1685 : return &A4q1g1685_eval;
case 1690 : return &A4q1g1690_eval;
case 1720 : return &A4q1g1720_eval;
case 1750 : return &A4q1g1750_eval;
case 1755 : return &A4q1g1755_eval;
case 1840 : return &A4q1g1840_eval;
case 1860 : return &A4q1g1860_eval;
case 1865 : return &A4q1g1865_eval;
case 1870 : return &A4q1g1870_eval;
case 1875 : return &A4q1g1875_eval;
case 1890 : return &A4q1g1890_eval;
case 1895 : return &A4q1g1895_eval;
case 1905 : return &A4q1g1905_eval;
case 1930 : return &A4q1g1930_eval;
case 1935 : return &A4q1g1935_eval;
case 1960 : return &A4q1g1960_eval;
case 1970 : return &A4q1g1970_eval;
case 2020 : return &A4q1g2020_eval;
case 2040 : return &A4q1g2040_eval;
case 2045 : return &A4q1g2045_eval;
case 2050 : return &A4q1g2050_eval;
case 2090 : return &A4q1g2090_eval;
case 2100 : return &A4q1g2100_eval;
case 2105 : return &A4q1g2105_eval;
case 2120 : return &A4q1g2120_eval;
case 2140 : return &A4q1g2140_eval;
case 2150 : return &A4q1g2150_eval;
case 2170 : return &A4q1g2170_eval;
case 2175 : return &A4q1g2175_eval;
case 2180 : return &A4q1g2180_eval;
case 2185 : return &A4q1g2185_eval;
case 2200 : return &A4q1g2200_eval;
case 2220 : return &A4q1g2220_eval;
case 2225 : return &A4q1g2225_eval;
case 2230 : return &A4q1g2230_eval;
case 2235 : return &A4q1g2235_eval;
case 2250 : return &A4q1g2250_eval;
case 2255 : return &A4q1g2255_eval;
case 2265 : return &A4q1g2265_eval;
case 2270 : return &A4q1g2270_eval;
case 2280 : return &A4q1g2280_eval;
case 2285 : return &A4q1g2285_eval;
case 2300 : return &A4q1g2300_eval;
case 2305 : return &A4q1g2305_eval;
case 2310 : return &A4q1g2310_eval;
case 2315 : return &A4q1g2315_eval;
case 2335 : return &A4q1g2335_eval;
case 2350 : return &A4q1g2350_eval;
case 2355 : return &A4q1g2355_eval;
case 2360 : return &A4q1g2360_eval;
case 2365 : return &A4q1g2365_eval;
case 2440 : return &A4q1g2440_eval;
case 2470 : return &A4q1g2470_eval;
case 2475 : return &A4q1g2475_eval;
case 2500 : return &A4q1g2500_eval;
case 2510 : return &A4q1g2510_eval;
case 2530 : return &A4q1g2530_eval;
case 2535 : return &A4q1g2535_eval;
case 2540 : return &A4q1g2540_eval;
case 2545 : return &A4q1g2545_eval;
case 2650 : return &A4q1g2650_eval;
case 2655 : return &A4q1g2655_eval;
case 2685 : return &A4q1g2685_eval;
case 2710 : return &A4q1g2710_eval;
case 2715 : return &A4q1g2715_eval;
case 2720 : return &A4q1g2720_eval;
case 2725 : return &A4q1g2725_eval;
case 2745 : return &A4q1g2745_eval;
case 2755 : return &A4q1g2755_eval;
case 2830 : return &A4q1g2830_eval;
case 2835 : return &A4q1g2835_eval;
case 2920 : return &A4q1g2920_eval;
case 2940 : return &A4q1g2940_eval;
case 2945 : return &A4q1g2945_eval;
case 2950 : return &A4q1g2950_eval;
case 2955 : return &A4q1g2955_eval;
case 2970 : return &A4q1g2970_eval;
case 2975 : return &A4q1g2975_eval;
case 2985 : return &A4q1g2985_eval;
case 3010 : return &A4q1g3010_eval;
case 3015 : return &A4q1g3015_eval;
case 3045 : return &A4q1g3045_eval;
case 3135 : return &A4q1g3135_eval;
case 3150 : return &A4q1g3150_eval;
case 3155 : return &A4q1g3155_eval;
case 3165 : return &A4q1g3165_eval;
case 3225 : return &A4q1g3225_eval;
case 3250 : return &A4q1g3250_eval;
case 3255 : return &A4q1g3255_eval;
case 3260 : return &A4q1g3260_eval;
case 3265 : return &A4q1g3265_eval;
case 3280 : return &A4q1g3280_eval;
case 3300 : return &A4q1g3300_eval;
case 3305 : return &A4q1g3305_eval;
case 3310 : return &A4q1g3310_eval;
case 3315 : return &A4q1g3315_eval;
case 3330 : return &A4q1g3330_eval;
case 3335 : return &A4q1g3335_eval;
case 3345 : return &A4q1g3345_eval;
case 3350 : return &A4q1g3350_eval;
case 3360 : return &A4q1g3360_eval;
case 3365 : return &A4q1g3365_eval;
case 3380 : return &A4q1g3380_eval;
case 3385 : return &A4q1g3385_eval;
case 3390 : return &A4q1g3390_eval;
case 3395 : return &A4q1g3395_eval;
case 3415 : return &A4q1g3415_eval;
case 3430 : return &A4q1g3430_eval;
case 3435 : return &A4q1g3435_eval;
case 3440 : return &A4q1g3440_eval;
case 3445 : return &A4q1g3445_eval;
case 3465 : return &A4q1g3465_eval;
case 3475 : return &A4q1g3475_eval;
case 3495 : return &A4q1g3495_eval;
case 3510 : return &A4q1g3510_eval;
case 3515 : return &A4q1g3515_eval;
case 3525 : return &A4q1g3525_eval;
case 3565 : return &A4q1g3565_eval;
case 3570 : return &A4q1g3570_eval;
case 3575 : return &A4q1g3575_eval;
case 3595 : return &A4q1g3595_eval;
case 3645 : return &A4q1g3645_eval;
case 3655 : return &A4q1g3655_eval;
case 3730 : return &A4q1g3730_eval;
case 3735 : return &A4q1g3735_eval;
case 3765 : return &A4q1g3765_eval;
case 3790 : return &A4q1g3790_eval;
case 3795 : return &A4q1g3795_eval;
case 3800 : return &A4q1g3800_eval;
case 3805 : return &A4q1g3805_eval;
case 3825 : return &A4q1g3825_eval;
case 3835 : return &A4q1g3835_eval;
case 3940 : return &A4q1g3940_eval;
case 3950 : return &A4q1g3950_eval;
case 3970 : return &A4q1g3970_eval;
case 3975 : return &A4q1g3975_eval;
case 3980 : return &A4q1g3980_eval;
case 3985 : return &A4q1g3985_eval;
case 4010 : return &A4q1g4010_eval;
case 4040 : return &A4q1g4040_eval;
case 4045 : return &A4q1g4045_eval;
case 4120 : return &A4q1g4120_eval;
case 4130 : return &A4q1g4130_eval;
case 4180 : return &A4q1g4180_eval;
case 4200 : return &A4q1g4200_eval;
case 4205 : return &A4q1g4205_eval;
case 4210 : return &A4q1g4210_eval;
case 4250 : return &A4q1g4250_eval;
case 4260 : return &A4q1g4260_eval;
case 4265 : return &A4q1g4265_eval;
case 4280 : return &A4q1g4280_eval;
case 4300 : return &A4q1g4300_eval;
case 4310 : return &A4q1g4310_eval;
case 4330 : return &A4q1g4330_eval;
case 4335 : return &A4q1g4335_eval;
case 4340 : return &A4q1g4340_eval;
case 4345 : return &A4q1g4345_eval;
case 4360 : return &A4q1g4360_eval;
case 4380 : return &A4q1g4380_eval;
case 4385 : return &A4q1g4385_eval;
case 4390 : return &A4q1g4390_eval;
case 4395 : return &A4q1g4395_eval;
case 4410 : return &A4q1g4410_eval;
case 4415 : return &A4q1g4415_eval;
case 4425 : return &A4q1g4425_eval;
case 4430 : return &A4q1g4430_eval;
case 4440 : return &A4q1g4440_eval;
case 4445 : return &A4q1g4445_eval;
case 4460 : return &A4q1g4460_eval;
case 4465 : return &A4q1g4465_eval;
case 4470 : return &A4q1g4470_eval;
case 4475 : return &A4q1g4475_eval;
case 4495 : return &A4q1g4495_eval;
case 4510 : return &A4q1g4510_eval;
case 4515 : return &A4q1g4515_eval;
case 4520 : return &A4q1g4520_eval;
case 4525 : return &A4q1g4525_eval;
case 4550 : return &A4q1g4550_eval;
case 4610 : return &A4q1g4610_eval;
case 4620 : return &A4q1g4620_eval;
case 4625 : return &A4q1g4625_eval;
case 4640 : return &A4q1g4640_eval;
case 4730 : return &A4q1g4730_eval;
case 4760 : return &A4q1g4760_eval;
case 4765 : return &A4q1g4765_eval;
case 4790 : return &A4q1g4790_eval;
case 4800 : return &A4q1g4800_eval;
case 4805 : return &A4q1g4805_eval;
case 4820 : return &A4q1g4820_eval;
case 4825 : return &A4q1g4825_eval;
case 4830 : return &A4q1g4830_eval;
case 4835 : return &A4q1g4835_eval;
case 4855 : return &A4q1g4855_eval;
case 4940 : return &A4q1g4940_eval;
case 4945 : return &A4q1g4945_eval;
case 5020 : return &A4q1g5020_eval;
case 5030 : return &A4q1g5030_eval;
case 5050 : return &A4q1g5050_eval;
case 5055 : return &A4q1g5055_eval;
case 5060 : return &A4q1g5060_eval;
case 5065 : return &A4q1g5065_eval;
case 5090 : return &A4q1g5090_eval;
case 5120 : return &A4q1g5120_eval;
case 5125 : return &A4q1g5125_eval;
case 5230 : return &A4q1g5230_eval;
case 5235 : return &A4q1g5235_eval;
case 5240 : return &A4q1g5240_eval;
case 5245 : return &A4q1g5245_eval;
case 5265 : return &A4q1g5265_eval;
case 5275 : return &A4q1g5275_eval;
case 5300 : return &A4q1g5300_eval;
case 5305 : return &A4q1g5305_eval;
case 5335 : return &A4q1g5335_eval;
case 5410 : return &A4q1g5410_eval;
case 5415 : return &A4q1g5415_eval;
case 5420 : return &A4q1g5420_eval;
case 5425 : return &A4q1g5425_eval;
case 5440 : return &A4q1g5440_eval;
case 5460 : return &A4q1g5460_eval;
case 5465 : return &A4q1g5465_eval;
case 5470 : return &A4q1g5470_eval;
case 5475 : return &A4q1g5475_eval;
case 5490 : return &A4q1g5490_eval;
case 5495 : return &A4q1g5495_eval;
case 5505 : return &A4q1g5505_eval;
case 5510 : return &A4q1g5510_eval;
case 5520 : return &A4q1g5520_eval;
case 5525 : return &A4q1g5525_eval;
case 5540 : return &A4q1g5540_eval;
case 5545 : return &A4q1g5545_eval;
case 5550 : return &A4q1g5550_eval;
case 5555 : return &A4q1g5555_eval;
case 5575 : return &A4q1g5575_eval;
case 5590 : return &A4q1g5590_eval;
case 5595 : return &A4q1g5595_eval;
case 5600 : return &A4q1g5600_eval;
case 5605 : return &A4q1g5605_eval;
case 5625 : return &A4q1g5625_eval;
case 5635 : return &A4q1g5635_eval;
case 5655 : return &A4q1g5655_eval;
case 5670 : return &A4q1g5670_eval;
case 5675 : return &A4q1g5675_eval;
case 5685 : return &A4q1g5685_eval;
case 5725 : return &A4q1g5725_eval;
case 5730 : return &A4q1g5730_eval;
case 5735 : return &A4q1g5735_eval;
case 5755 : return &A4q1g5755_eval;
case 5805 : return &A4q1g5805_eval;
case 5815 : return &A4q1g5815_eval;
case 5840 : return &A4q1g5840_eval;
case 5845 : return &A4q1g5845_eval;
case 5870 : return &A4q1g5870_eval;
case 5880 : return &A4q1g5880_eval;
case 5885 : return &A4q1g5885_eval;
case 5900 : return &A4q1g5900_eval;
case 5905 : return &A4q1g5905_eval;
case 5910 : return &A4q1g5910_eval;
case 5915 : return &A4q1g5915_eval;
case 5935 : return &A4q1g5935_eval;
case 6020 : return &A4q1g6020_eval;
case 6025 : return &A4q1g6025_eval;
case 6055 : return &A4q1g6055_eval;
case 6085 : return &A4q1g6085_eval;
case 6090 : return &A4q1g6090_eval;
case 6095 : return &A4q1g6095_eval;
case 6115 : return &A4q1g6115_eval;
case 6235 : return &A4q1g6235_eval;
case 6310 : return &A4q1g6310_eval;
case 6315 : return &A4q1g6315_eval;
case 6320 : return &A4q1g6320_eval;
case 6325 : return &A4q1g6325_eval;
case 6345 : return &A4q1g6345_eval;
case 6355 : return &A4q1g6355_eval;
case 6380 : return &A4q1g6380_eval;
case 6385 : return &A4q1g6385_eval;
case 6415 : return &A4q1g6415_eval;
case 6760 : return &A4q1g6760_eval;
case 6790 : return &A4q1g6790_eval;
case 6795 : return &A4q1g6795_eval;
case 6820 : return &A4q1g6820_eval;
case 6830 : return &A4q1g6830_eval;
case 6850 : return &A4q1g6850_eval;
case 6855 : return &A4q1g6855_eval;
case 6860 : return &A4q1g6860_eval;
case 6865 : return &A4q1g6865_eval;
case 6970 : return &A4q1g6970_eval;
case 6975 : return &A4q1g6975_eval;
case 7005 : return &A4q1g7005_eval;
case 7030 : return &A4q1g7030_eval;
case 7035 : return &A4q1g7035_eval;
case 7040 : return &A4q1g7040_eval;
case 7045 : return &A4q1g7045_eval;
case 7065 : return &A4q1g7065_eval;
case 7075 : return &A4q1g7075_eval;
case 7180 : return &A4q1g7180_eval;
case 7190 : return &A4q1g7190_eval;
case 7210 : return &A4q1g7210_eval;
case 7215 : return &A4q1g7215_eval;
case 7220 : return &A4q1g7220_eval;
case 7225 : return &A4q1g7225_eval;
case 7250 : return &A4q1g7250_eval;
case 7280 : return &A4q1g7280_eval;
case 7285 : return &A4q1g7285_eval;
case 7390 : return &A4q1g7390_eval;
case 7395 : return &A4q1g7395_eval;
case 7400 : return &A4q1g7400_eval;
case 7405 : return &A4q1g7405_eval;
case 7425 : return &A4q1g7425_eval;
case 7435 : return &A4q1g7435_eval;
case 7460 : return &A4q1g7460_eval;
case 7465 : return &A4q1g7465_eval;
case 7495 : return &A4q1g7495_eval;

default: return 0;
		throw BHerror("case missing for tree amplitude!");

}
}



template complex<R>  (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A4q1g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
6090: m qm qm qbp qbp
5910: m qm qp qbm qbp
4830: m qm qp qbp qbm
5730: m qm qbm qp qbp
3570: m qm qbm qbp qp
5550: m qm qbp qm qbp
4470: m qm qbp qp qbm
3390: m qm qbp qbm qp
2310: m qm qbp qbp qm
5880: m qp qm qbm qbp
4800: m qp qm qbp qbm
4620: m qp qp qbm qbm
5520: m qp qbm qm qbp
4440: m qp qbm qp qbm
3360: m qp qbm qbm qp
2280: m qp qbm qbp qm
4260: m qp qbp qm qbm
2100: m qp qbp qbm qm
5670: m qbm qm qp qbp
3510: m qbm qm qbp qp
5490: m qbm qp qm qbp
4410: m qbm qp qp qbm
3330: m qbm qp qbm qp
2250: m qbm qp qbp qm
3150: m qbm qbm qp qp
2970: m qbm qbp qm qp
1890: m qbm qbp qp qm
5460: m qbp qm qm qbp
4380: m qbp qm qp qbm
3300: m qbp qm qbm qp
2220: m qbp qm qbp qm
4200: m qbp qp qm qbm
2040: m qbp qp qbm qm
2940: m qbp qbm qm qp
1860: m qbp qbm qp qm
1680: m qbp qbp qm qm
6085: qm m qm qbp qbp
5905: qm m qp qbm qbp
4825: qm m qp qbp qbm
5725: qm m qbm qp qbp
3565: qm m qbm qbp qp
5545: qm m qbp qm qbp
4465: qm m qbp qp qbm
3385: qm m qbp qbm qp
2305: qm m qbp qbp qm
6055: qm qm m qbp qbp
5335: qm qm qbp m qbp
1015: qm qm qbp qbp m
7495: qm qm qbp qbp p
6415: qm qm qbp p qbp
6235: qm qm p qbp qbp
5845: qm qp m qbm qbp
4765: qm qp m qbp qbm
5305: qm qp qbm m qbp
985: qm qp qbm qbp m
7465: qm qp qbm qbp p
6385: qm qp qbm p qbp
4045: qm qp qbp m qbm
805: qm qp qbp qbm m
7285: qm qp qbp qbm p
5125: qm qp qbp p qbm
6025: qm qp p qbm qbp
4945: qm qp p qbp qbm
5635: qm qbm m qp qbp
3475: qm qbm m qbp qp
5275: qm qbm qp m qbp
955: qm qbm qp qbp m
7435: qm qbm qp qbp p
6355: qm qbm qp p qbp
2755: qm qbm qbp m qp
595: qm qbm qbp qp m
7075: qm qbm qbp qp p
3835: qm qbm qbp p qp
5815: qm qbm p qp qbp
3655: qm qbm p qbp qp
5425: qm qbp m qm qbp
4345: qm qbp m qp qbm
3265: qm qbp m qbm qp
2185: qm qbp m qbp qm
5245: qm qbp qm m qbp
925: qm qbp qm qbp m
7405: qm qbp qm qbp p
6325: qm qbp qm p qbp
3985: qm qbp qp m qbm
745: qm qbp qp qbm m
7225: qm qbp qp qbm p
5065: qm qbp qp p qbm
2725: qm qbp qbm m qp
565: qm qbp qbm qp m
7045: qm qbp qbm qp p
3805: qm qbp qbm p qp
1465: qm qbp qbp m qm
385: qm qbp qbp qm m
6865: qm qbp qbp qm p
2545: qm qbp qbp p qm
5605: qm qbp p qm qbp
4525: qm qbp p qp qbm
3445: qm qbp p qbm qp
2365: qm qbp p qbp qm
6115: qm p qm qbp qbp
5935: qm p qp qbm qbp
4855: qm p qp qbp qbm
5755: qm p qbm qp qbp
3595: qm p qbm qbp qp
5575: qm p qbp qm qbp
4495: qm p qbp qp qbm
3415: qm p qbp qbm qp
2335: qm p qbp qbp qm
5870: qp m qm qbm qbp
4790: qp m qm qbp qbm
4610: qp m qp qbm qbm
5510: qp m qbm qm qbp
4430: qp m qbm qp qbm
3350: qp m qbm qbm qp
2270: qp m qbm qbp qm
4250: qp m qbp qm qbm
2090: qp m qbp qbm qm
5840: qp qm m qbm qbp
4760: qp qm m qbp qbm
5300: qp qm qbm m qbp
980: qp qm qbm qbp m
7460: qp qm qbm qbp p
6380: qp qm qbm p qbp
4040: qp qm qbp m qbm
800: qp qm qbp qbm m
7280: qp qm qbp qbm p
5120: qp qm qbp p qbm
6020: qp qm p qbm qbp
4940: qp qm p qbp qbm
4550: qp qp m qbm qbm
4010: qp qp qbm m qbm
770: qp qp qbm qbm m
7250: qp qp qbm qbm p
5090: qp qp qbm p qbm
4730: qp qp p qbm qbm
5420: qp qbm m qm qbp
4340: qp qbm m qp qbm
3260: qp qbm m qbm qp
2180: qp qbm m qbp qm
5240: qp qbm qm m qbp
920: qp qbm qm qbp m
7400: qp qbm qm qbp p
6320: qp qbm qm p qbp
3980: qp qbm qp m qbm
740: qp qbm qp qbm m
7220: qp qbm qp qbm p
5060: qp qbm qp p qbm
2720: qp qbm qbm m qp
560: qp qbm qbm qp m
7040: qp qbm qbm qp p
3800: qp qbm qbm p qp
1460: qp qbm qbp m qm
380: qp qbm qbp qm m
6860: qp qbm qbp qm p
2540: qp qbm qbp p qm
5600: qp qbm p qm qbp
4520: qp qbm p qp qbm
3440: qp qbm p qbm qp
2360: qp qbm p qbp qm
4130: qp qbp m qm qbm
1970: qp qbp m qbm qm
3950: qp qbp qm m qbm
710: qp qbp qm qbm m
7190: qp qbp qm qbm p
5030: qp qbp qm p qbm
1430: qp qbp qbm m qm
350: qp qbp qbm qm m
6830: qp qbp qbm qm p
2510: qp qbp qbm p qm
4310: qp qbp p qm qbm
2150: qp qbp p qbm qm
5900: qp p qm qbm qbp
4820: qp p qm qbp qbm
4640: qp p qp qbm qbm
5540: qp p qbm qm qbp
4460: qp p qbm qp qbm
3380: qp p qbm qbm qp
2300: qp p qbm qbp qm
4280: qp p qbp qm qbm
2120: qp p qbp qbm qm
5655: qbm m qm qp qbp
3495: qbm m qm qbp qp
5475: qbm m qp qm qbp
4395: qbm m qp qp qbm
3315: qbm m qp qbm qp
2235: qbm m qp qbp qm
3135: qbm m qbm qp qp
2955: qbm m qbp qm qp
1875: qbm m qbp qp qm
5625: qbm qm m qp qbp
3465: qbm qm m qbp qp
5265: qbm qm qp m qbp
945: qbm qm qp qbp m
7425: qbm qm qp qbp p
6345: qbm qm qp p qbp
2745: qbm qm qbp m qp
585: qbm qm qbp qp m
7065: qbm qm qbp qp p
3825: qbm qm qbp p qp
5805: qbm qm p qp qbp
3645: qbm qm p qbp qp
5415: qbm qp m qm qbp
4335: qbm qp m qp qbm
3255: qbm qp m qbm qp
2175: qbm qp m qbp qm
5235: qbm qp qm m qbp
915: qbm qp qm qbp m
7395: qbm qp qm qbp p
6315: qbm qp qm p qbp
3975: qbm qp qp m qbm
735: qbm qp qp qbm m
7215: qbm qp qp qbm p
5055: qbm qp qp p qbm
2715: qbm qp qbm m qp
555: qbm qp qbm qp m
7035: qbm qp qbm qp p
3795: qbm qp qbm p qp
1455: qbm qp qbp m qm
375: qbm qp qbp qm m
6855: qbm qp qbp qm p
2535: qbm qp qbp p qm
5595: qbm qp p qm qbp
4515: qbm qp p qp qbm
3435: qbm qp p qbm qp
2355: qbm qp p qbp qm
3045: qbm qbm m qp qp
2685: qbm qbm qp m qp
525: qbm qbm qp qp m
7005: qbm qbm qp qp p
3765: qbm qbm qp p qp
3225: qbm qbm p qp qp
2835: qbm qbp m qm qp
1755: qbm qbp m qp qm
2655: qbm qbp qm m qp
495: qbm qbp qm qp m
6975: qbm qbp qm qp p
3735: qbm qbp qm p qp
1395: qbm qbp qp m qm
315: qbm qbp qp qm m
6795: qbm qbp qp qm p
2475: qbm qbp qp p qm
3015: qbm qbp p qm qp
1935: qbm qbp p qp qm
5685: qbm p qm qp qbp
3525: qbm p qm qbp qp
5505: qbm p qp qm qbp
4425: qbm p qp qp qbm
3345: qbm p qp qbm qp
2265: qbm p qp qbp qm
3165: qbm p qbm qp qp
2985: qbm p qbp qm qp
1905: qbm p qbp qp qm
5440: qbp m qm qm qbp
4360: qbp m qm qp qbm
3280: qbp m qm qbm qp
2200: qbp m qm qbp qm
4180: qbp m qp qm qbm
2020: qbp m qp qbm qm
2920: qbp m qbm qm qp
1840: qbp m qbm qp qm
1660: qbp m qbp qm qm
5410: qbp qm m qm qbp
4330: qbp qm m qp qbm
3250: qbp qm m qbm qp
2170: qbp qm m qbp qm
5230: qbp qm qm m qbp
910: qbp qm qm qbp m
7390: qbp qm qm qbp p
6310: qbp qm qm p qbp
3970: qbp qm qp m qbm
730: qbp qm qp qbm m
7210: qbp qm qp qbm p
5050: qbp qm qp p qbm
2710: qbp qm qbm m qp
550: qbp qm qbm qp m
7030: qbp qm qbm qp p
3790: qbp qm qbm p qp
1450: qbp qm qbp m qm
370: qbp qm qbp qm m
6850: qbp qm qbp qm p
2530: qbp qm qbp p qm
5590: qbp qm p qm qbp
4510: qbp qm p qp qbm
3430: qbp qm p qbm qp
2350: qbp qm p qbp qm
4120: qbp qp m qm qbm
1960: qbp qp m qbm qm
3940: qbp qp qm m qbm
700: qbp qp qm qbm m
7180: qbp qp qm qbm p
5020: qbp qp qm p qbm
1420: qbp qp qbm m qm
340: qbp qp qbm qm m
6820: qbp qp qbm qm p
2500: qbp qp qbm p qm
4300: qbp qp p qm qbm
2140: qbp qp p qbm qm
2830: qbp qbm m qm qp
1750: qbp qbm m qp qm
2650: qbp qbm qm m qp
490: qbp qbm qm qp m
6970: qbp qbm qm qp p
3730: qbp qbm qm p qp
1390: qbp qbm qp m qm
310: qbp qbm qp qm m
6790: qbp qbm qp qm p
2470: qbp qbm qp p qm
3010: qbp qbm p qm qp
1930: qbp qbm p qp qm
1540: qbp qbp m qm qm
1360: qbp qbp qm m qm
280: qbp qbp qm qm m
6760: qbp qbp qm qm p
2440: qbp qbp qm p qm
1720: qbp qbp p qm qm
5470: qbp p qm qm qbp
4390: qbp p qm qp qbm
3310: qbp p qm qbm qp
2230: qbp p qm qbp qm
4210: qbp p qp qm qbm
2050: qbp p qp qbm qm
2950: qbp p qbm qm qp
1870: qbp p qbm qp qm
1690: qbp p qbp qm qm
6095: p qm qm qbp qbp
5915: p qm qp qbm qbp
4835: p qm qp qbp qbm
5735: p qm qbm qp qbp
3575: p qm qbm qbp qp
5555: p qm qbp qm qbp
4475: p qm qbp qp qbm
3395: p qm qbp qbm qp
2315: p qm qbp qbp qm
5885: p qp qm qbm qbp
4805: p qp qm qbp qbm
4625: p qp qp qbm qbm
5525: p qp qbm qm qbp
4445: p qp qbm qp qbm
3365: p qp qbm qbm qp
2285: p qp qbm qbp qm
4265: p qp qbp qm qbm
2105: p qp qbp qbm qm
5675: p qbm qm qp qbp
3515: p qbm qm qbp qp
5495: p qbm qp qm qbp
4415: p qbm qp qp qbm
3335: p qbm qp qbm qp
2255: p qbm qp qbp qm
3155: p qbm qbm qp qp
2975: p qbm qbp qm qp
1895: p qbm qbp qp qm
5465: p qbp qm qm qbp
4385: p qbp qm qp qbm
3305: p qbp qm qbm qp
2225: p qbp qm qbp qm
4205: p qbp qp qm qbm
2045: p qbp qp qbm qm
2945: p qbp qbm qm qp
1865: p qbp qbm qp qm
1685: p qbp qbp qm qm
*/

