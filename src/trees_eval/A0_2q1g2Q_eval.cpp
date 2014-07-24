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

namespace BH {


template <class T> complex<T> A2q1g2Qzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q1g2Q317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q1052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q1162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1397_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1457_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q1472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q1502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1757_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1865_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q1870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q1877_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1895_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q1902_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1905_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q1912_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1930_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q1932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q1935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2105_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2140_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2192_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2255_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2262_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2360_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2452_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2470_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2472_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2500_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q2535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q2657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2662_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2747_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2767_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2782_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q2797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2837_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2982_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q2985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q2992_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q3010_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3012_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q3015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3305_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3415_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3430_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q3467_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3487_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3515_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3522_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3575_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q3637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3642_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3645_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q3655_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q3682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3697_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q3805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q3817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q3825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q3835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q4205_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4265_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4390_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4475_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q4940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q4945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q5252_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q5287_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5372_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q5377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q5417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5505_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5582_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5592_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5600_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5647_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5682_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q5685_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q5755_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5802_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5815_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q5915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q5935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6020_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6025_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6272_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6277_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6302_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6320_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6337_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6342_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6345_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6355_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q6532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2)*ep.spb(4,2))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6542_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6562_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6577_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6632_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q6712_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6722_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6772_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6790_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q6792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q6820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q6830_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6852_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6855_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6860_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6922_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6937_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q6952_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6970_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q6972_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2Q6975_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q7030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2Q7045_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7062_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7065_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7075_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7190_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(2,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q1g2Q7280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7285_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2Q7425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2Q7465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
 ); }


template <class T> complex<T>  (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 317 : return &A2q1g2Q317_eval;
case 322 : return &A2q1g2Q322_eval;
case 377 : return &A2q1g2Q377_eval;
case 392 : return &A2q1g2Q392_eval;
case 412 : return &A2q1g2Q412_eval;
case 422 : return &A2q1g2Q422_eval;
case 497 : return &A2q1g2Q497_eval;
case 502 : return &A2q1g2Q502_eval;
case 587 : return &A2q1g2Q587_eval;
case 607 : return &A2q1g2Q607_eval;
case 622 : return &A2q1g2Q622_eval;
case 637 : return &A2q1g2Q637_eval;
case 917 : return &A2q1g2Q917_eval;
case 932 : return &A2q1g2Q932_eval;
case 947 : return &A2q1g2Q947_eval;
case 967 : return &A2q1g2Q967_eval;
case 1052 : return &A2q1g2Q1052_eval;
case 1057 : return &A2q1g2Q1057_eval;
case 1132 : return &A2q1g2Q1132_eval;
case 1142 : return &A2q1g2Q1142_eval;
case 1162 : return &A2q1g2Q1162_eval;
case 1177 : return &A2q1g2Q1177_eval;
case 1232 : return &A2q1g2Q1232_eval;
case 1237 : return &A2q1g2Q1237_eval;
case 1397 : return &A2q1g2Q1397_eval;
case 1402 : return &A2q1g2Q1402_eval;
case 1457 : return &A2q1g2Q1457_eval;
case 1472 : return &A2q1g2Q1472_eval;
case 1492 : return &A2q1g2Q1492_eval;
case 1502 : return &A2q1g2Q1502_eval;
case 1757 : return &A2q1g2Q1757_eval;
case 1762 : return &A2q1g2Q1762_eval;
case 1865 : return &A2q1g2Q1865_eval;
case 1870 : return &A2q1g2Q1870_eval;
case 1877 : return &A2q1g2Q1877_eval;
case 1895 : return &A2q1g2Q1895_eval;
case 1902 : return &A2q1g2Q1902_eval;
case 1905 : return &A2q1g2Q1905_eval;
case 1912 : return &A2q1g2Q1912_eval;
case 1930 : return &A2q1g2Q1930_eval;
case 1932 : return &A2q1g2Q1932_eval;
case 1935 : return &A2q1g2Q1935_eval;
case 2045 : return &A2q1g2Q2045_eval;
case 2050 : return &A2q1g2Q2050_eval;
case 2105 : return &A2q1g2Q2105_eval;
case 2120 : return &A2q1g2Q2120_eval;
case 2140 : return &A2q1g2Q2140_eval;
case 2150 : return &A2q1g2Q2150_eval;
case 2177 : return &A2q1g2Q2177_eval;
case 2192 : return &A2q1g2Q2192_eval;
case 2237 : return &A2q1g2Q2237_eval;
case 2255 : return &A2q1g2Q2255_eval;
case 2262 : return &A2q1g2Q2262_eval;
case 2265 : return &A2q1g2Q2265_eval;
case 2285 : return &A2q1g2Q2285_eval;
case 2300 : return &A2q1g2Q2300_eval;
case 2342 : return &A2q1g2Q2342_eval;
case 2352 : return &A2q1g2Q2352_eval;
case 2355 : return &A2q1g2Q2355_eval;
case 2360 : return &A2q1g2Q2360_eval;
case 2392 : return &A2q1g2Q2392_eval;
case 2402 : return &A2q1g2Q2402_eval;
case 2452 : return &A2q1g2Q2452_eval;
case 2470 : return &A2q1g2Q2470_eval;
case 2472 : return &A2q1g2Q2472_eval;
case 2475 : return &A2q1g2Q2475_eval;
case 2500 : return &A2q1g2Q2500_eval;
case 2510 : return &A2q1g2Q2510_eval;
case 2522 : return &A2q1g2Q2522_eval;
case 2532 : return &A2q1g2Q2532_eval;
case 2535 : return &A2q1g2Q2535_eval;
case 2540 : return &A2q1g2Q2540_eval;
case 2657 : return &A2q1g2Q2657_eval;
case 2662 : return &A2q1g2Q2662_eval;
case 2747 : return &A2q1g2Q2747_eval;
case 2767 : return &A2q1g2Q2767_eval;
case 2782 : return &A2q1g2Q2782_eval;
case 2797 : return &A2q1g2Q2797_eval;
case 2837 : return &A2q1g2Q2837_eval;
case 2842 : return &A2q1g2Q2842_eval;
case 2945 : return &A2q1g2Q2945_eval;
case 2950 : return &A2q1g2Q2950_eval;
case 2957 : return &A2q1g2Q2957_eval;
case 2975 : return &A2q1g2Q2975_eval;
case 2982 : return &A2q1g2Q2982_eval;
case 2985 : return &A2q1g2Q2985_eval;
case 2992 : return &A2q1g2Q2992_eval;
case 3010 : return &A2q1g2Q3010_eval;
case 3012 : return &A2q1g2Q3012_eval;
case 3015 : return &A2q1g2Q3015_eval;
case 3305 : return &A2q1g2Q3305_eval;
case 3310 : return &A2q1g2Q3310_eval;
case 3395 : return &A2q1g2Q3395_eval;
case 3415 : return &A2q1g2Q3415_eval;
case 3430 : return &A2q1g2Q3430_eval;
case 3445 : return &A2q1g2Q3445_eval;
case 3467 : return &A2q1g2Q3467_eval;
case 3487 : return &A2q1g2Q3487_eval;
case 3497 : return &A2q1g2Q3497_eval;
case 3515 : return &A2q1g2Q3515_eval;
case 3522 : return &A2q1g2Q3522_eval;
case 3525 : return &A2q1g2Q3525_eval;
case 3575 : return &A2q1g2Q3575_eval;
case 3595 : return &A2q1g2Q3595_eval;
case 3637 : return &A2q1g2Q3637_eval;
case 3642 : return &A2q1g2Q3642_eval;
case 3645 : return &A2q1g2Q3645_eval;
case 3655 : return &A2q1g2Q3655_eval;
case 3682 : return &A2q1g2Q3682_eval;
case 3697 : return &A2q1g2Q3697_eval;
case 3712 : return &A2q1g2Q3712_eval;
case 3730 : return &A2q1g2Q3730_eval;
case 3732 : return &A2q1g2Q3732_eval;
case 3735 : return &A2q1g2Q3735_eval;
case 3790 : return &A2q1g2Q3790_eval;
case 3805 : return &A2q1g2Q3805_eval;
case 3817 : return &A2q1g2Q3817_eval;
case 3822 : return &A2q1g2Q3822_eval;
case 3825 : return &A2q1g2Q3825_eval;
case 3835 : return &A2q1g2Q3835_eval;
case 4205 : return &A2q1g2Q4205_eval;
case 4210 : return &A2q1g2Q4210_eval;
case 4265 : return &A2q1g2Q4265_eval;
case 4280 : return &A2q1g2Q4280_eval;
case 4300 : return &A2q1g2Q4300_eval;
case 4310 : return &A2q1g2Q4310_eval;
case 4385 : return &A2q1g2Q4385_eval;
case 4390 : return &A2q1g2Q4390_eval;
case 4475 : return &A2q1g2Q4475_eval;
case 4495 : return &A2q1g2Q4495_eval;
case 4510 : return &A2q1g2Q4510_eval;
case 4525 : return &A2q1g2Q4525_eval;
case 4805 : return &A2q1g2Q4805_eval;
case 4820 : return &A2q1g2Q4820_eval;
case 4835 : return &A2q1g2Q4835_eval;
case 4855 : return &A2q1g2Q4855_eval;
case 4940 : return &A2q1g2Q4940_eval;
case 4945 : return &A2q1g2Q4945_eval;
case 5020 : return &A2q1g2Q5020_eval;
case 5030 : return &A2q1g2Q5030_eval;
case 5050 : return &A2q1g2Q5050_eval;
case 5065 : return &A2q1g2Q5065_eval;
case 5120 : return &A2q1g2Q5120_eval;
case 5125 : return &A2q1g2Q5125_eval;
case 5237 : return &A2q1g2Q5237_eval;
case 5252 : return &A2q1g2Q5252_eval;
case 5267 : return &A2q1g2Q5267_eval;
case 5287 : return &A2q1g2Q5287_eval;
case 5372 : return &A2q1g2Q5372_eval;
case 5377 : return &A2q1g2Q5377_eval;
case 5417 : return &A2q1g2Q5417_eval;
case 5432 : return &A2q1g2Q5432_eval;
case 5477 : return &A2q1g2Q5477_eval;
case 5495 : return &A2q1g2Q5495_eval;
case 5502 : return &A2q1g2Q5502_eval;
case 5505 : return &A2q1g2Q5505_eval;
case 5525 : return &A2q1g2Q5525_eval;
case 5540 : return &A2q1g2Q5540_eval;
case 5582 : return &A2q1g2Q5582_eval;
case 5592 : return &A2q1g2Q5592_eval;
case 5595 : return &A2q1g2Q5595_eval;
case 5600 : return &A2q1g2Q5600_eval;
case 5627 : return &A2q1g2Q5627_eval;
case 5647 : return &A2q1g2Q5647_eval;
case 5657 : return &A2q1g2Q5657_eval;
case 5675 : return &A2q1g2Q5675_eval;
case 5682 : return &A2q1g2Q5682_eval;
case 5685 : return &A2q1g2Q5685_eval;
case 5735 : return &A2q1g2Q5735_eval;
case 5755 : return &A2q1g2Q5755_eval;
case 5797 : return &A2q1g2Q5797_eval;
case 5802 : return &A2q1g2Q5802_eval;
case 5805 : return &A2q1g2Q5805_eval;
case 5815 : return &A2q1g2Q5815_eval;
case 5885 : return &A2q1g2Q5885_eval;
case 5900 : return &A2q1g2Q5900_eval;
case 5915 : return &A2q1g2Q5915_eval;
case 5935 : return &A2q1g2Q5935_eval;
case 6020 : return &A2q1g2Q6020_eval;
case 6025 : return &A2q1g2Q6025_eval;
case 6272 : return &A2q1g2Q6272_eval;
case 6277 : return &A2q1g2Q6277_eval;
case 6302 : return &A2q1g2Q6302_eval;
case 6312 : return &A2q1g2Q6312_eval;
case 6315 : return &A2q1g2Q6315_eval;
case 6320 : return &A2q1g2Q6320_eval;
case 6337 : return &A2q1g2Q6337_eval;
case 6342 : return &A2q1g2Q6342_eval;
case 6345 : return &A2q1g2Q6345_eval;
case 6355 : return &A2q1g2Q6355_eval;
case 6380 : return &A2q1g2Q6380_eval;
case 6385 : return &A2q1g2Q6385_eval;
case 6532 : return &A2q1g2Q6532_eval;
case 6542 : return &A2q1g2Q6542_eval;
case 6562 : return &A2q1g2Q6562_eval;
case 6577 : return &A2q1g2Q6577_eval;
case 6632 : return &A2q1g2Q6632_eval;
case 6637 : return &A2q1g2Q6637_eval;
case 6712 : return &A2q1g2Q6712_eval;
case 6722 : return &A2q1g2Q6722_eval;
case 6772 : return &A2q1g2Q6772_eval;
case 6790 : return &A2q1g2Q6790_eval;
case 6792 : return &A2q1g2Q6792_eval;
case 6795 : return &A2q1g2Q6795_eval;
case 6820 : return &A2q1g2Q6820_eval;
case 6830 : return &A2q1g2Q6830_eval;
case 6842 : return &A2q1g2Q6842_eval;
case 6852 : return &A2q1g2Q6852_eval;
case 6855 : return &A2q1g2Q6855_eval;
case 6860 : return &A2q1g2Q6860_eval;
case 6922 : return &A2q1g2Q6922_eval;
case 6937 : return &A2q1g2Q6937_eval;
case 6952 : return &A2q1g2Q6952_eval;
case 6970 : return &A2q1g2Q6970_eval;
case 6972 : return &A2q1g2Q6972_eval;
case 6975 : return &A2q1g2Q6975_eval;
case 7030 : return &A2q1g2Q7030_eval;
case 7045 : return &A2q1g2Q7045_eval;
case 7057 : return &A2q1g2Q7057_eval;
case 7062 : return &A2q1g2Q7062_eval;
case 7065 : return &A2q1g2Q7065_eval;
case 7075 : return &A2q1g2Q7075_eval;
case 7180 : return &A2q1g2Q7180_eval;
case 7190 : return &A2q1g2Q7190_eval;
case 7210 : return &A2q1g2Q7210_eval;
case 7225 : return &A2q1g2Q7225_eval;
case 7280 : return &A2q1g2Q7280_eval;
case 7285 : return &A2q1g2Q7285_eval;
case 7352 : return &A2q1g2Q7352_eval;
case 7357 : return &A2q1g2Q7357_eval;
case 7382 : return &A2q1g2Q7382_eval;
case 7392 : return &A2q1g2Q7392_eval;
case 7395 : return &A2q1g2Q7395_eval;
case 7400 : return &A2q1g2Q7400_eval;
case 7417 : return &A2q1g2Q7417_eval;
case 7422 : return &A2q1g2Q7422_eval;
case 7425 : return &A2q1g2Q7425_eval;
case 7435 : return &A2q1g2Q7435_eval;
case 7460 : return &A2q1g2Q7460_eval;
case 7465 : return &A2q1g2Q7465_eval;

default: return 0;

}
}




template complex<R>  (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q1g2Q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
7422: m qm qp Qm Qp
6342: m qm qp Qp Qm
7062: m qm Qm qp Qp
3822: m qm Qm Qp qp
5802: m qm Qp qp Qm
3642: m qm Qp Qm qp
7392: m qp qm Qm Qp
6312: m qp qm Qp Qm
6852: m qp Qm qm Qp
2532: m qp Qm Qp qm
5592: m qp Qp qm Qm
2352: m qp Qp Qm qm
6972: m Qm qm qp Qp
3732: m Qm qm Qp qp
6792: m Qm qp qm Qp
2472: m Qm qp Qp qm
3012: m Qm Qp qm qp
1932: m Qm Qp qp qm
5682: m Qp qm qp Qm
3522: m Qp qm Qm qp
5502: m Qp qp qm Qm
2262: m Qp qp Qm qm
2982: m Qp Qm qm qp
1902: m Qp Qm qp qm
7417: qm m qp Qm Qp
6337: qm m qp Qp Qm
7057: qm m Qm qp Qp
3817: qm m Qm Qp qp
5797: qm m Qp qp Qm
3637: qm m Qp Qm qp
7357: qm qp m Qm Qp
6277: qm qp m Qp Qm
7465: qm qp p Qm Qp
6385: qm qp p Qp Qm
6637: qm qp Qm m Qp
7285: qm qp Qm p Qp
1237: qm qp Qm Qp m
5125: qm qp Qm Qp p
5377: qm qp Qp m Qm
6025: qm qp Qp p Qm
1057: qm qp Qp Qm m
4945: qm qp Qp Qm p
7435: qm p qp Qm Qp
6355: qm p qp Qp Qm
7075: qm p Qm qp Qp
3835: qm p Qm Qp qp
5815: qm p Qp qp Qm
3655: qm p Qp Qm qp
6937: qm Qm m qp Qp
3697: qm Qm m Qp qp
6577: qm Qm qp m Qp
7225: qm Qm qp p Qp
1177: qm Qm qp Qp m
5065: qm Qm qp Qp p
7045: qm Qm p qp Qp
3805: qm Qm p Qp qp
2797: qm Qm Qp m qp
637: qm Qm Qp qp m
4525: qm Qm Qp qp p
3445: qm Qm Qp p qp
5647: qm Qp m qp Qm
3487: qm Qp m Qm qp
5287: qm Qp qp m Qm
5935: qm Qp qp p Qm
967: qm Qp qp Qm m
4855: qm Qp qp Qm p
5755: qm Qp p qp Qm
3595: qm Qp p Qm qp
2767: qm Qp Qm m qp
607: qm Qp Qm qp m
4495: qm Qp Qm qp p
3415: qm Qp Qm p qp
7382: qp m qm Qm Qp
6302: qp m qm Qp Qm
6842: qp m Qm qm Qp
2522: qp m Qm Qp qm
5582: qp m Qp qm Qm
2342: qp m Qp Qm qm
7352: qp qm m Qm Qp
6272: qp qm m Qp Qm
7460: qp qm p Qm Qp
6380: qp qm p Qp Qm
6632: qp qm Qm m Qp
7280: qp qm Qm p Qp
1232: qp qm Qm Qp m
5120: qp qm Qm Qp p
5372: qp qm Qp m Qm
6020: qp qm Qp p Qm
1052: qp qm Qp Qm m
4940: qp qm Qp Qm p
7400: qp p qm Qm Qp
6320: qp p qm Qp Qm
6860: qp p Qm qm Qp
2540: qp p Qm Qp qm
5600: qp p Qp qm Qm
2360: qp p Qp Qm qm
6722: qp Qm m qm Qp
2402: qp Qm m Qp qm
6542: qp Qm qm m Qp
7190: qp Qm qm p Qp
1142: qp Qm qm Qp m
5030: qp Qm qm Qp p
6830: qp Qm p qm Qp
2510: qp Qm p Qp qm
1502: qp Qm Qp m qm
422: qp Qm Qp qm m
4310: qp Qm Qp qm p
2150: qp Qm Qp p qm
5432: qp Qp m qm Qm
2192: qp Qp m Qm qm
5252: qp Qp qm m Qm
5900: qp Qp qm p Qm
932: qp Qp qm Qm m
4820: qp Qp qm Qm p
5540: qp Qp p qm Qm
2300: qp Qp p Qm qm
1472: qp Qp Qm m qm
392: qp Qp Qm qm m
4280: qp Qp Qm qm p
2120: qp Qp Qm p qm
7425: p qm qp Qm Qp
6345: p qm qp Qp Qm
7065: p qm Qm qp Qp
3825: p qm Qm Qp qp
5805: p qm Qp qp Qm
3645: p qm Qp Qm qp
7395: p qp qm Qm Qp
6315: p qp qm Qp Qm
6855: p qp Qm qm Qp
2535: p qp Qm Qp qm
5595: p qp Qp qm Qm
2355: p qp Qp Qm qm
6975: p Qm qm qp Qp
3735: p Qm qm Qp qp
6795: p Qm qp qm Qp
2475: p Qm qp Qp qm
3015: p Qm Qp qm qp
1935: p Qm Qp qp qm
5685: p Qp qm qp Qm
3525: p Qp qm Qm qp
5505: p Qp qp qm Qm
2265: p Qp qp Qm qm
2985: p Qp Qm qm qp
1905: p Qp Qm qp qm
6952: Qm m qm qp Qp
3712: Qm m qm Qp qp
6772: Qm m qp qm Qp
2452: Qm m qp Qp qm
2992: Qm m Qp qm qp
1912: Qm m Qp qp qm
6922: Qm qm m qp Qp
3682: Qm qm m Qp qp
6562: Qm qm qp m Qp
7210: Qm qm qp p Qp
1162: Qm qm qp Qp m
5050: Qm qm qp Qp p
7030: Qm qm p qp Qp
3790: Qm qm p Qp qp
2782: Qm qm Qp m qp
622: Qm qm Qp qp m
4510: Qm qm Qp qp p
3430: Qm qm Qp p qp
6712: Qm qp m qm Qp
2392: Qm qp m Qp qm
6532: Qm qp qm m Qp
7180: Qm qp qm p Qp
1132: Qm qp qm Qp m
5020: Qm qp qm Qp p
6820: Qm qp p qm Qp
2500: Qm qp p Qp qm
1492: Qm qp Qp m qm
412: Qm qp Qp qm m
4300: Qm qp Qp qm p
2140: Qm qp Qp p qm
6970: Qm p qm qp Qp
3730: Qm p qm Qp qp
6790: Qm p qp qm Qp
2470: Qm p qp Qp qm
3010: Qm p Qp qm qp
1930: Qm p Qp qp qm
2842: Qm Qp m qm qp
1762: Qm Qp m qp qm
2662: Qm Qp qm m qp
502: Qm Qp qm qp m
4390: Qm Qp qm qp p
3310: Qm Qp qm p qp
1402: Qm Qp qp m qm
322: Qm Qp qp qm m
4210: Qm Qp qp qm p
2050: Qm Qp qp p qm
2950: Qm Qp p qm qp
1870: Qm Qp p qp qm
5657: Qp m qm qp Qm
3497: Qp m qm Qm qp
5477: Qp m qp qm Qm
2237: Qp m qp Qm qm
2957: Qp m Qm qm qp
1877: Qp m Qm qp qm
5627: Qp qm m qp Qm
3467: Qp qm m Qm qp
5267: Qp qm qp m Qm
5915: Qp qm qp p Qm
947: Qp qm qp Qm m
4835: Qp qm qp Qm p
5735: Qp qm p qp Qm
3575: Qp qm p Qm qp
2747: Qp qm Qm m qp
587: Qp qm Qm qp m
4475: Qp qm Qm qp p
3395: Qp qm Qm p qp
5417: Qp qp m qm Qm
2177: Qp qp m Qm qm
5237: Qp qp qm m Qm
5885: Qp qp qm p Qm
917: Qp qp qm Qm m
4805: Qp qp qm Qm p
5525: Qp qp p qm Qm
2285: Qp qp p Qm qm
1457: Qp qp Qm m qm
377: Qp qp Qm qm m
4265: Qp qp Qm qm p
2105: Qp qp Qm p qm
5675: Qp p qm qp Qm
3515: Qp p qm Qm qp
5495: Qp p qp qm Qm
2255: Qp p qp Qm qm
2975: Qp p Qm qm qp
1895: Qp p Qm qp qm
2837: Qp Qm m qm qp
1757: Qp Qm m qp qm
2657: Qp Qm qm m qp
497: Qp Qm qm qp m
4385: Qp Qm qm qp p
3305: Qp Qm qm p qp
1397: Qp Qm qp m qm
317: Qp Qm qp qm m
4205: Qp Qm qp qm p
2045: Qp Qm qp p qm
2945: Qp Qm p qm qp
1865: Qp Qm p qp qm
*/

