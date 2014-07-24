/* The 2 quark amplitudes use a base 4 code to enumerate them
   m,qm,qp,p corresponds to 0,1,2,3.
   For example, (m,qm,p,qp) = 0 + 1*4 + 3*16 + 2*64
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


template <class T> complex<T> A2q3gzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q3g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   complex<T>(0,0)
); }

template <class T> complex<T> A2q3g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g18_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g24_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g27_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(2,0))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g30_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2))/(ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g33_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g36_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g39_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g45_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g54_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g57_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g66_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g72_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g75_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2)*ep.spb(3,0))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g78_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2)*ep.spb(3,1))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g96_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g99_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g108_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g111_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q3g114_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q3g126_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g129_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g141_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g144_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g156_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g159_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g189_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g201_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g216_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g219_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g231_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g246_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g249_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g258_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g264_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g267_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g270_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(1,0),2)*ep.spb(4,1))/
 (ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g288_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g291_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g300_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,1),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g303_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2))
); }

template <class T> complex<T> A2q3g306_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,0),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g312_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q3g318_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q3g384_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g387_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g396_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),3)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g399_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g432_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,0)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g444_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g447_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g450_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g456_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g459_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g462_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g480_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g483_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g492_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g498_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g504_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/(ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g507_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g510_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g513_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g516_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g519_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),3))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g528_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g531_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g540_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g543_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2))
); }

template <class T> complex<T> A2q3g561_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3)*ep.spb(2,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g564_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g567_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g573_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q3g576_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g579_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g588_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g591_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g624_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g636_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g639_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g705_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(3,0))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q3g708_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2)*ep.spb(3,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q3g711_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g717_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3)*ep.spa(2,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g720_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g723_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(1,4))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g732_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g753_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(1,4))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g756_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2))/(ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g759_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g765_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g774_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g777_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g786_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,3))
); }

template <class T> complex<T> A2q3g792_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g795_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g798_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g801_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g804_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q3g807_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g813_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g834_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,0),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q3g840_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,1),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g843_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g846_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g864_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),3))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q3g867_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g876_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g879_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g882_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g888_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g891_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g894_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g897_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(2,1)*ep.spb(3,2))
); }

template <class T> complex<T> A2q3g900_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(4,1))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q3g903_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g909_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g912_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2)*ep.spb(4,2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*
ep.spb(4,0))
); }

template <class T> complex<T> A2q3g915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g924_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3)*ep.spa(0,3))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(1,2)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g927_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(1,3))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g948_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*ep.spa(0,3))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g951_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g957_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g966_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,4)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g969_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g978_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,4)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g984_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g987_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g990_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g993_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q3g996_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*ep.spa(0,2))/
 (ep.spa(0,4)*ep.spa(1,2)*ep.spa(2,3)*
ep.spa(3,4))
); }

template <class T> complex<T> A2q3g999_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g1005_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g1014_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q3g1017_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
 ); }


template <class T> complex<T>  (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 6 : return &A2q3g6_eval;
case 9 : return &A2q3g9_eval;
case 18 : return &A2q3g18_eval;
case 24 : return &A2q3g24_eval;
case 27 : return &A2q3g27_eval;
case 30 : return &A2q3g30_eval;
case 33 : return &A2q3g33_eval;
case 36 : return &A2q3g36_eval;
case 39 : return &A2q3g39_eval;
case 45 : return &A2q3g45_eval;
case 54 : return &A2q3g54_eval;
case 57 : return &A2q3g57_eval;
case 66 : return &A2q3g66_eval;
case 72 : return &A2q3g72_eval;
case 75 : return &A2q3g75_eval;
case 78 : return &A2q3g78_eval;
case 96 : return &A2q3g96_eval;
case 99 : return &A2q3g99_eval;
case 108 : return &A2q3g108_eval;
case 111 : return &A2q3g111_eval;
case 114 : return &A2q3g114_eval;
case 120 : return &A2q3g120_eval;
case 123 : return &A2q3g123_eval;
case 126 : return &A2q3g126_eval;
case 129 : return &A2q3g129_eval;
case 132 : return &A2q3g132_eval;
case 135 : return &A2q3g135_eval;
case 141 : return &A2q3g141_eval;
case 144 : return &A2q3g144_eval;
case 147 : return &A2q3g147_eval;
case 156 : return &A2q3g156_eval;
case 159 : return &A2q3g159_eval;
case 177 : return &A2q3g177_eval;
case 180 : return &A2q3g180_eval;
case 183 : return &A2q3g183_eval;
case 189 : return &A2q3g189_eval;
case 198 : return &A2q3g198_eval;
case 201 : return &A2q3g201_eval;
case 210 : return &A2q3g210_eval;
case 216 : return &A2q3g216_eval;
case 219 : return &A2q3g219_eval;
case 222 : return &A2q3g222_eval;
case 225 : return &A2q3g225_eval;
case 228 : return &A2q3g228_eval;
case 231 : return &A2q3g231_eval;
case 237 : return &A2q3g237_eval;
case 246 : return &A2q3g246_eval;
case 249 : return &A2q3g249_eval;
case 258 : return &A2q3g258_eval;
case 264 : return &A2q3g264_eval;
case 267 : return &A2q3g267_eval;
case 270 : return &A2q3g270_eval;
case 288 : return &A2q3g288_eval;
case 291 : return &A2q3g291_eval;
case 300 : return &A2q3g300_eval;
case 303 : return &A2q3g303_eval;
case 306 : return &A2q3g306_eval;
case 312 : return &A2q3g312_eval;
case 315 : return &A2q3g315_eval;
case 318 : return &A2q3g318_eval;
case 384 : return &A2q3g384_eval;
case 387 : return &A2q3g387_eval;
case 396 : return &A2q3g396_eval;
case 399 : return &A2q3g399_eval;
case 432 : return &A2q3g432_eval;
case 435 : return &A2q3g435_eval;
case 444 : return &A2q3g444_eval;
case 447 : return &A2q3g447_eval;
case 450 : return &A2q3g450_eval;
case 456 : return &A2q3g456_eval;
case 459 : return &A2q3g459_eval;
case 462 : return &A2q3g462_eval;
case 480 : return &A2q3g480_eval;
case 483 : return &A2q3g483_eval;
case 492 : return &A2q3g492_eval;
case 495 : return &A2q3g495_eval;
case 498 : return &A2q3g498_eval;
case 504 : return &A2q3g504_eval;
case 507 : return &A2q3g507_eval;
case 510 : return &A2q3g510_eval;
case 513 : return &A2q3g513_eval;
case 516 : return &A2q3g516_eval;
case 519 : return &A2q3g519_eval;
case 525 : return &A2q3g525_eval;
case 528 : return &A2q3g528_eval;
case 531 : return &A2q3g531_eval;
case 540 : return &A2q3g540_eval;
case 543 : return &A2q3g543_eval;
case 561 : return &A2q3g561_eval;
case 564 : return &A2q3g564_eval;
case 567 : return &A2q3g567_eval;
case 573 : return &A2q3g573_eval;
case 576 : return &A2q3g576_eval;
case 579 : return &A2q3g579_eval;
case 588 : return &A2q3g588_eval;
case 591 : return &A2q3g591_eval;
case 624 : return &A2q3g624_eval;
case 627 : return &A2q3g627_eval;
case 636 : return &A2q3g636_eval;
case 639 : return &A2q3g639_eval;
case 705 : return &A2q3g705_eval;
case 708 : return &A2q3g708_eval;
case 711 : return &A2q3g711_eval;
case 717 : return &A2q3g717_eval;
case 720 : return &A2q3g720_eval;
case 723 : return &A2q3g723_eval;
case 732 : return &A2q3g732_eval;
case 735 : return &A2q3g735_eval;
case 753 : return &A2q3g753_eval;
case 756 : return &A2q3g756_eval;
case 759 : return &A2q3g759_eval;
case 765 : return &A2q3g765_eval;
case 774 : return &A2q3g774_eval;
case 777 : return &A2q3g777_eval;
case 786 : return &A2q3g786_eval;
case 792 : return &A2q3g792_eval;
case 795 : return &A2q3g795_eval;
case 798 : return &A2q3g798_eval;
case 801 : return &A2q3g801_eval;
case 804 : return &A2q3g804_eval;
case 807 : return &A2q3g807_eval;
case 813 : return &A2q3g813_eval;
case 822 : return &A2q3g822_eval;
case 825 : return &A2q3g825_eval;
case 834 : return &A2q3g834_eval;
case 840 : return &A2q3g840_eval;
case 843 : return &A2q3g843_eval;
case 846 : return &A2q3g846_eval;
case 864 : return &A2q3g864_eval;
case 867 : return &A2q3g867_eval;
case 876 : return &A2q3g876_eval;
case 879 : return &A2q3g879_eval;
case 882 : return &A2q3g882_eval;
case 888 : return &A2q3g888_eval;
case 891 : return &A2q3g891_eval;
case 894 : return &A2q3g894_eval;
case 897 : return &A2q3g897_eval;
case 900 : return &A2q3g900_eval;
case 903 : return &A2q3g903_eval;
case 909 : return &A2q3g909_eval;
case 912 : return &A2q3g912_eval;
case 915 : return &A2q3g915_eval;
case 924 : return &A2q3g924_eval;
case 927 : return &A2q3g927_eval;
case 945 : return &A2q3g945_eval;
case 948 : return &A2q3g948_eval;
case 951 : return &A2q3g951_eval;
case 957 : return &A2q3g957_eval;
case 966 : return &A2q3g966_eval;
case 969 : return &A2q3g969_eval;
case 978 : return &A2q3g978_eval;
case 984 : return &A2q3g984_eval;
case 987 : return &A2q3g987_eval;
case 990 : return &A2q3g990_eval;
case 993 : return &A2q3g993_eval;
case 996 : return &A2q3g996_eval;
case 999 : return &A2q3g999_eval;
case 1005 : return &A2q3g1005_eval;
case 1014 : return &A2q3g1014_eval;
case 1017 : return &A2q3g1017_eval;

default: _WARNING3("Unknown pointer amplitude (*A2q3g_Tree_Ptr(int hc)) - case:",hc , " - throw BH error.");
		throw BHerror("case missing for tree amplitude!");

}
}



template complex<R>  (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q3g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
576: m m m qm qp384: m m m qp qm528: m m qm m qp144: m m qm qp m912: m m qm qp p720: m m qm p qp288: m m qp m qm96: m m qp qm m864: m m qp qm p480: m m qp p qm624: m m p qm qp432: m m p qp qm516: m qm m m qp132: m qm m qp m900: m qm m qp p708: m qm m p qp36: m qm qp m m804: m qm qp m p228: m qm qp p m996: m qm qp p p564: m qm p m qp180: m qm p qp m948: m qm p qp p756: m qm p p qp264: m qp m m qm72: m qp m qm m840: m qp m qm p456: m qp m p qm24: m qp qm m m792: m qp qm m p216: m qp qm p m984: m qp qm p p312: m qp p m qm120: m qp p qm m888: m qp p qm p504: m qp p p qm588: m p m qm qp396: m p m qp qm540: m p qm m qp156: m p qm qp m924: m p qm qp p732: m p qm p qp300: m p qp m qm108: m p qp qm m876: m p qp qm p492: m p qp p qm636: m p p qm qp444: m p p qp qm513: qm m m m qp129: qm m m qp m897: qm m m qp p705: qm m m p qp33: qm m qp m m801: qm m qp m p225: qm m qp p m993: qm m qp p p561: qm m p m qp177: qm m p qp m945: qm m p qp p753: qm m p p qp9: qm qp m m m777: qm qp m m p201: qm qp m p m969: qm qp m p p57: qm qp p m m825: qm qp p m p249: qm qp p p m1017: qm qp p p p525: qm p m m qp141: qm p m qp m909: qm p m qp p717: qm p m p qp45: qm p qp m m813: qm p qp m p237: qm p qp p m1005: qm p qp p p573: qm p p m qp189: qm p p qp m957: qm p p qp p765: qm p p p qp258: qp m m m qm66: qp m m qm m834: qp m m qm p450: qp m m p qm18: qp m qm m m786: qp m qm m p210: qp m qm p m978: qp m qm p p306: qp m p m qm114: qp m p qm m882: qp m p qm p498: qp m p p qm6: qp qm m m m774: qp qm m m p198: qp qm m p m966: qp qm m p p54: qp qm p m m822: qp qm p m p246: qp qm p p m1014: qp qm p p p270: qp p m m qm78: qp p m qm m846: qp p m qm p462: qp p m p qm30: qp p qm m m798: qp p qm m p222: qp p qm p m990: qp p qm p p318: qp p p m qm126: qp p p qm m894: qp p p qm p510: qp p p p qm579: p m m qm qp387: p m m qp qm531: p m qm m qp147: p m qm qp m915: p m qm qp p723: p m qm p qp291: p m qp m qm99: p m qp qm m867: p m qp qm p483: p m qp p qm627: p m p qm qp435: p m p qp qm519: p qm m m qp135: p qm m qp m903: p qm m qp p711: p qm m p qp39: p qm qp m m807: p qm qp m p231: p qm qp p m999: p qm qp p p567: p qm p m qp183: p qm p qp m951: p qm p qp p759: p qm p p qp267: p qp m m qm75: p qp m qm m843: p qp m qm p459: p qp m p qm27: p qp qm m m795: p qp qm m p219: p qp qm p m987: p qp qm p p315: p qp p m qm123: p qp p qm m891: p qp p qm p507: p qp p p qm591: p p m qm qp399: p p m qp qm543: p p qm m qp159: p p qm qp m927: p p qm qp p735: p p qm p qp303: p p qp m qm111: p p qp qm m879: p p qp qm p495: p p qp p qm639: p p p qm qp447: p p p qp qm*/

