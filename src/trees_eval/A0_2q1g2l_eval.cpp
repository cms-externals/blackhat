/* The 2 quark 2 lepton amplitudes use a base 6 code to enumerate them
   m,qm,qp,p,lm,lp corresponds to 0,1,2,3,5,6.
   ep and em are left out of the code.
   For example, (qp,qm,m,p,lm,lp) = 2+ 1*6 + 0*36 + 3*216 + 4*1296 + 5*7776
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"



using namespace std;

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i,j,k) ep.s(i-1,j-1, k-1)

namespace BH {


template <class T> complex<T> A2q1g2lzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q1g2l322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/
 (ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,1))
); }

template <class T> complex<T> A2q1g2l1237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/
 (ep.spb(3,2)*ep.spb(4,0)*ep.spb(4,1))
); }

template <class T> complex<T> A2q1g2l1402_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l1502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l1762_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,2))
); }

template <class T> complex<T> A2q1g2l1870_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(2,4))
); }

template <class T> complex<T> A2q1g2l1932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,2),2))/
 (ep.spb(2,1)*ep.spb(3,0)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l1935_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,4),2))/
 (ep.spa(0,3)*ep.spa(0,4)*ep.spa(1,2))
); }

template <class T> complex<T> A2q1g2l2050_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,4),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l2150_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,3)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l2532_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l2535_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l2662_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l2797_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,2),2))/(ep.spb(2,1)*
ep.spb(3,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l2842_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,2))
); }

template <class T> complex<T> A2q1g2l2950_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(2,4))
); }

template <class T> complex<T> A2q1g2l3310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,1)*ep.spa(2,3)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l3445_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l3817_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,1))
); }

template <class T> complex<T> A2q1g2l3822_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(4,3),2))/(ep.spb(1,0)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l3825_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l3835_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,1)*
ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l4210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(2,4)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l4525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l5237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(3,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l5417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l5477_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,1)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l5495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(1,3))
); }

template <class T> complex<T> A2q1g2l5525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(3,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l5627_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,2)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l5657_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2))/(ep.spb(2,1)*
ep.spb(3,1)*ep.spb(4,0))
); }

template <class T> complex<T> A2q1g2l5675_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(1,3))
); }

template <class T> complex<T> A2q1g2l5735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,4),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l5885_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,4),2))/(ep.spa(0,4)*
ep.spa(1,3)*ep.spa(2,3))
); }

template <class T> complex<T> A2q1g2l7352_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(2,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7357_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(2,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7382_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,0),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,1),2))/
 (ep.spb(1,0)*ep.spb(2,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7395_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(0,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l7400_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l7417_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(4,2),2))/
 (ep.spb(1,0)*ep.spb(2,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q1g2l7425_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/
 (ep.spa(0,1)*ep.spa(0,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l7435_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l7460_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/
 (ep.spa(0,2)*ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q1g2l7465_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,2)*ep.spa(1,2)*ep.spa(3,4))
 ); }


template <class T> complex<T>  (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 322 : return &A2q1g2l322_eval;
case 422 : return &A2q1g2l422_eval;
case 637 : return &A2q1g2l637_eval;
case 1232 : return &A2q1g2l1232_eval;
case 1237 : return &A2q1g2l1237_eval;
case 1402 : return &A2q1g2l1402_eval;
case 1502 : return &A2q1g2l1502_eval;
case 1762 : return &A2q1g2l1762_eval;
case 1870 : return &A2q1g2l1870_eval;
case 1932 : return &A2q1g2l1932_eval;
case 1935 : return &A2q1g2l1935_eval;
case 2050 : return &A2q1g2l2050_eval;
case 2150 : return &A2q1g2l2150_eval;
case 2532 : return &A2q1g2l2532_eval;
case 2535 : return &A2q1g2l2535_eval;
case 2662 : return &A2q1g2l2662_eval;
case 2797 : return &A2q1g2l2797_eval;
case 2842 : return &A2q1g2l2842_eval;
case 2950 : return &A2q1g2l2950_eval;
case 3310 : return &A2q1g2l3310_eval;
case 3445 : return &A2q1g2l3445_eval;
case 3817 : return &A2q1g2l3817_eval;
case 3822 : return &A2q1g2l3822_eval;
case 3825 : return &A2q1g2l3825_eval;
case 3835 : return &A2q1g2l3835_eval;
case 4210 : return &A2q1g2l4210_eval;
case 4310 : return &A2q1g2l4310_eval;
case 4525 : return &A2q1g2l4525_eval;
case 5120 : return &A2q1g2l5120_eval;
case 5125 : return &A2q1g2l5125_eval;
case 5237 : return &A2q1g2l5237_eval;
case 5417 : return &A2q1g2l5417_eval;
case 5477 : return &A2q1g2l5477_eval;
case 5495 : return &A2q1g2l5495_eval;
case 5525 : return &A2q1g2l5525_eval;
case 5627 : return &A2q1g2l5627_eval;
case 5657 : return &A2q1g2l5657_eval;
case 5675 : return &A2q1g2l5675_eval;
case 5735 : return &A2q1g2l5735_eval;
case 5885 : return &A2q1g2l5885_eval;
case 7352 : return &A2q1g2l7352_eval;
case 7357 : return &A2q1g2l7357_eval;
case 7382 : return &A2q1g2l7382_eval;
case 7392 : return &A2q1g2l7392_eval;
case 7395 : return &A2q1g2l7395_eval;
case 7400 : return &A2q1g2l7400_eval;
case 7417 : return &A2q1g2l7417_eval;
case 7422 : return &A2q1g2l7422_eval;
case 7425 : return &A2q1g2l7425_eval;
case 7435 : return &A2q1g2l7435_eval;
case 7460 : return &A2q1g2l7460_eval;
case 7465 : return &A2q1g2l7465_eval;

default: return 0;

}
}


template complex<R>  (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q1g2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
7422: m qm qp em ep
3822: m qm em ep qp
7392: m qp qm em ep
2532: m qp em ep qm
1932: m em ep qp qm
7417: qm m qp em ep
3817: qm m em ep qp
7357: qm qp m em ep
7465: qm qp p em ep
1237: qm qp em ep m
5125: qm qp em ep p
7435: qm p qp em ep
3835: qm p em ep qp
2797: qm em ep m qp
637: qm em ep qp m
4525: qm em ep qp p
3445: qm em ep p qp
7382: qp m qm em ep
7352: qp qm m em ep
7460: qp qm p em ep
1232: qp qm em ep m
5120: qp qm em ep p
7400: qp p qm em ep
1502: qp em ep m qm
422: qp em ep qm m
4310: qp em ep qm p
2150: qp em ep p qm
7425: p qm qp em ep
3825: p qm em ep qp
7395: p qp qm em ep
2535: p qp em ep qm
1935: p em ep qp qm
2842: em ep m qm qp
1762: em ep m qp qm
2662: em ep qm m qp
3310: em ep qm p qp
1402: em ep qp m qm
322: em ep qp qm m
4210: em ep qp qm p
2050: em ep qp p qm
2950: em ep p qm qp
1870: em ep p qp qm
5657: ep m qm qp em
5477: ep m qp qm em
5627: ep qm m qp em
5735: ep qm p qp em
5417: ep qp m qm em
5237: ep qp qm m em
5885: ep qp qm p em
5525: ep qp p qm em
5675: ep p qm qp em
5495: ep p qp qm em
*/

