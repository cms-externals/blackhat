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


template <class T> complex<T> A4qzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A4q280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q315_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q340_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q350_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q370_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(2,3))-(complex<T>(0,1)*pow(ep.spb(2,0),2))/
(ep.spb(2,1)*ep.spb(3,0))
); }

template <class T> complex<T> A4q375_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q380_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q385_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q490_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A4q525_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A4q550_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q555_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
 ep.spa(1,2))-(complex<T>(0,1)*pow(ep.spb(3,1),2))/
(ep.spb(1,0)*ep.spb(3,2))
); }

template <class T> complex<T> A4q560_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q565_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q585_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A4q595_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q700_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q710_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q730_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q735_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q740_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
 ep.spa(2,3))+(complex<T>(0,1)*pow(ep.spb(2,0),2))/
(ep.spb(2,1)*ep.spb(3,0))
); }

template <class T> complex<T> A4q745_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q770_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q800_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A4q805_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q910_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q915_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A4q920_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A4q925_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
 ep.spa(1,2))+(complex<T>(0,1)*pow(ep.spb(3,1),2))/
(ep.spb(1,0)*ep.spb(3,2))
); }

template <class T> complex<T> A4q945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q955_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A4q980_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A4q985_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A4q1015_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2))
 ); }


template <class T> complex<T>  (*A4q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 280 : return &A4q280_eval;
case 310 : return &A4q310_eval;
case 315 : return &A4q315_eval;
case 340 : return &A4q340_eval;
case 350 : return &A4q350_eval;
case 370 : return &A4q370_eval;
case 375 : return &A4q375_eval;
case 380 : return &A4q380_eval;
case 385 : return &A4q385_eval;
case 490 : return &A4q490_eval;
case 495 : return &A4q495_eval;
case 525 : return &A4q525_eval;
case 550 : return &A4q550_eval;
case 555 : return &A4q555_eval;
case 560 : return &A4q560_eval;
case 565 : return &A4q565_eval;
case 585 : return &A4q585_eval;
case 595 : return &A4q595_eval;
case 700 : return &A4q700_eval;
case 710 : return &A4q710_eval;
case 730 : return &A4q730_eval;
case 735 : return &A4q735_eval;
case 740 : return &A4q740_eval;
case 745 : return &A4q745_eval;
case 770 : return &A4q770_eval;
case 800 : return &A4q800_eval;
case 805 : return &A4q805_eval;
case 910 : return &A4q910_eval;
case 915 : return &A4q915_eval;
case 920 : return &A4q920_eval;
case 925 : return &A4q925_eval;
case 945 : return &A4q945_eval;
case 955 : return &A4q955_eval;
case 980 : return &A4q980_eval;
case 985 : return &A4q985_eval;
case 1015 : return &A4q1015_eval;

default: return 0;
		throw BHerror("case missing for tree amplitude!");

}
}


template complex<R>  (*A4q_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A4q_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A4q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A4q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
1015: qm qm qbp qbp
985: qm qp qbm qbp
805: qm qp qbp qbm
955: qm qbm qp qbp
595: qm qbm qbp qp
925: qm qbp qm qbp
745: qm qbp qp qbm
565: qm qbp qbm qp
385: qm qbp qbp qm
980: qp qm qbm qbp
800: qp qm qbp qbm
770: qp qp qbm qbm
920: qp qbm qm qbp
740: qp qbm qp qbm
560: qp qbm qbm qp
380: qp qbm qbp qm
710: qp qbp qm qbm
350: qp qbp qbm qm
945: qbm qm qp qbp
585: qbm qm qbp qp
915: qbm qp qm qbp
735: qbm qp qp qbm
555: qbm qp qbm qp
375: qbm qp qbp qm
525: qbm qbm qp qp
495: qbm qbp qm qp
315: qbm qbp qp qm
910: qbp qm qm qbp
730: qbp qm qp qbm
550: qbp qm qbm qp
370: qbp qm qbp qm
700: qbp qp qm qbm
340: qbp qp qbm qm
490: qbp qbm qm qp
310: qbp qbm qp qm
280: qbp qbp qm qm
*/

