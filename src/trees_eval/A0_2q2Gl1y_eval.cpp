/* based on A0_2q2Q1y_eval */

/* The 2q 2Q Ny amplitudes use a base 6 code to enumerate them
   m,qm,qp,p,Qm,Qp corresponds to 0,1,2,3,5.
   For example, (qp,qm,m,p,Qm,Qp) = 2+ 1*6 + 0*36 + 3*216 + 4*1296 + 5*7776
  See bottom of file for table of switch numbers and values
*/

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"
#include "BH_error.h"
#include <iostream>

using namespace std;

#define SPA(i,j) ep.spa(i-1,j-1)
#define SPB(i,j) ep.spb(i-1,j-1)
#define S(i,j) ep.s(i-1,j-1)
#define SS(i1,i2,i3) ep.s(ind[i1-1],ind[i2-1], ind[i3-1])
#define SSS(i1,i2,i3,i4) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1])
#define SSSS(i1,i2,i3,i4,i5) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1], ind[i5-1])
#define SSSSS(i1,i2,i3,i4,i5,i6) ep.s(ind[i1-1],ind[i2-1], ind[i3-1], ind[i4-1], ind[i5-1], ind[i6-1], )

namespace BH {


template <class T> complex<T> A2q2G1yzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }



template <class T> complex<T> A2q2Q1y1052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,1)*ep.spb(4,2))-
 (complex<T>(0,1)*pow(ep.spb(2,0),2)*ep.spb(3,0))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,0)*ep.spb(4,3))+
 (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(1,0)*
 ep.spb(4,2)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2Q1y1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,0),2)*ep.spb(2,1))/
(ep.spb(1,0)*ep.spb(3,2)*
 ep.spb(4,1)*ep.spb(4,2))-
 (complex<T>(0,1)*pow(ep.spb(3,0),3))/(ep.spb(1,0)*
 ep.spb(3,2)*ep.spb(4,0)*
 ep.spb(4,3))+(complex<T>(0,1)*pow(ep.spb(3,0),2))/
(ep.spb(1,0)*ep.spb(4,2)*
 ep.spb(4,3))
); }



template <class T> complex<T> A_GGqqy_pmpmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(4,2),2))/
 (SPA(3,5)*SPA(4,5)*SPA(1,2))
); }

template <class T> complex<T> A_GGqqy_pmmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(3,2),2))/
 (SPA(3,5)*SPA(4,5)*SPA(1,2))
); }

template <class T> complex<T> A_GGqqy_mppmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(4,1),2))/
 (SPA(3,5)*SPA(4,5)*SPA(1,2))
); }

template <class T> complex<T> A_GGqqy_mpmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(3,1),2))/
 (SPA(3,5)*SPA(4,5)*SPA(1,2))
 ); }





template <class T> complex<T> A_GGqqy_mpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(4,2),2))/
 (SPB(3,5)*SPB(4,5)*SPB(1,2))
); }

template <class T> complex<T> A_GGqqy_mppmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(3,2),2))/
 (SPB(3,5)*SPB(4,5)*SPB(1,2))
); }

template <class T> complex<T> A_GGqqy_pmmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(4,1),2))/
 (SPB(3,5)*SPB(4,5)*SPB(1,2))
); }

template <class T> complex<T> A_GGqqy_pmpmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(3,1),2))/
 (SPB(3,5)*SPB(4,5)*SPB(1,2))
 ); }

template <class T> complex<T> A_qqGGy_pmpmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(2,4),2))/
 (SPA(1,5)*SPA(2,5)*SPA(3,4))
); }

template <class T> complex<T> A_qqGGy_mppmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(1,4),2))/
 (SPA(1,5)*SPA(2,5)*SPA(3,4))
); }

template <class T> complex<T> A_qqGGy_pmmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(2,3),2))/
 (SPA(1,5)*SPA(2,5)*SPA(3,4))
); }

template <class T> complex<T> A_qqGGy_mpmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(1,3),2))/
 (SPA(1,5)*SPA(2,5)*SPA(3,4))
 ); }


template <class T> complex<T> A_qqGGy_mpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(2,4),2))/
 (SPB(1,5)*SPB(2,5)*SPB(3,4))
); }

template <class T> complex<T> A_qqGGy_pmmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(1,4),2))/
 (SPB(1,5)*SPB(2,5)*SPB(3,4))
); }

template <class T> complex<T> A_qqGGy_mppmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(2,3),2))/
 (SPB(1,5)*SPB(2,5)*SPB(3,4))
); }

template <class T> complex<T> A_qqGGy_pmpmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(1,3),2))/
 (SPB(1,5)*SPB(2,5)*SPB(3,4))
 ); }

template <class T> complex<T> A_qGGqy_mpmpm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(4,2),2))/(SPB(3,2)*
SPB(5,1)*SPB(5,4))
); }

template <class T> complex<T>  A_qGGqy_mmppm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(4,3),2))/(SPB(3,2)*
SPB(5,1)*SPB(5,4))
); }

template <class T> complex<T>  A_qGGqy_ppmmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,1)*pow(SPB(2,1),2))/(SPB(3,2)*
SPB(5,1)*SPB(5,4))
); }

template <class T> complex<T>  A_qGGqy_pmpmm_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(SPB(3,1),2))/(SPB(3,2)*
SPB(5,1)*SPB(5,4))
); }


template <class T> complex<T> A_qGGqy_pmpmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(4,2),2))/(SPA(3,2)*
SPA(5,1)*SPA(5,4))
); }

template <class T> complex<T>  A_qGGqy_ppmmp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(4,3),2))/(SPA(3,2)*
SPA(5,1)*SPA(5,4))
); }

template <class T> complex<T>  A_qGGqy_mmppp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,-1)*pow(SPA(2,1),2))/(SPA(3,2)*
SPA(5,1)*SPA(5,4))
); }

template <class T> complex<T>  A_qGGqy_mpmpp_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(SPA(3,1),2))/(SPA(3,2)*
SPA(5,1)*SPA(5,4))
); }





template <class T> complex<T>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {

   case 30344: return &A_qqGGy_mpmpp_eval;
   case 29896: return &A_qqGGy_mppmp_eval;
   case 30337: return &A_qqGGy_pmmpp_eval;
   case 29889: return &A_qqGGy_pmpmp_eval;

   case 25793: return &A_qqGGy_pmpmm_eval;
   case 26241: return &A_qqGGy_pmmpm_eval;
   case 25800: return &A_qqGGy_mppmm_eval;
   case 26248: return &A_qqGGy_mpmpm_eval;

   case 25296: return &A_qGGqy_mmppm_eval;
   case 25240: return &A_qGGqy_mpmpm_eval;
   case 24785: return &A_qGGqy_pmpmm_eval;
   case 24729: return &A_qGGqy_ppmmm_eval;

   case 28825: return &A_qGGqy_ppmmp_eval;
   case 28881: return &A_qGGqy_pmpmp_eval;
   case 29336: return &A_qGGqy_mpmpp_eval;
   case 29392: return &A_qGGqy_mmppp_eval;

   case 29210: return &A_GGqqy_mpmpp_eval;
   case 28762: return &A_GGqqy_mppmp_eval;
   case 29203: return &A_GGqqy_pmmpp_eval;
   case 28755: return &A_GGqqy_pmpmp_eval;

   case 25114: return &A_GGqqy_mpmpm_eval;
   case 24666: return &A_GGqqy_mppmm_eval;
   case 25107: return &A_GGqqy_pmmpm_eval;
   case 24659: return &A_GGqqy_pmpmm_eval;

   default: return 0;

}
}

template complex<R>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);


#if BH_USE_GMP
template complex<RGMP>  (*A2q2G1y_Tree_Ptr_eval(long hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
5125: qm qp Qm Qp gap
4945: qm qp Qp Qm gap
637: qm Qm Qp qp gam
607: qm Qp Qm qp gam
4495: qm Qp Qm qp gap
1232: qp qm Qm Qp gam
5120: qp qm Qm Qp gap
1052: qp qm Qp Qm gam
4940: qp qm Qp Qm gap
1142: qp Qm qm Qp gam
5030: qp Qm qm Qp gap
422: qp Qm Qp qm gam
4310: qp Qm Qp qm gap
932: qp Qp qm Qm gam
4820: qp Qp qm Qm gap
392: qp Qp Qm qm gam
4280: qp Qp Qm qm gap
*/
