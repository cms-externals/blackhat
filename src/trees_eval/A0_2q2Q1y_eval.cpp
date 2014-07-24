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


template <class T> complex<T> A2q2Q1yzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q2Q1y392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2Q1y422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2Q1y607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,1),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2Q1y637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(4,0)*ep.spb(4,3))
); }

template <class T> complex<T> A2q2Q1y932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

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

template <class T> complex<T> A2q2Q1y1142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
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

template <class T> complex<T> A2q2Q1y4280_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2Q1y4310_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2Q1y4495_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,4)*
ep.spa(1,2)*ep.spa(3,4))
); }

template <class T> complex<T> A2q2Q1y4820_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q1y4940_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q1y4945_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q1y5030_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q1y5120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q1y5125_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),2))/
 (ep.spa(0,4)*ep.spa(1,4)*ep.spa(2,3))
 ); }


template <class T> complex<T>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 392 : return &A2q2Q1y392_eval;
case 422 : return &A2q2Q1y422_eval;
case 607 : return &A2q2Q1y607_eval;
case 637 : return &A2q2Q1y637_eval;
case 932 : return &A2q2Q1y932_eval;
case 1052 : return &A2q2Q1y1052_eval;
case 1142 : return &A2q2Q1y1142_eval;
case 1232 : return &A2q2Q1y1232_eval;
case 4280 : return &A2q2Q1y4280_eval;
case 4310 : return &A2q2Q1y4310_eval;
case 4495 : return &A2q2Q1y4495_eval;
case 4820 : return &A2q2Q1y4820_eval;
case 4940 : return &A2q2Q1y4940_eval;
case 4945 : return &A2q2Q1y4945_eval;
case 5030 : return &A2q2Q1y5030_eval;
case 5120 : return &A2q2Q1y5120_eval;
case 5125 : return &A2q2Q1y5125_eval;

default: return 0;

}
}


template complex<R>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2Q1y_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
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

