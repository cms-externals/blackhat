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


template <class T> complex<T> A2q2Qzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q2Q317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   (complex<T>(0,-1)*pow(ep.spa(1,3),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q377_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(1,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A2q2Q412_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A2q2Q497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,0),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A2q2Q502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/(ep.spb(1,0)*
ep.spb(3,2))
); }

template <class T> complex<T> A2q2Q587_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2Q622_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2Q917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2Q932_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(2,0),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A2q2Q967_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q1052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),2))/
 (ep.spa(0,1)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q1057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(2,1),2))/
 (ep.spb(1,0)*ep.spb(3,2))
); }

template <class T> complex<T> A2q2Q1132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),2))/(ep.spa(0,3)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2Q1142_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q1162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spb(3,2),2))/(ep.spb(2,1)*
ep.spb(3,0))
); }

template <class T> complex<T> A2q2Q1177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2Q1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2))/
 (ep.spa(0,1)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2Q1237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spb(3,1),2))/
 (ep.spb(1,0)*ep.spb(3,2))
 ); }


template <class T> complex<T>  (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 317 : return &A2q2Q317_eval;
case 322 : return &A2q2Q322_eval;
case 377 : return &A2q2Q377_eval;
case 392 : return &A2q2Q392_eval;
case 412 : return &A2q2Q412_eval;
case 422 : return &A2q2Q422_eval;
case 497 : return &A2q2Q497_eval;
case 502 : return &A2q2Q502_eval;
case 587 : return &A2q2Q587_eval;
case 607 : return &A2q2Q607_eval;
case 622 : return &A2q2Q622_eval;
case 637 : return &A2q2Q637_eval;
case 917 : return &A2q2Q917_eval;
case 932 : return &A2q2Q932_eval;
case 947 : return &A2q2Q947_eval;
case 967 : return &A2q2Q967_eval;
case 1052 : return &A2q2Q1052_eval;
case 1057 : return &A2q2Q1057_eval;
case 1132 : return &A2q2Q1132_eval;
case 1142 : return &A2q2Q1142_eval;
case 1162 : return &A2q2Q1162_eval;
case 1177 : return &A2q2Q1177_eval;
case 1232 : return &A2q2Q1232_eval;
case 1237 : return &A2q2Q1237_eval;

default: return 0;;

}
}


template complex<R>  (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2Q_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
1237: qmqpQmQp
1057: qmqpQpQm
1177: qmQmqpQp
637: qmQmQpqp
967: qmQpqpQm
607: qmQpQmqp
1232: qpqmQmQp
1052: qpqmQpQm
1142: qpQmqmQp
422: qpQmQpqm
932: qpQpqmQm
392: qpQpQmqm
1162: QmqmqpQp
622: QmqmQpqp
1132: QmqpqmQp
412: QmqpQpqm
502: QmQpqmqp
322: QmQpqpqm
947: QpqmqpQm
587: QpqmQmqp
917: QpqpqmQm
377: QpqpQmqm
497: QpQmqmqp
317: QpQmqpqm
*/

