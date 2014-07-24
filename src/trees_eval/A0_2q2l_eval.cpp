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

namespace BH {

template <class T> complex<T> A2q2l322_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPA(1,4),2))/(SPA(1,2)*SPA(3,4))
                      );}

template <class T> complex<T> A2q2l422_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPA(2,4),2))/(SPA(1,4)*SPA(2,3))
                      );}

template <class T> complex<T> A2q2l502_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPA(1,3),2))/(SPA(1,2)*
SPA(3,4))
                      );}

template <class T> complex<T> A2q2l637_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPA(1,2),2))/(SPA(1,4)*SPA(2,3)));
}

template <class T> complex<T> A2q2l917_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPA(3,4),2))/(SPA(1,4)*SPA(2,3)));
}

template <class T> complex<T> A2q2l947_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPA(2,4),2))/(SPA(1,4)*SPA(2,3))
                      );}

template <class T> complex<T> A2q2l1232_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPA(2,3),2))/ (SPA(1,2)*SPA(3,4))
                      );}

template <class T> complex<T> A2q2l1237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPA(1,3),2))/ (SPA(1,2)*SPA(3,4))
                      );}


template <class T> complex<T> A2q2l497_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPB(4,1),2))/(SPB(2,1)*SPB(4,3))
                      );}

template <class T> complex<T> A2q2l607_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPB(4,2),2))/(SPB(4,1)*SPB(3,2))
                      );}

template <class T> complex<T> A2q2l317_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPB(3,1),2))/(SPB(2,1)*
SPB(4,3))
                      );}

template <class T> complex<T> A2q2l392_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPB(2,1),2))/(SPB(4,1)*SPB(3,2)));
}

template <class T> complex<T> A2q2l1162_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPB(4,3),2))/(SPB(4,1)*SPB(3,2)));
}

template <class T> complex<T> A2q2l1132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,1.)*pow(SPB(4,2),2))/(SPB(4,1)*SPB(3,2))
                      );}

template <class T> complex<T> A2q2l1057_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPB(3,2),2))/ (SPB(2,1)*SPB(4,3))
                      );}

template <class T> complex<T> A2q2l1052_eval(const eval_param<T>& ep, const mass_param_coll& masses){
	return( (complex<T>(0.,-1.)*pow(SPB(3,1),2))/ (SPB(2,1)*SPB(4,3))
                      );}



#define _CASE_A2q2l(K) case K : return A2q2l ## K ## _eval (ep,masses)


#define _CASE_A2q2l_Ptr(K) case K : return A2q2l ## K ## _eval

template <class T> complex<T>  (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
   _CASE_A2q2l_Ptr(322);
   _CASE_A2q2l_Ptr(422);
   _CASE_A2q2l_Ptr(502);
   _CASE_A2q2l_Ptr(637);
   _CASE_A2q2l_Ptr(917);
   _CASE_A2q2l_Ptr(947);
   _CASE_A2q2l_Ptr(1232);
   _CASE_A2q2l_Ptr(1237);

   _CASE_A2q2l_Ptr(497);
   _CASE_A2q2l_Ptr(607);
   _CASE_A2q2l_Ptr(317);
   _CASE_A2q2l_Ptr(392);
   _CASE_A2q2l_Ptr(1162);
   _CASE_A2q2l_Ptr(1132);
   _CASE_A2q2l_Ptr(1057);
   _CASE_A2q2l_Ptr(1052);

default: return 0;

}
}

template complex<R>  (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2l_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
1237: qmqpemep
637: qmemepqp
1232: qpqmemep
422: qpemepqm
502: emepqmqp
322: emepqpqm
947: epqmqpem
917: epqpqmem
*/

