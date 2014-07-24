/* The 2 quark amplitudes use a base 4 code to enumerate them
   m,qm,qp,p corresponds to 0,1,2,3.
   For example, (m,qm,p,qp) = 0 + 1*4 + 3*16 + 2*64
See bottom of file for table of switch numbers and values
*/



#include "amplitudes_tree_eval.h"
 #include "eval_param.h"
using namespace std;

namespace BH {


template <class T> complex<T> A2q1g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spa(1,2),2))/ep.spa(0,1)
                      ); }

template <class T> complex<T> A2q1g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spa(0,2),2))/ep.spa(0,1)
                      ); }

template <class T> complex<T> A2q1g18_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spa(1,2),2))/ep.spa(0,2)
                      ); }

template <class T> complex<T> A2q1g24_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spa(0,2),2))/ep.spa(1,2)
                      ); }

template <class T> complex<T> A2q1g27_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spb(1,0),2))/ep.spb(2,1)
                      ); }

template <class T> complex<T> A2q1g30_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spb(1,0),2))/ep.spb(2,0)
                      ); }

template <class T> complex<T> A2q1g33_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spa(0,1),2))/ep.spa(0,2)
                      ); }

template <class T> complex<T> A2q1g36_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spa(0,1),2))/ep.spa(1,2)
                      ); }

template <class T> complex<T> A2q1g39_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spb(2,0),2))/ep.spb(2,1)
                      ); }

template <class T> complex<T> A2q1g45_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,-1.)*pow(ep.spb(2,1),2))/ep.spb(2,0)
                      ); }

template <class T> complex<T> A2q1g54_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spb(2,0),2))/ep.spb(1,0)
                      ); }

template <class T> complex<T> A2q1g57_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return( (complex<T>(0.,1.)*pow(ep.spb(2,1),2))/ep.spb(1,0)
                      ); }



template <class T> complex<T> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	 switch (hc) {
	   case 6 : return &A2q1g6_eval;
	   case 9 : return &A2q1g9_eval;
	   case 18: return &A2q1g18_eval;
	   case 24: return &A2q1g24_eval;
	   case 27: return &A2q1g27_eval;
	   case 30: return &A2q1g30_eval;
	   case 33: return &A2q1g33_eval;
	   case 36: return &A2q1g36_eval;
	   case 39: return &A2q1g39_eval;
	   case 45: return &A2q1g45_eval;
	   case 54: return &A2q1g54_eval;
	   case 57: return &A2q1g57_eval;

	 }
	 // should never get there
	 return 0;
 }


 template complex<R> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&) ;
 template complex<RHP> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&) ;
 template complex<RVHP> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*A2q1g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&) ;
#endif

 }


/* *************** table of switch values ************* */

/*
36: mqmqp
24: mqpqm
33: qmmqp
9: qmqpm
57: qmqpp
45: qmpqp
18: qpmqm
6: qpqmm
54: qpqmp
30: qppqm
39: pqmqp
27: pqpqm
*/

