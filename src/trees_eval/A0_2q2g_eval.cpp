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


template <class T> complex<T> A2q2gzero_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return complex<T>(0,0); }


template <class T> complex<T> A2q2g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(   complex<T>(0,0)
); }

template <class T> complex<T> A2q2g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g18_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g24_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g27_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,3)*ep.spa(1,2))
); }

template <class T> complex<T> A2q2g30_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2g33_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g36_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g39_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,3)*ep.spa(1,2))
); }

template <class T> complex<T> A2q2g45_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(0,1)*
ep.spa(1,2))
); }

template <class T> complex<T> A2q2g54_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g57_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g66_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g72_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g75_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2))/(ep.spa(0,1)*
ep.spa(0,3))
); }

template <class T> complex<T> A2q2g78_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(2,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,3)*ep.spa(1,2))
); }

template <class T> complex<T> A2q2g96_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g99_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,3),3))/(ep.spa(0,1)*
ep.spa(0,3)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g108_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g111_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g114_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,3),3))/(ep.spa(0,3)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g120_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,3),2))/(ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2g123_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g126_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g129_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g132_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g135_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,1)*
ep.spa(0,3))
); }

template <class T> complex<T> A2q2g141_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,3)*ep.spa(1,2))
); }

template <class T> complex<T> A2q2g144_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g147_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(1,3))/
 (ep.spa(0,1)*ep.spa(0,3)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g156_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g159_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g177_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2)*ep.spa(1,3))/
 (ep.spa(0,3)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g180_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2))/(ep.spa(1,2)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2g183_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g189_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g198_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(1,2),2)*ep.spa(0,2))/
 (ep.spa(0,1)*ep.spa(0,3)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g201_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,1)*
ep.spa(0,3)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g210_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(1,2),2))/(ep.spa(0,3)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2g216_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,2),3))/(ep.spa(0,3)*
ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g219_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g222_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g225_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,1)*pow(ep.spa(0,1),2))/(ep.spa(0,3)*
ep.spa(2,3))
); }

template <class T> complex<T> A2q2g228_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  (complex<T>(0,-1)*pow(ep.spa(0,1),2)*ep.spa(0,2))/
 (ep.spa(0,3)*ep.spa(1,2)*ep.spa(2,3))
); }

template <class T> complex<T> A2q2g231_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g237_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g246_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
); }

template <class T> complex<T> A2q2g249_eval(const eval_param<T>& ep, const mass_param_coll& masses){
 return(  complex<T>(0,0)
 ); }


template <class T> complex<T>  (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {

   switch (hc) {
case 6 : return &A2q2g6_eval;
case 9 : return &A2q2g9_eval;
case 18 : return &A2q2g18_eval;
case 24 : return &A2q2g24_eval;
case 27 : return &A2q2g27_eval;
case 30 : return &A2q2g30_eval;
case 33 : return &A2q2g33_eval;
case 36 : return &A2q2g36_eval;
case 39 : return &A2q2g39_eval;
case 45 : return &A2q2g45_eval;
case 54 : return &A2q2g54_eval;
case 57 : return &A2q2g57_eval;
case 66 : return &A2q2g66_eval;
case 72 : return &A2q2g72_eval;
case 75 : return &A2q2g75_eval;
case 78 : return &A2q2g78_eval;
case 96 : return &A2q2g96_eval;
case 99 : return &A2q2g99_eval;
case 108 : return &A2q2g108_eval;
case 111 : return &A2q2g111_eval;
case 114 : return &A2q2g114_eval;
case 120 : return &A2q2g120_eval;
case 123 : return &A2q2g123_eval;
case 126 : return &A2q2g126_eval;
case 129 : return &A2q2g129_eval;
case 132 : return &A2q2g132_eval;
case 135 : return &A2q2g135_eval;
case 141 : return &A2q2g141_eval;
case 144 : return &A2q2g144_eval;
case 147 : return &A2q2g147_eval;
case 156 : return &A2q2g156_eval;
case 159 : return &A2q2g159_eval;
case 177 : return &A2q2g177_eval;
case 180 : return &A2q2g180_eval;
case 183 : return &A2q2g183_eval;
case 189 : return &A2q2g189_eval;
case 198 : return &A2q2g198_eval;
case 201 : return &A2q2g201_eval;
case 210 : return &A2q2g210_eval;
case 216 : return &A2q2g216_eval;
case 219 : return &A2q2g219_eval;
case 222 : return &A2q2g222_eval;
case 225 : return &A2q2g225_eval;
case 228 : return &A2q2g228_eval;
case 231 : return &A2q2g231_eval;
case 237 : return &A2q2g237_eval;
case 246 : return &A2q2g246_eval;
case 249 : return &A2q2g249_eval;

default: _WARNING3("Unknown pointer amplitude (*A2q2g_Tree_Ptr(int hc)) - case:",hc , " - throw BH error.");
		throw BHerror("case missing for tree amplitude!");

}
}


template complex<R>  (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<R>& ep, const mass_param_coll& masses);
template complex<RHP>  (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<RHP>& ep, const mass_param_coll& masses);
template complex<RVHP>  (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>& ep, const mass_param_coll& masses);

#if BH_USE_GMP

template complex<RGMP>  (*A2q2g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>& ep, const mass_param_coll& masses);
#endif
}


/* *************** table of switch values ************* */

/*
144: mmqmqp
96: mmqpqm
132: mqmmqp
36: mqmqpm
228: mqmqpp
180: mqmpqp
72: mqpmqm
24: mqpqmm
216: mqpqmp
120: mqppqm
156: mpqmqp
108: mpqpqm
129: qmmmqp
33: qmmqpm
225: qmmqpp
177: qmmpqp
9: qmqpmm
201: qmqpmp
57: qmqppm
249: qmqppp
141: qmpmqp
45: qmpqpm
237: qmpqpp
189: qmppqp
66: qpmmqm
18: qpmqmm
210: qpmqmp
114: qpmpqm
6: qpqmmm
198: qpqmmp
54: qpqmpm
246: qpqmpp
78: qppmqm
30: qppqmm
222: qppqmp
126: qpppqm
147: pmqmqp
99: pmqpqm
135: pqmmqp
39: pqmqpm
231: pqmqpp
183: pqmpqp
75: pqpmqm
27: pqpqmm
219: pqpqmp
123: pqppqm
159: ppqmqp
111: ppqpqm
*/

