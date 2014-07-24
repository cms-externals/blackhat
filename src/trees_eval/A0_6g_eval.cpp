

#include <complex>
#include <vector>
#include "amplitudes_tree_eval.h"
#include "BH_typedefs.h"
#include "eval_param.h"

using namespace std;

namespace BH {

template <class T> complex<T> A6g0_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g1_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g2_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g3_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(1,0),3))/(ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g4_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g5_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(2,0),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g6_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(2,1),3))/(ep.spb(1,0)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g7_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(5,0,2) + ep.spab(5,1,2),3)/(ep.s(2,3,4)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,0,4) + ep.spab(1,5,4))*ep.spb(3,2)*ep.spb(4,3))) - pow(ep.spab(3,4,0) + ep.spab(3,5,0),3)/(ep.s(0,4,5)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,0,4) + ep.spab(1,5,4))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g8_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g9_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,0),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g10_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,1),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g11_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*((pow(ep.spa(4,5),3)*pow(ep.spb(3,1),4))/(ep.s(1,2,3)*ep.spa(0,5)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(2,1)*ep.spb(3,2)) - (pow(ep.spa(2,4),4)*pow(ep.spb(1,0),3))/(ep.s(0,1,5)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,0,5) + ep.spab(2,1,5))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(5,0)) + pow(ep.spab(2,0,3) + ep.spab(2,1,3),4)/(ep.s(0,1,2)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(2,0,5) + ep.spab(2,1,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g12_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(3,2),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g13_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*((pow(ep.spa(4,5),3)*pow(ep.spb(2,0),4))/(ep.s(0,1,2)*ep.spa(3,4)*(ep.spab(3,1,0) + ep.spab(3,2,0))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(1,0)*ep.spb(2,1)) - (pow(ep.spa(1,5),4)*pow(ep.spb(3,2),3))/(ep.s(0,1,5)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(4,3)) + pow(ep.spab(1,4,0) + ep.spab(1,5,0),4)/(ep.s(1,2,3)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g14_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(4,0,1) + ep.spab(4,5,1),3)/(ep.s(0,1,5)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,0,5) + ep.spab(2,1,5))*ep.spb(1,0)*ep.spb(5,0))) + pow(ep.spab(0,1,3) + ep.spab(0,2,3),3)/(ep.s(3,4,5)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(2,0,5) + ep.spab(2,1,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g15_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(4,5),3))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)); }

template <class T> complex<T> A6g16_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g17_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,0),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g18_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,1),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g19_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*((pow(ep.spa(3,5),4)*pow(ep.spb(1,0),3))/(ep.s(3,4,5)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,4,0) + ep.spab(3,5,0))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(2,1)) + pow(ep.spab(5,2,4) + ep.spab(5,3,4),4)/(ep.s(0,1,5)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(3,2)*ep.spb(4,3)) - (pow(ep.spa(2,3),3)*pow(ep.spb(4,0),4))/(ep.s(0,4,5)*ep.spa(1,2)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(3,4,0) + ep.spab(3,5,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g20_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,2),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g21_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*((pow(ep.spa(3,5),4)*pow(ep.spb(2,0),4))/(ep.s(0,1,2)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,1,0) + ep.spab(3,2,0))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(1,0)*ep.spb(2,1)) + (pow(ep.spa(1,5),4)*pow(ep.spb(4,2),4))/(ep.s(2,3,4)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(3,2)*ep.spb(4,3)) - (pow(ep.spa(1,3),4)*pow(ep.spb(4,0),4))/(ep.s(1,2,3)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g22_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*((pow(ep.spa(3,5),4)*pow(ep.spb(2,1),3))/(ep.s(0,1,2)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,1,0) + ep.spab(3,2,0))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(1,0)) + (pow(ep.spa(0,5),3)*pow(ep.spb(4,2),4))/(ep.s(2,3,4)*ep.spa(0,1)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(3,2)*ep.spb(4,3)) - pow(ep.spab(3,1,4) + ep.spab(3,2,4),4)/(ep.s(1,2,3)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g23_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(3,5),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g24_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(4,3),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g25_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(pow(ep.spab(5,3,0) + ep.spab(5,4,0),4)/(ep.s(3,4,5)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,4,0) + ep.spab(3,5,0))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(1,0)*ep.spb(2,1)) - (pow(ep.spa(1,5),4)*pow(ep.spb(4,3),3))/(ep.s(2,3,4)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(3,2)) + (pow(ep.spa(1,2),3)*pow(ep.spb(4,0),4))/(ep.s(0,4,5)*ep.spa(2,3)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(3,4,0) + ep.spab(3,5,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g26_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-((pow(ep.spa(0,5),3)*pow(ep.spb(3,1),4))/(ep.s(1,2,3)*ep.spa(4,5)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(2,1)*ep.spb(3,2))) + pow(ep.spab(2,0,1) + ep.spab(2,5,1),4)/(ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,0,5) + ep.spab(2,1,5))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(1,0)*ep.spb(5,0)) + (pow(ep.spa(0,2),4)*pow(ep.spb(4,3),3))/(ep.s(0,1,2)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(2,0,5) + ep.spab(2,1,5))*ep.spb(5,4))); }

template <class T> complex<T> A6g27_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,5),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g28_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(pow(ep.spab(5,0,2) + ep.spab(5,1,2),3)/(ep.s(0,1,2)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(1,0)*ep.spb(2,1)) - pow(ep.spab(1,2,4) + ep.spab(1,3,4),3)/(ep.s(0,4,5)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(5,0)*ep.spb(5,4))); }

template <class T> complex<T> A6g29_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,5),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g30_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,5),3))/(ep.spa(0,1)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g31_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g32_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g33_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(5,0),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,4)); }

template <class T> complex<T> A6g34_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(5,1),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g35_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(4,0,1) + ep.spab(4,5,1),3)/(ep.s(1,2,3)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,4,3) + ep.spab(0,5,3))*ep.spb(2,1)*ep.spb(3,2))) + pow(ep.spab(2,3,5) + ep.spab(2,4,5),3)/(ep.s(3,4,5)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,4,3) + ep.spab(0,5,3))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g36_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(5,2),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g37_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-((pow(ep.spa(3,4),3)*pow(ep.spb(2,0),4))/(ep.s(0,1,2)*ep.spa(4,5)*(ep.spab(3,4,0) + ep.spab(3,5,0))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(1,0)*ep.spb(2,1))) - pow(ep.spab(1,0,2) + ep.spab(1,5,2),4)/(ep.s(0,1,5)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(5,0,2) + ep.spab(5,1,2))*ep.spb(3,2)*ep.spb(4,3)) - (pow(ep.spa(1,3),4)*pow(ep.spb(5,0),3))/(ep.s(0,4,5)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,0,4) + ep.spab(1,5,4))*(ep.spab(3,4,0) + ep.spab(3,5,0))*ep.spb(5,4))); }

template <class T> complex<T> A6g38_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-((pow(ep.spa(0,4),4)*pow(ep.spb(2,1),3))/(ep.s(0,4,5)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(3,2))) - (pow(ep.spa(3,4),3)*pow(ep.spb(5,1),4))/(ep.s(0,1,5)*ep.spa(2,3)*(ep.spab(2,0,5) + ep.spab(2,1,5))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(1,0)*ep.spb(5,0)) - pow(ep.spab(0,3,5) + ep.spab(0,4,5),4)/(ep.s(0,1,2)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(2,0,5) + ep.spab(2,1,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g39_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(3,4),3))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(4,5)); }

template <class T> complex<T> A6g40_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(5,3),4))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)*ep.spb(5,4)); }

template <class T> complex<T> A6g41_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(pow(ep.spab(4,1,3) + ep.spab(4,2,3),4)/(ep.s(0,4,5)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(2,1)*ep.spb(3,2)) - (pow(ep.spa(2,4),4)*pow(ep.spb(5,0),3))/(ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,3,5) + ep.spab(2,4,5))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(1,0)) + (pow(ep.spa(1,2),3)*pow(ep.spb(5,3),4))/(ep.s(3,4,5)*ep.spa(0,1)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(2,3,5) + ep.spab(2,4,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g42_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-((pow(ep.spa(0,4),4)*pow(ep.spb(3,1),4))/(ep.s(1,2,3)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(2,1)*ep.spb(3,2))) - (pow(ep.spa(2,4),4)*pow(ep.spb(5,1),4))/(ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,3,5) + ep.spab(2,4,5))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(1,0)*ep.spb(5,0)) - (pow(ep.spa(0,2),4)*pow(ep.spb(5,3),4))/(ep.s(3,4,5)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(2,3,5) + ep.spab(2,4,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g43_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,4),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g44_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-((pow(ep.spa(0,4),4)*pow(ep.spb(3,2),3))/(ep.s(1,2,3)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(2,1))) - pow(ep.spab(4,2,5) + ep.spab(4,3,5),4)/(ep.s(2,3,4)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(2,3,5) + ep.spab(2,4,5))*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(1,0)*ep.spb(5,0)) - (pow(ep.spa(0,1),3)*pow(ep.spb(5,3),4))/(ep.s(3,4,5)*ep.spa(1,2)*(ep.spab(0,1,3) + ep.spab(0,2,3))*(ep.spab(2,3,5) + ep.spab(2,4,5))*ep.spb(4,3)*ep.spb(5,4))); }

template <class T> complex<T> A6g45_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,4),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g46_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,4),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g47_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g48_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spb(5,4),3))/(ep.spb(1,0)*ep.spb(2,1)*ep.spb(3,2)*ep.spb(4,3)*ep.spb(5,0)); }

template <class T> complex<T> A6g49_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(pow(ep.spab(3,4,0) + ep.spab(3,5,0),3)/(ep.s(0,1,2)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(1,0)*ep.spb(2,1)) - pow(ep.spab(1,2,4) + ep.spab(1,3,4),3)/(ep.s(2,3,4)*ep.spa(0,1)*ep.spa(0,5)*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(3,2)*ep.spb(4,3))); }

template <class T> complex<T> A6g50_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(0,4,1) + ep.spab(0,5,1),4)/(ep.s(0,4,5)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(2,1)*ep.spb(3,2))) + (pow(ep.spa(2,3),3)*pow(ep.spb(5,1),4))/(ep.s(0,1,5)*ep.spa(3,4)*(ep.spab(2,3,5) + ep.spab(2,4,5))*(ep.spab(4,0,1) + ep.spab(4,5,1))*ep.spb(1,0)*ep.spb(5,0)) + (pow(ep.spa(0,2),4)*pow(ep.spb(5,4),3))/(ep.s(3,4,5)*ep.spa(0,1)*ep.spa(1,2)*(ep.spab(0,4,3) + ep.spab(0,5,3))*(ep.spab(2,3,5) + ep.spab(2,4,5))*ep.spb(4,3))); }

template <class T> complex<T> A6g51_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(2,3),3))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g52_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(3,0,2) + ep.spab(3,1,2),4)/(ep.s(3,4,5)*ep.spa(3,4)*ep.spa(4,5)*(ep.spab(3,1,0) + ep.spab(3,2,0))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(1,0)*ep.spb(2,1))) - (pow(ep.spa(0,1),3)*pow(ep.spb(4,2),4))/(ep.s(2,3,4)*ep.spa(0,5)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(5,3,2) + ep.spab(5,4,2))*ep.spb(3,2)*ep.spb(4,3)) - (pow(ep.spa(1,3),4)*pow(ep.spb(5,4),3))/(ep.s(1,2,3)*ep.spa(1,2)*ep.spa(2,3)*(ep.spab(1,2,4) + ep.spab(1,3,4))*(ep.spab(3,1,0) + ep.spab(3,2,0))*ep.spb(5,0))); }

template <class T> complex<T> A6g53_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,3),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g54_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,3),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }


template <class T> complex<T> A6g56_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0,1)*(-(pow(ep.spab(0,1,3) + ep.spab(0,2,3),3)/(ep.s(1,2,3)*ep.spa(0,5)*ep.spa(4,5)*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(2,1)*ep.spb(3,2))) - pow(ep.spab(2,3,5) + ep.spab(2,4,5),3)/(ep.s(0,1,5)*ep.spa(2,3)*ep.spa(3,4)*(ep.spab(4,2,1) + ep.spab(4,3,1))*ep.spb(1,0)*ep.spb(5,0))); }

template <class T> complex<T> A6g57_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(1,2),3))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g58_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,2),4))/(ep.spa(0,1)*ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g60_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return (complex<T>(0,-1)*pow(ep.spa(0,1),3))/(ep.spa(0,5)*ep.spa(1,2)*ep.spa(2,3)*ep.spa(3,4)*ep.spa(4,5)); }

template <class T> complex<T> A6g55_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g59_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }


template <class T> complex<T> A6g61_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g62_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> A6g63_eval(const eval_param<T>& ep, const mass_param_coll& masses){ return complex<T>(0.,0.); }

template <class T> complex<T> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<T>&, const mass_param_coll&) {
	switch (hc) {
	case 0: return &A6g0_eval;
case 1: return &A6g1_eval;
case 2: return &A6g2_eval;
case 3: return &A6g3_eval;
case 4: return &A6g4_eval;
case 5: return &A6g5_eval;
case 6: return &A6g6_eval;
case 7: return &A6g7_eval;
case 8: return &A6g8_eval;
case 9: return &A6g9_eval;
case 10: return &A6g10_eval;
case 11: return &A6g11_eval;
case 12: return &A6g12_eval;
case 13: return &A6g13_eval;
case 14: return &A6g14_eval;
case 15: return &A6g15_eval;
case 16: return &A6g16_eval;
case 17: return &A6g17_eval;
case 18: return &A6g18_eval;
case 19: return &A6g19_eval;
case 20: return &A6g20_eval;
case 21: return &A6g21_eval;
case 22: return &A6g22_eval;
case 23: return &A6g23_eval;
case 24: return &A6g24_eval;
case 25: return &A6g25_eval;
case 26: return &A6g26_eval;
case 27: return &A6g27_eval;
case 28: return &A6g28_eval;
case 29: return &A6g29_eval;
case 30: return &A6g30_eval;
case 31: return &A6g31_eval;
case 32: return &A6g32_eval;
case 33: return &A6g33_eval;
case 34: return &A6g34_eval;
case 35: return &A6g35_eval;
case 36: return &A6g36_eval;
case 37: return &A6g37_eval;
case 38: return &A6g38_eval;
case 39: return &A6g39_eval;
case 40: return &A6g40_eval;
case 41: return &A6g41_eval;
case 42: return &A6g42_eval;
case 43: return &A6g43_eval;
case 44: return &A6g44_eval;
case 45: return &A6g45_eval;
case 46: return &A6g46_eval;
case 47: return &A6g47_eval;
case 48: return &A6g48_eval;
case 49: return &A6g49_eval;
case 50: return &A6g50_eval;
case 51: return &A6g51_eval;
case 52: return &A6g52_eval;
case 53: return &A6g53_eval;
case 54: return &A6g54_eval;
case 56: return &A6g56_eval;
case 57: return &A6g57_eval;
case 58: return &A6g58_eval;
case 60: return &A6g60_eval;
case 55: return &A6g55_eval;
case 59: return &A6g59_eval;
case 61: return &A6g61_eval;
case 62: return &A6g62_eval;
case 63: return &A6g63_eval;
    }
}



template complex<R> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<R>&, const mass_param_coll&) ;
 template complex<RHP> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<RHP>&, const mass_param_coll&) ;
 template complex<RVHP> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<RVHP>&, const mass_param_coll&) ;

#if BH_USE_GMP

 template complex<RGMP> (*A6g_Tree_Ptr_eval(int hc))(const eval_param<RGMP>&, const mass_param_coll&) ;
#endif

}

