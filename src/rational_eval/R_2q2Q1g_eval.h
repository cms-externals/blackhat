#ifndef R_2q2Q1g_EVAL_H_
#define R_2q2Q1g_EVAL_H_

#include <vector>
#include <complex>

namespace BH{

class process;
template <class T> class eval_param; class mass_param_coll;


template <class T> std::complex<T> (*R2q2Q1g_LLT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q1g_LRT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q1g_nfLLT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);


}
#endif
