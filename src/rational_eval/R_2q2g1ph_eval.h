#ifndef R_2q2g1ph_EVAL_H_
#define R_2q2g1ph_EVAL_H_

#include <vector>
#include <complex>

namespace BH{

class process;
template <class T> class eval_param; class mass_param_coll;


template <class T> std::complex<T> (*R2q2g1ph_LT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2g1ph_nfLT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2g1ph_RT_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);


}
#endif
