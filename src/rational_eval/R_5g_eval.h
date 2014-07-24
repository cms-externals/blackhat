#ifndef R_5g_EVAL_H_
#define R_5g_EVAL_H_

#include <vector>
#include <complex>

namespace BH {

class process;
template <class T> class eval_param; class mass_param_coll;

template <class T> std::complex<T> (*R5g_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&) ;

template <class T> std::complex<T> (*R5g_nf_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&) ;



}
#endif /* R_5g_EVAL_H_ */
