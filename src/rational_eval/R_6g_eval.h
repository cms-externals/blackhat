#ifndef R_6g_EVAL_H_
#define R_6g_EVAL_H_

#include <vector>
#include <complex>

namespace BH {

class process;
template <class T> class eval_param; class mass_param_coll;

 template <class T> std::complex<T> (*R6g_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&) ;
 template <class T> std::complex<T> (*R6g_Ptr_nf_eval(int hc))(const eval_param<T>&,const mass_param_coll&) ;
 

}
#endif /* R_6g_EVAL_H_ */
