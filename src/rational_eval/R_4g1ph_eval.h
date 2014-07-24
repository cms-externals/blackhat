#ifndef R_4g1ph_EVAL_H_
#define R_4g1ph_EVAL_H_

#include <vector>
#include <complex>

namespace BH{

class process;
template <class T> class eval_param; class mass_param_coll;


template <class T> std::complex<T> (*R4g1ph_G_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);
template <class T> std::complex<T> (*R4g1ph_nf_Ptr_eval(int hc))(const eval_param<T>&,const mass_param_coll&);


}
#endif
