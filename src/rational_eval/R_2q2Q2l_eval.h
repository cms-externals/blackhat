#ifndef R_2q2Q2l_EVAL_H_
#define R_2q2Q2l_EVAL_H_

#include <vector>
#include <complex>

namespace BH {

class process;
template <class T> class eval_param; class mass_param_coll;

template <class T> std::complex<T> (*R2q2Q2l_L_Ptr_eval( int hc))(const eval_param<T>& ,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q2l_nf_Ptr_eval( int hc))(const eval_param<T>& ,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q2l_nf_top_Ptr_eval( int hc))(const eval_param<T>& ,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q2l_sl_Ptr_eval( int hc))(const eval_param<T>& ,const mass_param_coll&);
template <class T> std::complex<T> (*R2q2Q2l_AX_Ptr_eval( int hc))(const eval_param<T>& ,const mass_param_coll&);

}
#endif /* R_2q2Q2l_EVAL_H_ */
