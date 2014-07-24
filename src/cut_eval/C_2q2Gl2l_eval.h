#ifndef C_2q2G2l_EVAL_H_
#define C_2q2G2l_EVAL_H_


#include <vector>

namespace BH {

template <class T> class SeriesC;
template <class T> class eval_param ;
class process;

template <class T> SeriesC<T> ( *C2q2G2l_L_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> ( *C2q2G2l_nf_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> ( *C2q2G2l_nf_top_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> ( *C2q2G2l_sl_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> ( *C2q2G2l_AX_Ptr_eval(int hc))(const eval_param<T>&, const T&);

}




#endif /* C_2q2G2l_EVAL_H_ */
