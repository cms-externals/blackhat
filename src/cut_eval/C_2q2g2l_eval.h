#ifndef C_2q2g2l_EVAL_H_
#define C_2q2g2l_EVAL_H_


#include <vector>

namespace BH {

template <class T> class SeriesC;
template <class T> class eval_param ;
class process;

template <class T> SeriesC<T> (*C2q2g2l_L_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_SLC_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_nf_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_nf_top_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_VECT_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_AX_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> (*C2q2g2l_AXSL_Ptr_eval(int hc))(const eval_param<T>&, const T&);

}




#endif /* C_2q2g2l_EVAL_H_ */
