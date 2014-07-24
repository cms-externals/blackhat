#ifndef C_2q3g_EVAL_H_
#define C_2q3g_EVAL_H_


#include <vector>

namespace BH {

template <class T> class SeriesC;
template <class T> class eval_param ;
class process;

template <class T> SeriesC<T> ( *C2q3g_L_Ptr_eval(int hc))(const eval_param<T>&, const T&);
template <class T> SeriesC<T> ( *C2q3g_nf_Ptr_eval(int hc))(const eval_param<T>&, const T&);

}




#endif /* C_2q3g_EVAL_H_ */
