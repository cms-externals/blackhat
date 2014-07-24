#ifndef C_2q2l_EVAL_H_
#define C_2q2l_EVAL_H_


#include <vector>

namespace BH {

template <class T> class SeriesC;
template <class T> class eval_param ;
class process;

template <class T> SeriesC<T> (*C2q2l_L_Ptr_eval(int hc))(const eval_param<T>& ,const  T& mu);

}
#endif /* C_2q2l_EVAL_H_ */
