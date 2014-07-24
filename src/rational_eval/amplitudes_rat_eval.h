#ifndef AMPLITUDES_RAT_EVAL_H_
#define AMPLITUDES_RAT_EVAL_H_

#include <vector>
#include <complex>
#include "BH_typedefs.h"

#define _RECURSION_TESTING_MODE 0 // if set to 1, only the rational vertices that can not be constructed recursively are set known


namespace BH {

template <class T> class eval_param; class mass_param_coll;
template <class T> class SeriesC;
class process;

template <class T> std::complex<T> (*R_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&,const mass_param_coll&) ;

}

#include "R_2g2q_eval.h"
#include "R_2q1g1y_eval.h"
#include "R_2q1g2l_eval.h"
#include "R_2q2g1y_eval.h"
#include "R_2q2g2l_eval.h"
#include "R_2q2l_eval.h"
#include "R_2q2Q1g_eval.h"
#include "R_2q2Q1y_eval.h"
#include "R_2q2Q2l_eval.h"
#include "R_2q2Q_eval.h"
#include "R_2q3g2l_eval.h"
#include "R_2q3g_eval.h"
#include "R_2q3g_r_eval.h"
#include "R_2q4g_eval.h"
#include "R_4g_eval.h"
#include "R_5g_eval.h"
#include "R_6g_eval.h"
#include "R_7g_eval.h"
#include "R_8g_eval.h"
#include "R_9g_eval.h"
#include "R_4q_eval.h"
#include "R_2q2G1g2l_eval.h"
#include "R_2q2g_eval.h"

#include "R_3g1ph_eval.h"
#include "R_4g1ph_eval.h"
#include "R_2q2g1ph_eval.h"
#include "R_2q2Q1ph_eval.h"

#endif /*AMPLITUDES_RAT_EVAL_H_*/
