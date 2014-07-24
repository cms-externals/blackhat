/*
 * amplitudes_cut_eval.h
 *
 *  Created on: Aug 21, 2008
 *      Author: daniel
 */

#ifndef AMPLITUDES_CUT_EVAL_H_
#define AMPLITUDES_CUT_EVAL_H_


#include <vector>
#include "BH_typedefs.h"


#define _FAST_COMPILE 0

namespace BH {

template <class T> class SeriesC;
template <class T> class eval_param ;
class process;

template <class T> SeriesC<T> (*C_Ptr_eval(const process& pro,color_structure cs))(const eval_param<T>&, const T&) ;
}
#include "C_2q3g_eval.h"
#include "C_2q2g1y_eval.h"


#include "C_2q2l_eval.h"
#include "C_2q1g2l_eval.h"
#include "C_2q2g2l_eval.h"



#include "C_2q2Q_eval.h"

#include "C_2q2Q1g_eval.h"


#include "C_2q2Q2l_eval.h"
#include "C_2q2Q1y_eval.h"

#include "C_2q2Gl2l_eval.h"


#include "C_2q1g1y_eval.h"
#include "C_2q2g1y_eval.h"
#include "C_2q3g2l_eval.h"






#endif /* AMPLITUDES_CUT_EVAL_H_ */
