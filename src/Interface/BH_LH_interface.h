/*
 * LH_interface.h
 *
 *  Created on: Jun 16, 2009
 *      Author: daniel
 */

#ifndef LH_INTERFACE_H_
#define LH_INTERFACE_H_

#define LINK_WITH_FORTRAN 1

#include <vector>


extern "C" {
extern void OLP_EvalSubProcess(int* Label,double* Momenta,double *mu,double *parameters,double *result);
extern void OLP_Start(const char* filename,int *status);
}


namespace BH {

namespace LesHouches {


void EvalSubprocess(int Label,double* Momenta,double mu,double alpha_s,double alpha_ew,double* result);
int Init(const char* filename);

}


}
#endif /* LH_INTERFACE_H_ */
