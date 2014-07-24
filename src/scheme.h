/*
 * scheme.h
 *
 *  Created on: Aug 13, 2008
 *      Author: daniel
 */

#ifndef SCHEME_H_
#define SCHEME_H_

#include "particles.h"

namespace BH {

enum scheme { FDH, CRD, HV };

template <class T> T scheme_shift(const process& pro,scheme sc);

}

#endif  /*SCHEME_H_*/
