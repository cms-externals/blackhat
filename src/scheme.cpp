/*
 * scheme.cpp
 *
 *  Created on: Aug 12, 2008
 *      Author: daniel
 */

#include "particles.h"
#include "BH_typedefs.h"
#include "scheme.h"
#include "process.h"

namespace BH {

template <class T> T scheme_shift(const process& pro,scheme sc){
	switch (pro.pcode()){
	case 220: switch (sc) {
		case FDH: return T(0);
		case HV : return T(-1)/T(2);
	}
	case 221: switch (sc) {
		case FDH: return T(0);
		case HV : return T(-1)/T(2);
	}
	default: return T(0);
	}
}

template R scheme_shift<R>(const process& pro,scheme sc);
template RHP scheme_shift<RHP>(const process& pro,scheme sc);
template RVHP scheme_shift<RVHP>(const process& pro,scheme sc);

}
