#ifndef RAT_EXT_H_
#define RAT_EXT_H_

#include <complex>
#include <vector>
#include "BH_typedefs.h"
#include "spinor.h"
#ifndef BH_PUBLIC
#include "cut_part_D_Dims.h"
#endif
#include "cut_part_factory.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#if BH_USE_GMP
#include "gmp_r.h"
#endif



namespace BH {

namespace ratext {

/*
 *
 *
 * Functions for setting the distance to zero in a type safe way
 *
 *
 */

template <class T> inline T DeltaZero();

template<> inline R DeltaZero(){return R(1e-13);}
template<> inline RHP DeltaZero(){return RHP(1e-29);}
template<> inline RVHP DeltaZero(){return RVHP(1e-61);}
#if BH_USE_GMP
template<> inline RGMP DeltaZero(){return RGMP(exp10(RGMP(-20*((RGMP::get_current_precision()+63)/64))));}
#endif

/*
 *
 *
 * A function to return the maximum number of digits a particular type has to work with
 *
 *
 */

template <class T> inline T MaxDigits();

template<> inline R MaxDigits(){return R(16);}
template<> inline RHP MaxDigits(){return RHP(32);}
template<> inline RVHP MaxDigits(){return RVHP(64);}

#if BH_USE_GMP
template<> inline RGMP MaxDigits(){return RGMP(RGMP::get_current_nbr_digits());}
#endif


/*
 *
 *
 * Useful structures
 *
 *
 */

template <class T> struct daughter_info {
	typedef T daughter_type;
};

}

}

#endif /*RAT_EXT_H_*/
