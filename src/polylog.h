/* polylog.h */

/*  David A. Kosower, November 19, 2007  */

/* Dilogarithm and Clausen function */

#ifndef PolyLogDefined
#define PolyLogDefined 1

#if BH_USE_GMP
#include "gmp_r.h"
#endif

namespace BH {

template <class T> inline T pi();

template<> inline R pi(){return M_PI;}
template<> inline RHP pi(){return dd_real::_pi;}
template<> inline RVHP pi(){return qd_real::_pi;}
template<> inline C pi(){return C(M_PI,0.);}
template<> inline CHP pi(){return CHP(dd_real::_pi,0.);}
template<> inline CVHP pi(){return CVHP(qd_real::_pi,0.);}

#if BH_USE_GMP
template<> inline RGMP pi(){return get_pi();}
template<> inline std::complex<RGMP> pi(){return std::complex<RGMP>(get_pi(),RGMP(0));}
#endif

R ReLi2(R);
C Li2(C);
R Cl2(R);

RHP ReLi2(RHP);
RHP Cl2(RHP);
CHP Li2(CHP);

RVHP ReLi2(RVHP);
RVHP Cl2(RVHP);
CVHP Li2(CVHP);

#if BH_USE_GMP
RGMP ReLi2(const RGMP&);
RGMP Cl2(const RGMP&);
std::complex<RGMP> Li2(const std::complex<RGMP>&);
#endif



}
#endif /* PolyLogDefined */
