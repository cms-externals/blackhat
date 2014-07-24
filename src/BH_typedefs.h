/*!\file BH_typedefs.h
  \brief Header for typedef definitions used in BH
*/
#ifndef BH_TYPEDEFS_H_
#define BH_TYPEDEFS_H_

#include <complex>
#include "qd/qd_real.h"

namespace BH {


//! double precision real
typedef double R;
//! high precision real (2 x double)
typedef dd_real RHP;
//! higher precision real (4 x double)
typedef qd_real RVHP;
//! double precision complex
typedef std::complex<R> C;
//! high precision complex (2 x double)
typedef std::complex<RHP> CHP;
//! higher precision complex (4 x double)
typedef std::complex<RVHP> CVHP;


enum color_structure { nf, 
	LT, 
	RT, 
	leading_color,
	sub_leading_color,
	slc_q,
	slc_G,
	glue,
	nS,
	VECT,
	AX,
	AXSL,
	nf_top,
	LLT,
	RLT,
	LRT,
	RRT,
	nfLT,
	nfRT,
	nfLLT,
	nfRLT,
	nfLRT,
	nfRRT,
	LLLT,
	RRRT,
	RLLT,
	LRLT,
	LLRT,
	RRLT,
	RLRT,
	LRRT,
	nfLLLT,
	nfRRRT,
	nfRLLT,
	nfLRLT,
	nfLLRT,
	nfRRLT,
	nfRLRT,
	nfLRRT,
	unspecified};

enum QCDorder {lo=1,nlo};

std::ostream& operator<<(std::ostream& s, color_structure cs);

}


#endif /*BH_TYPEDEFS_H_*/
