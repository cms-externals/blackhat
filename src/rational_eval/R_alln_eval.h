// last update 9/17/2008 - CFB

#ifndef R_ALLN_EVAL_H_
#define R_ALLN_EVAL_H_

#include <complex>
#include "qd/qd_real.h"
#include "spinor.h"
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "particles.h"
#include "eval_param.h"
#include "qd_suppl.h"

namespace BH 
{

// all-plus amplitude
template <class T> std::complex<T> Rallp(const eval_param<T>& ep,const mass_param_coll& mpc);

// all-minus amplitude
template <class T> std::complex<T> Rallm(const eval_param<T>& ep,const mass_param_coll& mpc);

// mpppppp -- note 1st mpcex = minus
template <class T> std::complex<T> Rallmp(const eval_param<T>& ep,const mass_param_coll& mpc);

// pmmmmmm -- note 1st mpcex = plus
template <class T> std::complex<T> Rallpm(const eval_param<T>& ep,const mass_param_coll& mpc);

// m_fpppj_fppp L -- note 1st mpcex = minus
template <class T> std::complex<T> Rallmp_L(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!

// m_fpppj_fppp L -- note 1st mpcex = minus
template <class T> std::complex<T> Rallmp_Ls(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!


// m_fpppj_fppp s -- note 1st mpcex = minus
template <class T> std::complex<T> Rallmp_s(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!

// p_fmmmj_fmmm L -- note 1st mpcex = plus
template <class T> std::complex<T> Rallpm_L(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!

// p_fmmmj_fmmm L -- note 1st mpcex = plus
template <class T> std::complex<T> Rallpm_Ls(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!

// p_fmmmj_fmmm s -- note 1st mpcex = plus
template <class T> std::complex<T> Rallpm_s(const eval_param<T>& ep,const mass_param_coll& mpc,
		const size_t j);
// NOTE: j labels the fermion mpcex, i.e. the momentum label is the jth entry in mpc!!!!


// auxiliary stuff ---------------------

// angle bracket denominator
template <class T> std::complex<T> denomang(const eval_param<T>& ep,const mass_param_coll& mpc);

// square bracket denominator
template <class T> std::complex<T> denomsqu(const eval_param<T>& ep,const mass_param_coll& mpc);

}

#endif /*R_ALLN_EVAL_H_*/
