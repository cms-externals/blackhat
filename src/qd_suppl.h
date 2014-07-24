#ifndef QD_SUPPL_H_
#define QD_SUPPL_H_

#define _USE_SUN_CC __SUNPRO_CC
#define _USE_GCC __GNUC__
#define _USE_PGCC __PGI

#include <complex>
#include "qd/qd_real.h"

namespace BH {



class CPU_FIX {
	unsigned int old_cw;
public:
	CPU_FIX(){ fpu_fix_start(&old_cw);};
	~CPU_FIX(){fpu_fix_end(&old_cw);};
};


std::complex<dd_real> pow(std::complex<dd_real>  __x,int __n) ;
std::complex<qd_real> pow(std::complex<qd_real>  __x,int __n) ;
//CHP my_sqrt(const CHP __z);
//C my_sqrt(const C __z);
}

// specialization for complex<qd_real>

#if _USE_GCC
namespace std {
template <> std::complex<dd_real> sqrt(const std::complex<dd_real>& __z);
template <> std::complex<qd_real> sqrt(const std::complex<qd_real>& __z);
}
#endif

#if _USE_SUN_CC || _USE_PGCC


std::complex<dd_real> sqrt(const std::complex<dd_real>& __z)
//{
//	 dd_real __x = __z.real();
//	 dd_real __y = __z.imag();
//
//	  if (__x ==dd_real(0.))
//	    {
//		 dd_real __t = sqrt(abs(__y) /dd_real(2.));
//	      return std::complex<dd_real>(__t, __y <dd_real(0.) ? -__t : __t);
//	    }
//	  else
//	    {
//		 dd_real __t = sqrt(2 * (std::abs(__z) + abs(__x)));
//		 dd_real __u = __t /dd_real(2.);
//	      return __x >dd_real(0.)
//	        ? std::complex<dd_real>(__u, __y / __t)
//	        : std::complex<dd_real>(abs(__y) / __t, __y <dd_real(0.) ? -__u : __u);
//	    }
//	}
;
std::complex<qd_real> sqrt(const std::complex<qd_real>& __z)
//{
//	 qd_real __x = __z.real();
//	 qd_real __y = __z.imag();
//
//	  if (__x ==qd_real(0.))
//	    {
//		 qd_real __t = sqrt(abs(__y) /qd_real(2.));
//	      return std::complex<qd_real>(__t, __y <qd_real(0.) ? -__t : __t);
//	    }
//	  else
//	    {
//		 qd_real __t = sqrt(2 * (std::abs(__z) + abs(__x)));
//		 qd_real __u = __t /qd_real(2.);
//	      return __x >qd_real(0.)
//	        ? std::complex<qd_real>(__u, __y / __t)
//	        : std::complex<qd_real>(abs(__y) / __t, __y <qd_real(0.) ? -__u : __u);
//	    }
//	}
;

namespace std {

template <> inline dd_real abs(const std::complex<dd_real>& __z)  {
    dd_real __x = __z.real();
    dd_real __y = __z.imag();
    const dd_real __s = std::max(abs(__x), abs(__y));
    if (__s == dd_real())  // well ...
      return __s;
    __x /= __s;
    __y /= __s;
    return __s * ::sqrt(__x * __x + __y * __y);
  }
;
template <> inline qd_real abs(const std::complex<qd_real>& __z)  {
    qd_real __x = __z.real();
    qd_real __y = __z.imag();
    const qd_real __s = std::max(abs(__x), abs(__y));
    if (__s == qd_real())  // well ...
      return __s;
    __x /= __s;
    __y /= __s;
    return __s * ::sqrt(__x * __x + __y * __y);
  }
;
}

//inline dd_real BH_abs(const std::complex<dd_real>& __z)
//  {
//    dd_real __x = __z.real();
//    dd_real __y = __z.imag();
//    const dd_real __s = std::max(abs(__x), abs(__y));
//    if (__s == dd_real())  // well ...
//      return __s;
//    __x /= __s;
//    __y /= __s;
//    return __s * sqrt(__x * __x + __y * __y);
//  }
//
//inline qd_real BH_abs(const std::complex<qd_real>& __z)
//  {
//    qd_real __x = __z.real();
//    qd_real __y = __z.imag();
//    const qd_real __s = std::max(abs(__x), abs(__y));
//    if (__s == qd_real())  // well ...
//      return __s;
//    __x /= __s;
//    __y /= __s;
//    return __s * sqrt(__x * __x + __y * __y);
//  }

inline dd_real arg(const std::complex<dd_real>& __z) {
return  atan2(__z.imag(), __z.real());
}
inline qd_real arg(const std::complex<qd_real>& __z) {
return  atan2(__z.imag(), __z.real());
}

inline std::complex<dd_real>
  polar(const dd_real& __rho, const dd_real& __theta)
  { return std::complex<dd_real>(__rho * cos(__theta), __rho * sin(__theta)); }

inline std::complex<qd_real>
    polar(const qd_real& __rho, const qd_real& __theta)
    { return std::complex<qd_real>(__rho * cos(__theta), __rho * sin(__theta)); }



inline std::complex<dd_real>
  exp(const std::complex<dd_real>& __z)
  { return polar(exp(__z.real()), __z.imag()); }

inline std::complex<qd_real> exp(const std::complex<qd_real>& __z)
    { return polar(exp(__z.real()), __z.imag()); }



inline std::complex<dd_real>  log(const std::complex<dd_real>& __z)
  { return std::complex<dd_real>(log(std::abs(__z)), arg(__z)); }
inline std::complex<qd_real>  log(const std::complex<qd_real>& __z)
  { return std::complex<qd_real>(log(std::abs(__z)), arg(__z)); }


#endif


#if _USE_PGCC
// 26.2.8/5 log(__z): Reurns the natural complex logaritm of __z.
//                    The branch cut is along the negative axis.
template<typename _Tp>
  inline complex<_Tp>
  log(const complex<_Tp>& __z)
  { return complex<_Tp>(log(std::abs(__z)), std::arg(__z)); }

#endif


namespace BH {
inline double to_double(const double& r){return r;}
inline std::complex<double> to_double(const std::complex<double>& c){return c;}
inline std::complex<double> to_double(const std::complex<dd_real>& c){return std::complex<double>(to_double(c.real()),to_double(c.imag()));}
inline std::complex<double> to_double(const std::complex<qd_real>& c){return std::complex<double>(to_double(c.real()),to_double(c.imag()));}


inline std::complex<dd_real> to_HP(const std::complex<qd_real>& c){return std::complex<dd_real>(to_dd_real(c.real()),to_dd_real(c.imag()));}
inline std::complex<dd_real> to_HP(const std::complex<dd_real>& c){return c;}
inline std::complex<dd_real> to_HP(const std::complex<double>& c){return std::complex<dd_real>(c.real(),c.imag());}

inline std::complex<qd_real> to_VHP(const std::complex<qd_real>& c){return c;}
inline std::complex<qd_real> to_VHP(const std::complex<dd_real>& c){return std::complex<qd_real>(c.real(),c.imag());}
inline std::complex<qd_real> to_VHP(const std::complex<double>& c){return std::complex<qd_real>(c.real(),c.imag());}


inline dd_real to_HP(const qd_real& c){return to_dd_real(c);}
}
#endif /*QD_SUPPL_H_*/
