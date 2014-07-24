#include "qd_suppl.h"
#include "BH_typedefs.h"
namespace BH {

std::complex<BH::RHP> pow(std::complex<BH::RHP> __x,int __n)    {
	if ( __n<0 ) {return pow(CHP(1.,0.)/__x,-__n);}
	std::complex<BH::RHP> __y;
	if (__n % 2 ) {__y=__x;} else {__y=CHP(1.,0.);}

    while (__n >>= 1)
      {
        __x = __x * __x;
        if (__n % 2)
          __y = __y * __x;
      }

    return __y;
  }
std::complex<BH::RVHP> pow(std::complex<BH::RVHP> __x,int __n)    {
	if ( __n<0 ) {return pow(CVHP(1.,0.)/__x,-__n);}
	std::complex<BH::RVHP> __y;
	if (__n % 2 ) {__y=__x;} else {__y=CVHP(1.,0.);}

    while (__n >>= 1)
      {
        __x = __x * __x;
        if (__n % 2)
          __y = __y * __x;
      }

    return __y;
  }
}

// specialization for complex<qd_real>

namespace std {
#if _USE_GCC
template<> BH::CHP sqrt<BH::RHP>(const BH::CHP& __z)
{
  BH::RHP __x = __z.real();
  BH::RHP __y = __z.imag();

  if (__x == BH::RHP(0))
    {
	  BH::RHP __t = sqrt(abs(__y) / BH::RHP(2));
      return BH::CHP(__t, __y < BH::RHP(0) ? -__t : __t);
    }
  else
    {
	  BH::RHP __t = sqrt(BH::RHP(2) * (std::abs(__z) + abs(__x)));
	  BH::RHP __u = __t / BH::RHP(2);
      return __x > BH::RHP(0)
        ? BH::CHP(__u, __y / __t)
        : BH::CHP(abs(__y) / __t, __y < BH::RHP(0) ? -__u : __u);
    }
}

template<> BH::CVHP sqrt<BH::RVHP>(const BH::CVHP& __z)
{
  BH::RVHP __x = __z.real();
  BH::RVHP __y = __z.imag();

  if (__x == BH::RVHP(0.))
    {
	  BH::RVHP __t = sqrt(abs(__y) / BH::RVHP(2));
      return BH::CVHP(__t, __y < BH::RVHP(0) ? -__t : __t);
    }
  else
    {
	  BH::RVHP __t = sqrt(BH::RVHP(2) * (std::abs(__z) + abs(__x)));
	  BH::RVHP __u = __t / BH::RVHP(2);
      return __x > BH::RVHP(0)
        ? BH::CVHP(__u, __y / __t)
        : BH::CVHP(abs(__y) / __t, __y < BH::RVHP(0) ? -__u : __u);
    }
}
#endif
}
#if _USE_SUN_CC || _USE_PGCC
BH::CHP sqrt(const BH::CHP& __z)
{
  BH::RHP __x = __z.real();
  BH::RHP __y = __z.imag();

  if (__x == BH::RHP(0.))
    {
	  BH::RHP __t = sqrt(abs(__y) / BH::RHP(2.));
      return BH::CHP(__t, __y < BH::RHP(0.) ? -__t : __t);
    }
  else
    {
	  BH::RHP __t = sqrt(2 * ( std::abs(__z) + abs(__x)));
	  BH::RHP __u = __t / BH::RHP(2.);
      return __x > BH::RHP(0.)
        ? BH::CHP(__u, __y / __t)
        : BH::CHP(abs(__y) / __t, __y < BH::RHP(0.) ? -__u : __u);
    }
}

BH::CVHP sqrt(const BH::CVHP& __z)
{
  BH::RVHP __x = __z.real();
  BH::RVHP __y = __z.imag();

  if (__x == BH::RVHP(0.))
    {
	  BH::RVHP __t = sqrt(abs(__y) / BH::RVHP(2.));
      return BH::CVHP(__t, __y < BH::RVHP(0.) ? -__t : __t);
    }
  else
    {
	  BH::RVHP __t = sqrt(2 * ( std::abs(__z) + abs(__x)));
	  BH::RVHP __u = __t / BH::RVHP(2.);
      return __x > BH::RVHP(0.)
        ? BH::CVHP(__u, __y / __t)
        : BH::CVHP(abs(__y) / __t, __y < BH::RVHP(0.) ? -__u : __u);
    }
}

BH::RHP sqrt(const BH::RHP& __z);
BH::RVHP sqrt(const BH::RVHP& __z);

#endif

