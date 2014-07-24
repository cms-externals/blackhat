/*  zero.cc  */

/*  David A. Kosower, June 3, 2008  */

/*  Determine whether a number is close enough to zero to be set to zero.
 */

#ifndef ZeroDefined
using namespace BH;

// Is it smaller in absolute value than zeroThreshold?
// Be careful that zT^2 is representable!
const double zeroThreshold = 4e-14;
const RHP zeroThresholdHP = RHP("1e-30");  // What should this value be???

inline bool IsZero(const complex<R>& value)
{return(real(value)*real(value)+imag(value)*imag(value) < 
        zeroThreshold*zeroThreshold);}

inline bool IsZero(const complex<RHP>& value)
{return(real(value)*real(value)+imag(value)*imag(value) < 
        zeroThresholdHP*zeroThresholdHP);}

inline bool IsZero(const R& value)
{return(value < zeroThreshold);}

inline bool IsZero(const RHP& value)
{return(value < zeroThresholdHP);}

#define ZeroDefined
#endif // ZeroDefined
