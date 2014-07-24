/* polylog.cc */

/*  David A. Kosower, November 19, 2007  */

/*  Implementation of real and complex polylogarithms.  The approach is
    based on the one described in 't Hooft and Veltman [NPB153:365 (1979)], 
    and also documented in hep-ph/0502165v2 [Hameren, Vollinga and Weinzierl].

*/

//#include "c++standard.h"


#include "BH_typedefs.h"
#define Re(v) real(v)

#include "polylog.h"
#include <math.h>

#include "polylog.hpp"


namespace BH {

// Gives error < 10^-15 over key interval [0,1/2]
#define Li2TermLimit 7
// Gives |error| < 10^-15 in key region for complex Li_2
#define CLi2TermLimit 8
// Gives error < 10^-15 over key interval [0,2 Pi/3],
// overall < 10^-14 for Cl2 in [-12,12]
#define Cl2TermLimit 14

const double PiSquaredOver3 = sq(M_PI)/3.;
const double PiSquaredOver6 = sq(M_PI)/6.;

#define C std::complex<double>

// B_{2 index}
const double Bernoulli2[] = {1,
1./6., -1./30., 1./42., -1./30., 5./66., -691./2730., 7./6., -3617./510., 
43867./798.,-174611./330., 854513./138., -236364091./2730., 8553103./6., 
-23749461029./870., 8615841276005./14322., -7709321041217./510., 
2577687858367./6., -26315271553053477373./1919190., 2929993913841559./6.,
-261082718496449122051./13530.};
const double B0 = 1;
const double B1over2 = -1./4.;

double ReLi2(double x)
{double added = 0, factor = 1;
 // Map the argument into the range [0,1/2]
 if (x >= 2)
   // Re[PolyLog[2, x]] ->  Pi^2/3 - Log[x]^2/2 - PolyLog[2, x^(-1)]}
    {added = PiSquaredOver3 - 0.5*sq(log(x));
     factor = -1;
     x = 1/x;}
 else if (x > 1)
   // Re[PolyLog[2, x]] -> 
   //           Pi^2/6 - Log[x-1]*Log[x] + Log[x]^2/2 + PolyLog[2, (x-1)/x]
    {double lnx = log(x);
     added = PiSquaredOver6 +(0.5*lnx - log(x-1))*lnx;
     factor = 1;
     x = (x-1)/x;}
 else if (x > 0.5)
   // PolyLog[2,x] -> -PolyLog[2,1-x] + Pi^2/6 - Log[x] Log[1-x]
    {added = PiSquaredOver6 - log(x)*log(1-x);
     factor = -1;
     x = 1-x;}
 else if (x > 0) {}
 else if (x >= -1)
   // PolyLog[2, y] ->  -Log[1 - y]^2/2 - PolyLog[2, y/(-1 + y)]
    {added = -0.5*sq(log(1-x));
     factor = -1;
     x = x/(x-1);}
 else 
   // PolyLog[2, y] ->  -Pi^2/6 + Log[1-y]^2/2 - Log[1-y]*Log[-y] 
   //                   + PolyLog[2, 1/(1-y)]
    {double ln1x = log(1-x);
     added = -PiSquaredOver6 + (0.5*ln1x-log(-x))*ln1x;
     factor = 1;
     x = 1/(1-x);}
 // Compute PolyLog[2,x], x now in [0,1/2], by Bernoulli series
 double z = -log(1-x);

 double li2 = (B0+B1over2*z)*z;
 double term = z, zsq = z*z;
 int limit = Li2TermLimit;
 if (x < 0.2) limit = 4;
 for (int j = 1;  j <= limit;  j += 1)
    {term *= zsq/(double)((2*j+1)*2*j);
    li2 += Bernoulli2[j]*term;}

 return(factor*li2+added);
}

// Complex version of the above
C Li2(C z)
  /* Transformations may not be for optimal regions yet; following
     't Hooft & Veltman, put the argument inside the circle |z|=1,
     and make the real part < 0.5 */
{C added = 0;
 double factor = 1;
 if (Re(z*conj(z)) > 1)
    {added = -PiSquaredOver6-0.5*sq(log(-z));
     factor = -1;
     z = 1./z;}
 if (Re(z) > 0.5)
    {added += factor*(PiSquaredOver6-log(z)*log(1.-z));
     factor *= -1;
     z = 1.-z;}

 // Compute PolyLog[2,z], z now in unit disc with Re z < 1/2, by Bernoulli series
 C w = -log(1.-z);

 C li2 = (B0+B1over2*w)*w;
 C term = w, wsq = w*w;
 int limit = CLi2TermLimit;
 if (Re(w*conj(w)) < 0.05) limit = 4;
 for (int j = 1;  j <= limit;  j += 1)
    {term *= wsq/(double)((2*j+1)*2*j);
    li2 += Bernoulli2[j]*term;}

 return(factor*li2+added);
}

// Clausen function, from hep-ph/0502165v2 [Hameren, Vollinga and Weinzierl]
double Cl2(double x)
{double factor = 1;
 if (x < 0) {x = -x;  factor = -1;}
 const double TwoPi = 2*M_PI;
 const double TwoPiOver3 = 2*M_PI/3;
 const double PiOver3 = M_PI/3;
 while (x > TwoPi) x -= TwoPi;
 // shift to range [0,2 Pi/3] (larger upper limit prevents too many recursions)
 if (x > TwoPiOver3) return(2*factor*(Cl2(0.5*x)-Cl2(M_PI-0.5*x)));
 // Compute by Bernoulli series
 double cl2 = x*(1-log(x));
 double term = -x, mxsq = -x*x;
 int limit = Cl2TermLimit;
 if (x < PiOver3) limit = 8;
 for (int j = 1;  j <= limit;  j += 1)
    {term *= mxsq/(double)((2*j+1)*2*j);
     cl2 += Bernoulli2[j]*term/(double)(2*j);}
 
 return(factor*cl2);
}
}

#include "polylog_HP.hpp"
#include "polylog_VHP.hpp"

#if BH_USE_GMP
namespace BH {
RGMP ReLi2(const RGMP& x){
	return ReLi2_tpl(x);
};
std::complex<RGMP> Li2(const std::complex<RGMP>& x){
	return Li2_tpl(x);
};
}
#endif
