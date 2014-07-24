/* integrals.h */

/*  David A. Kosower, November 15, 2007  */

/*  Definitions of one-loop integrals.  The integrals have a factor
    of c_\Gamma removed, so that those with 1/e^2 divergences start
    with precisely n/e^2 as their leading term.

 */

#ifndef IntegralsDefined
#define IntegralsDefined 1


#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "mom_conf.h"
#include "Series.h"

#if BH_USE_GMP
#include "gmp_r.h"
#endif


#define MomentumConfiguration mom_conf
#define MomentumSet mom_conf

namespace BH {

class part;

// Define a vector whose square is the mu^2 scale
template <class T> int DefineMu(momentum_configuration<T>& k, T mu);

// Box wrappers
template <class T > std::complex<T> Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3,
      const std::vector<int>& corner4);

template <class T > std::complex<T> Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3,
      const std::vector<int>& corner4,
      const part& pa);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3,
            const std::vector<int>& corner4);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3,
            const std::vector<int>& corner4,
            const part&);

// Triangle wrappers
template <class T > std::complex<T> Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3);

template <class T > std::complex<T> Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3,
      const part& pa);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3,
            const part& pa);

// Bubble wrappers -- forms both with one mass argument (all that's really
// needed) and two for compatibility with above
template <class T > std::complex<T> Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2);


template <class T > Series<std::complex<T> > Int(int order, // which order in epsilon?
      momentum_configuration<T>& k,
      int mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1);

template <class T > Series<std::complex<T> > Int(momentum_configuration<T>& k,
            int mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1);

#if 0
//moved to integrals.cpp to allow easier testing
template<class T> inline T epsilon();
template <> inline R epsilon<R>(){return(1.0e-13);};
template <> inline RHP epsilon<RHP>(){return(RHP("1.0e-52"));};
template <> inline RVHP epsilon<RVHP>(){return(RVHP("1.0e-82"));}

#if BH_USE_GMP
template <> inline RGMP epsilon<RGMP>(){ return exp10(-RGMP(RGMP::get_current_nbr_digits()-3));}
#endif

#endif


}
#endif /* IntegralsDefined */
