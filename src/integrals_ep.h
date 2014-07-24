/* integrals_ep.h */

/*  Based on David A. Kosower, November 15, 2007  */
/*  Darren Forde, May 22, 2009  */

/*  Definitions of one-loop integrals.  The integrals have a factor
    of c_\Gamma removed, so that those with 1/e^2 divergences start
    with precisely n/e^2 as their leading term.

 */

#ifndef IntegralsDefined_ep
#define IntegralsDefined_ep 1


#include <vector>
#include <complex>
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "mom_conf.h"
#include "eval_param.h"
#include "Series.h"

namespace BH {

class boxD;
class triangleD;
// Define a vector whose square is the mu^2 scale
template <class T> int DefineMu(const eval_param<T>& k, T mu);

// Box wrappers
template <class T > std::complex<T> IntM(int order, // which order in ksilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3,
      const std::vector<int>& corner4, const std::vector<int>& masses);

template <class T > Series<std::complex<T> > IntM(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3,
            const std::vector<int>& corner4, const std::vector<int>& masses);

template <class T > std::complex<T> Int(int order, // which order in ksilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3,
      const std::vector<int>& corner4);

template <class T > Series<std::complex<T> > Int(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3,
            const std::vector<int>& corner4);

// Triangle wrappers
template <class T > std::complex<T> IntM(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3, const std::vector<int>& masses);

template <class T > Series<std::complex<T> > IntM(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3, const std::vector<int>& masses);

template <class T > std::complex<T> Int(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2,
      const std::vector<int>& corner3);

template <class T > Series<std::complex<T> > Int(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2,
            const std::vector<int>& corner3);

// Bubble wrappers -- forms both with one mass argument (all that's really
// needed) and two for compatibility with above
template <class T > std::complex<T> IntM(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2, const std::vector<int>& masses);

template <class T > Series<std::complex<T> > IntM(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2, const std::vector<int>& masses);

template <class T > std::complex<T> Int(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1,
      const std::vector<int>& corner2);

template <class T > Series<std::complex<T> > Int(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1,
            const std::vector<int>& corner2);


template <class T > Series<std::complex<T> > IntM(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1, const std::vector<int>& masses);

template <class T > Series<std::complex<T> > IntM(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1, const std::vector<int>& masses);

template <class T > Series<std::complex<T> > Int(int order, // which order in epsilon?
      const eval_param<T>& k,
      const T& mu, // index for vector whose square is mu^2
      const std::vector<int>& corner1);

template <class T > Series<std::complex<T> > Int(const eval_param<T>& k,
            const T& mu, // index for vector whose square is mu^2
            const std::vector<int>& corner1);

}
#endif /* IntegralsDefined_ep */
