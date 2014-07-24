/*  util.h  */

/*  David A. Kosower, June 2, 2008  */

/* Utility routines for counting particle types, extracting helicities, etc.;
   extracted from tree1.cc
 */

#define BaseIndex 0 // Must match value in tree1.cc

#define FIXZERO _OLD_PHASE_CONVENTION  // Force -0 to 0 in momenta, needed for old spinor phase convention

#define Re(v) real(v)
#define Im(v) imag(v)

#include <map>
#include <vector>
#include <particles.h>
#include "spinor.h"
#include "mom_conf.h"
using namespace std;
namespace BH {

template<class T> inline T Max(T x,T y) {return (x > y ? x : y);}
inline bool IsEven(int i) {return not(i&1);}
inline bool IsOdd(int i) {return(i&1);}
template<class T>
static inline T square(T x) {return(x*x);}


/* Counts the number of fermion legs of each different flavor in the
   given range of arguments within "id".  The returned vector has length
   equal to the highest flavor number in "id" plus one (the 0 component
   is not used)
*/
  vector<int> FermionCount(const vector<particle_ID>& id, int start, int end);

  // With an added leg included in the count (e.g. an internal line)
  vector<int> FermionCount(const vector<particle_ID>& id, int start, int end,
                           const particle_ID& addedLeg);

/* Determines whether the number of fermion legs of each different
   flavor is even or not.  Return vector as in FermionCount above.
*/
  vector<bool> FermionParity(const vector<particle_ID>& id, int start, int end);

// 6/3/08: Scalar equivalents
vector<int> ScalarCount(const vector<particle_ID>& id, int start, int end);
vector<bool> ScalarParity(const vector<particle_ID>& id, int start, int end);

// Helicity ignored here
  particle_ID FlavoredQuarkID(int flavor);
  particle_ID FlavoredScalarID(int flavor);
  particle_ID NParticleID(int helicity,const particle_ID& base);
vector<particle_ID> NParticleID(const vector<int>& helicity,
                               const vector<particle_ID>& base);

vector<particle_ID> NParticleID(const vector<int>& helicity,
                               const vector<particle_ID>& base,int n);
// returns a characteristic integer
  int ParticleCode(particle_ID id);
  vector<int> ParticleCode(const vector<particle_ID>& id);
  vector<int> Helicities(const vector<particle_ID>& id);

// 7/17/08
namespace Tree {
/* 8/4/08: for tracking masses; because the indices can be quite large
   (and not contiguous), use a map rather than a vector */
map<int,int> MassIndexCount(const vector<int>& massValue, int start, int end);

template <class T> inline void FixZero(momentum<complex<T> >& v) {
 #if FIXZERO
  complex<T> e = v.E(), x = v.X(), y = v.Y(), z = v.Z();
    // Make sure that "-0" doesn't show up as imaginary components, so
    // we always stay on one side of the branch cut; == 0 will pick up
    // both "0" and "-0"
  if (Im(e) == 0) e = complex<T>(Re(e),0);
  if (Im(x) == 0) x = complex<T>(Re(x),0);
  if (Im(y) == 0) y = complex<T>(Re(y),0);
  if (Im(z) == 0) z = complex<T>(Re(z),0);
  v = momentum<complex<T> >(e,x,y,z);
 #endif
}

template<class T> inline void
//  FixZero(momentum<complex<T> >& v) {
  FixZero(Cmom<T>& v) {
 #if FIXZERO
  complex<T> e = v.E(), x = v.X(), y = v.Y(), z = v.Z();
    // Make sure that "-0" doesn't show up as imaginary components, so
    // we always stay on one side of the branch cut; == 0 will pick up
    // both "0" and "-0"
  if (Im(e) == 0) e = complex<T>(Re(e),0);
  if (Im(x) == 0) x = complex<T>(Re(x),0);
  if (Im(y) == 0) y = complex<T>(Re(y),0);
  if (Im(z) == 0) z = complex<T>(Re(z),0);
  v = Cmom<T>(e,x,y,z);
 #endif
}

vector<int> Join(const vector<int>& v1, const vector<int>& v2);
vector<int> Join(const vector<int>& v1, const vector<int>& v2,
                 const vector<int>& v3);
vector<int> Join(const vector<int>& v1, const vector<int>& v2,
                 const vector<int>& v3, const vector<int>& v4);

template<class T> int
  MomentumSum(momentum_configuration<T>& k,
              const vector<int>& v, int start, int end,
              const vector<int>& extraK = empty);
template <class T> momentum<complex<T> > GenerateMomentum(const T& dummy);

template<class T> /* static */ int
  FlatSum(momentum_configuration<T>& k,
          int ref /* index of reference momentum */,
          const vector<int>& v, int start, int end,
          const vector<int>& extraK = empty);
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1);
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1,
                  int ev, const vector<int>& extraK = empty);
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1, int s2, int e2,
                  const vector<int>& extraK1 = empty,
                  const vector<int>& extraK2 = empty);
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1, int s2, int e2,
                  int s3, int e3);

template<class T> int
  Negative(momentum_configuration<T>& k, int i);

}

}
