/* Tree1.cc */

/*  David A. Kosower,  September 14, 2007  */

/* Implementation of light-cone gauge Berends-Giele recursion relations
   using BlackHat momentum library.

   3/22/08: Updated for new version of BlackHat

   5/22/08: Reference momentum is no longer an argument to J; it is cached,
            and created if no "ref" label exists in the momentum
            configuration.

*/

#define BRENDSGIELE_H_
#define BERENDSGIELE_IMPL_H_
//#include "BlackHat.h"
#include <vector>
#include <map>
#include <ctime>
#include "spinor.h"
#include "mom_conf.h"
#include "Tree.h"
#include "Tree_impl.h"
#include "util.h"
#define Re(v) real(v)
#define Im(v) imag(v)

#define RefTag "ref"

//#include "polylog.h"
#include <math.h>

#define BaseIndex 0 // [0] index not used in arg, helicity, or id.

#define isnt !=
#define is ==

#define FIXZERO _OLD_PHASE_CONVENTION  // Force -0 to 0 in momenta, needed for old spinor phase convention

using namespace std;

namespace BH {
#undef MomentumConfiguration
#define MomentumConfiguration momentum_configuration<R>




#if 0
template<class T> inline T Max(T x,T y) {return (x > y ? x : y);}
inline bool IsOdd(int i) {return(i&1);}
#endif

// Returns the number of elements from start to end, understood in
// a cyclic sense.
inline int CountCyclic(const vector<int>& v,int start,int end)
{if (start <= end) return (end-start+1);
 else return(v.size()-start+end+1);}

// Returns the given index mod n
inline int IndexCyclic(const vector<int>& v,int index)
{while (index >= v.size()) index -= v.size();
 while (index < 0) index += v.size();
 return index;}

// Sums elements of the vector from the given start to the given end
// position, these boundaries understood in a cyclic sense.
inline int SumCyclic(const vector<int>& v,int start,int end)
{int sum = 0;
 if (start <= end) for (int j = start;  j <= end;  j += 1)  sum += v[j];
 else {for (int j = start;  j < v.size();  j += 1)  sum += v[j];
       for (int j = BaseIndex;  j <= end;  j += 1)  sum += v[j];}
 return sum;
}

#define GGG VType3(Gluon,Gluon,Gluon)
#define FFG VType3(Quark,Quark,Gluon)
#define FGF VType3(Quark,Gluon,Quark)
#define GFF VType3(Gluon,Quark,Quark)

// 6/4/08
#define SSG VType3(Scalar,Scalar,Gluon)
#define SGS VType3(Scalar,Gluon,Scalar)
#define GSS VType3(Gluon,Scalar,Scalar)

#define GGGG VType4(Gluon,Gluon,Gluon,Gluon)
#define FFGG VType4(Quark,Quark,Gluon,Gluon)
#define FGFG VType4(Quark,Gluon,Quark,Gluon)
#define FGGF VType4(Quark,Gluon,Gluon,Quark)
#define GFGF VType4(Gluon,Quark,Gluon,Quark)
#define GGFF VType4(Gluon,Gluon,Quark,Quark)
#define GFFG VType4(Gluon,Quark,Quark,Gluon)
#define FFFF VType4(Quark,Quark,Quark,Quark)

// 6/4/08
#define SSGG VType4(Scalar,Scalar,Gluon,Gluon)
#define SGSG VType4(Scalar,Gluon,Scalar,Gluon)
#define SGGS VType4(Scalar,Gluon,Gluon,Scalar)
#define GSGS VType4(Gluon,Scalar,Gluon,Scalar)
#define GGSS VType4(Gluon,Gluon,Scalar,Scalar)
#define GSSG VType4(Gluon,Scalar,Scalar,Gluon)
#define SSSS VType4(Scalar,Scalar,Scalar,Scalar)


#define PPP HType3(1,1,1)
#define MPP HType3(-1,1,1)
#define PMP HType3(1,-1,1)
#define PPM HType3(1,1,-1)
#define MMP HType3(-1,-1,1)
#define MPM HType3(-1,1,-1)
#define PMM HType3(1,-1,-1)
#define MMM HType3(-1,-1,-1)

#define PPPP HType4(1,1,1,1)

#define PPPM HType4(1,1,1,-1)
#define PPMP HType4(1,1,-1,1)
#define PMPP HType4(1,-1,1,1)
#define MPPP HType4(-1,1,1,1)

#define PPMM HType4(1,1,-1,-1)
#define PMPM HType4(1,-1,1,-1)
#define PMMP HType4(1,-1,-1,1)
#define MPMP HType4(-1,1,-1,1)
#define MMPP HType4(-1,-1,1,1)
#define MPPM HType4(-1,1,1,-1)

#define MMMP HType4(-1,-1,-1,1)
#define MMPM HType4(-1,-1,1,-1)
#define MPMM HType4(-1,1,-1,-1)
#define PMMM HType4(1,-1,-1,-1)

#define MMMM HType4(-1,-1,-1,-1)

namespace Tree {
#if 0
inline void FixZero(Cmom<R>& v) {
 #if FIXZERO
  C e = v.E(), x = v.X(), y = v.Y(), z = v.Z();
    // Make sure that "-0" doesn't show up as imaginary components, so
    // we always stay on one side of the branch cut; == 0 will pick up
    // both "0" and "-0"
  if (Im(e) == 0) e = C(Re(e),0);
  if (Im(x) == 0) x = C(Re(x),0);
  if (Im(y) == 0) y = C(Re(y),0);
  if (Im(z) == 0) z = C(Re(z),0);
  v = Cmom<R>(e,x,y,z);
 #endif
}

template<class T> inline void
  FixZero(momentum<complex<T> >& v) {
 #if FIXZERO
  complex<T> e = v.E(), x = v.X(), y = v.Y(), z = v.Z();
    // Make sure that "-0" doesn't show up as imaginary components, so
    // we always stay on one side of the branch cut; == 0 will pick up
    // both "0" and "-0"
  if (Im(e) == 0) e = C(Re(e),0);
  if (Im(x) == 0) x = C(Re(x),0);
  if (Im(y) == 0) y = C(Re(y),0);
  if (Im(z) == 0) z = C(Re(z),0);
  v = momentum<complex<T> >(e,x,y,z);
 #endif
}
#endif

#if 0
// Computes the sum of the momenta v[start..end], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
template<class T> int
  MomentumSum(momentum_configuration<T>& k,
                const vector<int>& v, int start, int end)
{/* For massless complex momenta, we must be careful to preserve the
    separate lambda and lambda-tilde spinors, which will in general not
    be identical to those produced by the "insert" below, even if k[sum]
    is identical to the original momentum.  In addition, this is faster... */
 if (start is end) return(v[start]);
 string key = GenKey("ms",start,end,v);
// cout << "In MS " << start << ", " << end << endl;
 size_t index;
 if (not k.get_label(key,index)) {
    momentum<complex<T> > sum;  // Initialized to 0; don't use Cmom, to avoid recomputing
    // lambda & lambda-tilde at every +=
    if (start <= end)
       for (int j = start;  j <= end;  j += 1)
{
         sum += k.mom(v[j]);
         //         cout << v[j] << ": " << k.mom(v[j]) << endl;
}
    else {for (int j = start;  j < v.size();  j += 1)  sum += k.mom(v[j]);
       for (int j = BaseIndex;  j <= end;  j += 1)
         sum += k.mom(v[j]);}
 #if FIXZERO
    FixZero(sum);
 #endif
    //    cout << "sum: " << sum << endl;
    index = k.insert(sum);
    k.put_label(key,index);}
 // cout << "sum [" << index << "]: " << k.p(index) << endl;
 return index;
}
#endif

template <class T> momentum<complex<T> > GenerateMomentum(const T& dummy);

#include "zero.h"

#if 0
// Computes the flatted sum of the momenta v[start..end], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
template<class T> /* static */ int FlatSum(momentum_configuration<T>& k,
            int ref /* index of reference momentum */,
            const vector<int>& v, int start, int end)
{// cout << "In FS " << start << ", " << end << endl;
 /* For massless complex momenta, we must be careful to preserve the
    separate lambda and lambda-tilde spinors, which will in general not
    be identical to those produced by the "insert" below, even if k[sum]
    is identical to the original momentum. */
 if (start is end and IsZero(k.m2(v[start]))) return(v[start]);
 int sum = MomentumSum(k,v,start,end);
 string key = GenKey("fs",start,end,ref,v);
 size_t index;
 if (not k.get_label(key,index)) {
    Cmom<T> flat = k[sum] - (k.m2(sum)/(T(2)*(k[sum]*k[ref]))) * k[ref];
 #if FIXZERO
    FixZero(flat);
 #endif
    index = k.insert(flat);
#if 0
    if (!IsZero(Im(k[sum].E())))
       {T dummy; momentum<complex<T> > refMom2 = GenerateMomentum(dummy);
       int ref2 = k.insert(refMom2);
       momentum<complex<T> > refMom3 = GenerateMomentum(dummy);
       int ref3 = k.insert(refMom3);
        C value1 = k.spab(ref,sum,ref2);
        C value2 = k.spab(ref,index,ref2);
        C value3 = k.spa(ref,index)*k.spb(index,ref2);
        cout << "CFS: " << value2/value1 << " "
             <<value3/value1 << endl;
        value1 = k.spab(ref3,sum,ref2);
        value2 = k.spab(ref3,index,ref2);
        value3 = k.spa(ref3,index)*k.spb(index,ref2);
        cout << "CFS a: " << value2/value1 << " "
             <<value3/value1 << endl;
        }
#endif
    k.put_label(key,index);}
 return index;
}
#endif

#if 0
// Computes the flatted sum of the negative of the
//   momenta v[start1..end1], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1)
{// cout << "In NFS " << s1 << ", " << e1 << "; " << s2 << ", " << e2 << endl;

 string key = GenKey("nf",s1,e1,ref,v);
 size_t index;
 if (not k.get_label(key,index)) {
   /* For massless complex momenta, we must be careful to preserve the
      separate lambda and lambda-tilde spinors, which will in general not
      be identical to those produced by the "insert" after flattening below,
      even if k[sum] is identical to the original momentum. */
   if (s1 is e1 and IsZero(k.m2(v[s1])))
      {Cmom<T> negative = -k[v[s1]]; // This operation does transfer spinors properly
      index = k.insert(negative);}
   else {
     int sum = MomentumSum(k,v,s1,e1);
     // The spinor products in BlackHat alas have their branch cuts along the
     // real axis (with the old phase convention), which means we have to be
     // VERY careful to ensure that
     // Negative(FlatSum(...)) is identical -- down to the sign of
     // (infinitesimal) imaginary parts to NegativeFlatSum(...) -- hence
     // do this in two steps
     momentum<complex<T> > flat = k[sum] - (k.m2(sum)/(2.*(k[sum]*k[ref]))) * k[ref];
     flat = -flat;
 #if FIXZERO
     FixZero(flat);
 #endif
     index = k.insert(flat);
     k.put_label(key,index);}}
 return index;
}
#endif

#if 0
// Computes the flatted sum of the negative of the
//   momenta v[start1..end1] + v[start2..end2], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1, int s2, int e2)
{// cout << "In NFS " << s1 << ", " << e1 << "; " << s2 << ", " << e2 << endl;
 int sum1 = MomentumSum(k,v,s1,e1);
 int sum2 = MomentumSum(k,v,s2,e2);

 // string key = GenKey("nf",s1,e1,s2,e2,ref,v);
 // until appropriate routine is available:
 string key = GenKey("nf",MakeVector(s1,e1,s2,e2,ref),v);
 size_t index;
 if (not k.get_label(key,index)) {
   // Can't use k[sum1] notation here, ugh
   momentum<complex<T> > ksum = k.mom(sum1)+k.mom(sum2);
   int sum = k.insert(ksum);
   momentum<complex<T> > flat =
     ksum - (k.m2(sum)/(T(2)*(k[sum]*k[ref]))) * k.mom(ref);
#if 0
    cout << "NF sum: " << ksum << endl;
    cout << "NF sum^2: " << k.m2(sum) << endl;
    cout << "NF ref [" << ref << "]: " << k[ref] << endl;
    cout << "NF den: " << (2.*(k[sum]*k[ref])) << endl;
    cout << "NF den alt: " << (k.s(sum,ref)-k.m2(sum)) << endl;
    cout << "NF den alt 2: " << (2.*(k[sum].E()*k[ref].E()
                                     -k[sum].X()*k[ref].X()
                                     -k[sum].Y()*k[ref].Y()
                                     -k[sum].Z()*k[ref].Z())) << endl;
#endif
    flat = -flat;
    //    cout << "NFS 1: " << flat << endl;
 #if FIXZERO
    FixZero(flat);
 #endif
#if 0
    if (Im(flat.E()) == 0) {cout << "zero" << endl;
    flat = Cmom<R>(C(Re(flat.E()),0),C(Re(flat.X()),0),C(Re(flat.Y()),0),
                C(Re(flat.Z()),0));}
    //    if (Im(flat.E()) == 0) flat.E() *= copysign(1.0,flat.E());
#endif
    //    cout << "NFS 2: " << flat << endl;
    index = k.insert(flat);
    k.put_label(key,index);}
 // cout << "NFS 2 [" << index << "]: " << k[index] << endl;
 return index;
}
#endif

#if 0
template<class T> int
  NegativeFlatSum(momentum_configuration<T>& k,
                  int ref /* index of reference momentum */,
                  const vector<int>& v, int s1, int e1, int s2, int e2,
                  int s3, int e3)
{
  // cout << "In NFS " << s1 << ", " << e1 << "; " << s2 << ", " << e2
  //      << s3 << ", " << e3 << endl;
 int sum1 = MomentumSum(k,v,s1,e1);
 int sum2 = MomentumSum(k,v,s2,e2);
 int sum3 = MomentumSum(k,v,s3,e3);

 // string key = GenKey("nf",s1,e1,s2,e2,s2,e3,ref,v);
 // until appropriate routine is available:
 string key = GenKey("nf",MakeVector(s1,e1,s2,e2,s3,e3,ref),v);
 size_t index;
 if (not k.get_label(key,index)) {
   // Can't use k[sum1] notation here, ugh
    momentum<complex<T> > ksum = k.mom(sum1)+k.mom(sum2)+k.mom(sum3);
    int sum = k.insert(ksum);
    momentum<complex<T> > flat = ksum - (k.m2(sum)/(T(2)*(k[sum]*k[ref]))) * k.mom(ref);
    flat = -flat;
 #if FIXZERO
    FixZero(flat);
 #endif
    index = k.insert(flat);
    k.put_label(key,index);}
 return index;
}
#endif

#if 0
template<class T> int
  Negative(momentum_configuration<T>& k, int i)
{string key = GenKey("neg",i);
 size_t index;
 if (not k.get_label(key,index)) {
   Cmom<T> kneg = -k[i]; // Note this transfers spinors properly for complex momenta
    //    cout << "N 1: " << kneg << endl;
 #if FIXZERO
    FixZero(kneg);
 #endif
    //    cout << "N 2: " << kneg << endl;
    index = k.insert(kneg);
    k.put_label(key,index);}
 return index;
}
#endif

// Vertices for one off-shell leg (with no explicit momentum argument, its
// momentum is given by momentum conservation), and some number of
// arguments with given indices, helicities, and particle IDs

#define VgggMPP(k1,k2,k3) \
          (halfI) * k.spa(k1,ref) * k.spb(k2,k3) \
             * (k.s(k1,ref)-k.s(k2,ref)-k.s(k3,ref))/ \
          (k.spa(k2,ref)*k.spa(k3,ref)*k.spb(k1,ref))

#define VgggPMM(k1,k2,k3) \
          (-halfI) * k.spb(k1,ref) * k.spa(k2,k3) \
              * (k.s(k1,ref)-k.s(k2,ref)-k.s(k3,ref))/ \
          (k.spb(k2,ref) * k.spb(k3,ref) * k.spa(k1,ref))

// Light-cone three-gluon vertex
// V_3(-K_{s1..e1}-K_{s2..e2},K_{s1..e1},K_{s2..e2})
// Lightcone means that it's the tree-gluon vertex with flatted momenta
template<class T> complex<T> inline
  Vggg(momentum_configuration<T>& k,
       int ref /* index of reference momentum */,
       int helicity0,
       const vector<int>& arg /* indices of momenta from which legs are taken */,
       int s1, int e1, int helicity1,
       int s2, int e2, int helicity2)
{//string key = GenKey("Vggg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vggg",MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                                               ref),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy
 static complex<T> halfI(0,0.5); // keep macros below happy

 // cout << "ref [" << ref << "]: " << k[ref] << endl;
#define $_DEBUG1_ 0
#define $_DEBUG2_ 0
 if ( !k.get_value(key,result) || true ) {
    int k0 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2);
    int k1 = FlatSum(k,ref,arg,s1,e1);
    int k2 = FlatSum(k,ref,arg,s2,e2);
#if 0
    cout << "k0 [" <<s1 << ":" << e1 <<"," << s2 <<":" << e2 << "] (" << k0 << "): "
         << k[k0] << endl;
    cout << "k1 [" <<s1 << ":" << e1 << "]: ";
    cout << k[k1] << endl;
    cout << "k2 [" <<s2 << ":" << e2 << "]: ";
    cout << k[k2] << endl;
#endif
    switch (HelicityType(helicity0,helicity1,helicity2)) {
    case MPP:
#if $_DEBUG1_
      cout << "-++" << endl;
#endif
      result = VgggMPP(k0,k1,k2);  break;
    case PMP:
#if $_DEBUG1_
      cout << "+-+" << endl;
#endif
   #if 0
      cout << "l_0: " << k.L(k0) << endl;
      cout << "l_q: " << k.L(ref) << endl;
      cout << "lt_0: " << k.Lt(k0) << endl;
      cout << "lt_q: " << k.Lt(ref) << endl;
      cout << "[0q]: " << k.spb(k0,ref) << endl;
      cout << "[1q]/[0q]: " << k.spb(k1,ref)/k.spb(k0,ref) << endl;
      cout << "[2q]/[0q]: " << k.spb(k2,ref)/k.spb(k0,ref) << endl;
      cout << "<12>: " << k.spa(k1,k2) << endl;
      cout << "<0q>: " << k.spa(k0,ref) << endl;
      cout << "<1q>/<0q>: " << k.spa(k1,ref)/k.spa(k0,ref) << endl;
      cout << "<2q>/<0q>: " << k.spa(k2,ref)/k.spa(k0,ref) << endl;
      cout << "q^2: " << k.m2(ref) << endl;
      cout << "(0q): " << k.s(k0,ref) << endl;
      cout << "(1q)/(0q): " << k.s(k1,ref)/k.s(k0,ref) << endl;
      cout << "(2q)/(0q): " << k.s(k2,ref)/k.s(k0,ref) << endl;
      cout << "delta/(0q): " << (k.s(k0,ref)-k.s(k1,ref)-k.s(k2,ref))/k.s(k0,ref) << endl;
      cout << "v: " << (VgggPMM(k0,k1,k2)) << endl;
   #endif
      result = VgggMPP(k1,k2,k0);  break;
    case PPM:
#if $_DEBUG1_
      cout << "++-" << endl;
#endif
      result = VgggMPP(k2,k0,k1);  break;
    case PMM:
#if $_DEBUG2_
      cout << "+--" << endl;
#endif
   #if 0
      cout << "[0q]: " << k.spb(k0,ref) << endl;
      cout << "[1q]/[0q]: " << k.spb(k1,ref)/k.spb(k0,ref) << endl;
      cout << "[2q]/[0q]: " << k.spb(k2,ref)/k.spb(k0,ref) << endl;
      cout << "<12>: " << k.spa(k1,k2) << endl;
      cout << "<0q>: " << k.spa(k0,ref) << endl;
      cout << "<1q>/<0q>: " << k.spa(k1,ref)/k.spa(k0,ref) << endl;
      cout << "<2q>/<0q>: " << k.spa(k2,ref)/k.spa(k0,ref) << endl;
      cout << "q^2: " << k.m2(ref) << endl;
      cout << "(0q): " << k.s(k0,ref) << endl;
      cout << "(1q)/(0q): " << k.s(k1,ref)/k.s(k0,ref) << endl;
      cout << "(2q)/(0q): " << k.s(k2,ref)/k.s(k0,ref) << endl;
      cout << "delta/(0q): " << (k.s(k0,ref)-k.s(k1,ref)-k.s(k2,ref))/k.s(k0,ref) << endl;
      cout << "v: " << (VgggPMM(k0,k1,k2)) << endl;
   #endif
      result = VgggPMM(k0,k1,k2);  break;
    case MPM:
#if $_DEBUG2_
      cout << "-+-" << endl;
#endif
      result = VgggPMM(k1,k2,k0);  break;
    case MMP:
   #if 0
      cout << "k1 b: " << k.p(k1) << endl;
      cout << "k2 b: " << k.p(k2) << endl;
      cout << k.p(k0) << endl;
      cout << "[2q]: " << k.spb(k2,ref) << endl;
      cout << "<01>: " << k.spa(k0,k1) << endl;
      cout << "[1q]: " << k.spb(k1,ref) << endl;
      cout << "[0q]: " << k.spb(k0,ref) << endl;
      cout << "<2q>: " << k.spa(k2,ref) << endl;
      cout << "(0q): " << k.s(k0,ref) << endl;
      cout << "(1q): " << k.s(k1,ref) << endl;
      cout << "(2q): " << k.s(k2,ref) << endl;
      cout << "delta: " << (k.s(k2,ref)-k.s(k0,ref)-k.s(k1,ref)) << endl;
      cout << "v: " << (VgggPMM(k2,k0,k1)) << endl;
   #endif
#if $_DEBUG2_
      cout << "--+" << endl;
#endif
      result = VgggPMM(k2,k0,k1);  break;
    case PPP:
    case MMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vggg]";
    }
    k.put_value(key,result);}
 // cout << "Vggg: " << result << endl;
 return result;
}

// In Blackhat, spa(q,x)/spa(q,-x) is 1 and spb(q,x)/spb(q,-x) is -1
// as opposed to the original Mathematica implementation with both being i.
// do we need to compensate somewhere???
#define VffgMPP(k1,k2n,k3) \
 -I * k.spa(ref,k1) * k.spb(k2n,k3)/(k.spa(ref,k3))

#define VffgPMP(k1n,k2,k3) \
  I * k.spb(k1n,k3) * k.spa(ref,k2)/(k.spa(ref,k3))

#define VffgMPM(k1,k2n,k3) \
  -I * k.spa(k1,k3) * k.spb(ref,k2n)/(k.spb(k3,ref))

#define VffgPMM(k1n,k2,k3) \
  I * k.spb(ref,k1n) * k.spa(k2,k3)/(k.spb(k3,ref))


// Light-cone fermion-fermion-gluon vertex
// Convention: fermion arrow always flows form - helicity to + helicity
// V_3(-K_{s1..e1}-K_{s2..e2},K_{s1..e1},K_{s2..e2}) or cyclic perm
// Lightcone means that it's the tree vertex with flatted momenta
// To allow it to be used for ffg, gff, and fgf, we indicate
// a rotation of arguments rightwards by 0, 1, or 2.
template<class T> complex<T> inline
  Vffg(momentum_configuration<T>& k,
       int ref /* index of reference momentum */,
       int helicity0,
       const vector<int>& arg /* indices of momenta from which legs are taken */,
       int s1, int e1, int helicity1,
       int s2, int e2, int helicity2,
       int rotateRight)
{//string key = GenKey("Vffg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vffg",MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                                               ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k0 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2);
    int k1 = FlatSum(k,ref,arg,s1,e1);
    int k2 = FlatSum(k,ref,arg,s2,e2);
    int k0n, k1n;
    if (rotateRight is 2) {// fgf
      int ks = k2;  k2 = k1;  k1 = k0;  k0 = ks;
      int helicityS = helicity2;  helicity2 = helicity1;
      helicity1 = helicity0; helicity0 = helicityS;
    } else if (rotateRight is 1) {// gff
      int ks = k2;  k2 = k0;  k0 = k1;  k1 = ks;
      int helicityS = helicity2;  helicity2 = helicity0;
      helicity0 = helicity1; helicity1 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2)) {
    case MPP:
      k1n = Negative(k,k1);
      //      cout << "ffg -++" << endl;
      result = VffgMPP(k0,k1n,k2);  break;
    case PMP:
      k0n = Negative(k,k0);
      //      cout << "ffg +-+" << endl;
      result = VffgPMP(k0n,k1,k2);  break;
    case PMM:
      k0n = Negative(k,k0);
      //      cout << "ffg +--" << endl;
      result = VffgPMM(k0n,k1,k2);  break;
    case MPM:
      k1n = Negative(k,k1);
      //      cout << "ffg -+-" << endl;
      result = VffgMPM(k0,k1n,k2);  break;
    case PPM:
    case MMP:
    case PPP:
    case MMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vffg]";
    }

    // 5/20/08 New phase convention for spinors requires additional phase
    // here in order to get conventional phase for two-quark amplitudes
    result = -result;
    k.put_value(key,result);}
 return result;
}

// Note that it's OK to use either flatted k1 or unflatted k1
#define VssgMPP(k0,k1,k2) \
  I*k.spab(ref,k1,k2)/k.spa(ref,k2)
#define VssgPMP(k0,k1,k2) VssgMPP(k0,k1,k2)

#define VssgMPM(k0,k1,k2) \
  -I*k.spab(k2,k1,ref)/k.spb(ref,k2)
#define VssgPMM(k0,k1,k2) VssgMPM(k0,k1,k2)

#if 0
// By convention, we take scalar "helicity" label to be conserved even
// in the massive case
#define VssgMMP(k0,k1,k2) VssgMPP(k0,k1,k2)
#define VssgPPP(k0,k1,k2) VssgMPP(k0,k1,k2)
#define VssgMMM(k0,k1,k2) VssgMPM(k0,k1,k2)
#define VssgPPM(k0,k1,k2) VssgMPM(k0,k1,k2)
#endif

// Light-cone massive-scalar gluon three-point vertex
template<class T> complex<T> inline
  Vssg(momentum_configuration<T>& k,
              int ref /* index of reference momentum */,
              int helicity0,
              const vector<int>& arg /* indices of momenta from which legs are taken */,
              int s1, int e1, int helicity1,
              int s2, int e2, int helicity2,
              int rotateRight)
{//string key = GenKey("Vssg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vssg",MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                                               ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k0 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2);
    int k1 = FlatSum(k,ref,arg,s1,e1);
    int k2 = FlatSum(k,ref,arg,s2,e2);
    if (rotateRight is 2) {// sgs
      int ks = k2;  k2 = k1;  k1 = k0;  k0 = ks;
      int helicityS = helicity2;  helicity2 = helicity1;
      helicity1 = helicity0; helicity0 = helicityS;
    } else if (rotateRight is 1) {// gss
      int ks = k2;  k2 = k0;  k0 = k1;  k1 = ks;
      int helicityS = helicity2;  helicity2 = helicity0;
      helicity0 = helicity1; helicity1 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2)) {
    case MPP:
      //      cout << "ssg -++" << endl;
      result = VssgMPP(k0,k1,k2);  break;
    case PMP:
      //      cout << "ssg +-+" << endl;
      result = VssgPMP(k0n,k1,k2);  break;
    case PMM:
      //      cout << "ssg +--" << endl;
      result = VssgPMM(k0n,k1,k2);  break;
    case MPM:
      //      cout << "ssg -+-" << endl;
      result = VssgMPM(k0,k1,k2);  break;
#if 0
      // Needed for massive case:
    case PPM:
      result = VssgPPM(k0,k1,k2);  break;
    case MMP:
      result = VssgMMP(k0,k1,k2);  break;
    case PPP:
      result = VssgPPP(k0,k1,k2);  break;
    case MMM:
      result = VssgMMM(k0,k1,k2);  break;
#else
    case PPM:
    case MMP:
    case PPP:
    case MMM:
      result = complex<T>(0,0);  break;
#endif
    default:
      throw "Illegal helicity configuration [Vssg]";
    }

    // Phase convention
    result = -result;
    k.put_value(key,result);}
 return result;
}

#define VMMPP(k1,k2,k3,k4) \
  -((k.s(k1,ref)*k.s(k3,ref) + k.s(k2,ref)*k.s(k4,ref)) \
          *k.spa(k1,ref)*k.spa(k2,ref)*k.spb(k3,ref)*k.spb(k4,ref))/ \
  (square(k.s(k2,ref) + k.s(k3,ref))*k.spa(k3,ref)*k.spa(k4,ref) \
   *k.spb(k1,ref)*k.spb(k2,ref))
#define VMPMP(k1,k2,k3,k4) \
 (- ((k.s(k2,ref)*k.s(k3,ref) + k.s(k1,ref)*k.s(k4,ref)) \
       *k.spa(k1,ref)*k.spa(k3,ref)*k.spb(k2,ref)*k.spb(k4,ref))/ \
  ((k.s(k1,ref) + k.s(k2,ref))*(k.s(k3,ref) + k.s(k4,ref))*k.spa(k2,ref)* \
   k.spa(k4,ref)*k.spb(k1,ref)*k.spb(k3,ref)) - \
 ((k.s(k1,ref)*k.s(k2,ref) + k.s(k3,ref)*k.s(k4,ref))*k.spa(k1,ref) \
    *k.spa(k3,ref)*k.spb(k2,ref)* k.spb(k4,ref))/ \
   ((k.s(k2,ref) + k.s(k3,ref))*(k.s(k1,ref) + k.s(k4,ref))*k.spa(k2,ref)* \
    k.spa(k4,ref)*k.spb(k1,ref)*k.spb(k3,ref)) )

// Light-cone four-gluon vertex
// V_4(-K_{s1..e1}-K_{s2..e2}-K_{s3..e3},K_{s1..e1},K_{s2..e2},K_{s3..e3})
// Lightcone means that it's the tree-gluon vertex with flatted momenta
template<class T> complex<T> inline
  Vgggg(momentum_configuration<T>& k,
        int ref /* index of reference momentum */,
        int helicity0,
        const vector<int>& arg /* indices of momenta from which legs
                                  are taken */,
        int s1, int e1, int helicity1,
        int s2, int e2, int helicity2,
        int s3, int e3, int helicity3)
{//string key = GenKey("Vgggg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vgggg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref),arg);
 complex<T> result;

 if ( !k.get_value(key,result)||true ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case PPMM:
      result = VMMPP(k3,k4,k1,k2);  break;
    case MPPM:
      result = VMMPP(k4,k1,k2,k3);  break;
    case MMPP:
      result = VMMPP(k1,k2,k3,k4);  break;
    case PMMP:
      result = VMMPP(k2,k3,k4,k1);  break;
    case PMPM:
      result = VMPMP(k2,k3,k4,k1);  break;
    case MPMP:
      result = VMPMP(k1,k2,k3,k4);  break;
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vgggg]";
    }
    k.put_value(key,result);}
 return result;
}

// V4[{k1_,F,P},{k2_,F,M},{k3_,G,M},{k4_,G,P}] :=
#define VffggPMMP(k1n,k2,k3,k4) \
  I * k.spab(k2,ref,k1n) * (k.s(ref,k3)-k.s(ref,k4))* \
                 k.spa(ref,k3) * k.spb(ref,k4)/ \
  (k.spb(ref,k3) * k.spa(ref,k4) * square(k.s(ref,k1n)-k.s(ref,k2)))
// V4[{k1_,F,P},{k2_,F,M},{k3_,G,P},{k4_,G,M}] :=
#define VffggPMPM(k1n,k2,k3,k4) \
  -I * (k.spb(ref,k1n) * k.spab(k4,ref,k3)* k.spa(ref,k2)/ \
        (k.spb(k4,ref) * k.spa(ref,k3) * \
          (k.s(k4,ref)-k.s(k1n,ref))) \
   -k.spab(k2,ref,k1n) * (k.s(ref,k3)-k.s(ref,k4)) * \
        k.spa(ref,k4) * k.spb(ref,k3)/ \
     (k.spb(ref,k4) * k.spa(ref,k3) * square(k.s(ref,k1n)-k.s(ref,k2))) )

// V4[{k1_,F,M},{k2_,F,P},{k3_,G,P},{k4_,G,M}] :=
#define VffggMPPM(k1,k2n,k3,k4) \
  -I * k.spab(k1,ref,k2n) * (k.s(ref,k3)-k.s(ref,k4)) * \
                 k.spa(ref,k4) * k.spb(ref,k3)/ \
         (k.spb(ref,k4) * k.spa(ref,k3) * square(k.s(ref,k1)-k.s(ref,k2n)))

// V4[{k1_,F,M},{k2_,F,P},{k3_,G,M},{k4_,G,P}] :=
#define VffggMPMP(k1,k2n,k3,k4) \
  I * (k.spb(ref,k2n) * k.spab(k3,ref,k4) * k.spa(ref,k1)/ \
      (k.spb(k3,ref) * k.spa(ref,k4) * (k.s(k1,ref)+k.s(k4,ref))) \
     -k.spab(k1,ref,k2n) * (k.s(ref,k3)-k.s(ref,k4)) * \
           k.spa(ref,k3) * k.spb(ref,k4)/ \
      (k.spb(ref,k3) * k.spa(ref,k4) * square(k.s(ref,k1)-k.s(ref,k2n))))

//V4[{k1_,F,P},{k2_,G,P},{k3_,F,M},{k4_,G,M}] :=
#define VfgfgPPMM(k1n,k2,k3,k4) \
  I * k.spb(ref,k1n) * k.spab(k4,ref,k2) * k.spa(ref,k3)/ \
         (k.spb(k4,ref) * k.spa(ref,k2) * (k.s(k4,ref)-k.s(k1n,ref)))
//V4[{k1_,F,P},{k2_,G,M},{k3_,F,M},{k4_,G,P}] :=
#define VfgfgPMMP(k1n,k2,k3,k4) \
  I * k.spb(ref,k1n) * k.spab(k2,ref,k4) * k.spa(ref,k3)/ \
         (k.spb(k2,ref) * k.spa(ref,k4) * (k.s(k2,ref)-k.s(k1n,ref)))
//V4[{k1_,F,M},{k2_,G,M},{k3_,F,P},{k4_,G,P}] :=
#define VfgfgMMPP(k1,k2,k3n,k4) \
  -I * k.spb(ref,k3n) * k.spab(k2,ref,k4) * k.spa(ref,k1)/ \
         (k.spb(k2,ref) * k.spa(ref,k4) * (k.s(k1,ref)+k.s(k4,ref)))
//V4[{k1_,F,M},{k2_,G,P},{k3_,F,P},{k4_,G,M}] :=
#define VfgfgMPPM(k1,k2,k3n,k4) \
  -I * k.spb(ref,k3n) * k.spab(k4,ref,k2) * k.spa(ref,k1)/ \
         (k.spb(k4,ref) * k.spa(ref,k2) * (k.s(k1,ref)+k.s(k2,ref)))

// Light-cone two-gluon two-fermion vertex ffgg
// V_4(-K_{s1..e1}-K_{s2..e2}-K_{s3..e3},K_{s1..e1},K_{s2..e2},K_{s3..e3})
// Lightcone means that it's the tree vertex with flatted momenta
// To allow it to be used for ffgg, gffg, ggff, and fggf, we indicate
// a rotation of arguments rightwards by 0, 1, 2, or 3.
template<class T> complex<T> inline
  Vffgg(momentum_configuration<T>& k,
               int ref /* index of reference momentum */,
               int helicity0,
               const vector<int>& arg /* indices of momenta from which legs
                                         are taken */,
               int s1, int e1, int helicity1,
               int s2, int e2, int helicity2,
               int s3, int e3, int helicity3,
               int rotateRight)
{//string key = GenKey("Vffgg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vffgg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy
 // const C fudge = 0.5 * I;

 if ( !k.get_value(key,result) ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
    int k1n, k2n;
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    if (rotateRight is 1) {// gffg
      int ks = k4;  k4 = k1;  k1 = k2;  k2 = k3;  k3 = ks;
      int helicityS = helicity3;  helicity3 = helicity0;
      helicity0 = helicity1; helicity1 = helicity2; helicity2 = helicityS;
    }
    else if (rotateRight is 2) {// ggff
      int ks = k4;  k4 = k2;  k2 = ks;  ks = k3;  k3 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity1;
      helicity1 = helicityS;
      helicityS = helicity2;  helicity2 = helicity0;  helicity0 = helicityS;
    }
    else if (rotateRight is 3) {// fggf
      int ks = k4;  k4 = k3;  k3 = k2;  k2 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity2;
      helicity2 = helicity1; helicity1 = helicity0; helicity0 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case PMMP:
      k1n = Negative(k,k1);
      result = VffggPMMP(k1n,k2,k3,k4);  break;
    case PMPM:
      k1n = Negative(k,k1);
      result = VffggPMPM(k1n,k2,k3,k4);  break;
    case MPPM:
      k2n = Negative(k,k2);
      result = VffggMPPM(k1,k2n,k3,k4);  break;
    case MPMP:
      k2n = Negative(k,k2);
      result = VffggMPMP(k1,k2n,k3,k4);  break;
    case PPMM:
    case MMPP:
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vffgg]";
    }
    // Extra phase seems necessary -- why??
    // 5/20/08 New phase convention requires I rather than -I to get
    // conventional overall phase for two-quark amplitudes
    result *= I;
    k.put_value(key,result);}
 return result;
}

// Light-cone two-gluon two-fermion vertex fgfg
// V_4(-K_{s1..e1}-K_{s2..e2}-K_{s3..e3},K_{s1..e1},K_{s2..e2},K_{s3..e3})
// Lightcone means that it's the tree vertex with flatted momenta
// To allow it to be used for gfgf as well, we indicate
// a rotation of arguments rightwards by 0 or 1.
template<class T> complex<T> inline
  Vfgfg(momentum_configuration<T>& k,
        int ref /* index of reference momentum */,
        int helicity0,
        const vector<int>& arg /* indices of momenta from which legs
                                  are taken */,
        int s1, int e1, int helicity1,
        int s2, int e2, int helicity2,
        int s3, int e3, int helicity3,
        int rotateRight)
{//string key = GenKey("Vfgfg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vfgfg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
    int k1n, k3n;
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    if (rotateRight is 1) {// gfgf
      int ks = k4;  k4 = k3;  k3 = k2;  k2 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity2;
      helicity2 = helicity1; helicity1 = helicity0; helicity0 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case PPMM:
      k1n = Negative(k,k1);
      result = VfgfgPPMM(k1n,k2,k3,k4);  break;
    case PMMP:
      k1n = Negative(k,k1);
      result = VfgfgPMMP(k1n,k2,k3,k4);  break;
    case MMPP:
      k3n = Negative(k,k3);
      result = VfgfgMMPP(k1,k2,k3n,k4);  break;
    case MPPM:
      k3n = Negative(k,k3);
      result = VfgfgMPPM(k1,k2,k3n,k4);  break;
    case PMPM:
    case MPMP:
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vfgfg]";
    }
    // Extra phase seems necessary -- why??
    // 5/20/08 New phase convention requires I rather than -I to get
    // conventional overall phase for two-quark amplitudes
    result *= I;
    k.put_value(key,result);}
 return result;
}

//V4[{k1_,F,M},{k2_,F,P},{k3_,F',M},{k4_,F',P}] :=
#define VffhhMPMP(k1,k2n,k3,k4n) \
  twoI * k.spab(k1,ref,k2n) * k.spab(k3,ref,k4n)/ \
         ((k.s(k1,ref)-k.s(k2n,ref))*(k.s(k3,ref)-k.s(k4n,ref)))
// Where does this sign come from??? but it is needed
#define VffhhMPPM(k1,k2n,k3n,k4) -VffhhMPMP(k1,k2n,k4,k3n)
#define VffhhPMMP(k1n,k2,k3,k4n) -VffhhMPMP(k2,k1n,k3,k4n)

#define VffhhPMPM(k1n,k2,k3n,k4) VffhhMPMP(k2,k1n,k4,k3n)

// Light-cone four-fermion vertex fff'f'
// V_4(-K_{s1..e1}-K_{s2..e2}-K_{s3..e3},K_{s1..e1},K_{s2..e2},K_{s3..e3})
// Lightcone means that it's the tree vertex with flatted momenta
// To allow it to be used for f'fff' as well, we indicate
// a rotation of arguments rightwards by 0 or 1.
template<class T> complex<T> inline
  Vffhh(momentum_configuration<T>& k,
               int ref /* index of reference momentum */,
               int helicity0,
               const vector<int>& arg /* indices of momenta from which legs
                                         are taken */,
               int s1, int e1, int helicity1,
               int s2, int e2, int helicity2,
               int s3, int e3, int helicity3,
               int rotateRight)
{//string key = GenKey("Vffhh",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vffhh",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy
 static complex<T> twoI(0,2); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
    int k1n, k2n, k3n, k4n;
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    if (rotateRight is 1) {// f'fff'
      int ks = k4;  k4 = k3;  k3 = k2;  k2 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity2;
      helicity2 = helicity1; helicity1 = helicity0; helicity0 = helicityS;
    }
#if 0
    cout << "ffhh: " << helicity0 << " " << helicity1 << " " << helicity2 <<
         " " << helicity3 << endl;
#endif
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case MPMP:
      k2n = Negative(k,k2);
      k4n = Negative(k,k4);
      result = VffhhMPMP(k1,k2n,k3,k4n);  break;
    case PMMP:
      k1n = Negative(k,k1);
      k4n = Negative(k,k4);
      result = VffhhPMMP(k1n,k2,k3,k4n);  break;
    case MPPM:
      k2n = Negative(k,k2);
      k3n = Negative(k,k3);
      result = VffhhMPPM(k1,k2n,k3n,k4);  break;
    case PMPM:
      k1n = Negative(k,k1);
      k3n = Negative(k,k3);
      result = VffhhPMPM(k1n,k2,k3n,k4);  break;
    case MMPP:
    case PPMM:
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vfgfg]";
    }
    // Extra phase seems necessary -- why??
    result *= -I;
    k.put_value(key,result);}
 return result;
}

// k1 & k2 are massive; must use flatted k1 & k2 for the following formulae
#define VssggPMMP(k1,k2,k3,k4) \
 I*k.spa(ref,k3)*k.spb(ref,k4)*(k.s(k2,ref)*k.s(k3,ref)+k.s(k1,ref)*k.s(k4,ref))/ \
   (k.spb(ref,k3)*k.spa(ref,k4)*square(k.s(ref,k3)+k.s(ref,k4)))
#define VssggPMPM(k1,k2,k3,k4) \
 I*k.spa(ref,k4)*k.spb(ref,k3)*(k.s(k2,ref)*k.s(k3,ref)+k.s(k1,ref)*k.s(k4,ref))/ \
   (k.spb(ref,k4)*k.spa(ref,k3)*square(k.s(ref,k3)+k.s(ref,k4)))
#define VssggMPPM(k1,k2,k3,k4) VssggPMPM(k1,k2,k3,k4)
#define VssggMPMP(k1,k2,k3,k4) VssggPMMP(k1,k2,k3,k4)


// Light-cone massive-scalar gluon four-vertex ssgg
// To allow it to be used for ssgg, gssg, ggss, and sggs, we indicate
// a rotation of arguments rightwards by 0, 1, 2, or 3.
template<class T> complex<T> inline
  Vssgg(momentum_configuration<T>& k,
        int ref /* index of reference momentum */,
        int helicity0,
        const vector<int>& arg /* indices of momenta from which legs
                                  are taken */,
        int s1, int e1, int helicity1,
        int s2, int e2, int helicity2,
        int s3, int e3, int helicity3,
        int rotateRight)
{//string key = GenKey("Vssgg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vssgg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    if (rotateRight is 1) {// gssg
      int ks = k4;  k4 = k1;  k1 = k2;  k2 = k3;  k3 = ks;
      int helicityS = helicity3;  helicity3 = helicity0;
      helicity0 = helicity1; helicity1 = helicity2; helicity2 = helicityS;
    }
    else if (rotateRight is 2) {// ggss
      int ks = k4;  k4 = k2;  k2 = ks;  ks = k3;  k3 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity1;
      helicity1 = helicityS;
      helicityS = helicity2;  helicity2 = helicity0;  helicity0 = helicityS;
    }
    else if (rotateRight is 3) {// sggs
      int ks = k4;  k4 = k3;  k3 = k2;  k2 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity2;
      helicity2 = helicity1; helicity1 = helicity0; helicity0 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case PMMP:
      result = VssggPMMP(k1,k2,k3,k4);  break;
    case PMPM:
      result = VssggPMPM(k1,k2,k3,k4);  break;
    case MPPM:
      result = VssggMPPM(k1,k2,k3,k4);  break;
    case MPMP:
      result = VssggMPMP(k1,k2,k3,k4);  break;
    case PPMM:
    case MMPP:
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vssgg]";
    }
    // Phase convention
    result *= I;
    k.put_value(key,result);}
 return result;
}

#define VsgsgPPMM(k1,k2,k3,k4) \
   I*k.spa(ref,k4)*k.spb(ref,k2)/(k.spa(ref,k2)*k.spb(ref,k4))
#define VsgsgPMMP(k1,k2,k3,k4) \
   I*k.spa(ref,k2)*k.spb(ref,k4)/(k.spa(ref,k4)*k.spb(ref,k2))
#define VsgsgMMPP(k1,k2,k3,k4) VsgsgPMMP(k1,k2,k3,k4)
#define VsgsgMPPM(k1,k2,k3,k4) VsgsgPPMM(k1,k2,k3,k4)

// Light-cone two-gluon two-massive-scalar vertex sgsg
// V_4(-K_{s1..e1}-K_{s2..e2}-K_{s3..e3},K_{s1..e1},K_{s2..e2},K_{s3..e3})
// Lightcone means that it's the tree vertex with flatted momenta
// To allow it to be used for gsgs as well, we indicate
// a rotation of arguments rightwards by 0 or 1.
template<class T> complex<T> inline
  Vsgsg(momentum_configuration<T>& k,
        int ref /* index of reference momentum */,
        int helicity0,
        const vector<int>& arg /* indices of momenta from which legs
                                  are taken */,
        int s1, int e1, int helicity1,
        int s2, int e2, int helicity2,
        int s3, int e3, int helicity3,
        int rotateRight)
{//string key = GenKey("Vsgsg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vsgsg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref,rotateRight),arg);
 complex<T> result;
 static complex<T> I(0,1); // keep macros below happy

 if ( !k.get_value(key,result) ) {
    int k1 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2,s3,e3);
    int k2 = FlatSum(k,ref,arg,s1,e1);
    int k3 = FlatSum(k,ref,arg,s2,e2);
    int k4 = FlatSum(k,ref,arg,s3,e3);
 #if 0
    cout << "k_0: " << k.p(k1) << endl;
    cout << s1 << ":" << e1 << ",";
    cout << s2 << ":" << e2 << ",";
    cout << s3 << ":" << e3 << endl;
    cout << "q1f: " << k.spa(ref,k1) << "[" << abs(k.spa(ref,k1)) << "]" << endl;
    cout << "q2f: " << k.spa(ref,k2) << "[" << abs(k.spa(ref,k2)) << "]" << endl;
    cout << "q3f: " << k.spa(ref,k3) << "[" << abs(k.spa(ref,k3)) << "]" << endl;
    cout << "q4f: " << k.spa(ref,k4) << "[" << abs(k.spa(ref,k4)) << "]" << endl;
   #endif
    if (rotateRight is 1) {// gsgs
      int ks = k4;  k4 = k3;  k3 = k2;  k2 = k1;  k1 = ks;
      int helicityS = helicity3;  helicity3 = helicity2;
      helicity2 = helicity1; helicity1 = helicity0; helicity0 = helicityS;
    }
    switch (HelicityType(helicity0,helicity1,helicity2,helicity3)) {
    case PPMM:
#if 0
      cout << "V: " << (VsgsgPPMM(k1,k2,k3,k4)) << endl;
      cout << "Va': " << (VsgsgPMMP(k1,k2,k3,k4)) << endl;
#endif
      result = VsgsgPPMM(k1,k2,k3,k4);  break;
    case PMMP:
#if 0
      cout << "V: " << (VsgsgPMMP(k1,k2,k3,k4)) << endl;
      cout << "Va': " << (VsgsgPPMM(k1,k2,k3,k4)) << endl;
#endif
      result = VsgsgPMMP(k1,k2,k3,k4);  break;
    case MMPP:
#if 0
      cout << "V: " << (VsgsgMMPP(k1,k2,k3,k4)) << endl;
      cout << "Va': " << (VsgsgMPPM(k1,k2,k3,k4)) << endl;
#endif
      result = VsgsgMMPP(k1,k2,k3,k4);  break;
    case MPPM:
#if 0
      cout << "V: " << (VsgsgMPPM(k1,k2,k3,k4)) << endl;
      cout << "Va': " << (VsgsgMMPP(k1,k2,k3,k4)) << endl;
#endif
      result = VsgsgMPPM(k1,k2,k3,k4);  break;
    case PMPM:
    case MPMP:
    case PPPP:
    case PPPM:
    case PPMP:
    case PMPP:
    case MPPP:
    case PMMM:
    case MPMM:
    case MMPM:
    case MMMP:
    case MMMM:
      result = complex<T>(0,0);  break;
    default:
      throw "Illegal helicity configuration [Vsgsg]";
    }
    // Phase convention
    result *= I;
    k.put_value(key,result);}
 return result;
}



// All vertices: three-point
template<class T> complex<T> static
  Vertex(momentum_configuration<T>& k,
         int ref /* index of reference momentum */,
         int helicity0,
         int id0,
         const vector<int>& arg,
         int s1, int e1, int helicity1, int id1,
         int s2, int e2, int helicity2, int id2)
{
#if 0
cout << hex << "VT: " << TypeName(id0) << TypeName(id1) << TypeName(id2)
      << VertexType(id0,id1,id2)
     << " [" << s1 <<":"<<e1<<","<<s2<<":"<<e2<<"]" << endl << dec;
#endif
// cout << hex << "VT: " << id0 << " " << id1 << " " << id2 << endl;
switch (VertexType(id0,id1,id2)) {
 case GGG:
    return(Vggg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2));
 case FFG:
    return(Vffg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,0));
 case FGF:
    return(Vffg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,2));
 case GFF:
    return(Vffg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,1));
 // 6/4/08
 case SSG:
    return(Vssg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,0));
 case SGS:
    return(Vssg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,2));
 case GSS:
    return(Vssg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,1));
 default:
   throw "Illegal vertex type [Vertex]";
}
}

// And four-point (includes contributions from shuffling around light-cone
// terms)
template<class T> complex<T>
  Vertex(momentum_configuration<T>& k,
         int ref /* index of reference momentum */,
         int helicity0,
         int id0,
         const vector<int>& arg,
         int s1, int e1, int helicity1, int id1,
         int s2, int e2, int helicity2, int id2,
         int s3, int e3, int helicity3, int id3)

{switch (VertexType(id0,id1,id2,id3)) {
 case GGGG:
    return(Vgggg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3));
 case FFGG:
    return(Vffgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,0));
 case FGFG:
    return(Vfgfg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,0));
 case FGGF:
    return(Vffgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,3));
 case GFFG:
    return(Vffgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,1));
 case GFGF:
    return(Vfgfg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,1));
 case GGFF:
    return(Vffgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,2));
 case FFFF:
    // Can we assume here that the two fermion lines have distinct flavors?
    if (id0 is id1 and id2 is id3)
     return(Vffhh(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                  s3,e3,helicity3,0));
    else if (id0 is id3 and id1 is id2)
     return(Vffhh(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                  s3,e3,helicity3,1));
 // 6/4/08
 case SSGG:
    return(Vssgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,0));
 case SGSG:
    return(Vsgsg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,0));
 case SGGS:
    return(Vssgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,3));
 case GSSG:
    return(Vssgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,1));
 case GSGS:
    return(Vsgsg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,1));
 case GGSS:
    return(Vssgg(k,ref,helicity0,arg,s1,e1,helicity1,s2,e2,helicity2,
                 s3,e3,helicity3,2));
 case SSSS: // Still missing -- only one scalar line allowed for the moment
 default:
   throw "Illegal vertex type [Vertex]";
}
}

void PrintVector(vector<ParticleID> v)
{cout << "{"; for (int i = BaseIndex;  i < v.size();  i += 1)
   {cout << v[i]; if (i < v.size()-1) cout << " ";}
 cout << "}";}

void PrintVector(vector<bool> v)
{cout << "{"; for (int i = BaseIndex;  i < v.size();  i += 1)
   {cout << v[i]; if (i < v.size()-1) cout << " ";}
 cout << "}";}

void PrintVector2(vector<int> v)
{cout << "{"; for (int i = BaseIndex;  i < v.size();  i += 1)
   {cout << v[i]; if (i < v.size()-1) cout << " ";}
 cout << "}";}

#if 0 // Now in count.cc
/* Counts the number of fermion legs of each different flavor in the
   given range of arguments within "id".  The returned vector has length
   equal to the highest flavor number in "id" plus one (the 0 component
   is not used)
*/
vector<int> FermionCount(const vector<ParticleID>& id, int start, int end)
{int max = 0;
 for (int j = BaseIndex;  j < id.size();  j += 1)
    if (IsFermion(id[j])) max = Max(max,(int)FlavorOf(id[j]));
 vector<int> count(max+1,0); // Initialized to 0
 if (start <= end)
   {for (int j = start;  j <= end;  j += 1)
     if (IsFermion(id[j])) count[FlavorOf(id[j])] += 1;}
 else {for (int j = start;  j < id.size();  j += 1)
      if (IsFermion(id[j])) count[FlavorOf(id[j])] += 1;
    for (int j = BaseIndex;  j <= end;  j += 1)
      if (IsFermion(id[j])) count[FlavorOf(id[j])] += 1;}

 return count;
}

/* Determines whether the number of fermion legs of each different
   flavor is even or not.  Return vector as in FermionCount above.
*/
vector<bool> FermionParity(const vector<ParticleID>& id, int start, int end)
{vector<int> count = FermionCount(id,start,end);
#if 0
 cout << "FP: "; PrintVector(id); cout << "[" << start << ":" << end << "]: ";
 PrintVector(count); cout << endl;
#endif
 vector<bool> parity(count.size());
 for (int j = 1;  j < count.size();  j += 1) parity[j] = IsOdd(count[j]);
 return parity;
}
#endif

void print(const vector<int>& v,int s, int e) {
  for (int i = s;  i <= e; i += 1)
    cout << v[i] << " ";
}

template <class T> /* volatile breaks T = RHP */ inline T randomR(const T& dummy){
	srand(time(NULL));
	return T(double(rand())/double(RAND_MAX));
}

// Generate null momentum; use of "unsafe" random-number generator
// should be OK here, as it's only used to generate a reference vector
template <class T> momentum<complex<T> > GenerateMomentum(const T& dummy) {
  T x = randomR(dummy);
  T y = randomR(dummy);
  T z = randomR(dummy);
  momentum<T> bare(0.,x,y,z);
  //	  momentum<R> bare(0.,3.,2.,1.);
  //  momentum<R> bare(0.,0.911981314988006,-0.98243072252415, -0.282967015045739);

 // 5/31/08 Wrong answers because of incorrect lambda generation with
 // momentum<T>, use momentum<complex<T> > instead
 return momentum<complex<T> >(sqrt(-bare*bare),bare.X(),bare.Y(),bare.Z());
}

#if 0
// Helicity ignored here
ParticleID FlavoredQuarkID(int flavor) {
  return ParticleID(quark,1,flavor,false);
}

ParticleID NParticleID(int helicity,const ParticleID& base) {
  return ParticleID(base.type(),helicity,base.flavor(),false);
}

vector<ParticleID> NParticleID(const vector<int>& helicity,
                               const vector<ParticleID>& base) {
  vector<ParticleID> result(base.size());
  for (int j = 0;  j < base.size();  j += 1)
    result[j] = NParticleID(helicity[j],base[j]);
  return result;
}

vector<ParticleID> NParticleID(const vector<int>& helicity,
                               const vector<ParticleID>& base,int n) {
  size_t size = n < base.size() ? n : base.size();
  vector<ParticleID> result(size);
  for (int j = 0;  j < size;  j += 1)
    result[j] = NParticleID(helicity[j],base[j]);
  return result;
}

// returns a characteristic integer
int ParticleCode(ParticleID id) {
  if (id.is_a(gluon)) return Gluon;
  else if (id.is_a(quark)) return(FlavoredQuark(id.flavor()));
}

vector<int> ParticleCode(const vector<ParticleID>& id) {
  vector<int> result(id.size());
  for (int j = 0;  j < id.size();  j += 1)
    result[j] = ParticleCode(id[j]);
  return result;
}

vector<int> Helicities(const vector<ParticleID>& id) {
  vector<int> result(id.size());
  for (int j = 0;  j < id.size();  j += 1)
    result[j] = id[j].helicity();
  return result;
}
#endif

bool HasQuarks(const vector<ParticleID>& leg, int start, int end)
{if (start <= end) for (int j = start;  j <= end;  j += 1)
  if (IsQuark(leg[j])) return(true);
  else {for (int j = start;  j < leg.size();  j += 1)
    if (IsQuark(leg[j])) return(true);
  for (int j = BaseIndex;  j <= end;  j += 1)
    if (IsQuark(leg[j])) return(true);}
 return(false);
}


/* Scope out fermion and scalar flavors; only one fermion flavor can be
   imbalanced, and it must match that of the off-shell leg.  Because
   we don't track the flavor structure of Yukawa couplings, any number
   of scalars are OK with a fermion line; in the absence of a fermion
   line, only one scalar flavor can be imbalanced, and it must match that
   of the off-shell leg.  Leptons, if any, must be balanced, and
   must match the flavor of one of the quark pairs. */
bool FlavorsOK(const ParticleID& offshell /* of the current's offshell leg */,
               const vector<ParticleID>& leg /* of remaining legs */,
               int start, int end /* indices into the vectors */,
               const vector<int>& coupleTo /* quark flavor */ = empty,
               int offshellMass = defaultMass,
               const vector<int>& massValue = empty)
{
  vector<bool> flavorParity = FermionParity(leg,start,end);
  bool alreadyImbalanced = false;
  int flavorO = FlavorOf(offshell);
  for (int f = 0;  f < flavorParity.size();  f += 1)
    if (flavorParity[f])
      if (alreadyImbalanced or (IsQuark(offshell) and f isnt flavorO))
        return(false);
      else alreadyImbalanced = true;

  if (IsQuark(offshell) and not alreadyImbalanced) return(false);
  // Is there something for the vectors to couple to?
  if (coupleTo.size() > 0)
     {vector<int> flavorCount = FermionCount(leg,start,end);
      for (int v = 0;  v < coupleTo.size();  v += 1)
        if (flavorCount[coupleTo[v]] is 0) return(false);}

  /* Because the indices can be quite large (and not contiguous),
     use a map rather than a vector */
  map<int,int> massCount = MassIndexCount(massValue,start,end);
  bool massFixed = false;

  /* Each mass index appearing must appear exactly twice,
     unless it's the offshell leg's mass index */
  for (map<int,int>::iterator j = massCount.begin();
       j isnt massCount.end();  ++j)
     if (j->second is 1)
       {if (massFixed or j->first isnt offshellMass) return(false);
        else massFixed = true;}
     else if ((j->second isnt 0 and j->second isnt 2)
              or j->first is offshellMass) return(false);

  if (flavorParity.size() > 0) return(true);
  flavorParity = ScalarParity(leg,start,end);
  for (int s = 0;  s < flavorParity.size();  s += 1)
    if (flavorParity[s])
      if (alreadyImbalanced or (IsScalar(offshell) and s isnt flavorO))
         return(false);
      else alreadyImbalanced = true;
  if (IsScalar(offshell) and not alreadyImbalanced) return(false);
  return(not alreadyImbalanced);
}

// Determine scalar imbalance; if no lone imbalance, return false,
// otherwise set "sImbalance" appropriately
bool ScalarImbalance(const vector<ParticleID>& leg,
                            int start, int end,
                            int scalarFlavor /* to resolve ambiguous cases */,
                            int& sImbalance)
{bool hasQuarks = HasQuarks(leg,start,end);
 sImbalance = 0;
 if (hasQuarks and scalarFlavor > 0)
    {sImbalance = scalarFlavor;  return true;}
 if (not hasQuarks)
    {// Need to known about scalar imbalance directly
      vector<bool> scalarParity = ScalarParity(leg,start,end);
      for (int s = 0;  s < scalarParity.size();  s += 1)
        if (scalarParity[s])
          if (sImbalance) return(false);
          else sImbalance = s; /* flavor-to-be of internal leg */}
 return true;
}

/* Classify internal leg in three-point vertex for current below:
     if the current's offshell leg is a gluon -> any lone imbalance (fermion
       or scalar) or none is allowed, and the internal leg has that ID;
     if the offshell leg is a fermion -> the same fermion, any scalar, or a
       gluon are allowed, and the internal leg has that ID.  If there is
       a fermion line in the subsidiary current, any scalar is allowed as
       the internal leg, even one not appearing in the current, for
       which we make use of the "scalarFlavor" hint
      (0 means a gluon in each case)
     if the offshell leg is a scalar -> the same scalar, any fermion, or a
       gluon are allowed, and the internal leg has that ID;
     in these cases, we set id1 and id2, and return true; otherwise false.

 */
bool Classify3(const ParticleID& offshell /* of the current's offshell leg */,
                      const vector<ParticleID>& leg /* of remaining legs */,
                      int start1, int end1 /* indices into the vectors of subsidiary current 1*/,
                      int start2, int end2 /* indices into the vectors of subsidiary current 2*/,
                      ParticleID& id1, ParticleID& id2 /* resulting ids */,
               int& internal1, int& internal2 /* resulting mass indices */,
                      int scalarFlavor /* of first (or second) subsidiary current,
                                          if needed to resolve ambiguities
                                          between gluons & scalars */,
               int offshellMass = defaultMass,
               const vector<int>& massValue = empty)
{int imbalance = 0;
 vector<bool> flavorParity = FermionParity(leg,start1,end1);
 int flavorO = FlavorOf(offshell);
 for (int f = 0;  f < flavorParity.size();  f += 1)
   if (flavorParity[f])
     if (imbalance) return(false);
     else imbalance = f;  // flavor-to-be of internal leg

 if (IsQuark(offshell) and imbalance isnt 0 and imbalance isnt flavorO)
   return(false);

 map<int,int> massCount = MassIndexCount(massValue,start1,end1);
 internal1 = defaultMass;  internal2 = defaultMass;
 bool massFixed = false;
 /* Each mass index appearing must appear exactly once or twice;
    if once, it's also the mass index for the internal line,
    which may not be a gluon */
 for (map<int,int>::iterator j = massCount.begin();
      j isnt massCount.end();  ++j)
    if (j->second is 1)
      if (massFixed) return(false);
      else {internal1 = j->first;
        if (offshellMass isnt j->first) internal2 = internal1;
         massFixed = true;}
    else if (j->second isnt 0 and j->second isnt 2) return(false);
  if (internal1 < 0) internal2 = offshellMass;

 /* No need to check the second subsidiary current, because the overall
    current's arguments have been checked & this means the second is OK;
    but we do need to see if there are any fermions in either current */
 if (imbalance is 0)
    {int sImbalance1;
    if (not ScalarImbalance(leg,start1,end1,scalarFlavor,sImbalance1))
      return false;
    if (IsQuark(offshell))
       {id1 = sImbalance1 > 0 ? FlavoredScalarID(sImbalance1) : GluonID;
       id2 = offshell;}
    else if (IsGluon(offshell))
       id1 = id2 = sImbalance1 > 0 ? FlavoredScalarID(sImbalance1) : GluonID;
    else {// offshell is Scalar; which way should the scalar line go?
      id1 = sImbalance1 > 0 ? offshell : GluonID;
      id2 = IsScalar(id1) ? GluonID : offshell;}}
 else /* imbalance != 0 */ {id1 = FlavoredQuarkID(imbalance);
   if (not IsQuark(offshell))
     id2 = id1;
   else {// Need to know about scalar imbalances, leg 2 could be gluon or scalar
     int sImbalance2;
     if (not ScalarImbalance(leg,start2,end2,scalarFlavor,sImbalance2))
       return false;
     id2 = sImbalance2 is 0 ? GluonID : FlavoredScalarID(sImbalance2);}}

 if (IsGluon(id1) and internal1 >= 0) return(false);
 if (IsGluon(id2) and internal2 >= 0) return(false);
 return(true);
}

/* Classify internal leg in four-point vertex for current below:

     if the current's offshell leg is a gluon -> any lone imbalance or
       none in the first daughter current is allowed; the same lone
       imbalanced flavor should flow into the second or third daughter
       current;
     if the offshell leg is a fermion, it can flow into any daughter current;
       if it flows into the first or third, any lone imbalance or none
       is allowed in the other two.

 */
inline bool Classify4(const ParticleID& offshell /* of the current's offshell leg */,
                      const vector<ParticleID>& leg /* of remaining legs */,
                      int start1, int end1 /* indices into the vectors of subsidiary current 1*/,
                      int start2, int end2 /* indices into the vectors of subsidiary current 2*/,
                      int start3, int end3 /* indices into the vectors of subsidiary current 3*/,
                      ParticleID& id1, ParticleID& id2, ParticleID& id3 /* resulting ids */,
                      int& internal1, int& internal2, int& internal3
                      /* resulting mass indices */,
                      int scalarFlavor1 /* of first subsidiary current,
                                          if needed to resolve ambiguities
                                          between gluons & scalars */,
                      int scalarFlavor2 /* of second subsidiary current,
                                          if needed to resolve ambiguities
                                          between gluons & scalars */,
                      int scalarFlavor3 /* of second subsidiary current,
                                          if needed to resolve ambiguities
                                          between gluons & scalars */,
                      int offshellMass = defaultMass,
                      const vector<int>& massValue = empty)
{int imbalance1 = 0, imbalance2 = 0;
 int flavorO = FlavorOf(offshell);
 vector<bool> flavorParity = FermionParity(leg,start1,end1);
 for (int f = 0;  f < flavorParity.size();  f += 1)
   if (flavorParity[f])
     if (imbalance1) return(false);
     else imbalance1 = f;  // flavor-to-be of internal leg

 flavorParity = FermionParity(leg,start2,end2);
 for (int f = 0;  f < flavorParity.size();  f += 1)
   if (flavorParity[f])
     if (imbalance2) return(false);
     else imbalance2 = f;  // flavor-to-be of internal leg

 if (IsQuark(offshell))
    {if (imbalance1 isnt 0 and imbalance1 isnt flavorO
         and imbalance1 isnt imbalance2) return(false);
    if (imbalance1 is 0 and imbalance2 isnt 0 and
        imbalance2 isnt flavorO) return(false);}
 else {if (imbalance1 isnt 0 and imbalance2 isnt 0
           and imbalance1 isnt imbalance2) return(false);}

 map<int,int> massCount1 = MassIndexCount(massValue,start1,end1);
 internal1 = defaultMass;  internal2 = defaultMass;  internal3 = defaultMass;
 bool massFixed = false;
 /* Each mass index appearing must appear exactly once or twice;
    if once, it's also the mass index for the internal line,
    which may not be a gluon */
 for (map<int,int>::iterator j = massCount1.begin();
      j isnt massCount1.end();  ++j)
    if (j->second is 1)
      if (massFixed) return(false);
      else {internal1 = j->first;  massFixed = true;}
    else if (j->second isnt 0 and j->second isnt 2) return(false);
 map<int,int> massCount2 = MassIndexCount(massValue,start2,end2);
 massFixed = false;
 for (map<int,int>::iterator j = massCount2.begin();
      j isnt massCount2.end();  ++j)
   if (j->second is 1)
     if (massFixed) return(false);
     else {internal2 = j->first;  massFixed = true;}
   else if (j->second isnt 0 and j->second isnt 2) return(false);
 if ((internal1 < 0 and internal2 < 0)
     or (internal1 is internal2)) internal3 = offshellMass;
 else if (internal1 < 0 and internal2 isnt offshellMass)
   internal3 = internal2;
 else if (internal2 < 0 and internal1 isnt offshellMass)
   internal3 = internal1;
 else if (internal1 is offshellMass) internal3 = internal2;
 else if (internal2 is offshellMass) internal3 = internal1;

 if (imbalance1 is 0 and imbalance2 is 0)
    {int sImbalance1, sImbalance2;
    ScalarImbalance(leg,start1,end1,scalarFlavor1,sImbalance1);
    ScalarImbalance(leg,start2,end2,scalarFlavor2,sImbalance2);
    // Note that a ffss vertex only arises from collapsing a (ffg gss)
    // system into a contact interaction -- thus the scalars must have
    // the same flavor
    if (IsQuark(offshell))
       {if (sImbalance1 is sImbalance2)
          id1 = id2 = sImbalance1 > 0 ? FlavoredScalarID(scalarFlavor1)
            : GluonID;
       else return(false);
       id3 = offshell;}
    else {// Need to know about scalar imbalances
      int sImbalance3;
      ScalarImbalance(leg,start3,end3,scalarFlavor3,sImbalance3);
      flavorO = IsScalar(offshell) ? FlavorOf(offshell) : 0;
      if (sImbalance1 is flavorO /* and sImbalance2 is sImbalance3 is
                                    then implicitly true */)
         {id1 = offshell;
         id2 = id3 = sImbalance2 > 0 ? FlavoredScalarID(sImbalance2) : GluonID;}
      else if (sImbalance2 is flavorO /* and then 1,3 are implicitly gluons
                                         or identical scalars */)
         {id2 = offshell;
         if (sImbalance1 is sImbalance3)
           id1 = id3 = sImbalance1 > 0 ? FlavoredScalarID(sImbalance1) : GluonID;
         else return false;}
      else if (sImbalance3 is flavorO)
         {id3 = offshell;
         id1 = id2 = sImbalance2 > 0 ? FlavoredScalarID(sImbalance2) : GluonID;}
      else return false;}}
 else if (imbalance1 is 0)
    {if (IsQuark(offshell) /* and flavorO is imbalance2 -- implicitly true
                              because overall arguments have been checked */)
       {int sImbalance1;
       ScalarImbalance(leg,start1,end1,scalarFlavor1,sImbalance1);
       int sImbalance3;
       ScalarImbalance(leg,start3,end3,scalarFlavor3,sImbalance3);
       id1 = sImbalance1 > 0 ? FlavoredScalarID(sImbalance1) : GluonID;
       id3 = sImbalance3 > 0 ? FlavoredScalarID(sImbalance3) : GluonID;
       id2 = offshell;}
    else {/* fermion line goes from 2->3; (implicitly) flavorO is sImbalance1 */
      id1 = offshell;
      id2 = id3 = FlavoredQuarkID(imbalance2);}}
 else if (imbalance2 is 0)
    {if (IsQuark(offshell) /* and implicitly flavorO is imbalance1 */)
       {id1 = offshell;
       int sImbalance2;
       ScalarImbalance(leg,start2,end2,scalarFlavor2,sImbalance2);
       int sImbalance3;
       ScalarImbalance(leg,start3,end3,scalarFlavor3,sImbalance3);
       id2 = sImbalance2 > 0 ? FlavoredScalarID(sImbalance2) : GluonID;
       id3 = sImbalance3 > 0 ? FlavoredScalarID(sImbalance3) : GluonID;}
    else {// fermion line goes from 1->3; others must be gluons
      id1 = id3 = FlavoredQuarkID(imbalance1);
      id2 = GluonID;}}
 else if (imbalance1 is imbalance2)
    {id1 = id2 = FlavoredQuarkID(imbalance1);
    id3 = offshell; /* Note that fermion lines are assumed to have
                       distinct flavors by the time we get here */}
 else {// Four-quark, overall arguments have been checked so this is the
   // only possibility
   id1 = offshell;  id2 = id3 = FlavoredQuarkID(imbalance2);}
 if (IsGluon(id1) and internal1 >= 0) return(false);
 if (IsGluon(id2) and internal2 >= 0) return(false);
 if (IsGluon(id3) and internal3 >= 0) return(false);
 return(true);
}

#define Complex C
#define Power pow
C KnownGluonMPPM(momentum_configuration<R>& k,int ref0,
                 int i1, int i2, int i3, int i4)
{ size_t ref /* index of reference momentum */;
 const string refKey(RefTag);
 if (ref0 >= 0) ref = ref0;
 else if (not k.get_label(refKey,ref)) {
   R dummy;
   // inserting momentum<R> doesn't seem to yield correct lambdas
   momentum<complex<R> > refMom = GenerateMomentum(dummy);
   ref = k.insert(refMom);
   k.put_label(refKey,ref);
 }

int ksp12 = FlatSum(k,ref,MakeVector(i1,i2),0,1);
 int ksp14 = FlatSum(k,ref,MakeVector(i1,i4),0,1);
 int ksp23 = FlatSum(k,ref,MakeVector(i2,i3),0,1);
 int ksp34 = FlatSum(k,ref,MakeVector(i3,i4),0,1);
 cout << "[KG: " << ref << " " << i1 << " " << i2 << " " << i3 << " "
      << i4 ;
 cout << k.spa(i1,ref) << " ";
 cout << k.spa(i2,ref) << " ";
 cout << k.spa(i3,ref) << " ";
 cout << k.spa(i4,ref) << " ";
#if 0
 cout << k.spa(ksp12,ref)/k.spa(ksp34,ref) << endl;
 cout << k.spb(ksp12,ref)/k.spb(ksp34,ref) << endl;
 cout << k.spa(i4,ksp12)*k.spb(i2,ksp34)/(k.spb(i2,ksp12)*k.spa(i4,ksp34));
 cout << k.spa(ksp14,ref)/k.spa(ksp23,ref) << endl;
 cout << k.spb(ksp14,ref)/k.spb(ksp23,ref) << endl;
 cout << k.spa(i4,ksp14)*k.spb(i2,ksp23)/(k.spb(i2,ksp14)*k.spa(i4,ksp23));
#endif
 cout << ";";
 cout << (k.spa(i1,ksp34)*k.spa(ref,ksp34)/
         (k.spa(i1,ksp12)*k.spa(ref,ksp12))) << " ";
 cout << (k.spb(i1,ksp34)*k.spb(ref,ksp34)/
         (k.spb(i1,ksp12)*k.spb(ref,ksp12))) << " ";
 cout << (k.spa(i1,ksp34)*k.spb(ref,ksp34)/
         (k.spa(i1,ksp12)*k.spb(ref,ksp12))) << " ";
 cout << "]" << endl;
 C result1, result2, result3, result4;
 C result2a,result3a,result4a;
result1 =
(Complex(0,-1)*(k.s(i1,ref)*k.s(i3,ref) + k.s(i2,ref)*k.s(i4,ref))*
 k.spa(i1,ref)*k.spa(i4,ref)*k.spb(i2,ref)*k.spb(i3,ref))/
(Power(k.s(i1,ref) + k.s(i2,ref),2)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*
 k.spb(i4,ref));
result2 =
+ (Complex(0,0.25)*(k.s(i1,ref) - k.s(i2,ref) - k.s(i3,ref) - k.s(i4,ref))*
 (-k.s(i1,ref)-k.s(i2,ref)+k.s(i3,ref)-k.s(i4,ref))*k.spa(i1,ref)*k.spb(i3,ref)
   *k.spa(i4,ksp12)*k.spb(i2,ksp34))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spb(ksp12,ref)*k.spa(ksp34,ref));
result2a =
+ (Complex(0,0.25)*(2.0*k.s(i1,ref))*(2.0*k.s(i3,ref))
    *k.spa(i1,ref)*k.spb(i3,ref)*k.spa(i4,ksp12)*k.spb(i2,ksp34))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spb(ksp12,ref)*k.spa(ksp34,ref));
result3 =
+ (Complex(0,0.25)* (-k.s(i1,ref) + k.s(i2,ref) + k.s(i3,ref) - k.s(i4,ref))*
(k.s(i1,ref)-k.s(i2,ref)-k.s(i3,ref)+k.s(i4,ref))*k.spa(i4,i1)*k.spa(ksp14,ref)*
                    k.spb(i2,i3)*k.spb(ksp23,ref))/
(k.s(i2,i3)*k.spa(i2,ref)*k.spa(i3,ref)*k.spa(ksp23,ref)*k.spb(i1,ref)*k.spb(i4,ref)*
 k.spb(ksp14,ref));
result3a =
- (Complex(0,0.25)*4.0* Power(k.s(i2,ref) + k.s(i3,ref),2)
   *k.spa(i4,i1)*k.spa(ksp14,ref)*k.spb(i2,i3)*k.spb(ksp23,ref))/
(k.s(i2,i3)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spa(ksp23,ref)*k.spb(ksp14,ref));
result4 =
+ (Complex(0,0.25)* (-k.s(i1,ref) + k.s(i2,ref) - k.s(i3,ref) - k.s(i4,ref))*
 (-k.s(i1,ref)-k.s(i2,ref)-k.s(i3,ref) + k.s(i4,ref))*k.spa(i4,ref)*k.spb(i2,ref)
  *k.spa(ksp34,i1)*k.spb(ksp12,i3))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spa(ksp12,ref)*k.spb(ksp34,ref));
result4a =
+ (Complex(0,0.25)* (2.0*k.s(i2,ref))*(2.0*k.s(i4,ref))
     *k.spa(i4,ref)*k.spb(i2,ref)*k.spa(ksp34,i1)*k.spb(ksp12,i3))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spa(ksp12,ref)*k.spb(ksp34,ref));
C result24b =
   (Complex(0,1)*(k.s(i1,ref)*k.s(i3,ref)*k.spa(i1,ref)*k.spa(i4,ksp12)*
         k.spb(i2,ksp12)*k.spb(i3,ref) +
        k.s(i2,ref)*k.s(i4,ref)*k.spa(i1,ksp12)*k.spa(i4,ref)*k.spb(i2,ref)*
         k.spb(i3,ksp12)))/
    (k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spa(ref,ksp12)*k.spb(i1,ref)*
     k.spb(i4,ref)*k.spb(ref,ksp12));
 cout << "{";
 cout << result1 << " ";
 cout << result2 << " ";
 cout << "(" << result2a << ") ";
 cout << result3 << " ";
 cout << "(" << result3a << ") ";
 cout << result4 << " ";
 cout << "(" << result4a << "); ";
 // Here, there is a minus sign; in checkrec1.m, there's no sign...
 // issue with spinor phases...
 cout << (result2+result4) << ", " << result24b << " ";
 cout << "}" << endl;

C result =
(Complex(0,-1)*(k.s(i1,ref)*k.s(i3,ref) + k.s(i2,ref)*k.s(i4,ref))*
 k.spa(i1,ref)*k.spa(i4,ref)*k.spb(i2,ref)*k.spb(i3,ref))/
(Power(k.s(i1,ref) + k.s(i2,ref),2)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*
 k.spb(i4,ref))
+ (Complex(0,0.25)*(k.s(i1,ref) - k.s(i2,ref) - k.s(i3,ref) - k.s(i4,ref))*
 (-k.s(i1,ref)-k.s(i2,ref)+k.s(i3,ref)-k.s(i4,ref))*k.spa(i1,ref)*k.spa(i4,ksp12)*
                k.spb(i2,ksp34)*k.spb(i3,ref))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spa(ksp34,ref)*k.spb(i1,ref)*k.spb(i4,ref)*
 k.spb(ksp12,ref))
+ (Complex(0,0.25)* (-k.s(i1,ref) + k.s(i2,ref) + k.s(i3,ref) - k.s(i4,ref))*
(k.s(i1,ref)-k.s(i2,ref)-k.s(i3,ref)+k.s(i4,ref))*k.spa(i4,i1)*k.spa(ksp14,ref)*
                    k.spb(i2,i3)*k.spb(ksp23,ref))/
(k.s(i2,i3)*k.spa(i2,ref)*k.spa(i3,ref)*k.spa(ksp23,ref)*k.spb(i1,ref)*k.spb(i4,ref)*
 k.spb(ksp14,ref))
+ (Complex(0,0.25)* (-k.s(i1,ref) + k.s(i2,ref) - k.s(i3,ref) - k.s(i4,ref))*
 (-k.s(i1,ref)-k.s(i2,ref)-k.s(i3,ref) + k.s(i4,ref))*k.spa(i4,ref)*k.spb(i2,ref)
  *k.spa(ksp34,i1)*k.spb(ksp12,i3))/
(k.s(i3,i4)*k.spa(i2,ref)*k.spa(i3,ref)*k.spb(i1,ref)*k.spb(i4,ref)
 *k.spa(ksp12,ref)*k.spb(ksp34,ref));

 return result;
}


void DumpHelicities(int h0, int h1, int h2) ;
void DumpHelicities(int h0, int h1, int h2, int h3) ;
void DumpTaggedHelicities(int id0, int h0,
                          int id1, int h1,
                          int id2, int h2,
                          int id3, int h3);

/* Amputated light-cone current for pure QCD tree amplitudes (the propagator
   for the off-shell leg is NOT included).
   The arguments are as follows:
      a momentum configuration;
      the helicity of the offshell leg, -1 or +1;
      the particle id of the offshell leg;
      a vector of indices into the momentum configuration;
      a vector of the corresponding particles' helicities;
      a vector of the corresponding particles' ids;
      the start and end positions within the vector.  These wrap
      around if start > end.

   Each fermion line must have a different flavor, and both ends
   must have the same particle id (but opposite helicities).  (These
   flavors will ideally have been renamed from the original ones so
   that they run from 1..# fermion lines, but this is not crucial
   to obtaining correct results.)

   The "start" and "end" arguments refer to the position in the vector,
   NOT indices within the momenumt configuration.  They are present in
   order to avoid constantly creating and destroying vectors of indices
   as we descend the recursion.  All legs assumed to be colored (not
   checked).

   With the new phase convention for spinor products introduced 5/20/08,
   the fermion amplitudes are adjusted to come out with the "canonical"
   (susy) phase conventions.

   Massive scalars introduced 6/3/08; they are "charged" (complex) and
   hence carry "helicity" ("charge").  We don't keep track of the
   flavor structure of Yukawa couplings, whose bare structure is
   "inherited" from N=4.  There are no four-scalar couplings except
   those produced by the light-cone reshuffling of terms.  These
   couplings are intended solely to allow scalars to replace gluons in
   the D-dimensional parts of loop amplitudes.

   8/4/08: For masses, an additional argument, giving the indices of
   four-vectors for each entry (ignored for gluons) whose square is the
   mass of the corresponding leg.  Each distinct index tracks a "mass
   line" (which may implicitly be the D-dimensional component of a
   D-vector), which can flow from scalar to scalar, fermion to fermion,
   or scalar to fermion (or vice versa).  Each index must appear exactly
   twice.

*/
template<class T> complex<T>
  J(momentum_configuration<T>& k,
    int ref0 /* Supply negative number to use default one */,
    ParticleID offshell /* helicity, type etc */,
#if 0
    int helicityO /* helicity of offshell leg */,
    int idO /* ParticleID of offshell leg */,
#endif
    const vector<int>& arg /* momentum labels of colored arguments */,
#if 0
    const vector<int>& helicity,
    const vector<ParticleID>& id,
#endif
    const vector<ParticleID>& leg,
    int start, int end /* indices into the vectors */,
    int offshellMass = defaultMass,
    const vector<int>& massValue = empty)
#if 0
    const vector<int>& leptons = empty
                     /* of leptons coupling via photon to the current,
                           if any; must be an even number.  Each lepton
                           pair couples to quarks of the same flavor */,
    const vector<ParticleID>& leptonIDs = emptyID)
#endif
{ // Get reference momentum
 size_t ref /* index of reference momentum */;
 const string refKey(RefTag);
 if (ref0 >= 0) ref = ref0;
 else if (not k.get_label(refKey,ref)) {
   T dummy;
   // inserting momentum<R> doesn't seem to yield correct lambdas
   momentum<complex<T> > refMom = GenerateMomentum(dummy);
   ref = k.insert(refMom);
   k.put_label(refKey,ref);
 }

#if 0
 cout << start << "," << end << endl;
#endif

   int helicityO = offshell.helicity();
   int idcodeO = ParticleCode(offshell);
#if 0
   cout << "id O: " << idcodeO << endl;
#endif

   vector<int> helicity = Helicities(leg);
   vector<int> idcode = ParticleCode(leg);

#if 0
   cout << "ids: ";
   Tree::PrintVector2(idcode);
   cout << endl;
#endif

 //string key = GenKey("J",helicityO,idO,start,end,arg,helicity,id);
 // until the right routine is available:
   string key = GenKey("J",start,end,MakeVector(ref,helicityO,idcodeO,offshellMass),
                       Join(massValue,arg),helicity,idcode);

 complex<T> result;

 if (not k.get_value(key,result)) {
    // Must compute it
    complex<T> term;

    // Number of arguments including off-shell leg
    int n = CountCyclic(arg,start,end)+1;

    //       cout << "---" << endl;
    //    cout << "J_" << n <<": " << idO << "; ";
    // PrintVector(id); cout << " [" << start << ":" << end << "] " << endl;

    // First look at special cases:
    // End-points of recursion: three-point currents are just vertices

    if (n is 3)
       {result = Vertex(k,ref,helicityO,idcodeO,arg,
                        start,start,helicity[start],idcode[start],
                        end,end,helicity[end],idcode[end]);
       //       cout << "--- " << result << endl;
        k.put_value(key,result);
        return result;}

    // All-plus vanishes, likewise all-minus (independent of number of leptons)
    int netHelicity = helicityO+SumCyclic(helicity,start,end);

    if (netHelicity == (n=CountCyclic(arg,start,end)+1)
        || netHelicity == -n)
       {result = complex<T>(0,0);  k.put_value(key,result);  return result;}
    // Small number of negative helicities, expand using MHV vertices
    // or perhaps ordinary on-shell recursion
    //else {}

    // General case: sum over currents joined to three-point and
    // four-point vertices

    if (not FlavorsOK(offshell,leg,start,end,empty,offshellMass,massValue))
      return(complex<T>(0,0));

    vector<bool> flavorParity;
    int flavorO = FlavorOf(offshell);
#if 0
    // Scope out fermion flavors; only one flavor can be imbalanced,
    // and it must match that of the off-shell leg
    vector<bool> flavorParity = FermionParity(leg,start,end);
    bool alreadyImbalanced = false;
    for (int f = 0;  f < flavorParity.size();  f += 1)
       if (flavorParity[f])
          if (alreadyImbalanced or f isnt flavorO) return(0);
          else alreadyImbalanced = true;
#endif

    result = complex<T>(0,0);

    // 5/20/08 New phase convention for spinors requires different phase
    // for two-particle "current"
    complex<T> externalLeg = complex<T>(0,1);

    // Three-point vertex terms
    for (int j1 = 1;  j1 <= n-2;  j1 += 1)
       {int midL, midR; // Indices of j1th & j1+1st legs "arg"
       // Figure out internal labels (common to all helicities)
       midL=IndexCyclic(arg,start+j1-1);
       midR=IndexCyclic(arg,start+j1);
       /* Particle type of offshell leg of first daughter current should
          be the lone imbalanced fermion flavor (argument) or gluon if its
          arguments are fermion-number balanced; as the overall imbalance
          has already been checked above, this ensures correct type
          for other daughter current */
#if 0
       int imbalance = 0;
       bool skip = false;
       flavorParity = FermionParity(leg,start,midL);
#if 0
       cout << "id[" << start << ":" << midL <<"]: ";
       PrintVector(id); cout << "; p: "; PrintVector(flavorParity);
       cout << endl;
#endif
       for (int f = 0;  f < flavorParity.size();  f += 1)
         if (flavorParity[f])
           if (imbalance) {skip = true;  break;}
           else imbalance = f;  // flavor-to-be of internal leg
       // offshell is gluon -> any lone imbalance or zero allowed
       // offshell is fermion -> same fermion or gluon allowed
       if (skip or (IsQuark(offshell) and imbalance isnt 0
                    and imbalance isnt flavorO)) continue;
       ParticleID id1 = imbalance is 0 ?
         GluonID : FlavoredQuarkID(imbalance);
       ParticleID id2 = imbalance is 0 ?
         offshell : (IsQuark(offshell) ? GluonID : id1);
#else
       ParticleID id1, id2;
       int internal1, internal2;  // Internal legs' mass indices
       if (not Classify3(offshell,leg,start,midL,midR,end,id1,id2,
                         internal1,internal2,0,offshellMass,massValue))
         continue;
#endif

#if 0
       cout << "n: " << n << " (" << j1 << "," << midL << ")" << endl;
       cout << "id: " << idcodeO << " " << ParticleCode(id1) << " " << ParticleCode(id2) << endl;
       cout << IsQuark(offshell) << "; " << imbalance << endl;
#endif

       // Sum over helicities: note that h_i is the helicity on the vertex,
        // the negative of the helicity on the current
        for (int h1 = -1;  h1 <= 1;  h1 += 2)
        for (int h2 = -1;  h2 <= 1;  h2 += 2)
           {complex<T> J1, J2;
            // If not an allowed configuration, skip it
            if ((h1+h2+helicityO == 3 or h1+h2+helicityO == -3)
                and massValue.size() is 0 and offshellMass < 0) continue;
            /* For allowed configurations, vertex is always non-vanishing;
               we may want to compute shorter current first, but just do
               them in order for now. */
            if (j1 is 1) // 2-pt current: do helicities match?
               if (helicity[start] isnt h1) continue; else J1 = externalLeg;
            else J1 = J(k,ref,NParticleID(-h1,id1),arg,leg,start,midL,
                        internal1,massValue);
            if (J1 is complex<T>(0,0)) continue;
            // Fermion phase convention: extra factor of -i for +;- configuration
            // Is this phase independent of the sign of the energy?
            //            if (j1 is 1 and IsQuark(id1) and h1 < 0) J1 = -J1;
            if (j1 is n-2) // 2-pt current: do helicities match?
               if (helicity[end] isnt h2) continue; else J2 = externalLeg;
            else J2 = J(k,ref,NParticleID(-h2,id2),arg,leg,midR,end,
                        internal2,massValue);
            if (J2 is complex<T>(0,0)) continue;
            //            if (j1 is n-2 and IsQuark(id2) and h2 < 0) J2 = -J2;
            int sum1 = MomentumSum(k,arg,start,midL);
            int sum2 = MomentumSum(k,arg,midR,end);
            // Need to get at scalar masses for massive case...
            complex<T> prop1 = k.m2(sum1);
            complex<T>  prop2 = k.m2(sum2);
            complex<T> props;
            if (j1 is 1) // Because n>3 here, j1 cannot be both 1 and n-2
               {// No massive quarks yet...
                 if (IsScalar(id2) and internal2 >= 0)
                   prop2 -= k.m2(internal2);
                 props = prop2;}
            else if (j1 is n-2)
               {if (IsScalar(id1) and internal1 >= 0)
                   prop1 -= k.m2(internal1);
                 props = prop1;}
            else {if (IsScalar(id1) and internal1 >= 0)
                 prop1 -= k.m2(internal1);
              if (IsScalar(id2) and internal2 >= 0)
                 prop2 -= k.m2(internal2);
               props = prop1 * prop2;}
            // Phase convention for fermions?
#if 0
            cout << "h_in: ";
            print(helicity,start,end);
            cout << endl;
            cout << "h: " << helicityO << " " << h1 << " " << h2 << endl;
            cout << "J1: " << J1 << "; J2: " << J2;
            cout << "; V: " << Vertex(k,ref,helicityO,idO,arg,
                                      start,midL,h1,id1,  midR,end,h2,id2);
            cout << endl;
#endif
#if 0
            if (n is 4) {
            cout << "V3 [";  DumpHelicities(helicityO,h1,h2);
            cout << "; " << start << ".." << midL << ";" << midR << ".." << end;
            cout << "]: " << Vertex(k,ref,helicityO,idcodeO,arg,
                          start,midL,h1,ParticleCode(id1),
                          midR,end,h2,ParticleCode(id2)) << endl;
            }
#endif
            term = Vertex(k,ref,helicityO,idcodeO,arg,
                          start,midL,h1,ParticleCode(id1),
                          midR,end,h2,ParticleCode(id2))
                      * J1 * J2/props;
            // 5/20/08 New phase convention for spinors requires additional sign
#if 0
            if (IsQuark(id1) and j1 isnt 1) term = -term;
            if (IsQuark(id2) and j1 isnt n-2) term = -term;

            if (IsQuark(id1) and j1 is 1) term = -term;
            if (IsQuark(id2) and j1 is n-2) term = -term;
#endif
#if QuarkFlip
            if (IsQuark(id1)) term = -term;
            if (IsQuark(id2)) term = -term;
#endif
#if 0
            cout << "n: " << n << "; " << j1 << "; ref: " << ref << endl;
            cout << "mv: " << massValue << endl;
            cout << "h: " << helicityO << " " << h1 << " " << h2 << endl;
            cout << "J1: " << J1 << "; J2: " << J2;
            cout << "; V: " << Vertex(k,ref,helicityO,idcodeO,arg,
                                      start,midL,h1,ParticleCode(id1),
                                      midR,end,h2,ParticleCode(id2));
            cout << "; props: " << props;
            cout << endl;
            cout << "id: " << idcodeO << " " << ParticleCode(id1) << " " << ParticleCode(id2) << endl;
#endif
#if 0
            if (n is 6) {
              cout << "j1: " << j1 << "; ref: " << ref << endl;
              cout << "J1: " << J1 << "; J2: " << J2 << "; ";
#if 0 // ***ss
            if (j1 is 2)
               {C J2e = -1.0/(k.spa(j1+1,j1+2)*(k.s(j1+2,j1+3)-k.m2(j1+3)))*
                  (k.m2(j1+3)*k.spb(j1+1,j1+2)
                   -k.spab(ref,j1+3,j1+2)*(k.s(j1+1,j1+2,j1+3)-k.m2(j1+3))/
                   k.spa(ref,j1+1));
               cout << "J2 e: " << J2e;}
            else if (j1 is 3)
               {C J2e = -1.0/(k.spa(ref,j1+1))*k.spab(ref,j1+2,j1+1);
               cout << "J2 e(" << (j1+1) << "," << (j1+2) << "): " << J2e;}
#else // s***s
            if (j1 is 4)
               {C J1e = -1.0/(k.spa(2,3)*k.spa(3,4)*(k.s(1,2)-k.m2(1)))*
                  (k.m2(1)*(k.spbb(4,2,1,2)+k.spbb(4,3,1,2))/
                   (k.s(1,2,3)-k.m2(1))
                   -(k.spab(ref,1,2)-k.m2(1)*k.spab(ref,3,2)/
                   (k.s(1,2,3)-k.m2(1)))
                   *(k.s(1,2,3,4)-k.m2(1))/k.spa(ref,4));
               cout << "J1 e(2,3,4): " << J1e;}
            else if (j1 is 3)
               {C J1e = -1.0/(k.spa(2,3)*(k.s(1,2)-k.m2(1)))*
                  (k.m2(1)*k.spb(3,2)
                   -k.spab(ref,1,2)*(k.s(1,2,3)-k.m2(1))/
                   k.spa(ref,3));
               cout << "J1 e(2,3): " << J1e;}
            else if (j1 is 2)
               {C J1e = -1.0/(k.spa(ref,2))*k.spab(ref,1,2);
               cout << "J1 e(" << (1) << "," << (2) << "): " << J1e;}
#endif
            cout << endl;
            cout << "term: " << term << endl;
            }
#endif
#if 0
            if (n is 4) {
            cout << "term3 [";  DumpHelicities(helicityO,h1,h2);
            cout << "; " << start << ".." << midL << ";" << midR << ".." << end;
            cout << "]: " << term << endl;
            if (j1 isnt 1 and (h1 is -1 or h1 is 1))
              cout << "J1 [" << h1 << "]: " << J1 << endl;
            if (/* j1 isnt n-2 and */(h2 is -1 or h2 is 1))
              cout << "J2 [" << h2 << "]: " << J2 << endl;
            cout << "props: " << props << endl;
            }
#endif
            result -= term;}}

    // Four-point vertex terms
    for (int j1 = 1;  j1 <= n-3;  j1 += 1)
    for (int j2 = j1+1;  j2 <= n-2;  j2 += 1)
       {
         int mid1L, mid1R, mid2L, mid2R;
         /* Particle type of offshell leg of first daughter current should
            be the lone imbalanced fermion flavor (argument) or gluon if its
            arguments are fermion-number balanced; as the overall imbalance
            has already been checked above, this ensures correct type
            for other daughter current */

         /* For allowed configurations, vertex is always non-vanishing
            so compute shortest current first. */
         mid1L=IndexCyclic(arg,start+j1-1);
         mid1R=IndexCyclic(arg,start+j1);
         mid2L=IndexCyclic(arg,start+j2-1);
         mid2R=IndexCyclic(arg,start+j2);
         /* Particle type of offshell leg of first daughter current should
            be the lone imbalanced fermion flavor (argument), the same lone
            imbalanced fermion flavor of the second daughter current, or
            gluon if its
            arguments are fermion-number balanced; as the overall imbalance
            has already been checked above, this ensures correct type
            for other daughter current */
#if 1
         ParticleID id1, id2, id3;
         int internal1, internal2, internal3; // Internal lines' mass indices
         if (not Classify4(offshell,leg,start,mid1L,mid1R,mid2L,mid2R,end,
                           id1,id2,id3,internal1,internal2,internal3,0,0,0,
                           offshellMass,massValue))
            continue;
#else
         int imbalance1 = 0, imbalance2 = 0;
         bool skip = false;
         flavorParity = FermionParity(leg,start,mid1L);
         for (int f = 0;  f < flavorParity.size();  f += 1)
           if (flavorParity[f])
             if (imbalance1) {skip = true;  break;}
             else imbalance1 = f;  // flavor-to-be of internal leg
         if (skip) continue;
         flavorParity = FermionParity(leg,mid1R,mid2L);
#if 0
         cout << "id(2): "; PrintVector(id);  cout << " " << mid1R << " " << mid2L << endl;
         cout << "fp2: "; PrintVector(flavorParity);  cout << endl;
         vector<int> fc = FermionCount(leg,mid1R,mid2L);
         cout << "fc2: "; PrintVector(fc);  cout << endl;
#endif
         for (int f = 0;  f < flavorParity.size();  f += 1)
           if (flavorParity[f])
             if (imbalance2) {skip = true;  break;}
             else imbalance2 = f;  // flavor-to-be of internal leg
         if (skip) continue;
         if (IsQuark(offshell))
            {if (imbalance1 isnt 0 and imbalance1 isnt flavorO
                 and imbalance1 isnt imbalance2) continue;
            if (imbalance1 is 0 and imbalance2 isnt 0 and
                imbalance2 isnt flavorO) continue;}
         else {if (imbalance1 isnt 0 and imbalance2 isnt 0
                   and imbalance1 isnt imbalance2) continue;}

         //            cout << "imb: " << imbalance1 << " " << imbalance2 << endl;
         ParticleID id1 = imbalance1 is 0 ?
           GluonID : FlavoredQuarkID(imbalance1);
         ParticleID id2 = imbalance2 is 0 ?
           GluonID : FlavoredQuarkID(imbalance2);
         ParticleID id3 =
           imbalance1 is 0 ?
           (imbalance2 is 0 ?
            offshell : (IsQuark(offshell) ? GluonID : id2))
           : (imbalance2 is 0 ?
              (IsQuark(offshell) ? GluonID : id1)
              : (imbalance1 is imbalance2 ? offshell : id2));
#endif
        // Sum over helicities
        for (int h1 = -1;  h1 <= 1;  h1 += 2)
        for (int h2 = -1;  h2 <= 1;  h2 += 2)
        for (int h3 = -1;  h3 <= 1;  h3 += 2)
           {complex<T> J1, J2, J3;
            // If not an allowed configuration, skip it
             if (h1+h2+h3+helicityO != 0 and massValue.size() is 0
                 and offshellMass < 0) continue;
            if (j1 is 1) // 2-pt current: do helicities match?
              if (helicity[start] isnt h1) continue; else J1 = externalLeg;
            else J1 = J(k,ref,NParticleID(-h1,id1),arg,leg,start,mid1L,
                        internal1,massValue);
            if (J1 is complex<T>(0,0)) continue;
            //            if (j1 is 1 and IsQuark(id1) and h1 < 0) J1 = -J1;
            if (j2 is j1+1) // 2-pt current: do helicities match?
              if (helicity[mid2L] isnt h2) continue; else J2 = externalLeg;
            else J2 = J(k,ref,NParticleID(-h2,id2),arg,leg,mid1R,mid2L,
                        internal2,massValue);
            if (J2 is complex<T>(0,0)) continue;
            //            if (j2 is j1+1 and IsQuark(id2) and h2 < 0) J2 = -J2;
            if (j2 is n-2) // 2-pt current: do helicities match?
              if (helicity[end] isnt h3) continue; else J3 = externalLeg;
            else J3 = J(k,ref,NParticleID(-h3,id3),arg,leg,mid2R,end,
                        internal3,massValue);
            if (J3 is complex<T>(0,0)) continue;
            //            if (j2 is n-2 and IsQuark(id3) and h3 < 0) J3 = -J3;
            int sum1 = MomentumSum(k,arg,start,mid1L);
            int sum2 = MomentumSum(k,arg,mid1R,mid2L);
            int sum3 = MomentumSum(k,arg,mid2R,end);
            complex<T> prop1 = k.m2(sum1);
            complex<T> prop2 = k.m2(sum2);
            complex<T> prop3 = k.m2(sum3);
            complex<T> props = complex<T>(1,0);
            if (j1 > 1)
               {if (IsScalar(id1) and internal1 >= 0)
                   prop1 -= k.m2(internal1);
               props *= prop1;}
            if (j2 > j1+1)
               {if (IsScalar(id2) and internal2 >= 0)
                 prop2 -= k.m2(internal2);
               props *= prop2;}
            if (n-1-j2 > 1)
               {if (IsScalar(id3) and internal3 >= 0)
                 prop3 -= k.m2(internal3);
               props *= prop3;}
#if 0
            cout << "h_in: ";
            print(helicity,start,end);
            cout << endl;
            cout << "h: " << helicityO << " " << h1 << " " << h2
                 << " " << h3 << endl;
            cout << "id: " << idO << " " <<id1 << " " << id2 << " " << id3 << endl;
            cout << "J1: " << J1 << "; J2: " << J2 << "; J3: " << J3;
            cout << "; V: " << Vertex(k,ref,helicityO,idO,arg,
                       start,mid1L,h1,id1,  mid1R,mid2L,h2,id2,
                                      mid2R,end,h3,id3);
            cout << endl;
#endif
            term = Vertex(k,ref,helicityO,idcodeO,arg,
                          start,mid1L,h1,ParticleCode(id1),
                          mid1R,mid2L,h2,ParticleCode(id2),
                       mid2R,end,h3,ParticleCode(id3)) * J1 * J2 * J3/props;
            // 5/20/08 New phase convention for spinors requires phase???
            term = -term;
            // 5/20/08 New phase convention for spinors requires additional sign too
#if 0
            if (IsQuark(id1) and j1 isnt 1) term = -term;
            if (IsQuark(id2) and j2 isnt j1+1) term = -term;
            if (IsQuark(id3) and j2 isnt n-2) term = -term;

            if (IsQuark(id1) and j1 is 1) term = -term;
            if (IsQuark(id2) and j2 is j1+1) term = -term;
            if (IsQuark(id3) and j2 is n-2) term = -term;
#endif
#if QuarkFlip
            if (IsQuark(id1)) term = -term;
            if (IsQuark(id2)) term = -term;
            if (IsQuark(id3)) term = -term;
#endif
#if 0
            if (n is 5)
            cout << "term 4: " << term << endl;
#endif
#define DebugSigns 0
#if DebugSigns
            if (n is 4) {
              cout << "term4 [";
              //DumpHelicities(helicityO,h1,h2,h3);
              DumpTaggedHelicities(idcodeO,helicityO,
                                   ParticleCode(id1),h1,
                                   ParticleCode(id2),h2,
                                   ParticleCode(id3),h3);
              cout << "; " << start << ".." << mid1L << ";" << mid1R << ".." << mid2L << ";" << mid2R << ".." << end;
            cout << "]: " << term << endl;
            }
#endif
            result += term;}}
#if QuarkFlip
    if (IsQuark(offshell)) result = -result;
#endif
    //    if (helicityO < 0 and IsQuark(idO)) result = -result;
    k.put_value(key,result);
#if 0
 if (n is 4 and helicity[0] is -1 and helicity[1] is 1 and helicity[2] is 1
     and helicityO is -1)
   cout << "LC: " << KnownGluonMPPM(k,ref,arg[0],arg[1],arg[2],arg[3]) << " " << result << endl;
#endif
 }
 //       cout << "---" << endl;
 return(result);
}

#if 0
template <class T> Cmom<T> operator*(complex<T> c,Cmom<T> p){
  if (c==0.) return Cmom<T>(0,0,0,0);
  if (c.imag()==0) return c.real()*p;
  return Cmom<T>(c*p.P(),sqrt(c)*p.L(),sqrt(c)*p.Lt());
};
#endif

// Wrapper routine
template<class T> complex<T> A0a(const process& p,
                                 momentum_configuration<R>& k,
                                 const vector<int>& ind)
{
}

//EXPLICIT INSTANTIATION

template class complex<R>
  J(momentum_configuration<R>& k,
    int ref0 /* index of reference momentum */,
    ParticleID offshell /* helicity, type etc */,
    const vector<int>& arg /* indices of arguments */,
    const vector<ParticleID>& leg,
    int start, int end /* indices into the vectors */,
    int offshellMass ,
    const vector<int>& massValue );

template class complex<RHP>
  J(momentum_configuration<RHP>& k,
    int ref0 /* index of reference momentum */,
    ParticleID offshell /* helicity, type etc */,
    const vector<int>& arg /* indices of arguments */,
    const vector<ParticleID>& leg,
    int start, int end /* indices into the vectors */,
    int offshellMass,
    const vector<int>& massValue);

}}

