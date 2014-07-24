/* Tree1.cc */

/*  David A. Kosower,  September 14, 2007  */

/* Implementation of light-cone gauge Berends-Giele recursion relations
   using BlackHat momentum library.

   3/22/08: Updated for new version of BlackHat

   5/22/08: Reference momentum is no longer an argument to J; it is cached,
            and created if no "ref" label exists in the momentum
            configuration.

*/

#include "mom_conf.h"
#include "BerendsGiele.h"
#include "BerendsGiele_impl.h"
#include "polylog.h"
#include <math.h>
#include <ctime>

#define Re(v) real(v)
#define Im(v) imag(v)

#define RefTag "ref"


#define BaseIndex 1 // [0] index not used in arg, helicity, or id.

#define isnt !=
#define is ==

#define FIXZERO 0  // Force -0 to 0 in momenta, needed for old spinor phase convention

#define ParticleID particle_ID
#undef IsGluon
#define IsGluon(id) id.is_a(gluon)
#undef IsQuark
#define IsQuark(id) id.is_a(quark)
#undef FlavorOf
#define FlavorOf(id) id.flavor()
#undef HelicityOf
#define HelicityOf(id) id.helicity()

using namespace std;

namespace BH {
namespace BerendsGiele {
#undef MomentumConfiguration
#define MomentumConfiguration momentum_configuration<R>

template<class T> inline T Max(T x,T y) {return (x > y ? x : y);}
inline bool IsOdd(int i) {return(i&1);}

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

#define GGGG VType4(Gluon,Gluon,Gluon,Gluon)
#define FFGG VType4(Quark,Quark,Gluon,Gluon)
#define FGFG VType4(Quark,Gluon,Quark,Gluon)
#define FGGF VType4(Quark,Gluon,Gluon,Quark)
#define GFGF VType4(Gluon,Quark,Gluon,Quark)
#define GGFF VType4(Gluon,Gluon,Quark,Quark)
#define GFFG VType4(Gluon,Quark,Quark,Gluon)
#define FFFF VType4(Quark,Quark,Quark,Quark)


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

// Computes the sum of the momenta v[start..end], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
int MomentumSum(MomentumConfiguration& k,
                const vector<int>& v, int start, int end)
{if (start is end) return(start);
 string key = GenKey("ms",start,end,v);
// cout << "In MS " << start << ", " << end << endl;
 size_t index;
 if (not k.get_label(key,index)) {
    Cmom<R> sum;  // Initialized to 0
    if (start <= end)
       for (int j = start;  j <= end;  j += 1)
         sum += k[v[j]];
    else {for (int j = start;  j < v.size();  j += 1)  sum += k[v[j]];
       for (int j = BaseIndex;  j <= end;  j += 1)
         sum += k[v[j]];}
    FixZero(sum);
    index = k.insert(sum);
    k.put_label(key,index);}
 return index;
}


// Computes the flatted sum of the momenta v[start..end], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
int FlatSum(MomentumConfiguration& k,
            int ref /* index of reference momentum */,
            const vector<int>& v, int start, int end)
{// cout << "In FS " << start << ", " << end << endl;
 int sum = MomentumSum(k,v,start,end);
 string key = GenKey("fs",start,end,ref,v);
 size_t index;
 if (not k.get_label(key,index)) {
    momC flat = k[sum] - (k.m2(sum)/(2.*(k[sum]*k[ref]))) * k[ref];
    FixZero(flat);
    index = k.insert(flat);
    k.put_label(key,index);}
 return index;
}

// Computes the flatted sum of the negative of the
//   momenta v[start1..end1], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
int NegativeFlatSum(MomentumConfiguration& k,
              int ref /* index of reference momentum */,
        const vector<int>& v, int s1, int e1)
{// cout << "In NFS " << s1 << ", " << e1 << "; " << s2 << ", " << e2 << endl;
 int sum = MomentumSum(k,v,s1,e1);

 string key = GenKey("nf",s1,e1,ref,v);
 size_t index;
 if (not k.get_label(key,index)) {
   // The spinor products in BlackHat alas have their branch cuts along the
   // real axis, which means we have to be VERY careful to ensure that
   // Negative(FlatSum(...)) is identical -- down to the sign of
   // (infinitesimal) imaginary parts to NegativeFlatSum(...) -- hence
   // do this in two steps
    momC flat = k[sum] - (k.m2(sum)/(2.*(k[sum]*k[ref]))) * k[ref];
    flat = -flat;
    FixZero(flat);
    index = k.insert(flat);
    k.put_label(key,index);}
 return index;
}

// Computes the flatted sum of the negative of the
//   momenta v[start1..end1] + v[start2..end2], if not already
// known, and creates a new momentum entry corresponding to it.
// "start" and "end" are understood in a cyclic sense
int NegativeFlatSum(MomentumConfiguration& k,
              int ref /* index of reference momentum */,
        const vector<int>& v, int s1, int e1, int s2, int e2)
{// cout << "In NFS " << s1 << ", " << e1 << "; " << s2 << ", " << e2 << endl;
 int sum1 = MomentumSum(k,v,s1,e1);
 int sum2 = MomentumSum(k,v,s2,e2);

 // string key = GenKey("nf",s1,e1,s2,e2,ref,v);
 // until appropriate routine is available:
 string key = GenKey("fs",MakeVector(s1,e1,s2,e2,ref),v);
 size_t index;
 if (not k.get_label(key,index)) {
    momC ksum = k[sum1]+k[sum2];
    int sum = k.insert(ksum);
    momC flat = ksum - (k.m2(sum)/(2.*(k[sum]*k[ref]))) * k[ref];
    flat = -flat;
    //    cout << "NFS 1: " << flat << endl;
    FixZero(flat);
#if 0
    if (Im(flat.E()) == 0) {cout << "zero" << endl;
    flat = Cmom<R>(C(Re(flat.E()),0),C(Re(flat.X()),0),C(Re(flat.Y()),0),
                C(Re(flat.Z()),0));}
    //    if (Im(flat.E()) == 0) flat.E() *= copysign(1.0,flat.E());
#endif
    //    cout << "NFS 2: " << flat << endl;
    index = k.insert(flat);
    k.put_label(key,index);}
 return index;
}

int NegativeFlatSum(MomentumConfiguration& k,
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
 string key = GenKey("fs",MakeVector(s1,e1,s2,e2,s3,e3,ref),v);
 size_t index;
 if (not k.get_label(key,index)) {
    momC ksum = k[sum1]+k[sum2]+k[sum3];
    int sum = k.insert(ksum);
    momC flat = ksum - (k.m2(sum)/(2.*(k[sum]*k[ref]))) * k[ref];
    flat = -flat;
    FixZero(flat);
    index = k.insert(flat);
    k.put_label(key,index);}
 return index;
}

int Negative(MomentumConfiguration& k, int i)
{string key = GenKey("neg",i);
 size_t index;
 if (not k.get_label(key,index)) {
    momC kneg = -k[i];
    //    cout << "N 1: " << kneg << endl;
    FixZero(kneg);
    //    cout << "N 2: " << kneg << endl;
    index = k.insert(kneg);
    k.put_label(key,index);}
 return index;
}

// Vertices for one off-shell leg (with no explicit momentum argument, its
// momentum is given by momentum conservation), and some number of
// arguments with given indices, helicities, and particle IDs



#define VgggMPP(k1,k2,k3) \
          (I/2.) * k.spa(k1,ref) * k.spb(k2,k3) \
             * (k.s(k1,ref)-k.s(k2,ref)-k.s(k3,ref))/ \
          (k.spa(k2,ref)*k.spa(k3,ref)*k.spb(k1,ref))

#define VgggPMM(k1,k2,k3) \
          (-I/2.) * k.spb(k1,ref) * k.spa(k2,k3) \
              * (k.s(k1,ref)-k.s(k2,ref)-k.s(k3,ref))/ \
          (k.spb(k2,ref) * k.spb(k3,ref) * k.spa(k1,ref))

// Light-cone three-gluon vertex
// V_3(-K_{s1..e1}-K_{s2..e2},K_{s1..e1},K_{s2..e2})
// Lightcone means that it's the tree-gluon vertex with flatted momenta
inline C Vggg(MomentumConfiguration& k,
              int ref /* index of reference momentum */,
              int helicity0,
              const vector<int>& arg /* indices of momenta from which legs are taken */,
              int s1, int e1, int helicity1,
              int s2, int e2, int helicity2)
{//string key = GenKey("Vggg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vggg",MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                                               ref),arg);
 C result;

 if ( !k.get_value(key,result) ) {
    int k0 = NegativeFlatSum(k,ref,arg,s1,e1,s2,e2);
    int k1 = FlatSum(k,ref,arg,s1,e1);
    int k2 = FlatSum(k,ref,arg,s2,e2);
    switch (HelicityType(helicity0,helicity1,helicity2)) {
    case MPP:
      result = VgggMPP(k0,k1,k2);  break;
    case PMP:
      result = VgggMPP(k1,k2,k0);  break;
    case PPM:
      result = VgggMPP(k2,k0,k1);  break;
    case PMM:
      result = VgggPMM(k0,k1,k2);  break;
    case MPM:
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
      result = VgggPMM(k2,k0,k1);  break;
    case PPP:
    case MMM:
      result = 0;  break;
    default:
      throw "Illegal helicity configuration [Vggg]";
    }
    k.put_value(key,result);}
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
#if 0
inline
#endif
inline C Vffg(MomentumConfiguration& k,
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
 C result;

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
      result = 0;  break;
    default:
      throw "Illegal helicity configuration [Vffg]";
    }

    // 5/20/08 New phase convention for spinors requires additional phase
    // here in order to get conventional phase for two-quark amplitudes
    result = -result;
    k.put_value(key,result);}
 return result;
}

static inline C square(C x) {return(x*x);}

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
inline C Vgggg(MomentumConfiguration& k,
              int ref /* index of reference momentum */,
              int helicity0,
              const vector<int>& arg /* indices of momenta from which legs
                                        are taken */,
              int s1, int e1, int helicity1,
              int s2, int e2, int helicity2,
              int s3, int e3, int helicity3)
{//string key = GenKey("Vggg",helicity0,s1,e1,helicity1,s2,e2,helicity2,ref,arg);

 // until appropriate routine is available:
 string key = GenKey("Vgggg",
          MakeVector(helicity0,s1,e1,helicity1,s2,e2,helicity2,
                               s3,e3,helicity3,ref),arg);
 C result;

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
      result = 0;  break;
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
inline C Vffgg(MomentumConfiguration& k,
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
 C result;
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
      result = 0;  break;
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
inline C Vfgfg(MomentumConfiguration& k,
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
 C result;

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
      result = 0;  break;
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
  2. * I * k.spab(k1,ref,k2n) * k.spab(k3,ref,k4n)/ \
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
inline C Vffhh(MomentumConfiguration& k,
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
 C result;

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
      result = 0;  break;
    default:
      throw "Illegal helicity configuration [Vfgfg]";
    }
    // Extra phase seems necessary -- why??
    result *= -I;
    k.put_value(key,result);}
 return result;
}

// All vertices: three-point
static C Vertex(MomentumConfiguration& k, int ref /* index of reference momentum */,
         int helicity0,
         int id0,
         const vector<int>& arg,
         int s1, int e1, int helicity1, int id1,
         int s2, int e2, int helicity2, int id2)
{
#if 0
cout << hex << "VT: " << TypeName(id0) << TypeName(id1) << TypeName(id2)
      << VertexType(id0,id1,id2)
      << " [" << s1 <<":"<<e1<<","<<s2<<":"<<e2<<"]" << endl;
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
 default:
   throw "Illegal vertex type [Vertex]";
}
}

// And four-point (includes contributions from shuffling around light-cone
// terms)
static C Vertex(MomentumConfiguration& k, int ref /* index of reference momentum */,
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

/* Counts the number of fermion legs of each different flavor in the
   given range of arguments within "id".  The returned vector has length
   equal to the highest flavor number in "id" plus one (the 0 component
   is not used)
*/
vector<int> FermionCount(const vector<ParticleID>& id, int start, int end)
{int max = 0;
 for (int j = BaseIndex;  j < id.size();  j += 1)
    if (IsQuark(id[j])) max = Max(max,(int)FlavorOf(id[j]));
 vector<int> count(max+1,0); // Initialized to 0
 if (start <= end)
   {for (int j = start;  j <= end;  j += 1)
     if (IsQuark(id[j])) count[FlavorOf(id[j])] += 1;}
 else {for (int j = start;  j < id.size();  j += 1)
      if (IsQuark(id[j])) count[FlavorOf(id[j])] += 1;
    for (int j = BaseIndex;  j <= end;  j += 1)
      if (IsQuark(id[j])) count[FlavorOf(id[j])] += 1;}

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

void print(const vector<int>& v,int s, int e) {
  for (int i = s;  i <= e; i += 1)
    cout << v[i] << " ";
}

//extern "C" void fillrandom01(double *x,int n);

// Generate null momentum
//static momC GenerateMomentum() {
//  R components[4];
//  components[0] = 0;
//  fillrandom01(components+1,3);
//  momC bare(components[0],components[1],components[2],components[3]);
//  return momC(sqrt(-bare*bare),bare.X(),bare.Y(),bare.Z());
//}

//-----------------------------------------
// DM: alternative function that does not require the C implemetation

inline R randomR(){
	srand(time(NULL));
	return double(rand())/double(RAND_MAX);
}

static momC GenerateMomentum() {
	  momentum<R> bare(0.,randomR(),randomR(),randomR());

	return momC(sqrt(-bare*bare),bare.X(),bare.Y(),bare.Z());
}
//-----------------------------------------


// Helicity ignored here
//const ParticleID GluonID(gluon,1,0,false);
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
   the fermion amplitudes seem to come out with a minus sign compared
   to the "canonical" (susy) phase conventions.
*/
template<class T> complex<T> J(momentum_configuration<T>& k,
    int ref0 /* Supply negative number to use default one */,
    ParticleID offshell /* helicity, type etc */,
#if 0
    int helicityO /* helicity of offshell leg */,
    int idO /* ParticleID of offshell leg */,
#endif
    vector<int> arg /* momentum labels of arguments */,
#if 0
    vector<int> helicity,
    vector<ParticleID> id,
#endif
    vector<ParticleID> leg,
    int start, int end /* indices into the vectors */)
{ // Get reference momentum
 size_t ref /* index of reference momentum */;
 const string refKey(RefTag);
 if (ref0 >= 0) ref = ref0;
 else if (not k.get_label(refKey,ref)) {
   momC refMom = GenerateMomentum();
   ref = k.insert(refMom);
   k.put_label(refKey,ref);
 }

 vector<int> helicity = Helicities(leg);
 vector<int> idcode = ParticleCode(leg);

   int helicityO = offshell.helicity();
   int idcodeO = ParticleCode(offshell);

 //string key = GenKey("J",helicityO,idO,start,end,arg,helicity,id);
 // until the right routine is available:
 string key = GenKey("J",start,end,MakeVector(ref,helicityO,idcodeO),
                         arg,helicity,idcode);

 C result;

 if (not k.get_value(key,result)) {
    // Must compute it
    C term;

    // Number of arguments including off-shell leg
    int n = CountCyclic(arg,start,end)+1;

    //       cout << "---" << endl;
    //    cout << "J_" << n <<": " << idO << "; ";
    // PrintVector(id); cout << " [" << start << ":" << end << "] " << endl;

    // First look at special cases:
    // End-points of recursion: three-point currents are just vertices

    if (n == 3)
       {result = Vertex(k,ref,helicityO,idcodeO,arg,
                                      start,start,helicity[start],idcode[start],
                                      end,end,helicity[end],idcode[end]);
       //       cout << "--- " << result << endl;
        k.put_value(key,result);
        return result;}

    // All-plus vanishes, likewise all-minus
    int netHelicity = helicityO+SumCyclic(helicity,start,end);

    if (netHelicity == (n=CountCyclic(arg,start,end)+1)
        || netHelicity == -n)
       {result = 0;  k.put_value(key,result);  return result;}
    // Small number of negative helicities, expand using MHV vertices
    // or perhaps ordinary on-shell recursion
    //else {}

    // General case: sum over currents joined to three-point and
    // four-point vertices

    // Scope out fermion flavors; only one flavor can be imbalanced,
    // and it must match that of the off-shell leg
    vector<bool> flavorParity = FermionParity(leg,start,end);
    bool alreadyImbalanced = false;
    int flavorO = FlavorOf(offshell);
    for (int f = 0;  f < flavorParity.size();  f += 1)
       if (flavorParity[f])
          if (alreadyImbalanced or f isnt flavorO) return(0);
          else alreadyImbalanced = true;

    // Just gluons for now -- no checking of IDs
    result = 0;

    // 5/20/08 New phase convention for spinors requires different phase
    // for two-particle "current"
    C externalLeg = I;

    // Three-point vertex terms
    for (int j1 = 1;  j1 <= n-2;  j1 += 1)
       {// Sum over helicities: note that h_i is the helicity on the vertex,
        // the negative of the helicity on the current
        for (int h1 = -1;  h1 <= 1;  h1 += 2)
        for (int h2 = -1;  h2 <= 1;  h2 += 2)
           {C J1, J2;
            int midL, midR; // Indices of j1th & j1+1st legs "arg"
            // If not an allowed configuration, skip it
            if (h1+h2+helicityO == 3 or h1+h2+helicityO == -3) continue;
            /* For allowed configurations, vertex is always non-vanishing;
               we may want to compute shorter current first, but just do
               them in order for now. */
            midL=IndexCyclic(arg,start+j1-1);
            /* Particle type of offshell leg of first daughter current should
               be the lone imbalanced fermion flavor (argument) or gluon if its
               arguments are fermion-number balanced; as the overall imbalance
               has already been checked above, this ensures correct type
               for other daughter current */
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
#if 0
            cout << "n: " << n << " (" << j1 << "," << midL << ")" << endl;
            cout << "id: " << idO << " " << id1 << " " << id2 << endl;
            cout << IsQuark(idO) << "; " << imbalance << endl;
#endif
            if (j1 is 1) // 2-pt current: do helicities match?
               if (helicity[start] isnt h1) continue; else J1 = externalLeg;
               else J1 = J(k,ref,NParticleID(-h1,id1),arg,leg,start,midL);
            if (J1 is 0.) continue;
            // Fermion phase convention: extra factor of -i for +;- configuration
            // Is this phase independent of the sign of the energy?
            //            if (j1 is 1 and IsQuark(id1) and h1 < 0) J1 = -J1;
            midR=IndexCyclic(arg,start+j1);
            if (j1 is n-2) // 2-pt current: do helicities match?
               if (helicity[end] isnt h2) continue; else J2 = externalLeg;
            else J2 = J(k,ref,NParticleID(-h2,id2),arg,leg,midR,end);
            if (J2 is 0.) continue;
            //            if (j1 is n-2 and IsQuark(id2) and h2 < 0) J2 = -J2;
            int sum1 = MomentumSum(k,arg,start,midL);
            int sum2 = MomentumSum(k,arg,midR,end);
            C prop1 = k.m2(sum1);
            C prop2 = k.m2(sum2);
            C props;
            if (j1 is 1) // Because n>3 here, j1 cannot be both 1 and n-2
               props = prop2;
            else if (j1 is n-2)
               props = prop1;
            else props = prop1 * prop2;
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
            if (IsQuark(id1)) term = -term;
            if (IsQuark(id2)) term = -term;
#if 0
            cout << "n: " << n << endl;
            cout << "h: " << helicityO << " " << h1 << " " << h2 << endl;
            cout << "id: " << idO << " " << id1 << " " << id2 << endl;
#endif
#if 0
            cout << "term: " << term << endl;
#endif
            result -= term;}}

    // Four-point vertex terms
    for (int j1 = 1;  j1 <= n-3;  j1 += 1)
    for (int j2 = j1+1;  j2 <= n-2;  j2 += 1)
       {// Sum over helicities
        for (int h1 = -1;  h1 <= 1;  h1 += 2)
        for (int h2 = -1;  h2 <= 1;  h2 += 2)
        for (int h3 = -1;  h3 <= 1;  h3 += 2)
           {C J1, J2, J3;
            int mid1L, mid1R, mid2L, mid2R;
            // If not an allowed configuration, skip it
            if (h1+h2+h3+helicityO != 0) continue;
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
            if (j1 is 1) // 2-pt current: do helicities match?
              if (helicity[start] isnt h1) continue; else J1 = externalLeg;
            else J1 = J(k,ref,NParticleID(-h1,id1),arg,leg,start,mid1L);
            if (J1 is 0.) continue;
            //            if (j1 is 1 and IsQuark(id1) and h1 < 0) J1 = -J1;
            if (j2 is j1+1) // 2-pt current: do helicities match?
              if (helicity[mid2L] isnt h2) continue; else J2 = externalLeg;
            else J2 = J(k,ref,NParticleID(-h2,id2),arg,leg,mid1R,mid2L);
            if (J2 is 0.) continue;
            //            if (j2 is j1+1 and IsQuark(id2) and h2 < 0) J2 = -J2;
            if (j2 is n-2) // 2-pt current: do helicities match?
              if (helicity[end] isnt h3) continue; else J3 = externalLeg;
            else J3 = J(k,ref,NParticleID(-h3,id3),arg,leg,mid2R,end);
            if (J3 is 0.) continue;
            //            if (j2 is n-2 and IsQuark(id3) and h3 < 0) J3 = -J3;
            int sum1 = MomentumSum(k,arg,start,mid1L);
            int sum2 = MomentumSum(k,arg,mid1R,mid2L);
            int sum3 = MomentumSum(k,arg,mid2R,end);
            C prop1 = k.m2(sum1);
            C prop2 = k.m2(sum2);
            C prop3 = k.m2(sum3);
            C props = 1;
            if (j1 > 1) props *= prop1;
            if (j2 > j1+1) props *= prop2;
            if (n-1-j2 > 1) props *= prop3;
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
            if (IsQuark(id1)) term = -term;
            if (IsQuark(id2)) term = -term;
            if (IsQuark(id3)) term = -term;
#if 0
            cout << "term 4: " << term << endl;
#endif
            result += term;}}
    if (IsQuark(offshell)) result = -result;
    //    if (helicityO < 0 and IsQuark(idO)) result = -result;
    k.put_value(key,result);
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

//EXPLICIT INSTANTIATION

template complex<R> J(momentum_configuration<R>& k,
    int ref0 /* index of reference momentum */,
    ParticleID offshell /* helicity, type etc */,
    vector<int> arg /* indices of arguments */,
    vector<ParticleID> leg,
    int start, int end /* indices into the vectors */);
}
}

//=============================================================================================================================================================

