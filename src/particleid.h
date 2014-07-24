/* ParticleID.h */

/* David A. Kosower, March 23, 2008 */

/* An internal particle-ID framework; to be entirely replaced by BH's
   once particle codes are implemented */

/* Extracted from tree1.cc 3/28/08
*/

#ifndef ParticleCDDefined
typedef int ParticleCD;

namespace BH {
/* Model for encoding of particle IDs: bitwise-or of following tags:

     Boson 0x0
     Fermion 0x1
     Vector 0x2
     Colored 0x4

     Flavor 0x10 .. 0xF0

 */
const ParticleCD BosonTag = 0x0;
const ParticleCD FermionTag = 0x1;
const ParticleCD VectorTag = 0x2;
const ParticleCD ColoredTag = 0x4;

const ParticleCD Gluon = BosonTag|VectorTag|ColoredTag;
const ParticleCD Quark = FermionTag|ColoredTag;
// Intended for D-dimensional components of gluons
const ParticleCD Scalar = BosonTag|ColoredTag;
const ParticleCD MaskFlavor = 0xF;
const ParticleCD PlaceFlavor = 0xF0;

inline ParticleCD FlavoredQuark(int flavor) {
  return(Quark | (flavor<<4&PlaceFlavor));}

inline ParticleCD FlavoredScalar(int flavor) {
  return(Scalar | (flavor<<4&PlaceFlavor));}

inline bool IsGluon(ParticleCD id) {return(id == Gluon);}
inline bool IsQuark(ParticleCD id) {return((id & MaskFlavor) == Quark);}
inline bool IsScalar(ParticleCD id) {return((id & MaskFlavor) == Scalar);}
inline int FlavorOf(ParticleCD id) {return(id >> 4);}
inline bool IsFermion(ParticleCD id) {return(id & FermionTag);}

inline ParticleCD AddFlavor(ParticleCD id, int flavor)
  {return (id | (flavor << 4));}

// Suppress fermion flavors
#define VType2(id1, id2) \
 ( (id1&MaskFlavor) << 8 | (id2&MaskFlavor) )
#define VType3(id1, id2, id3) \
 ( (id1&MaskFlavor) << 16 | (id2&MaskFlavor) << 8 | (id3&MaskFlavor) )
#define VType4(id1, id2, id3, id4) \
 ( (id1&MaskFlavor) << 24 | (id2&MaskFlavor) << 16 | (id3&MaskFlavor) << 8 \
   | (id4&MaskFlavor) )
inline ParticleCD VertexType(ParticleCD id1, ParticleCD id2)
{return(VType2(id1,id2));}

inline ParticleCD VertexType(ParticleCD id1, ParticleCD id2, ParticleCD id3)
{return(VType3(id1,id2,id3));}

// Suppress fermion flavors
inline ParticleCD VertexType(ParticleCD id1, ParticleCD id2, ParticleCD id3,
       ParticleCD id4)
{return(VType4(id1,id2,id3,id4));}


inline std::string TypeName(ParticleCD id) {
  if (IsQuark(id)) return("q"); else if (id == Gluon) return("g");
  return("?");
}

/* 7/17/08: Need to encode -3 & +3 (for 'L' and 'R' 'helicities',
   see cubic.tex) as well as -1 & + 1 */
const int MaskH = 0xFF;
#define HType2(h1,h2) \
  ( (h1&MaskH) << 8 | (h2&MaskH) )

#define HType3(h1,h2,h3) \
  ( (h1&MaskH) << 16 | (h2&MaskH) << 8 | (h3&MaskH) )

#define HType4(h1,h2,h3,h4) \
  ( (h1&MaskH) << 24 | (h2&MaskH) << 16 | (h3&MaskH) << 8 | (h4&MaskH) )

inline int HelicityType(int h1, int h2)
{return(HType2(h1,h2));}

inline int HelicityType(int h1, int h2, int h3)
{return(HType3(h1,h2,h3));}

inline int HelicityType(int h1, int h2, int h3, int h4)
{return(HType4(h1,h2,h3,h4));}

typedef particle_ID ParticleID;
#undef IsGluon
inline bool IsGluon(ParticleID id) {return id.is_a(gluon);}
#undef IsQuark
inline bool IsQuark(ParticleID id) {return id.is_a(quark);}
#undef FlavorOf
#define FlavorOf(id) id.flavor()
#undef HelicityOf
#define HelicityOf(id) id.helicity()

#undef IsScalar
inline bool IsScalar(ParticleID id) {return id.is_a(scalar_massive);}

inline bool IsFermion(ParticleID id) {return id.is_a(quark) or
                                        id.is_a(lepton);}
}

#define ParticleCDDefined 1
#endif /* ParticleCDDefined */


