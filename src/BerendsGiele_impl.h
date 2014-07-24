#ifndef BERENDSGIELE_IMPL_H_
#define BERENDSGIELE_IMPL_H_

#include "particleid.h"

namespace BH {

namespace BerendsGiele {



inline std::vector<int> MakeVector(int i1,int i2)
  {std::vector<int> v(2);  v[0] = i1;  v[1] = i2;  return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3)
  {std::vector<int> v(3);  v[0] = i1;  v[1] = i2;  v[2] = i3;  return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4)
  {std::vector<int> v(4);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5)
  {std::vector<int> v(5);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6)
  {std::vector<int> v(6);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7)
  {std::vector<int> v(7);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7,
                              int i8)
  {std::vector<int> v(8);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;  v[7] = i8;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7,
                              int i8,int i9)
  {std::vector<int> v(9);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;  v[7] = i8;  v[8] = i9;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7,
                              int i8,int i9,int i10)
  {std::vector<int> v(10);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;  v[7] = i8;  v[8] = i9;  v[9] = i10;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7,
                              int i8,int i9,int i10,int i11)
  {std::vector<int> v(11);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;  v[7] = i8;  v[8] = i9;  v[9] = i10;
   v[10] = i11;
   return v;}

inline std::vector<int> MakeVector(int i1,int i2,int i3,int i4,int i5,int i6,int i7,
                              int i8,int i9,int i10,int i11,int i12)
  {std::vector<int> v(12);  v[0] = i1;  v[1] = i2;  v[2] = i3;  v[3] = i4;
   v[4] = i5;  v[5] = i6; v[6] = i7;  v[7] = i8;  v[8] = i9;  v[9] = i10;
   v[10] = i11; v[11] = i12;
   return v;}


particle_ID NParticleID(int helicity,const particle_ID& base);

std::vector<particle_ID> NParticleID(const std::vector<int>& helicity,
                               const std::vector<particle_ID>& base);

const particle_ID GluonID(gluon,1,0,false);
const C I(0,1);

particle_ID FlavoredQuarkID(int flavor);

}
}
#endif /*BERENDSGIELE_IMPL_H_*/
