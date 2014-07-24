#ifndef BH_A0_H_
#define BH_A0_H_

#include "BH_typedefs.h"
#include <vector>
//!Tree level amplitude
/*!
     works so far for
     	3-6 gluons,
     	q-qb +2-4 gluons.

      \param p is the process
      \param ps is the momentum configuration of type mom_conf
      \param ind is a vector containing the integer indices of the momenta in the momentum configuration ps that are to be used to evaluate the amplitude.
      \return the complex value of the amplitude
      \sa A1
    */
namespace BH {

template <class T> class momentum_configuration;
class  process;

C A0(const process&,momentum_configuration<R>&,const std::vector<int>&);
C A0(long pc,const process& p,momentum_configuration<R>&,const std::vector<int>&);
CHP A0(const process&,momentum_configuration<RHP>&,const std::vector<int>&);
CHP A0(long,const process& p,momentum_configuration<RHP>&,const std::vector<int>&);
CVHP A0(const process&,momentum_configuration<RVHP>&,const std::vector<int>&);
CVHP A0(long,const process& p,momentum_configuration<RVHP>&,const std::vector<int>&);
template <class T> std::complex<T> A0_safe(const process&,momentum_configuration<T>&,const std::vector<int>&);
template <class T> std::complex<T> A0_safe(long,const process&,momentum_configuration<T>&,const std::vector<int>&);
C A1(const process&,momentum_configuration<R>&,const std::vector<int>&);
C A0(const process& p,momentum_configuration<R>& ps,int i1);
C A0(const process& p,momentum_configuration<R>& ps,int i1,int i2);
C A0(const process& p,momentum_configuration<R>& ps,int i1,int i2,int i3);
C A0(const process& p,momentum_configuration<R>& ps,int i1,int i2,int i3,int i4);
C A0(const process& p,momentum_configuration<R>& ps,int i1,int i2,int i3,int i4,int i5);
C A0(const process& p,momentum_configuration<R>& ps,int i1,int i2,int i3,int i4,int i5,int i6);
}

#endif /*BH_A0_H_*/
