#ifndef TREE_H_
#define TREE_H_

#include <vector>


namespace BH {

namespace Tree {

const std::vector<int> empty(0);
const std::vector<particle_ID> emptyID(0);

const int defaultAlgorithm = 0;

const int defaultMass = -1;

template<class T> std::complex<T>
  Aosrr(momentum_configuration<T>& k,
        const std::vector<int>& arg /* indices of arguments */,
        const std::vector<particle_ID>& leg /* helicities and particle ids */,
        const std::vector<int>& vectorK /* momenta */ = empty,
        const std::vector<int>& polarization = empty,
        const std::vector<int>& coupleTo /* quark flavor */ = empty,
        const std::vector<int>& massValue = empty);

template<class T> std::complex<T> A(momentum_configuration<T>& k,
                               const std::vector<int>& arg0 /* indices of arguments */,
                               const std::vector<particle_ID>& leg,
                               int algorithm = defaultAlgorithm,
                               const std::vector<int>& vectorK /* momenta */ = empty,
                               const std::vector<int>& polarization = empty,
                               const std::vector<int>& coupleTo /* quark flavor */ = empty,
                               const std::vector<int>& massValue = empty);

template<class T> std::complex<T> A(momentum_configuration<T>& k,
                               const std::vector<int>& arg0 /* indices of arguments */,
                               const process& p,
                               int algorithm = defaultAlgorithm,
                               //                               const std::vector<int>& vectorK /* momenta */ = empty,
                               //                               const std::vector<int>& polarization = empty,
                               //                               const std::vector<int>& coupleTo /* quark flavor */ = empty,
                               const std::vector<int>& massValue = empty);

template<class T> std::complex<T>
  J(momentum_configuration<T>& k,
    int ref0 /* index of reference momentum */,
    particle_ID offshell /* helicity, type etc */,
    const std::vector<int>& arg /* indices of arguments */,
    const std::vector<particle_ID>& leg,
    int start, int end /* indices into the vectors */,
    int offshellMass = defaultMass,
    const std::vector<int>& massValue = empty);

template<class T> std::complex<T>
  Jc(momentum_configuration<T>& k,
     int ref0 /* index of reference momentum */,
     particle_ID offshell /* helicity, type etc */,
     const std::vector<int>& arg /* indices of arguments */,
     const std::vector<particle_ID>& leg,
     int start, int end /* indices into the vectors */,
     int flow,
     const std::vector<int>& vectorK /* momenta */ = empty,
     const std::vector<int>& polarization = empty,
     const std::vector<int>& coupleTo /* quark flavor */ = empty,
     int offshellMass = defaultMass,
     const std::vector<int>& massValue = empty);


}

using Tree::J;

}


#endif /*TREE_H_*/
