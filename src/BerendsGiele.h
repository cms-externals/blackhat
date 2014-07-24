#ifndef BRENDSGIELE_H_
#define BRENDSGIELE_H_

namespace BH {

namespace BerendsGiele {



template<class T> std::complex<T> J(momentum_configuration<T>& k,
    int ref0 /* index of reference momentum */,
    particle_ID offshell /* helicity, type etc */,
    std::vector<int> arg /* indices of arguments */,
    std::vector<particle_ID> leg,
    int start, int end /* indices into the vectors */);


}

using BerendsGiele::J;

}


#endif /*BRENDSGIELE_H_*/
