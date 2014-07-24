#include <iostream>
#include <iomanip>

#include <stdio.h>
#include <complex>

#include "gmp_r.h"


#include "spinor.h"
#include "splib.hpp"
#include "mom_conf.h"
#include "mom_conf_inline.h"
#include "mom_conf.hpp"


namespace BH {

#define INSTANTIATE_SPLIB(TYPE) \
template class momentum<TYPE >;\
template class momentum<complex<TYPE > >;\
template class spinor<TYPE >;\
template class lambda<TYPE >;\
template class Cmom<TYPE >;\
template Cmom<TYPE > operator*(const TYPE  &,const Cmom<TYPE >&);\
template Cmom<TYPE > operator/(const Cmom<TYPE >&,const TYPE &);\
template Cmom<TYPE > operator/(const Cmom<TYPE >&,const complex<TYPE >&);\
template std::ostream& operator<<(std::ostream& s, const Cmom<TYPE >& p);\
template class smatrix<TYPE >;\
template std::ostream& operator<<(std::ostream& s, const momentum<TYPE >& p);\
template std::ostream& operator<<(std::ostream& s, const momentum<complex<TYPE > >& p);\
template std::ostream& operator<<(std::ostream& s, const smatrix<TYPE >& p);\
template std::ostream& operator<<(std::ostream& s, const spinor<TYPE >& p);\
template spinor<TYPE > la(const momentum<complex<TYPE > >&);\
template spinor<TYPE > lat<TYPE >(const momentum<complex<TYPE > >&);\
template complex<TYPE > spaa(const momentum<TYPE >&,const momentum<TYPE >&,const momentum<TYPE >&,const momentum<TYPE >&);\
template complex<TYPE > spaa(const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&);\
template complex<TYPE > spbb(const momentum<TYPE >&,const momentum<TYPE >&,const momentum<TYPE >&,const momentum<TYPE >&);\
template complex<TYPE > spbb(const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&,const momentum<complex<TYPE > >&);


INSTANTIATE_SPLIB(RGMP)
INSTANTIATE_SPLIB(RGMP_FP<200>)

template class momentum_configuration<RGMP_FP<200> >;
template class mom_conf_reader<RGMP_FP<200> >;

template class momentum_configuration<RGMP>;
template class sub_momentum_configuration<RGMP>;
template class mom_conf_reader<RGMP>;


}

