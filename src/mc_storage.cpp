#include "mc_storage.hpp"
#include "environment.hpp"
#include "BH_typedefs.h"
#include "spinor.h"

namespace BH {


namespace Tools {

int FSArray_settings::s_initialMaxChunkNumber=readFromEnvironment<int>("BH_FSARRAY_INITIAL_CHUNK_NUMBER",40);

template class FSArray<Cmom<R>,1000>;
template class FSArray<Cmom<RHP>,1000>;
template class FSArray<Cmom<RVHP>,1000>;

template class FSArray<std::complex<R>,1000>;
template class FSArray<std::complex<RHP>,1000>;
template class FSArray<std::complex<RVHP>,1000>;

template class FSArray<long volatile, 20>;

}

}

