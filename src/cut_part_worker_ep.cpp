
#include "cut_part.h"

#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include<vector>
#include<complex>

#include "BH_typedefs.h"
#include "amplitudes.h"
#include "integrals_ep.h"
#include "BH_A0.h"
#include "rational.h"
#include "cut_eval/amplitudes_cut_eval.h"
#include "process_utils.h"
#include "settings.h"
#include "cut_part_worker.h"
#include "standard_cut_part.h"

#define _VERBOSE 0 // 0 is silent, 1 gives some information

using namespace std;

namespace BH {

namespace cut {



namespace worker {




SeriesC<R> worker_cut_part::eval(const eval_param<R>& ep){
    return eval_fn(ep);
}
SeriesC<RHP> worker_cut_part::eval(const eval_param<RHP>& ep){
    return eval_fn(ep);
}
SeriesC<RVHP> worker_cut_part::eval(const eval_param<RVHP>& ep){
    return eval_fn(ep);
}

}
}

}
