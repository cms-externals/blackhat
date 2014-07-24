/*
 * standard_cut_part_ep.cpp
 *
 *  Created on: Jun 29, 2009
 *      Author: daniel
 */


#include "cut_part.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include<vector>
#include<complex>

#include "BH_typedefs.h"
#include "amplitudes.h"
#include "partitions.h"
#ifndef BH_PUBLIC
#include "options.h"
#endif
#include "integrals_ep.h"
#include "BH_A0.h"
#include "rational.h"
#include "cut_eval/amplitudes_cut_eval.h"
#include "process_utils.h"
#include "settings.h"
#include "cut_part_worker.h"
#ifndef BH_PUBLIC
#include "cut_part_normal.h"
#endif

using namespace std;


namespace BH {

namespace cut {


template <class CutDbase> vector<int> Indicesm1(CutDbase* cdb,int cor)
{vector<int> indices;
  for (size_t j = 1;  j <= cdb->corner_size(cor);  j += 1)
    indices.push_back(cdb->corner_ind(cor,j)-1);
 return(indices);
}

template <> vector<int> Indicesm1(worker::worker_cutD* cdb,int cor)
{
	return cdb->corner(cor);
}

}
}

#include "standard_cut_part_ep.hpp"



namespace BH {

namespace cut {


typedef standard_cut_part<worker::worker_boxDarren,worker::worker_triangleDarren,worker::worker_bubbleDarren> SCP_worker ;
#ifndef BH_PUBLIC
 typedef standard_cut_part<normal_boxDarren,normal_triangleDarren,normal_bubbleDarren> SCP_Darren ;
typedef standard_cut_part<boxD,triangleD,bubbleD> SCP_plain ;
typedef standard_cut_part<boxFHZ,triangleFHZ,bubbleFHZ> SCP_FHZ ;
typedef standard_cut_part<higgs_boxDarren,higgs_triangleDarren,higgs_bubbleDarren> SCP_higgs ;
#endif


template SeriesC<R> SCP_worker::eval_with_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_worker::eval_with_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_worker::eval_with_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_worker::eval_without_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_worker::eval_without_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_worker::eval_without_check(const eval_param<RVHP>& ep);
#ifndef BH_PUBLIC
template SeriesC<R> SCP_Darren::eval_with_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_Darren::eval_with_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_Darren::eval_with_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_Darren::eval_without_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_Darren::eval_without_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_Darren::eval_without_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_plain::eval_with_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_plain::eval_with_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_plain::eval_with_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_plain::eval_without_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_plain::eval_without_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_plain::eval_without_check(const eval_param<RVHP>& ep);



template SeriesC<R> SCP_FHZ::eval_with_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_FHZ::eval_with_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_FHZ::eval_with_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_FHZ::eval_without_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_FHZ::eval_without_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_FHZ::eval_without_check(const eval_param<RVHP>& ep);


template SeriesC<R> SCP_higgs::eval_with_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_higgs::eval_with_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_higgs::eval_with_check(const eval_param<RVHP>& ep);

template SeriesC<R> SCP_higgs::eval_without_check(const eval_param<R>& ep);
template SeriesC<RHP> SCP_higgs::eval_without_check(const eval_param<RHP>& ep);
template SeriesC<RVHP> SCP_higgs::eval_without_check(const eval_param<RVHP>& ep);
#endif

#if BH_USE_GMP

template SeriesC<RGMP> SCP_worker::eval_with_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_worker::eval_without_check(const eval_param<RGMP>& ep);
#ifndef BH_PUBLIC
template SeriesC<RGMP> SCP_Darren::eval_with_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_Darren::eval_without_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_plain::eval_with_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_plain::eval_without_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_FHZ::eval_with_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_FHZ::eval_without_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_higgs::eval_without_check(const eval_param<RGMP>& ep);
template SeriesC<RGMP> SCP_higgs::eval_with_check(const eval_param<RGMP>& ep);
#endif

#endif


}

}


