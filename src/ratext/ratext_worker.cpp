/*
 * ratext_worker.cpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */


#include "ratext/pentagon_ratext.hpp"
#include "ratext/box_ratext.hpp"
#include "ratext/triangle_ratext.hpp"
#include "ratext/bubble_ratext.hpp"
#include "ratext/rat_worker.h"
#include "rec_tree.h"


using namespace std;

namespace BH {

namespace ratext {

/*
 *
 *
 * Normal template definitions
 *
 *
 */

template <> struct daughter_info<box_Rat<rat_worker,Normal_RatBox_Specification<rat_worker> > > {
	typedef pentagon_Rat<rat_worker,Normal_RatPent_Specification<rat_worker> > daughter_type;
};
template <> struct daughter_info<triangle_Rat<rat_worker,Normal_RatTri_Specification<rat_worker> > > {
	typedef box_Rat<rat_worker,Normal_RatBox_Specification<rat_worker> > daughter_type;
};
template <> struct daughter_info<bubble_Rat<rat_worker,Normal_RatBub_Specification<rat_worker> > > {
	typedef triangle_Rat<rat_worker,Normal_RatTri_Specification<rat_worker> > daughter_type;
};


template class pentagon_Rat<rat_worker,Normal_RatPent_Specification<rat_worker> >;
template class box_Rat<rat_worker,Normal_RatBox_Specification<rat_worker> >;
template class triangle_Rat<rat_worker,Normal_RatTri_Specification<rat_worker> >;
template class bubble_Rat<rat_worker,Normal_RatBub_Specification<rat_worker> >;

//
// We need to explicitly instantiate each different triangle possibility
//

// For the standard renormalisable case we have

template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<R> const&, triangle_param<R,2,3>&);
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RHP> const&, triangle_param<RHP,2,3>&);
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RVHP> const&, triangle_param<RVHP,2,3>&);

template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<R> const&, triangle_param<R,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RHP> const&, triangle_param<RHP,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RVHP> const&, triangle_param<RVHP,2,3>&);

template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<R>&, std::vector<int> const&, triangle_param<R,2,3>&);
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RHP>&, std::vector<int> const&, triangle_param<RHP,2,3>&);
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RVHP>&, std::vector<int> const&, triangle_param<RVHP,2,3>&);

template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<R>&, std::vector<int> const&, triangle_param<R,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RHP>&, std::vector<int> const&, triangle_param<RHP,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RVHP>&, std::vector<int> const&, triangle_param<RVHP,2,3>&);

#if BH_USE_GMP
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RGMP> const&, triangle_param<RGMP,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RGMP> const&, triangle_param<RGMP,2,3>&);
template void triangle_Rat_plusminus<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RGMP>&, std::vector<int> const&, triangle_param<RGMP,2,3>&);
template void triangle_Rat_3mass<rat_worker,Normal_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RGMP>&, std::vector<int> const&, triangle_param<RGMP,2,3>&);
#endif

/*
 *
 *
 * Higgs template definitions
 *
 *
 */

template <> struct daughter_info<box_Rat<rat_worker,Higgs_RatBox_Specification<rat_worker> > > {
	typedef pentagon_Rat<rat_worker,Higgs_RatPent_Specification<rat_worker> > daughter_type;
};
template <> struct daughter_info<triangle_Rat<rat_worker,Higgs_RatTri_Specification<rat_worker> > > {
	typedef box_Rat<rat_worker,Higgs_RatBox_Specification<rat_worker> > daughter_type;
};
template <> struct daughter_info<bubble_Rat<rat_worker,Higgs_RatBub_Specification<rat_worker> > > {
	typedef triangle_Rat<rat_worker,Higgs_RatTri_Specification<rat_worker> > daughter_type;
};


template class pentagon_Rat<rat_worker,Higgs_RatPent_Specification<rat_worker> >;
template class box_Rat<rat_worker,Higgs_RatBox_Specification<rat_worker> >;
template class triangle_Rat<rat_worker,Higgs_RatTri_Specification<rat_worker> >;
template class bubble_Rat<rat_worker,Higgs_RatBub_Specification<rat_worker> >;

//
// We need to explicitly instantiate each different triangle possibility
//

// For the higgs case we have

template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<R> const&, triangle_param<R,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RHP> const&, triangle_param<RHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RVHP> const&, triangle_param<RVHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);

template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<R> const&, triangle_param<R,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RHP> const&, triangle_param<RHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RVHP> const&, triangle_param<RVHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);

template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<R>&, std::vector<int> const&, triangle_param<R,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RHP>&, std::vector<int> const&, triangle_param<RHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RVHP>&, std::vector<int> const&, triangle_param<RVHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);

template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<R>&, std::vector<int> const&, triangle_param<R,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RHP>&, std::vector<int> const&, triangle_param<RHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RVHP>&, std::vector<int> const&, triangle_param<RVHP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);

#if BH_USE_GMP
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RGMP> const&, triangle_param<RGMP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(eval_param<RGMP> const&, triangle_param<RGMP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_plusminus<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RGMP>&, std::vector<int> const&, triangle_param<RGMP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);
template void triangle_Rat_3mass<rat_worker,Higgs_RatTri_Specification<rat_worker> >::get_sub_terms_work(momentum_configuration<RGMP>&, std::vector<int> const&, triangle_param<RGMP,MUBUBPOINTS_HIGGS,YPOINTS_HIGGS>&);

#endif

}
}
