/*
 * bubble_specification.h
 *
 *  Created on: 12 Aug 2009
 *      Author: daniel
 */

#ifndef BUBBLE_SPECIFICATION_H_
#define BUBBLE_SPECIFICATION_H_

#include "momenta_evaluator.h"
#include "triangle_specification.h"

namespace BH {

namespace cut {
namespace Darren {
template <class T> class Darren_wrapper;
}
}

template <class cutDtype> struct DaughterType {
	typedef void type;
};

template <> struct DaughterType<cut::Darren::Darren_wrapper<BH::bubbleD> > {
	typedef cut::Darren::Darren_wrapper<BH::triangleD> type;
};

template <> struct DaughterType<BH::cut::worker::worker_cutD> {
	typedef BH::cut::worker::worker_cutD type;
};

template <class CutType> struct Normal_Bubble_Specification {
public:
	enum { CPOINTS=7 };
	enum { TPOINTSBUB=4 };
	enum { TPOINTSTRI=8 };
	enum { YPOINTS=2 };

	typedef cut::Darren::bubble_Darren_eval_points<TPOINTSBUB, YPOINTS> MomentaEvaluator;
	typedef cut::Darren::Normal_Corner_Tree_Strategy <MomentaEvaluator,CutType,CPOINTS>  CornerTreeStrategy;
	typedef cut::Darren::Normal_Bubble_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> BubbleCombiner;
	typedef Normal_Triangle_Specification<typename DaughterType<CutType>::type> TriangleSpecs;
};

template <class CutType,int CPOINTS_T, int TPOINTSBUB_T, int TPOINTSTRI_T, int YPOINTS_T> struct General_Bubble_Specification {
public:
	enum { CPOINTS = CPOINTS_T };
	enum { TPOINTSBUB = TPOINTSBUB_T };
	enum { TPOINTSTRI = TPOINTSTRI_T };
	enum { YPOINTS = YPOINTS_T };

	typedef cut::Darren::bubble_Darren_eval_points<TPOINTSBUB_T, YPOINTS_T> MomentaEvaluator;
	typedef cut::Darren::Normal_Corner_Tree_Strategy<MomentaEvaluator,CutType,CPOINTS_T> CornerTreeStrategy;
	typedef cut::Darren::General_Bubble_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> BubbleCombiner;
	typedef General_Triangle_Specification<typename DaughterType<CutType>::type,CPOINTS_T,TPOINTSBUB_T,TPOINTSTRI_T,YPOINTS_T> TriangleSpecs;
};

template <class CutType> struct Skeleton_Bubble_Specification {
public:
	enum { CPOINTS=7 };
	enum { TPOINTSBUB=4 };
	enum { TPOINTSTRI=8 };
	enum { YPOINTS=2 };

	typedef cut::Darren::bubble_Darren_eval_points<TPOINTSBUB, YPOINTS> MomentaEvaluator;
	typedef cut::Darren::Skeleton_Corner_Tree_Strategy <MomentaEvaluator,CutType,CPOINTS>  CornerTreeStrategy;
	typedef cut::Darren::Normal_Bubble_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> BubbleCombiner;
	typedef Normal_Triangle_Specification<CutType> TriangleSpecs;
};


template <class CutType> struct Normal_Higgs_Bubble_Specification : public  General_Bubble_Specification<CutType,9,5,9,4> {};
template <class CutType> struct Normal_test_Bubble_Specification : public  General_Bubble_Specification<CutType,7,6,10,5> {};
template <class CutType> struct Normal_3y_Bubble_Specification : public  General_Bubble_Specification<CutType,7,4,8,3> {};


}
#endif /* BUBBLE_SPECIFICATION_H_ */
