/*
 * triangle_specification.h
 *
 *  Created on: 9 Nov 2009
 *      Author: daniel
 */

#ifndef TRIANGLE_SPECIFICATION_H_
#define TRIANGLE_SPECIFICATION_H_


#include "momenta_evaluator.h"
#include "triangle_subtraction.h"

namespace BH {


template <class CutType> struct Normal_Triangle_Specification {
public:
	enum { CPOINTS=7 };
	enum { TPOINTSBUB=4 };
	enum { TPOINTSTRI=8 };
	enum { YPOINTS=2 };

	typedef cut::Darren::triangle_Darren_eval_points<CPOINTS, TPOINTSTRI> MomentaEvaluator;
	typedef cut::Darren::Normal_Corner_Tree_Strategy <MomentaEvaluator,CutType,CPOINTS>  CornerTreeStrategy;
//	typedef cut::Darren::Normal_Triangle_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> TriangleCombiner;
	typedef cut::Darren::Normal_Triangle_Subtraction<CutType,CPOINTS,TPOINTSBUB> Subtraction;

};

template <class CutType,int CPOINTS_T, int TPOINTSBUB_T, int TPOINTSTRI_T, int YPOINTS_T> struct General_Triangle_Specification {
public:
	enum { CPOINTS = CPOINTS_T };
	enum { TPOINTSBUB = TPOINTSBUB_T };
	enum { TPOINTSTRI = TPOINTSTRI_T };
	enum { YPOINTS = YPOINTS_T };

	typedef cut::Darren::triangle_Darren_eval_points<CPOINTS, TPOINTSTRI> MomentaEvaluator;
	typedef cut::Darren::Normal_Corner_Tree_Strategy<MomentaEvaluator,CutType,CPOINTS> CornerTreeStrategy;
//	typedef cut::Darren::General_Triangle_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> TriangleCombiner;
	typedef cut::Darren::General_Triangle_Subtraction<CutType,CPOINTS,TPOINTSBUB,YPOINTS> Subtraction;
};

template <class CutType> struct Skeleton_Triangle_Specification {
public:
	enum { CPOINTS=7 };
	enum { TPOINTSBUB=4 };
	enum { TPOINTSTRI=8 };
	enum { YPOINTS=2 };

	typedef cut::Darren::triangle_Darren_eval_points<CPOINTS, TPOINTSTRI> MomentaEvaluator;
//	typedef cut::Darren::Skeleton_Corner_Tree_Strategy <MomentaEvaluator,CutType,CPOINTS>  CornerTreeStrategy;
	typedef cut::Darren::Normal_Corner_Tree_Strategy <MomentaEvaluator,CutType,CPOINTS>  CornerTreeStrategy;
//	typedef cut::Darren::Normal_Triangle_Combiner<YPOINTS,TPOINTSBUB,MomentaEvaluator> TriangleCombiner;
	typedef cut::Darren::Normal_Triangle_Subtraction<CutType,CPOINTS,TPOINTSBUB> Subtraction;
};


template <class CutType> struct Normal_Higgs_Triangle_Specification : public  General_Triangle_Specification<CutType,9,5,9,4> {};
template <class CutType> struct Normal_test_Triangle_Specification : public  General_Triangle_Specification<CutType,7,6,10,5> {};

}

#endif /* TRIANGLE_SPECIFICATION_H_ */
