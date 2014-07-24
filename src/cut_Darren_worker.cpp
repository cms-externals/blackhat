/*
 * cut_Darren_worker.cpp
 *
 *  Created on: 1 Jul 2009
 *      Author: daniel
 */

#include "cut_Darren.h"
#include "cut_Darren.hpp"
#include "cut_Darren_ep.hpp"

#include "bubble_specification.h"


namespace BH {

namespace cut {

namespace Darren {
	
//
// General definitions and declarations
//
template <class TriangleSpecs> triangle_Darren<BH::cut::worker::worker_cutD,TriangleSpecs>* triangle_Darren_factory<BH::cut::worker::worker_cutD,TriangleSpecs>::new_triangle(std::istream& is)
{
	string title;
	is >> title;
	assert(title=="TType");
	int tri_type;
	is >> tri_type;
	switch (tri_type){
		case 1: case 2: return new triangle_Darren_plusminus<BH::cut::worker::worker_cutD,TriangleSpecs>(is);
		case 3: return new triangle_Darren_3mass<BH::cut::worker::worker_cutD,TriangleSpecs>(is);
	}
}
	

//
//The standard template
//

typedef cut::Darren::bubble_Darren_eval_points<BUBPOINTS_STD, BUBYPOINTS_STD> StandardMomentaEvaluator;
typedef cut::Darren::Normal_Corner_Tree_Strategy<StandardMomentaEvaluator,BH::cut::worker::worker_cutD,TRIPOINTS_STD> StandardCornerTreeStrategy;

	
typedef Normal_Triangle_Specification<BH::cut::worker::worker_cutD> NTS;
typedef Normal_Bubble_Specification<BH::cut::worker::worker_cutD> NBS;
//typedef General_Bubble_Specification<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,BUBPOINTS_STD,TRIPOINTS_STD,BUBYPOINTS_STD> NBS;
//typedef General_Triangle_Specification<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,BUBPOINTS_STD,TRIPOINTS_STD,BUBYPOINTS_STD> NTS;
	
typedef BH::cut::worker::worker_cutD WCD;

template <> template <> box_Darren<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,TRIPOINTS_STD>::box_Darren(std::istream& is): worker_cutD(is) {

	string title;
	is >> title;
	assert(title == "CDspecific");

	int kleg1,kleg2,kleg3,kleg4, masslessleg_type;

	is >> _k1leg;
	is >> _k2leg;
	is >> _k3leg;
	is >> _k4leg;
	is >> _masslessleg_type;
	is >> _massless_K1;



	init();
}
template <> template <> triangle_Darren<BH::cut::worker::worker_cutD,NTS>::triangle_Darren(std::istream& is): worker_cutD(is) {

	string title;
	is >> title;
	assert(title == "CDspecific");

	int kleg1,kleg2,kleg3, masslessleg_type;

	is >> _k1leg;
	is >> _k2leg;
	is >> _k3leg;
	is >> _masslessleg_type;

	// The default is for reverse to be 1 and it seems it is never changed!
	_reverse=1;


	init();
}
template <> template <> bubble_Darren<BH::cut::worker::worker_cutD,NBS>::bubble_Darren(std::istream& is): worker_cutD(is) {

	init();
}

template <> template <> triangle_Darren_3mass<BH::cut::worker::worker_cutD,NTS>::triangle_Darren_3mass(std::istream& is): triangle_Darren<BH::cut::worker::worker_cutD,NTS>(is){};
template <> template <> triangle_Darren_plusminus<BH::cut::worker::worker_cutD,NTS>::triangle_Darren_plusminus(std::istream& is): triangle_Darren<worker_cutD,NTS>(is){};

template class box_Darren<BH::cut::worker::worker_cutD,CTRIPOINTS_STD,TRIPOINTS_STD>;
template class triangle_Darren<BH::cut::worker::worker_cutD,NTS>;
template class bubble_Darren<BH::cut::worker::worker_cutD,NBS>;

template triangle_Darren<BH::cut::worker::worker_cutD,NTS>* triangle_Darren_factory<BH::cut::worker::worker_cutD,NTS>::new_triangle(std::istream& is);


//
// The higgs case
//

typedef cut::Darren::bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS> HiggsMomentaEvaluator;
typedef cut::Darren::Normal_Corner_Tree_Strategy<HiggsMomentaEvaluator,BH::cut::worker::worker_cutD,CTRIPOINTS_HIGGS> HiggsCornerTreeStrategy;
typedef General_Bubble_Specification<BH::cut::worker::worker_cutD,CTRIPOINTS_HIGGS,BUBPOINTS_HIGGS,TRIPOINTS_HIGGS,BUBYPOINTS_HIGGS> HiggsBubbleSpecs;
typedef General_Triangle_Specification<BH::cut::worker::worker_cutD,CTRIPOINTS_HIGGS,BUBPOINTS_HIGGS,TRIPOINTS_HIGGS,BUBYPOINTS_HIGGS> HiggsTriangleSpecs;

template <> template <> triangle_Darren<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>::triangle_Darren(std::istream& is): worker_cutD(is) {
	init();
}
template <> template <> bubble_Darren<BH::cut::worker::worker_cutD,HiggsBubbleSpecs>::bubble_Darren(std::istream& is): worker_cutD(is) {
	init();
}

template <> template <> triangle_Darren_3mass<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>::triangle_Darren_3mass(std::istream& is): triangle_Darren<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>(is){};;
template <> template <> triangle_Darren_plusminus<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>::triangle_Darren_plusminus(std::istream& is): triangle_Darren<worker_cutD,HiggsTriangleSpecs>(is){};

template class triangle_Darren<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>;
template class bubble_Darren<BH::cut::worker::worker_cutD,HiggsBubbleSpecs>;

template void BH::cut::Darren::box_Darren<BH::cut::worker::worker_cutD,CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::get_sub_terms(momentum_configuration<R>& mc, const vector<int>& ind, complex<R>* coeffsret, coeffparam<R,CTRIPOINTS_HIGGS>& tp);
template void BH::cut::Darren::box_Darren<BH::cut::worker::worker_cutD,CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::get_coeffs_fn(momentum_configuration<R>& mc, const vector<int>& ind, coeffparam<R,CTRIPOINTS_HIGGS>& tp);


template triangle_Darren<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>* triangle_Darren_factory<BH::cut::worker::worker_cutD,HiggsTriangleSpecs>::new_triangle(std::istream& is);



//
// A general test
//


typedef Normal_test_Triangle_Specification<BH::cut::worker::worker_cutD> NTTS;

template <> template <> triangle_Darren<BH::cut::worker::worker_cutD,NTTS>::triangle_Darren(std::istream& is): worker_cutD(is) {
	init();
}
	template <> template <> bubble_Darren<BH::cut::worker::worker_cutD,Normal_test_Bubble_Specification<BH::cut::worker::worker_cutD > >::bubble_Darren(std::istream& is): worker_cutD(is) {
	init();
}
	template class triangle_Darren<BH::cut::worker::worker_cutD,NTTS>;
	template class bubble_Darren<BH::cut::worker::worker_cutD,Normal_test_Bubble_Specification<BH::cut::worker::worker_cutD > >;


	template void BH::cut::Darren::box_Darren<BH::cut::worker::worker_cutD,TRICOEFFS,TRICOEFFSEP>::get_sub_terms(momentum_configuration<R>& mc, const vector<int>& ind, complex<R>* coeffsret, coeffparam<R,TRICOEFFS>& tp);
	template void BH::cut::Darren::box_Darren<BH::cut::worker::worker_cutD,TRICOEFFS,TRICOEFFSEP>::get_coeffs_fn(momentum_configuration<R>& mc, const vector<int>& ind, coeffparam<R,TRICOEFFS>& tp);

template <> template <> triangle_Darren_3mass<BH::cut::worker::worker_cutD,NTTS>::triangle_Darren_3mass(std::istream& is): triangle_Darren<BH::cut::worker::worker_cutD,NTTS>(is){};;
template <> template <> triangle_Darren_plusminus<BH::cut::worker::worker_cutD,NTTS>::triangle_Darren_plusminus(std::istream& is): triangle_Darren<worker_cutD,NTTS>(is){};

template triangle_Darren<BH::cut::worker::worker_cutD,NTTS>* triangle_Darren_factory<BH::cut::worker::worker_cutD,NTTS>::new_triangle(std::istream& is);


}
}
}
