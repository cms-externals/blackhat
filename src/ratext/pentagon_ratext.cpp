/*
 * pentagon_ratext.cpp
 *
 *  Created on: 21 Aug 2009
 *      Author: darrenforde
 */

#include "BH_utilities.h"
#include "rec_tree.h"
#include <iostream>
#include <cassert>
#include <string>
#include "pentagon_ratext.h"
#include "coeff_param.h"
#ifndef BH_PUBLIC
#include "ratext/basecutRat.h"
#endif
#include "ratext/rat_worker.h"

#define _VERBOSE 0

using namespace std;

namespace BH {

namespace ratext {

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> void pentagon_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init(){
}


//Explicitly Instantiate for all used cases

//
//The standard template
//

template void pentagon_rat_eval_points<3,2,2,7,3,3>::init();

template <class BoxType,class Specs,class T> struct initializer {
	void init(T,bool& k1massive,bool& k2massive);
};
#ifndef BH_PUBLIC

template <class Specs> struct initializer<basepentagonRat,Specs,const pentagonRat_comp*> {
	static void init(const pentagonRat_comp* t,bool& k1massive,bool& k2massive){
		if(t->corner_size(1)>1||t->get_process(1).p(2).mass_label()>0){
			k1massive=true;
		}
		else{
			k1massive=false;
		}
		if(t->corner_size(5)>1||t->get_process(5).p(2).mass_label()>0){
			k2massive=true;
		}
		else{
			k2massive=false;
		}

	};
};
#endif
template <class Specs> struct initializer<rat_worker,Specs,std::istream&> {
	static void init(std::istream& is,bool& k1massive,bool& k2massive){
		string title,value;
		is >> title;
		assert(title == "Pent_Ratspecific");

		is >> title;
		assert(title == "k1m");

		is >> value;
		if (value == "t" ){
			k1massive=true;
		} else if (value == "f" ){
			k1massive=false;
		} else {
			_WARNING3("Unexpected input ",value, " in box_Rat constructor.");
			abort();
		}

		is >> title;
		assert(title == "k2m");

		is >> value;
		if (value == "t" ){
			k2massive=true;
		} else if (value == "f" ){
			k2massive=false;
		} else {
			_WARNING3("Unexpected input ",value, " in box_Rat constructor.");
			abort();
		}
	};
};


template <class PentType,class PentSpecs> template <class T> pentagon_Rat<PentType,PentSpecs>::pentagon_Rat(T* t) : PentType(t), PentSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(5))
{
	initializer<PentType,PentSpecs,T*>::init(t,_k1massive,_k2massive);
	init();
}
template <class PentType,class PentSpecs> template <class T> pentagon_Rat<PentType,PentSpecs>::pentagon_Rat(T& t) : PentType(t), PentSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(5))
{
	initializer<PentType,PentSpecs,T&>::init(t,_k1massive,_k2massive);
	init();
}


#ifndef BH_PUBLIC
template BH::ratext::pentagon_Rat<BH::basepentagonRat, BH::ratext::Normal_RatPent_Specification<BH::basepentagonRat> >::pentagon_Rat(BH::pentagonRat_comp const*);
#endif
template BH::ratext::pentagon_Rat<BH::ratext::rat_worker, BH::ratext::Normal_RatPent_Specification<BH::ratext::rat_worker> >::pentagon_Rat(std::istream&);


//
//The Higgs template
//

template void pentagon_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::init();

#ifndef BH_PUBLIC
template BH::ratext::pentagon_Rat<BH::basepentagonRat, BH::ratext::Higgs_RatPent_Specification<BH::basepentagonRat> >::pentagon_Rat(BH::pentagonRat_comp const*);
#endif
template BH::ratext::pentagon_Rat<BH::ratext::rat_worker, BH::ratext::Higgs_RatPent_Specification<BH::ratext::rat_worker> >::pentagon_Rat(std::istream&);


}

}
