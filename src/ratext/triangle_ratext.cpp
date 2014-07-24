/*
 * triangle_ratext.cpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */



#include "BH_utilities.h"
#include "rec_tree.h"
#include <iostream>
#include <cassert>
#include <string>
#include "triangle_ratext.h"
#include "box_ratext.h"
#ifndef BH_PUBLIC
#include "ratext/basecutRat.h"
#endif
#include "ratext/rat_worker.h"
#include "polylog.h"  // for pi

#define _VERBOSE 0

#define _T_EVAL_RADIUS 1

using namespace std;

namespace BH {

namespace ratext {


template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> const bool tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points(){

	for(int i=0;i<CPOINTS_T ;i++){
#if BH_USE_GMP
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i]=RGMP(_T_EVAL_RADIUS)*exp(CGMP(0,1)*pi<RGMP>()/RGMP(CPOINTS_T))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(CPOINTS_T));
#endif
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]=RVHP(_T_EVAL_RADIUS)*exp(CVHP(0,1)*pi<RVHP>()/RVHP(CPOINTS_T))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(CPOINTS_T));
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_HP[i]=to_HP(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i]=to_double(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
	}
#if _VERBOSE==1
	std::cout << " tri t points: {";
	for(int i=0;i<CPOINTS_T ;i++){
		std::cout << tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i] << ", ";
	}
	std::cout << "}" << std::endl;
#endif

	for(int j=0;j<CPOINTS_T;j++){
		for(int i=0;i<CPOINTS_T ;i++){
#if BH_USE_GMP
			tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*CPOINTS_T ]=pow(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i],j-(CPOINTS_T-1)/2)/RGMP(CPOINTS_T);
#endif
			tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CPOINTS_T ]=pow(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i],j-(CPOINTS_T-1)/2)/RVHP(CPOINTS_T);
			tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[i+j*CPOINTS_T ]=to_HP(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CPOINTS_T ]);
			tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*CPOINTS_T ]=to_double(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CPOINTS_T ]);
		}
	}
#ifdef BH_USE_GMP
    s_GMP_precision=RGMP::get_current_precision();
#endif

	return true;
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  long int  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances=0;

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> void tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init(){

	//Setup the evaluation points if this has not been dont yet
	if(tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances==0){
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances++;
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points();

		//Set up the pointers to the rational functions
		_rat_int[0]=&Rat_Int_Tri_1<R>;
		_rat_int[1]=&Rat_Int_Tri_2<R>;
		_rat_int[2]=&Rat_Int_Tri_3<R>;
		_rat_int[3]=&Rat_Int_Tri_4<R>;
		_rat_int_HP[0]=&Rat_Int_Tri_1<RHP>;
		_rat_int_HP[1]=&Rat_Int_Tri_2<RHP>;
		_rat_int_HP[2]=&Rat_Int_Tri_3<RHP>;
		_rat_int_HP[3]=&Rat_Int_Tri_4<RHP>;
		_rat_int_VHP[0]=&Rat_Int_Tri_1<RVHP>;
		_rat_int_VHP[1]=&Rat_Int_Tri_2<RVHP>;
		_rat_int_VHP[2]=&Rat_Int_Tri_3<RVHP>;
		_rat_int_VHP[3]=&Rat_Int_Tri_4<RVHP>;

#if BH_USE_GMP
		_rat_int_GMP[0]=&Rat_Int_Tri_1<RGMP>;
		_rat_int_GMP[1]=&Rat_Int_Tri_2<RGMP>;
		_rat_int_GMP[2]=&Rat_Int_Tri_3<RGMP>;
		_rat_int_GMP[3]=&Rat_Int_Tri_4<RGMP>;
#endif
		//It also might be the case that there are no boxes in which case we must force the computation of the
		// m0points.
		box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init();
	}
}

	
// The rational integral results table with each entry multipied by 2 so the first entry is -1/2 etc. (So far only first two entrys are correct)
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  Rat_Int_Fac_Fn_Ptr  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int[4]={NULL,NULL,NULL,NULL};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  Rat_Int_Fac_Fn_Ptr_HP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_HP[4]={NULL,NULL,NULL,NULL};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  Rat_Int_Fac_Fn_Ptr_VHP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_VHP[4]={NULL,NULL,NULL,NULL};
#if BH_USE_GMP
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  Rat_Int_Fac_Fn_Ptr_GMP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_GMP[4]={NULL,NULL,NULL,NULL};
#endif
	
// A function for returning the appropriate rational integral function
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> inline const C  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const C& s1, const C& s2, const C& s3){
	return (*_rat_int[element])(s1,s2,s3);
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> inline const CHP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CHP& s1, const CHP& s2, const CHP& s3){
	return (*_rat_int_HP[element])(s1,s2,s3);
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> inline const CVHP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CVHP& s1, const CVHP& s2, const CVHP& s3){
	return (*_rat_int_VHP[element])(s1,s2,s3);
}
	
#if BH_USE_GMP
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> inline const CGMP  tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CGMP& s1, const CGMP& s2, const CGMP& s3){
	return (*_rat_int_GMP[element])(s1,s2,s3);
}

#endif
//
//Explicitly Instantiate or specialise for all used cases
//


//
//The standard template

template <> C tri_rat_eval_points<3,2,2,7,3,3>::_circpos[7]={C(0,0)};
template <> CHP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_HP[7]={CHP(0,0)};
template <> CVHP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_VHP[7]={CVHP(0,0)};

template <> C tri_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix[7*7]={C(0,0)};
template <> CHP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_HP[7*7]={CHP(0,0)};
template <> CVHP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_VHP[7*7]={CVHP(0,0)};

template void tri_rat_eval_points<3,2,2,7,3,3>::init();
template const C tri_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const C& s1, const C& s2, const C& s3);
template const CHP tri_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CHP& s1, const CHP& s2, const CHP& s3);
template const CVHP tri_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CVHP& s1, const CVHP& s2, const CVHP& s3);	

#if BH_USE_GMP
template <> CGMP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_GMP[7]={CGMP(0,0)};
template <> CGMP tri_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_GMP[7*7]={CGMP(0,0)};
template const CGMP tri_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CGMP& s1, const CGMP& s2, const CGMP& s3);
template <> int tri_rat_eval_points<3,2,2,7,3,3>::s_GMP_precision=0;
#endif

template <class TriangleType,class Specs,class T> struct initializer {
	void init(T,bool& k1massive,bool& k2massive,bool& k3massive);
};

#ifndef BH_PUBLIC
template <class Specs> struct initializer<basetriangleRat,Specs,const triangleRat_comp*> {
	static void init(const triangleRat_comp* t,bool& k1massive,bool& k2massive,bool& k3massive){
		if(t->corner_size(1)>1||t->get_process(1).p(2).mass_label()>0){
			k1massive=true;
		}
		else{
			k1massive=false;
		}
		if(t->corner_size(3)>1||t->get_process(3).p(2).mass_label()>0){
			k2massive=true;
		}
		else{
			k2massive=false;
		}
		if(t->corner_size(2)>1||t->get_process(2).p(2).mass_label()>0){
			k3massive=true;
		}
		else{
			k3massive=false;
		}
	};
};
#endif

template <class Specs> struct initializer<rat_worker,Specs,std::istream&> {
	static void init(std::istream& is,bool& k1massive,bool& k2massive,bool& k3massive){
		string title,value;
		is >> title;
		assert(title == "Triangle_Ratspecific");

		is >> title;
		assert(title == "k1m");

		is >> value;
		if (value == "t" ){
			k1massive=true;
		} else if (value == "f" ){
			k1massive=false;
		} else {
			_WARNING3("Unexpected input ",value, " in triangle_Rat constructor.");
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
			_WARNING3("Unexpected input ",value, " in triangle_Rat constructor.");
			abort();
		}

		is >> title;
		assert(title == "k3m");

		is >> value;
		if (value == "t" ){
			k3massive=true;
		} else if (value == "f" ){
			k3massive=false;
		} else {
			_WARNING3("Unexpected input ",value, " in triangle_Rat constructor.");
			abort();
		}
};
};


template <class TriType,class TriSpecs> template <class T> triangle_Rat<TriType,TriSpecs>::triangle_Rat(T* t) : TriType(t), TriSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(3))
{
	initializer<TriType,TriSpecs,T*>::init(t,_k1massive,_k2massive,_k3massive);
	init();
}
template <class TriType,class TriSpecs> template <class T> triangle_Rat<TriType,TriSpecs>::triangle_Rat(T& t) : TriType(t), TriSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(3))
{
	initializer<TriType,TriSpecs,T&>::init(t,_k1massive,_k2massive,_k3massive);
	init();
}

#ifndef BH_PUBLIC
template BH::ratext::triangle_Rat<BH::basetriangleRat, BH::ratext::Normal_RatTri_Specification<BH::basetriangleRat> >::triangle_Rat(BH::triangleRat_comp const*);
#endif
template BH::ratext::triangle_Rat<BH::ratext::rat_worker, BH::ratext::Normal_RatTri_Specification<BH::ratext::rat_worker> >::triangle_Rat(std::istream&);


//
//The Higgs template

template <> C tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos[CPOINTS_HIGGS]={C(0,0)};
template <> CHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_HP[CPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_VHP[CPOINTS_HIGGS]={CVHP(0,0)};

template <> C tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix[CPOINTS_HIGGS*CPOINTS_HIGGS]={C(0,0)};
template <> CHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_HP[CPOINTS_HIGGS*CPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_VHP[CPOINTS_HIGGS*CPOINTS_HIGGS]={CVHP(0,0)};

template void tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::init();
template const C tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const C& s1, const C& s2, const C& s3);
template const CHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CHP& s1, const CHP& s2, const CHP& s3);
template const CVHP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CVHP& s1, const CVHP& s2, const CVHP& s3);
	
#if BH_USE_GMP
template <> CGMP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_GMP[CPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_GMP[CPOINTS_HIGGS*CPOINTS_HIGGS]={CGMP(0,0)};
template const CGMP tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CGMP& s1, const CGMP& s2, const CGMP& s3);
template <> int tri_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::s_GMP_precision=0;
#endif

#ifndef BH_PUBLIC
template BH::ratext::triangle_Rat<BH::basetriangleRat, BH::ratext::Higgs_RatTri_Specification<BH::basetriangleRat> >::triangle_Rat(BH::triangleRat_comp const*);
#endif
template BH::ratext::triangle_Rat<BH::ratext::rat_worker, BH::ratext::Higgs_RatTri_Specification<BH::ratext::rat_worker> >::triangle_Rat(std::istream&);


} /* ratext */

} /* BH */
