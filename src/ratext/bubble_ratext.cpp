/*
 * bubble_ratext.cpp
 *
 *  Created on: 28 Aug 2009
 *      Author: darrenforde
 */


#include "BH_utilities.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include <iostream>
#include <cassert>
#include "bubble_ratext.h"
#include "box_ratext.h"
#include "polylog.h"  // for pi

#define _VERBOSE 0

#define _T_EVAL_RADIUS 1
#define _Y_EVAL_RADIUS 1

namespace BH {

namespace ratext {


template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> const bool bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points(){

	for(int i=0;i<CBUBPOINTS_T ;i++){
#if BH_USE_GMP
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i]=RGMP(_T_EVAL_RADIUS)*exp(CGMP(0,2)*pi<RGMP>()/RGMP(CBUBPOINTS_T+1))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(CBUBPOINTS_T));
#endif
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]=RVHP(_T_EVAL_RADIUS)*exp(CVHP(0,2)*pi<RVHP>()/RVHP(CBUBPOINTS_T+1))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(CBUBPOINTS_T));
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_HP[i]=to_HP(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i]=to_double(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
	}

#if _VERBOSE==1
	std::cout << " bub t points: {";
	for(int i=0;i<CBUBPOINTS_T ;i++){
		std::cout << bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i] << ", ";	
    }
    std::cout << "}" << std::endl;
#endif
    
	for(int j=0;j<CBUBPOINTS_T;j++){
		for(int i=0;i<CBUBPOINTS_T ;i++){
#if BH_USE_GMP
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*CBUBPOINTS_T]=pow(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i],j-(CBUBPOINTS_T-1)/2)/RGMP(CBUBPOINTS_T);
#endif
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CBUBPOINTS_T]=pow(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i],j-(CBUBPOINTS_T-1)/2)/RVHP(CBUBPOINTS_T);
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[i+j*CBUBPOINTS_T]=to_HP(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CBUBPOINTS_T]);
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*CBUBPOINTS_T]=to_double(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*CBUBPOINTS_T]);
		}
	}

#if _VERBOSE==1
	std::cout << " bub t matrix: {";
	for(int j=0;j<CBUBPOINTS_T;j++){
		for(int i=0;i<CBUBPOINTS_T ;i++){
		std::cout << bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*CBUBPOINTS_T] << ", ";	
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
#endif


	// Set up the points for y
	for(int i=0;i<YPOINTS_T;i++){
#if BH_USE_GMP
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_GMP[i]=RGMP(_Y_EVAL_RADIUS)*exp(CGMP(0,1)*pi<RGMP>()/RGMP(YPOINTS_T))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(YPOINTS_T));
#endif
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_VHP[i]=RVHP(_Y_EVAL_RADIUS)*exp(CVHP(0,1)*pi<RVHP>()/RVHP(YPOINTS_T))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(YPOINTS_T));
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_HP[i]=to_HP(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_VHP[i]);
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos[i]=to_double(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_VHP[i]);
	}
	for(int j=0;j<YPOINTS_T;j++){
		for(int i=0;i<YPOINTS_T;i++){
#if BH_USE_GMP
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix_GMP[i+j*YPOINTS_T]=pow(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_GMP[i],-j)/RGMP(YPOINTS_T);
#endif
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix_VHP[i+j*YPOINTS_T]=pow(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_VHP[i],-j)/RVHP(YPOINTS_T);
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix_HP[i+j*YPOINTS_T]=to_HP(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix_VHP[i+j*YPOINTS_T]);
			bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix[i+j*YPOINTS_T]=to_double(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix_VHP[i+j*YPOINTS_T]);
		}
	}

#if _VERBOSE==1
		std::cout << " bub y points: {";
		for(int i=0;i<YPOINTS_T ;i++){
			std::cout << bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos[i] << ", ";
		}
		std::cout << "}" << std::endl;

    std::cout << " bub y matrix: {";
    for(int j=0;j<YPOINTS_T;j++){
		for(int i=0;i<YPOINTS_T;i++){
       std::cout << bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_y_circpos_matrix[i+j*YPOINTS_T] << ", ";
		}
        std::cout << std::endl;
	}    
    std::cout << "}" << std::endl;
#endif

#ifdef BH_USE_GMP
    s_GMP_precision=RGMP::get_current_precision();
#endif

	return true;
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  long int  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances=0;

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> void bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init(){

	//Setup the evaluation points if this has not been done yet
	if(bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances==0){
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances++;
		bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points();

		//It also might be the case that there are no boxes/ and or triangles in which case we must force the computation of the
		// circle of points to evaluate on.
		tri_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init();
	}
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const C  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int[4]={C(0,1)/R(6),C(0,1)/R(60),C(0,1)/R(420),C(0,1)/R(2520)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_HP[4]={CHP(0,1)/RHP(6),CHP(0,1)/RHP(60),CHP(0,1)/RHP(420),CHP(0,1)/RHP(2520)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CVHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_VHP[4]={CVHP(0,1)/RVHP(6),CVHP(0,1)/RVHP(60),CVHP(0,1)/RVHP(420),CVHP(0,1)/RVHP(2520)};

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const C  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_extra_fac[4]={-C(3,0),C(2,0),C(1,0),C(1,0)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_extra_fac_HP[4]={-CHP(3,0),CHP(2,0),CHP(1,0),CHP(1,0)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CVHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_extra_fac_VHP[4]={-CVHP(3,0),CVHP(2,0),CVHP(1,0),CVHP(1,0)};

// A function for returning the approriate rational integral function
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const C  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const C& s1){
	return pow(s1,element+1)*_rat_int[element];
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CHP& s1){
	return pow(s1,element+1)*_rat_int_HP[element];
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CVHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CVHP& s1){
	return pow(s1,element+1)*_rat_int_VHP[element];
}


// A function for returning the approriate extra coefficient factor
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const C  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_extra_fac(int element, const C& s1){
	return R(1)/(_extra_fac[element]*s1);
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_extra_fac(int element, const CHP& s1){
	return RHP(1)/(_extra_fac_HP[element]*s1);
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CVHP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_extra_fac(int element, const CVHP& s1){
	return RVHP(1)/(_extra_fac_VHP[element]*s1);
}

#if BH_USE_GMP
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>   CGMP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_GMP[4]={CGMP(0,1)/RGMP(6),CGMP(0,1)/RGMP(60),CGMP(0,1)/RGMP(420),CGMP(0,1)/RGMP(2520)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>   CGMP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_extra_fac_GMP[4]={-CGMP(3,0),CGMP(2,0),CGMP(1,0),CGMP(1,0)};

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CGMP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CGMP& s1){
	return pow(s1,element+1)*_rat_int_GMP[element];
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  inline const CGMP  bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_extra_fac(int element, const CGMP& s1){
	if (s_GMP_precision< RGMP::get_current_precision()){
		recomuteGMP();
	}
	return RGMP(1)/(_extra_fac_GMP[element]*s1);
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  void bub_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::recomuteGMP() {

	BH_DEBUG_MESSAGE2("recomputing bubble points for precision ",RGMP::get_current_precision());get_points();
	_extra_fac_GMP[0]=-CGMP(3,0);
	_extra_fac_GMP[1]=CGMP(2,0);
	_extra_fac_GMP[2]=CGMP(1,0);
	_extra_fac_GMP[3]=CGMP(1,0);

	_rat_int_GMP[0]=CGMP(0,1)/RGMP(6);
	_rat_int_GMP[1]=CGMP(0,1)/RGMP(60);
	_rat_int_GMP[2]=CGMP(0,1)/RGMP(420);
	_rat_int_GMP[3]=CGMP(0,1)/RGMP(2520);
}


#endif
//Explicitly Instantiate for all used cases

//
//The standard template
//
template <> C bub_rat_eval_points<3,2,2,7,3,3>::_circpos[3]={C(0,0)};
template <> CHP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_HP[3]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_VHP[3]={CVHP(0,0)};

template <> C bub_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix[3*3]={C(0,0)};
template <> CHP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_HP[3*3]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_VHP[3*3]={CVHP(0,0)};

template <> C bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos[3]={C(0,0)};
template <> CHP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_HP[3]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_VHP[3]={CVHP(0,0)};

template <> C bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_matrix[3*3]={C(0,0)};
template <> CHP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_matrix_HP[3*3]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_matrix_VHP[3*3]={CVHP(0,0)};

template void bub_rat_eval_points<3,2,2,7,3,3>::init();
template const C bub_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const C& s1);
template const CHP bub_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CHP& s1);
template const CVHP bub_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CVHP& s1);

template const C bub_rat_eval_points<3,2,2,7,3,3>::get_extra_fac(int element, const C& s1);
template const CHP bub_rat_eval_points<3,2,2,7,3,3>::get_extra_fac(int element, const CHP& s1);
template const CVHP bub_rat_eval_points<3,2,2,7,3,3>::get_extra_fac(int element, const CVHP& s1);


#if BH_USE_GMP
template <> CGMP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_GMP[3]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_GMP[3*3]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_GMP[3]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<3,2,2,7,3,3>::_y_circpos_matrix_GMP[3*3]={CGMP(0,0)};
template const CGMP bub_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CGMP& s1);
template const CGMP bub_rat_eval_points<3,2,2,7,3,3>::get_extra_fac(int element, const CGMP& s1);
template <> int bub_rat_eval_points<3,2,2,7,3,3>::s_GMP_precision=0;



#endif

//
//The Higgs template
//
template <> C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos[CBUBPOINTS_HIGGS]={C(0,0)};
template <> CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_HP[CBUBPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_VHP[CBUBPOINTS_HIGGS]={CVHP(0,0)};

template <> C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix[CBUBPOINTS_HIGGS*CBUBPOINTS_HIGGS]={C(0,0)};
template <> CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_HP[CBUBPOINTS_HIGGS*CBUBPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_VHP[CBUBPOINTS_HIGGS*CBUBPOINTS_HIGGS]={CVHP(0,0)};

template <> C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos[YPOINTS_HIGGS]={C(0,0)};
template <> CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_HP[YPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_VHP[YPOINTS_HIGGS]={CVHP(0,0)};

template <> C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_matrix[YPOINTS_HIGGS*YPOINTS_HIGGS]={C(0,0)};
template <> CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_matrix_HP[YPOINTS_HIGGS*YPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_matrix_VHP[YPOINTS_HIGGS*YPOINTS_HIGGS]={CVHP(0,0)};

template void bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::init();
template const C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const C& s1);
template const CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CHP& s1);
template const CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CVHP& s1);

template const C bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_extra_fac(int element, const C& s1);
template const CHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_extra_fac(int element, const CHP& s1);
template const CVHP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_extra_fac(int element, const CVHP& s1);

#if BH_USE_GMP
template <> CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_GMP[CBUBPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_GMP[CBUBPOINTS_HIGGS*CBUBPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_GMP[YPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_y_circpos_matrix_GMP[YPOINTS_HIGGS*YPOINTS_HIGGS]={CGMP(0,0)};
template const CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CGMP& s1);
template const CGMP bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_extra_fac(int element, const CGMP& s1);
template <> int bub_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::s_GMP_precision=0;

#endif

}

}

