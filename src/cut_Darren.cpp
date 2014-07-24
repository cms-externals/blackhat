/*
 * cut_Darren.cpp
 *
 *  Created on: 6 Jul 2009
 *      Author: darrenforde
 */

#include "BH_utilities.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include <iostream>
#include <cassert>
#include "cut_worker.h"
#include "cut_Darren.h"
#include "coeff_param.h"
#include "polylog.h"  //  for pi

#define _VERBOSE 0

namespace BH {

namespace cut {

namespace Darren {
	
//Explicitly Instantiate for all used cases

//
//The standard template
//

//Definition for the fixed 2 point solution
template <> C bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y[BUBYPOINTS_STD]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_HP[BUBYPOINTS_STD]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_VHP[BUBYPOINTS_STD]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_matrix[BUBYPOINTS_STD*BUBYPOINTS_STD]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_matrix_HP[BUBYPOINTS_STD*BUBYPOINTS_STD]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_matrix_VHP[BUBYPOINTS_STD*BUBYPOINTS_STD]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos[BUBPOINTS_STD]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_HP[BUBPOINTS_STD]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_VHP[BUBPOINTS_STD]={CVHP(0,0)};

template <> C triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::coeffext[CTRIPOINTS_STD*TRIPOINTS_STD]={C(0,0)};
template <> CHP triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::coeffext_HP[CTRIPOINTS_STD*TRIPOINTS_STD]={CHP(0,0)};
template <> CVHP triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::coeffext_VHP[CTRIPOINTS_STD*TRIPOINTS_STD]={CVHP(0,0)};

template <> C box_Darren_eval_points<TRIPOINTS_STD>::circpos[TRIPOINTS_STD]={C(0,0)};
template <> CHP box_Darren_eval_points<TRIPOINTS_STD>::circpos_HP[TRIPOINTS_STD]={CHP(0,0)};
template <> CVHP box_Darren_eval_points<TRIPOINTS_STD>::circpos_VHP[TRIPOINTS_STD]={CVHP(0,0)};



#if BH_USE_GMP
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_GMP[BUBYPOINTS_STD]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_y_matrix_GMP[BUBYPOINTS_STD*BUBYPOINTS_STD]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::_circpos_GMP[BUBPOINTS_STD]={CGMP(0,0)};
template <> CGMP triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::coeffext_GMP[CTRIPOINTS_STD*TRIPOINTS_STD]={CGMP(0,0)};
template <> CGMP box_Darren_eval_points<TRIPOINTS_STD>::circpos_GMP[TRIPOINTS_STD]={CGMP(0,0)};
template <> int box_Darren_eval_points<TRIPOINTS_STD>::s_GMP_precision=0;
template <> int triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::s_GMP_precision=0;
template <> int bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::s_GMP_precision=0;

#endif


//
// The higgs case
//

template <> C bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y[BUBYPOINTS_HIGGS]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_HP[BUBYPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_VHP[BUBYPOINTS_HIGGS]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_matrix[BUBYPOINTS_HIGGS*BUBYPOINTS_HIGGS]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_matrix_HP[BUBYPOINTS_HIGGS*BUBYPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_matrix_VHP[BUBYPOINTS_HIGGS*BUBYPOINTS_HIGGS]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos[BUBPOINTS_HIGGS]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_HP[BUBPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_VHP[BUBPOINTS_HIGGS]={CVHP(0,0)};

template <> C triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::coeffext[CTRIPOINTS_HIGGS*TRIPOINTS_HIGGS]={C(0,0)};
template <> CHP triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::coeffext_HP[CTRIPOINTS_HIGGS*TRIPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::coeffext_VHP[CTRIPOINTS_HIGGS*TRIPOINTS_HIGGS]={CVHP(0,0)};

template <> C box_Darren_eval_points<TRIPOINTS_HIGGS>::circpos[TRIPOINTS_HIGGS]={C(0,0)};
template <> CHP box_Darren_eval_points<TRIPOINTS_HIGGS>::circpos_HP[TRIPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP box_Darren_eval_points<TRIPOINTS_HIGGS>::circpos_VHP[TRIPOINTS_HIGGS]={CVHP(0,0)};

#if BH_USE_GMP
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_GMP[BUBYPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_y_matrix_GMP[BUBYPOINTS_HIGGS*BUBYPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::_circpos_GMP[BUBPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::coeffext_GMP[CTRIPOINTS_HIGGS*TRIPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP box_Darren_eval_points<TRIPOINTS_HIGGS>::circpos_GMP[TRIPOINTS_HIGGS]={CGMP(0,0)};
template <> int box_Darren_eval_points<TRIPOINTS_HIGGS>::s_GMP_precision=0;
template <> int triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::s_GMP_precision=0;
template <> int bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::s_GMP_precision=0;
#endif


//
// A special test case
//

template <> C bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y[BUBCOEFFS]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_HP[BUBCOEFFS]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_VHP[BUBCOEFFS]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_matrix[BUBCOEFFS*BUBCOEFFS]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_matrix_HP[BUBCOEFFS*BUBCOEFFS]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_matrix_VHP[BUBCOEFFS*BUBCOEFFS]={CVHP(0,0)};

template <> C bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos[BUBCOEFFSEP]={C(0,0)};
template <> CHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_HP[BUBCOEFFSEP]={CHP(0,0)};
template <> CVHP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_VHP[BUBCOEFFSEP]={CVHP(0,0)};

template <> C triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::coeffext[TRICOEFFS*TRICOEFFSEP]={C(0,0)};
template <> CHP triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::coeffext_HP[TRICOEFFS*TRICOEFFSEP]={CHP(0,0)};
template <> CVHP triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::coeffext_VHP[TRICOEFFS*TRICOEFFSEP]={CVHP(0,0)};

template <> C box_Darren_eval_points<TRICOEFFSEP>::circpos[TRICOEFFSEP]={C(0,0)};
template <> CHP box_Darren_eval_points<TRICOEFFSEP>::circpos_HP[TRICOEFFSEP]={CHP(0,0)};
template <> CVHP box_Darren_eval_points<TRICOEFFSEP>::circpos_VHP[TRICOEFFSEP]={CVHP(0,0)};


#if BH_USE_GMP
template <> CGMP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_GMP[BUBCOEFFS]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_y_matrix_GMP[BUBCOEFFS*BUBCOEFFS]={CGMP(0,0)};
template <> CGMP bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::_circpos_GMP[BUBCOEFFSEP]={CGMP(0,0)};
template <> CGMP triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::coeffext_GMP[TRICOEFFS*TRICOEFFSEP]={CGMP(0,0)};
template <> CGMP box_Darren_eval_points<TRICOEFFSEP>::circpos_GMP[TRICOEFFSEP]={CGMP(0,0)};
template <> int bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::s_GMP_precision=0;
template <> int triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::s_GMP_precision=0;
template <> int box_Darren_eval_points<TRICOEFFSEP>::s_GMP_precision=0;
#endif

	//
	// Implementations
	//
	
template <> const bool bubble_Darren_eval_points<4,2>::get_y_points(){


	for(int i=0;i<4;i++){
		bubble_Darren_eval_points<4,2>::_circpos[i]=exp(C(0,2)*pi<R>()*R(i)/R(BUBPOINTS_STD));
		bubble_Darren_eval_points<4,2>::_circpos_HP[i]=exp(CHP(0,2)*pi<RHP>()*RHP(i)/RHP(BUBPOINTS_STD));
		bubble_Darren_eval_points<4,2>::_circpos_VHP[i]=exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(BUBPOINTS_STD));
#if BH_USE_GMP

		CGMP r1=pi<RGMP>()*RGMP(i);
		CGMP r2=CGMP(0,2)*pi<RGMP>()*RGMP(i);
		CGMP r3=CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(BUBPOINTS_STD);
		CGMP r4=exp(r3);
		BH_DEBUG_PRINT(r4);
		bubble_Darren_eval_points<4,2>::_circpos_GMP[i].real().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<4,2>::_circpos_GMP[i].imag().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<4,2>::_circpos_GMP[i]=r4;
		BH_DEBUG_MESSAGE(_circpos_GMP[i]);

		s_GMP_precision=RGMP::get_current_precision();
#endif
	}

	bubble_Darren_eval_points<4,2>::_circpos_y[0]=C(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y[1]=C(R(2)/R(3),0);
	bubble_Darren_eval_points<4,2>::_circpos_y_HP[0]=CHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_HP[1]=CHP(RHP(2)/RHP(3),0);
	bubble_Darren_eval_points<4,2>::_circpos_y_VHP[0]=CVHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_VHP[1]=CVHP(RVHP(2)/RVHP(3),0);
#if BH_USE_GMP
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[0].real().set_prec(RGMP::get_current_precision());
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[0].imag().set_prec(RGMP::get_current_precision());
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[1].real().set_prec(RGMP::get_current_precision());
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[1].imag().set_prec(RGMP::get_current_precision());
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[0]=CGMP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_GMP[1]=CGMP(RGMP(2)/RGMP(3),0);
#endif

	bubble_Darren_eval_points<4,2>::_circpos_y_matrix[0]=C(1,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix[1]=C(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix[2]=C(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix[3]=C(6,0);

	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_HP[0]=CHP(1,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_HP[1]=CHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_HP[2]=CHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_HP[3]=CHP(6,0);

	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_VHP[0]=CVHP(1,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_VHP[1]=CVHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_VHP[2]=CVHP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_VHP[3]=CVHP(6,0);
#if BH_USE_GMP
	for (int ii=0;ii<4;ii++){
		bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[ii].real().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[ii].real().set_prec(RGMP::get_current_precision());
	}
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[0]=CGMP(1,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[1]=CGMP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[2]=CGMP(0,0);
	bubble_Darren_eval_points<4,2>::_circpos_y_matrix_GMP[3]=CGMP(6,0);
#endif


	return true;
}


template <int TPOINTSBUB,int YPOINTS> const bool bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_y_points(){

#if _VERBOSE
	_PRINT(TPOINTSBUB);
	_PRINT(YPOINTS);
	cout << "bubble::circpos={";
#endif
	for(int i=0;i<TPOINTSBUB;i++){
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos[i]=exp(C(0,1)*pi<R>()/R(TPOINTSBUB+1))*exp(C(0,2)*pi<R>()*R(i)/R(TPOINTSBUB));
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_HP[i]=exp(CHP(0,1)*pi<RHP>()/RHP(TPOINTSBUB+1))*exp(CHP(0,2)*pi<RHP>()*RHP(i)/RHP(TPOINTSBUB));
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_VHP[i]=exp(CVHP(0,1)*pi<RVHP>()/RVHP(TPOINTSBUB+1))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(TPOINTSBUB));
#if BH_USE_GMP
		CGMP r1=exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(TPOINTSBUB));
		CGMP r2=exp(CGMP(0,1)*pi<RGMP>()/RGMP(TPOINTSBUB+1));
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_GMP[i].real().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_GMP[i].imag().set_prec(RGMP::get_current_precision());

		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_GMP[i]=exp(CGMP(0,1)*pi<RGMP>()/RGMP(TPOINTSBUB+1))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(TPOINTSBUB));
		s_GMP_precision=RGMP::get_current_precision();
#endif

#if _VERBOSE
		cout << bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos[i] << ",";
#endif
	}
#if _VERBOSE
	cout << "}" << endl;

	cout << "bubble::circpos_y={";
#endif
	for(int i=0;i<YPOINTS;i++){
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y[i]=exp(C(0,1)*pi<R>()/R(YPOINTS+1))*exp(C(0,2)*pi<R>()*R(i)/R(YPOINTS));
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_HP[i]=exp(CHP(0,1)*pi<RHP>()/RHP(YPOINTS+1))*exp(CHP(0,2)*pi<RHP>()*RHP(i)/RHP(YPOINTS));
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_VHP[i]=exp(CVHP(0,1)*pi<RVHP>()/RVHP(YPOINTS+1))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(YPOINTS));
#if BH_USE_GMP
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_GMP[i].real().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_GMP[i].imag().set_prec(RGMP::get_current_precision());
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_GMP[i]=exp(CGMP(0,1)*pi<RGMP>()/RGMP(YPOINTS+1))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(YPOINTS));
#endif

#if _VERBOSE
		cout << bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y[i] << ",";
#endif
	}
#if _VERBOSE
	cout << "}" << endl;

	cout << "bubble::circpos_y_matrix={";
#endif
	for(int j=0;j<YPOINTS;j++){
		for(int i=0;i<YPOINTS;i++){
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix[i+j*YPOINTS]=pow(bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y[i],-j)/R((j+1)*YPOINTS);
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix_HP[i+j*YPOINTS]=pow(bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_HP[i],-j)/RHP((j+1)*YPOINTS);
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix_VHP[i+j*YPOINTS]=pow(bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_VHP[i],-j)/RVHP((j+1)*YPOINTS);
#if BH_USE_GMP
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix_GMP[i+j*YPOINTS].real().set_prec(RGMP::get_current_precision());
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix_GMP[i+j*YPOINTS].imag().set_prec(RGMP::get_current_precision());
			bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix_GMP[i+j*YPOINTS]=pow(bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_GMP[i],-j)/RGMP((j+1)*YPOINTS);
#endif

#if _VERBOSE
			cout << bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_circpos_y_matrix[i+j*YPOINTS] << ",";
#endif
		}
#if _VERBOSE
		cout << endl;
#endif
	}
#if _VERBOSE
	cout << "}" << endl;
#endif

	return true;
}

template <int CPOINTS, int TPOINTSTRI> const bool triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::get_matrix_points(){

	// We construct the points that we evaluate in t, this is just a copy of the box version of this
	C coeffpts[TPOINTSTRI];
	CHP coeffpts_HP[TPOINTSTRI];
	CVHP coeffpts_VHP[TPOINTSTRI];
#if BH_USE_GMP
	CGMP coeffpts_GMP[TPOINTSTRI];
#endif
	for(int i=0;i<TPOINTSTRI;i++){
		coeffpts[i]=exp(C(0,1)*pi<R>()/R(TPOINTSTRI+1))*exp(C(0,2)*pi<R>()*R(i)/R(TPOINTSTRI));
		coeffpts_HP[i]=exp(CHP(0,1)*pi<RHP>()/RHP(TPOINTSTRI+1))*exp(CHP(0,2)*pi<RHP>()*RHP(i)/RHP(TPOINTSTRI));
		coeffpts_VHP[i]=exp(CVHP(0,1)*pi<RVHP>()/RVHP(TPOINTSTRI+1))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(TPOINTSTRI));
#if BH_USE_GMP
		coeffpts_GMP[i]=exp(CGMP(0,1)*pi<RGMP>()/RGMP(TPOINTSTRI+1))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(TPOINTSTRI));
		s_GMP_precision=RGMP::get_current_precision();
#endif
	}
#if _VERBOSE
	cout << "triangle::coeffext={";
#endif
	for(int j=0;j<CPOINTS;j++){
		for(int i=0;i<TPOINTSTRI;i++){
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext[i+j*TPOINTSTRI]=pow(coeffpts[i],j-(CPOINTS-1)/2)/R(TPOINTSTRI);
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext_HP[i+j*TPOINTSTRI]=pow(coeffpts_HP[i],j-(CPOINTS-1)/2)/RHP(TPOINTSTRI);
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext_VHP[i+j*TPOINTSTRI]=pow(coeffpts_VHP[i],j-(CPOINTS-1)/2)/RVHP(TPOINTSTRI);
#if BH_USE_GMP
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext_GMP[i+j*TPOINTSTRI].real().set_prec(RGMP::get_current_precision());
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext_GMP[i+j*TPOINTSTRI].imag().set_prec(RGMP::get_current_precision());
			triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext_GMP[i+j*TPOINTSTRI]=pow(coeffpts_GMP[i],j-(CPOINTS-1)/2)/RGMP(TPOINTSTRI);
#endif

#if _VERBOSE
			cout << triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::coeffext[i+j*TPOINTSTRI] << ",";
#endif
		}
#if _VERBOSE
		cout << endl;
#endif
	}
#if _VERBOSE
	cout << "}" << endl;
#endif

	return true;
}

template <int TPOINTSTRI> const bool box_Darren_eval_points<TPOINTSTRI>::get_t_points(){
#if _VERBOSE
	cout << "box::circpos={";
#endif
	for(int i=0;i<TPOINTSTRI;i++){
		box_Darren_eval_points<TPOINTSTRI>::circpos[i]=exp(C(0,1)*pi<R>()/R(TPOINTSTRI+1))*exp(C(0,2)*pi<R>()*R(i)/R(TPOINTSTRI));
		box_Darren_eval_points<TPOINTSTRI>::circpos_HP[i]=exp(CHP(0,1)*pi<RHP>()/RHP(TPOINTSTRI+1))*exp(CHP(0,2)*pi<RHP>()*RHP(i)/RHP(TPOINTSTRI));
		box_Darren_eval_points<TPOINTSTRI>::circpos_VHP[i]=exp(CVHP(0,1)*pi<RVHP>()/RVHP(TPOINTSTRI+1))*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(TPOINTSTRI));
#if BH_USE_GMP
		box_Darren_eval_points<TPOINTSTRI>::circpos_GMP[i].real().set_prec(RGMP::get_current_precision());
		box_Darren_eval_points<TPOINTSTRI>::circpos_GMP[i].imag().set_prec(RGMP::get_current_precision());
		box_Darren_eval_points<TPOINTSTRI>::circpos_GMP[i]=exp(CGMP(0,1)*pi<RGMP>()/RGMP(TPOINTSTRI+1))*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(TPOINTSTRI));
		s_GMP_precision=RGMP::get_current_precision();
#endif

#if _VERBOSE
		cout << box_Darren_eval_points<TPOINTSTRI>::circpos[i] << ",";
#endif
	}
#if _VERBOSE
	cout << "}" << endl;
#endif

	return true;
}

template <int TPOINTSBUB, int YPOINTS> long int bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_nbr_instances=0;
template <int CPOINTS, int TPOINTSTRI> long int triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::_nbr_instances=0;
template <int TPOINTSTRI>  long int  box_Darren_eval_points<TPOINTSTRI>::_nbr_instances=0;


template <int TPOINTSBUB, int YPOINTS> void  bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::init(){

	//Setup the evaluation points if this has not been dont yet
	if(bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_nbr_instances==0){
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::_nbr_instances++;
		bubble_Darren_eval_points<TPOINTSBUB,YPOINTS>::get_y_points();
	}
}

template <int CPOINTS, int TPOINTSTRI> void  triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::init(){

	//Setup the evaluation points if this has not been dont yet
	if(triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::_nbr_instances==0){
		triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::_nbr_instances++;
		triangle_Darren_eval_points<CPOINTS,TPOINTSTRI>::get_matrix_points();
	}
}

template <int TPOINTSTRI> void box_Darren_eval_points<TPOINTSTRI>::init() {

	//Setup the evaluation points if this has not been dont yet
	if(box_Darren_eval_points<TPOINTSTRI>::_nbr_instances==0){
		box_Darren_eval_points<TPOINTSTRI>::_nbr_instances++;
		box_Darren_eval_points<TPOINTSTRI>::get_t_points();
	}
}

	
	
//
// Explicit Instantiation of the init functions
//
	
template void bubble_Darren_eval_points<BUBPOINTS_STD,BUBYPOINTS_STD>::init();
template void triangle_Darren_eval_points<CTRIPOINTS_STD,TRIPOINTS_STD>::init();
template void box_Darren_eval_points<TRIPOINTS_STD>::init();	

template void bubble_Darren_eval_points<BUBPOINTS_HIGGS,BUBYPOINTS_HIGGS>::init();
template void triangle_Darren_eval_points<CTRIPOINTS_HIGGS,TRIPOINTS_HIGGS>::init();
template void box_Darren_eval_points<TRIPOINTS_HIGGS>::init();

template void bubble_Darren_eval_points<BUBCOEFFSEP,BUBCOEFFS>::init();
template void triangle_Darren_eval_points<TRICOEFFS,TRICOEFFSEP>::init();
template void box_Darren_eval_points<TRICOEFFSEP>::init();
	
	


}
}
}
