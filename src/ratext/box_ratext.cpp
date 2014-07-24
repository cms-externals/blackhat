/*
 * box_ratext.cpp
 *
 *  Created on: 27 Aug 2009
 *      Author: darrenforde
 */


#include "BH_utilities.h"
#include "rec_tree.h"
#include <iostream>
#include <cassert>
#include <string>
#include "box_ratext.h"
#include "triangle_ratext.h"
#ifndef BH_PUBLIC
#include "ratext/basecutRat.h"
#endif
#include "ratext/rat_worker.h"
#include "polylog.h"  // for pi

#define _VERBOSE 0
#define _MU_POINTS_RADIUS 1 // Sets the radius of the circle we use to evalaute the mu points

using namespace std;

namespace BH {

namespace ratext {


template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> const bool box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points(){

    // If we have the same number of points for box, triangle and bubble then we can use a simple form of the eval points
    if((MUPOINTS_T==MUTRIPOINTS_T)&&(MUTRIPOINTS_T==MUBUBPOINTS_T)){
        CVHP rotfac_VHP=RVHP(_MU_POINTS_RADIUS)*exp(CVHP(0,1)*pi<RVHP>()/RVHP(2*MUBUBPOINTS_T));
#if BH_USE_GMP
        CGMP rotfac_GMP=exp(CGMP(0,1)*pi<RGMP>()/RGMP(2*MUBUBPOINTS_T));
#endif
        for(int i=0;i<MUBUBPOINTS_T;i++){
#if BH_USE_GMP
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i].real().set_prec(RGMP::get_current_precision());
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i].imag().set_prec(RGMP::get_current_precision());
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i]=rotfac_GMP*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(MUBUBPOINTS_T));
#endif
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]=rotfac_VHP*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(MUBUBPOINTS_T));
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_HP[i]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
        }
        
        
        CVHP extpt_VHP=box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[MUPOINTS_T-1];
#if BH_USE_GMP
        CGMP extpt_GMP=box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[MUPOINTS_T-1];
#endif
        for(int j=0;j<MUPOINTS_T;j++){
            for(int i=0;i<MUPOINTS_T;i++){
#if BH_USE_GMP
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*MUPOINTS_T].real().set_prec(RGMP::get_current_precision());
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*MUPOINTS_T].imag().set_prec(RGMP::get_current_precision());
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*MUPOINTS_T]=RGMP(1)/(RGMP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i],j));
#endif
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]=RVHP(1)/(RVHP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i],j));
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[i+j*MUPOINTS_T]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]);
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*MUPOINTS_T]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]);
            }
        }
    }
    else{
        // We should only ever have at most one power of difference between the power of mu^2 in the box/triangle/bubble the number given by the power of mu in the bubble
        // and so we construct the points on the circle for this reduced number of mu^2 and then we add another point which we pick to be between two other points to
        // deal with the box and for every other power of l, the triangle
        CVHP rotfac_VHP=RVHP(_MU_POINTS_RADIUS)*exp(CVHP(0,1)*pi<RVHP>()/RVHP(2*MUBUBPOINTS_T));
#if BH_USE_GMP
        CGMP rotfac_GMP=exp(CGMP(0,1)*pi<RGMP>()/RGMP(2*MUBUBPOINTS_T));
#endif
        for(int i=0;i<MUBUBPOINTS_T;i++){
#if BH_USE_GMP
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i]=rotfac_GMP*exp(CGMP(0,2)*pi<RGMP>()*RGMP(i)/RGMP(MUBUBPOINTS_T));
#endif
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]=rotfac_VHP*exp(CVHP(0,2)*pi<RVHP>()*RVHP(i)/RVHP(MUBUBPOINTS_T));
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_HP[i]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[i]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i]);
        }
        // Now add the final point on the circle
#if BH_USE_GMP
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[MUPOINTS_T-1]=rotfac_GMP*exp(CGMP(0,1)*pi<RGMP>()/RGMP(MUBUBPOINTS_T));
#endif
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[MUPOINTS_T-1]=rotfac_VHP*exp(CVHP(0,1)*pi<RVHP>()/RVHP(MUBUBPOINTS_T));
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_HP[MUPOINTS_T-1]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[MUPOINTS_T-1]);
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[MUPOINTS_T-1]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[MUPOINTS_T-1]);

        // Construct the matrix to compute the m^2 coefficients, the m^0 coefficient will be wrong but we do not need it anyway so that is not a problem
        //  the m^0 coefficient will be the m^0_coeff-m^(2n)_coeff and so we simply add in the last coefficient to get it
        CVHP extpt_VHP=box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[MUPOINTS_T-1];
#if BH_USE_GMP
        CGMP extpt_GMP=box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[MUPOINTS_T-1];
#endif
        for(int j=0;j<MUPOINTS_T-1;j++){
            for(int i=0;i<MUPOINTS_T-1;i++){
#if BH_USE_GMP
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*MUPOINTS_T]=RGMP(1)/(RGMP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i],j));
#endif
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]=RVHP(1)/(RVHP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i],j));
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[i+j*MUPOINTS_T]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]);
                box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*MUPOINTS_T]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T]);
            }
#if BH_USE_GMP
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[MUPOINTS_T-1+j*MUPOINTS_T]=CGMP(0,0);
#endif
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[MUPOINTS_T-1+j*MUPOINTS_T]=CVHP(0,0);
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[MUPOINTS_T-1+j*MUPOINTS_T]=CHP(0,0);
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[MUPOINTS_T-1+j*MUPOINTS_T]=C(0,0);
        }

        // We must first construct the final point of the last line in the inversion matrix by summing all the other terms on that would be on that line if there were no higher powers of mu^2
#if BH_USE_GMP
        CGMP finalpt_GMP(0,0);
#endif
        CVHP finalpt_VHP(0,0);
        for(int j=0;j<MUPOINTS_T-1;j++){
            finalpt_VHP-=RVHP(1)/(RVHP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[j],MUPOINTS_T-1));
#if BH_USE_GMP
            finalpt_GMP-=RGMP(1)/(RGMP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[j],MUPOINTS_T-1));
#endif
        }
#if BH_USE_GMP
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]=finalpt_GMP/(RGMP(1)+finalpt_GMP*pow(extpt_GMP,MUPOINTS_T-1));
#endif
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]=finalpt_VHP/(RVHP(1)+finalpt_VHP*pow(extpt_VHP,MUPOINTS_T-1));
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]=to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]);
        box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]=to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[MUPOINTS_T-1+(MUPOINTS_T-1)*MUPOINTS_T]);

        // The remaining terms on the final line (avoiding the final point) can now be computed
        for(int i=0;i<MUPOINTS_T-1;i++){
#if BH_USE_GMP
            CGMP subfac_GMP(0,0);
#endif
            CVHP subfac_VHP(0,0);
            for(int j=1;j<MUPOINTS_T-1;j++){
                subfac_VHP+=pow(extpt_VHP,j)*box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+j*MUPOINTS_T];
#if BH_USE_GMP
                subfac_GMP+=pow(extpt_GMP,j)*box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+j*MUPOINTS_T];
#endif
            }
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+(MUPOINTS_T-1)*MUPOINTS_T]
            =(RVHP(1)/(RVHP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_VHP[i],MUPOINTS_T-1))-finalpt_VHP*subfac_VHP)/(RVHP(1)+finalpt_VHP*pow(extpt_VHP,MUPOINTS_T-1));
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_HP[i+(MUPOINTS_T-1)*MUPOINTS_T]
            =to_HP(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+(MUPOINTS_T-1)*MUPOINTS_T]);
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+(MUPOINTS_T-1)*MUPOINTS_T]
            =to_double(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_VHP[i+(MUPOINTS_T-1)*MUPOINTS_T]);
#if BH_USE_GMP
            box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix_GMP[i+(MUPOINTS_T-1)*MUPOINTS_T]
            =(RGMP(1)/(RGMP(MUBUBPOINTS_T)*pow(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_GMP[i],MUPOINTS_T-1))-finalpt_GMP*subfac_GMP)/(RGMP(1)+finalpt_GMP*pow(extpt_GMP,MUPOINTS_T-1));
#endif
        }
    }
    
    
#if _VERBOSE==1
    std::cout << "mu pts are : {";
    for(int j=0;j<MUPOINTS_T;j++){
        std::cout << box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[j] << ",";
    }
    std::cout << "}" << std::endl;
    
    std::cout << "mu matrix is : {";
    for(int j=0;j<MUPOINTS_T;j++){
        for(int i=0;i<MUPOINTS_T;i++){
            std::cout << box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[i+j*MUPOINTS_T] << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
#endif
    
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[0]=C(2.82842712474619,2.82842712474619);
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[1]=C(-2.82842712474619,-2.82842712474619);
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[2]=C(2.82842712474619,-2.82842712474619);
//
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[0]=C(0,0);
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[1]=C(0,0);
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[2]=C(0,0);
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[3]=C(1,0)/box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[0];
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[4]=C(1,0)/box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[1];
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[5]=C(1,0)/box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[2];
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[6]=-box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[0]/(R(2)*sqrt(R(MASS_IND))*R(MASS_IND));
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[7]=box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos[2]/(R(2)*sqrt(R(MASS_IND))*R(MASS_IND));
//	box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_circpos_matrix[8]=C(0,1)/R(MASS_IND);

#ifdef BH_USE_GMP
    s_GMP_precision=RGMP::get_current_precision();
#endif
	return true;
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  long int  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances=0;

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T> void box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::init(){
	//Setup the evaluation points if this has not been dont yet
	if(box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances==0){
		box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_nbr_instances++;
		box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_points();
	}
}

// Here we include the factor of 1/6 from the mu^4 term and the 1/2 from summing over the two solutions.
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const C  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int[4]={R(1)/C(0,12),C(1,0),C(1,0),C(1,0)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CHP  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_HP[4]={RHP(1)/CHP(0,12),CHP(1,0),CHP(1,0),CHP(1,0)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CVHP  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_VHP[4]={RVHP(1)/CVHP(0,12),CVHP(1,0),CVHP(1,0),CVHP(1,0)};

// A function for returning the appropriate rational integral function
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const C  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const C& s1){
	return _rat_int[element];
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CHP  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CHP& s1){
	return _rat_int_HP[element];
}

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CVHP  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CVHP& s1){
	return _rat_int_VHP[element];
}

#if BH_USE_GMP
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  CGMP  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::_rat_int_GMP[4]={RGMP(1)/CGMP(0,12),CGMP(1,0),CGMP(1,0),CGMP(1,0)};
template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  const CGMP box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::get_rat_integral(int element, const CGMP& s1){
	if (s_GMP_precision< RGMP::get_current_precision()){
		recomputeGMP();
	}
	return _rat_int_GMP[element];
}



#endif


//Explicitly Instantiate for all used cases

//
//The standard template
//
template <> C box_rat_eval_points<3,2,2,7,3,3>::_circpos[3]={C(0,0)};
template <> CHP box_rat_eval_points<3,2,2,7,3,3>::_circpos_HP[3]={CHP(0,0)};
template <> CVHP box_rat_eval_points<3,2,2,7,3,3>::_circpos_VHP[3]={CVHP(0,0)};

template <> C box_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix[3*3]={C(0,0)};
template <> CHP box_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_HP[3*3]={CHP(0,0)};
template <> CVHP box_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_VHP[3*3]={CVHP(0,0)};

template void box_rat_eval_points<3,2,2,7,3,3>::init();
template const C box_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const C& s1);
template const CHP box_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CHP& s1);
template const CVHP box_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CVHP& s1);

#if BH_USE_GMP
template <> CGMP box_rat_eval_points<3,2,2,7,3,3>::_circpos_GMP[3]={CGMP(0,0)};
template <> CGMP box_rat_eval_points<3,2,2,7,3,3>::_circpos_matrix_GMP[3*3]={CGMP(0,0)};
template const CGMP box_rat_eval_points<3,2,2,7,3,3>::get_rat_integral(int element, const CGMP& s1);
template <> int box_rat_eval_points<3,2,2,7,3,3>::s_GMP_precision=0;
#endif


template <class BoxType,class Specs,class T> struct initializer {
	void init(T,bool& k1massive,bool& k2massive);
};
#ifndef BH_PUBLIC
template <class Specs> struct initializer<baseboxRat,Specs,const boxRat_comp*> {
	static void init(const boxRat_comp* t,bool& k1massive,bool& k2massive){
		if(t->corner_size(1)>1||t->get_process(1).p(2).mass_label()>0){
			k1massive=true;
		}
		else{
			k1massive=false;
		}
		if(t->corner_size(4)>1||t->get_process(4).p(2).mass_label()>0){
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
		assert(title == "Box_Ratspecific");

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


template <class BoxType,class BoxSpecs> template <class T> box_Rat<BoxType,BoxSpecs>::box_Rat(T* t) : BoxType(t), BoxSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(4))
{
	initializer<BoxType,BoxSpecs,T*>::init(t,_k1massive,_k2massive);
	init();
}
template <class BoxType,class BoxSpecs> template <class T> box_Rat<BoxType,BoxSpecs>::box_Rat(T& t) : BoxType(t), BoxSpecs::CornerTreeStrategy(), _leg_masses(new mass_param_coll(4))
{
	initializer<BoxType,BoxSpecs,T&>::init(t,_k1massive,_k2massive);
	init();
}


#ifndef BH_PUBLIC
template BH::ratext::box_Rat<BH::baseboxRat, BH::ratext::Normal_RatBox_Specification<BH::baseboxRat> >::box_Rat(BH::boxRat_comp const*);
#endif
template BH::ratext::box_Rat<BH::ratext::rat_worker, BH::ratext::Normal_RatBox_Specification<BH::ratext::rat_worker> >::box_Rat(std::istream&);


//
//The Higgs template
//
template <> C box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos[MUBOXPOINTS_HIGGS]={C(0,0)};
template <> CHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_HP[MUBOXPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_VHP[MUBOXPOINTS_HIGGS]={CVHP(0,0)};

template <> C box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix[MUBOXPOINTS_HIGGS*MUBOXPOINTS_HIGGS]={C(0,0)};
template <> CHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_HP[MUBOXPOINTS_HIGGS*MUBOXPOINTS_HIGGS]={CHP(0,0)};
template <> CVHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_VHP[MUBOXPOINTS_HIGGS*MUBOXPOINTS_HIGGS]={CVHP(0,0)};

template void box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::init();
template const C box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const C& s1);
template const CHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CHP& s1);
template const CVHP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CVHP& s1);

#ifndef BH_PUBLIC
template BH::ratext::box_Rat<BH::baseboxRat, BH::ratext::Higgs_RatBox_Specification<BH::baseboxRat> >::box_Rat(BH::boxRat_comp const*);
#endif
template BH::ratext::box_Rat<BH::ratext::rat_worker, BH::ratext::Higgs_RatBox_Specification<BH::ratext::rat_worker> >::box_Rat(std::istream&);

#if BH_USE_GMP
template <> CGMP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_GMP[MUBOXPOINTS_HIGGS]={CGMP(0,0)};
template <> CGMP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::_circpos_matrix_GMP[MUBOXPOINTS_HIGGS*MUBOXPOINTS_HIGGS]={CGMP(0,0)};
template const CGMP box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::get_rat_integral(int element, const CGMP& s1);

template <> int box_rat_eval_points<MUBOXPOINTS_HIGGS,MUTRIPOINTS_HIGGS,MUBUBPOINTS_HIGGS,CPOINTS_HIGGS,CBUBPOINTS_HIGGS,YPOINTS_HIGGS>::s_GMP_precision=0;

template <int MUPOINTS_T, int MUTRIPOINTS_T, int MUBUBPOINTS_T, int CPOINTS_T, int CBUBPOINTS_T, int YPOINTS_T>  void  box_rat_eval_points<MUPOINTS_T,MUTRIPOINTS_T,MUBUBPOINTS_T,CPOINTS_T,CBUBPOINTS_T,YPOINTS_T>::recomputeGMP(){
	BH_DEBUG_MESSAGE2("recomputing mu_eval points for precision ",RGMP::get_current_precision());
	_rat_int_GMP[0]=RGMP(1)/CGMP(0,12);
	_rat_int_GMP[1]=CGMP(1,0);
	_rat_int_GMP[2]=CGMP(1,0);
	_rat_int_GMP[3]=CGMP(1,0);
	get_points();
};


#endif


}

}
