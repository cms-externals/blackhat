

/*
 * matrix elements for y -> 2q3g
 *
 *  Created on: Oct 21, 2008
 *
 */



#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include <ctime>
#include "BH_interface_impl.h"



#define _VERBOSE 0
#define _CLOCK 1

// various switches: photonZW:  0 for photon
// 				1 for photon+Z
// 				2 for Z->neutrinos only
// 				3 for W only
// #define _PHOTON_ONLY 1     // 0 photon only; 1 for all other   replaced by a setting
#define LC_includes_nf 1 //0 corresponde to W+3jet-paper setup. 1 is new setup for Z+3 and W+4 jets.

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {


partial_amplitude_cached* A_loop_2q_3g_1y_6_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const vector<ph_type> _ph_type, int color,QCDorder lo_or_nlo=nlo){
        partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

        int i1=ind.at(0);
        int i2=ind.at(1);
        int i3=ind.at(2);
        int i4=ind.at(3);
        int i5=ind.at(4);
        int i6=ind.at(5);
        

        ph_type h1=pro.p(1);
        ph_type h2=pro.p(2);
        ph_type h3=pro.p(3);
        ph_type h4=pro.p(4);
        ph_type h5=pro.p(5);
        ph_type h6=pro.p(6);
        

        process pro_1=process(h1,h5,h4,h3,h2,h6);
        vector<int> ind_1;
        ind_1.push_back(i1);
        ind_1.push_back(i5);
        ind_1.push_back(i4);
        ind_1.push_back(i3);
        ind_1.push_back(i2);
        ind_1.push_back(i6);


//--------------------- propagator correction
        prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i6,_ph_type);
        PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------
        multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
        multi_precision_fraction nm4_over_2(3,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(1,2), r5(-1,2*n_c*n_c);
        multi_precision_fraction r6(1,6), r7(-1,6*n_c*n_c);


if(color==1){
//Leading color
        PA->add(pro,leading_color,ind,1,1);
#if LC_includes_nf==1
	PA->add(pro,nf,ind,n_f,n_c);
#endif
        // Scheme shift
        PA->add_subtraction(pro,ind,r4,0);
        // renormalization
        PA->add_subtraction(pro,ind,r1*nm4_over_2,-1);
}
if(color==0){
//Full color
        PA->add(pro,leading_color,ind,1,1);
        PA->add(pro_1,sub_leading_color,ind_1,1,n_c*n_c); //SIGN
	PA->add(pro,nf,ind,n_f,n_c);
	// Scheme shift
        PA->add_subtraction(pro,ind,(r4+r5),0);
        // renormalization
        PA->add_subtraction(pro,ind,(r1+r2+r3)*nm4_over_2,-1);
}
if(color==2){
//Full minus Leading color
        PA->add(pro_1,sub_leading_color,ind_1,1,n_c*n_c); //SIGN
#if LC_includes_nf==0
	PA->add(pro,nf,ind,n_f,n_c);
#endif
        // Scheme shift
        PA->add_subtraction(pro,ind,r5,0);
        // renormalization
        PA->add_subtraction(pro,ind,(r2+r3)*nm4_over_2,-1);
}
        return PA;
}






// partial amplitudes not needed for leading color
// notice labels below correspond to A_7;3(1_q,4,5_qb,2,3) as in notes. eqn. 1.63

partial_amplitude_cached* A_loop_2q_3g_1y_6_3(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);


ph_type h1=pro.p(1);
ph_type h2=pro.p(2);
ph_type h3=pro.p(3);
ph_type h4=pro.p(4);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);


//permutations from eqn. 4.7 hep-ph/9409393
//permutations generated with mathematica
vector<int> ind_1, ind_2, ind_3, ind_4, ind_5, ind_6, ind_7, ind_8, ind_9, ind_10, ind_11, ind_12, ind_13, ind_14, ind_15, ind_16, ind_17, ind_18, ind_19, ind_20, ind_21, ind_22, ind_23, ind_24;

process pro_1=process(h1,h2,h3,h4,h5,h6);
ind_1.push_back(i1);
ind_1.push_back(i2);
ind_1.push_back(i3);
ind_1.push_back(i4);
ind_1.push_back(i5);
ind_1.push_back(i6);

process pro_2=process(h1,h2,h4,h3,h5,h6);
ind_2.push_back(i1);
ind_2.push_back(i2);
ind_2.push_back(i4);
ind_2.push_back(i3);
ind_2.push_back(i5);
ind_2.push_back(i6);

process pro_3=process(h1,h3,h2,h4,h5,h6);
ind_3.push_back(i1);
ind_3.push_back(i3);
ind_3.push_back(i2);
ind_3.push_back(i4);
ind_3.push_back(i5);
ind_3.push_back(i6);

process pro_4=process(h1,h3,h4,h2,h5,h6);
ind_4.push_back(i1);
ind_4.push_back(i3);
ind_4.push_back(i4);
ind_4.push_back(i2);
ind_4.push_back(i5);
ind_4.push_back(i6);

process pro_5=process(h1,h4,h2,h3,h5,h6);
ind_5.push_back(i1);
ind_5.push_back(i4);
ind_5.push_back(i2);
ind_5.push_back(i3);
ind_5.push_back(i5);
ind_5.push_back(i6);

process pro_6=process(h1,h4,h3,h2,h5,h6);
ind_6.push_back(i1);
ind_6.push_back(i4);
ind_6.push_back(i3);
ind_6.push_back(i2);
ind_6.push_back(i5);
ind_6.push_back(i6);

process pro_7=process(h1,h3,h4,h5,h2,h6);
ind_7.push_back(i1);
ind_7.push_back(i3);
ind_7.push_back(i4);
ind_7.push_back(i5);
ind_7.push_back(i2);
ind_7.push_back(i6);

process pro_8=process(h1,h4,h3,h5,h2,h6);
ind_8.push_back(i1);
ind_8.push_back(i4);
ind_8.push_back(i3);
ind_8.push_back(i5);
ind_8.push_back(i2);
ind_8.push_back(i6);

process pro_9=process(h1,h4,h5,h2,h3,h6);
ind_9.push_back(i1);
ind_9.push_back(i4);
ind_9.push_back(i5);
ind_9.push_back(i2);
ind_9.push_back(i3);
ind_9.push_back(i6);

process pro_10=process(h1,h2,h4,h5,h3,h6);
ind_10.push_back(i1);
ind_10.push_back(i2);
ind_10.push_back(i4);
ind_10.push_back(i5);
ind_10.push_back(i3);
ind_10.push_back(i6);

process pro_11=process(h1,h4,h2,h5,h3,h6);
ind_11.push_back(i1);
ind_11.push_back(i4);
ind_11.push_back(i2);
ind_11.push_back(i5);
ind_11.push_back(i3);
ind_11.push_back(i6);

process pro_12=process(h1,h4,h5,h3,h2,h6);
ind_12.push_back(i1);
ind_12.push_back(i4);
ind_12.push_back(i5);
ind_12.push_back(i3);
ind_12.push_back(i2);
ind_12.push_back(i6);

//---nf terms
/* give tadpoles or bubbles on external legs
process pro_13=process(h1,h5,h4,h3,h2,h6);
ind_13.push_back(i1);
ind_13.push_back(i5);
ind_13.push_back(i4);
ind_13.push_back(i3);
ind_13.push_back(i2);
ind_13.push_back(i6);

process pro_14=process(h1,h5,h3,h4,h2,h6);
ind_14.push_back(i1);
ind_14.push_back(i5);
ind_14.push_back(i3);
ind_14.push_back(i4);
ind_14.push_back(i2);
ind_14.push_back(i6);

process pro_15=process(h1,h5,h4,h2,h3,h6);
ind_15.push_back(i1);
ind_15.push_back(i5);
ind_15.push_back(i4);
ind_15.push_back(i2);
ind_15.push_back(i3);
ind_15.push_back(i6);

process pro_16=process(h1,h5,h2,h4,h3,h6);
ind_16.push_back(i1);
ind_16.push_back(i5);
ind_16.push_back(i2);
ind_16.push_back(i4);
ind_16.push_back(i3);
ind_16.push_back(i6);

process pro_17=process(h1,h5,h3,h2,h4,h6);
ind_17.push_back(i1);
ind_17.push_back(i5);
ind_17.push_back(i3);
ind_17.push_back(i2);
ind_17.push_back(i4);
ind_17.push_back(i6);

process pro_18=process(h1,h5,h2,h3,h4,h6);
ind_18.push_back(i1);
ind_18.push_back(i5);
ind_18.push_back(i2);
ind_18.push_back(i3);
ind_18.push_back(i4);
ind_18.push_back(i6);

process pro_19=process(h1,h2,h5,h4,h3,h6);
ind_19.push_back(i1);
ind_19.push_back(i2);
ind_19.push_back(i5);
ind_19.push_back(i4);
ind_19.push_back(i3);
ind_19.push_back(i6);

process pro_20=process(h1,h2,h5,h3,h4,h6);
ind_20.push_back(i1);
ind_20.push_back(i2);
ind_20.push_back(i5);
ind_20.push_back(i3);
ind_20.push_back(i4);
ind_20.push_back(i6);
*/
process pro_21=process(h1,h3,h2,h5,h4,h6);
ind_21.push_back(i1);
ind_21.push_back(i3);
ind_21.push_back(i2);
ind_21.push_back(i5);
ind_21.push_back(i4);
ind_21.push_back(i6);

process pro_22=process(h1,h3,h5,h4,h2,h6);
ind_22.push_back(i1);
ind_22.push_back(i3);
ind_22.push_back(i5);
ind_22.push_back(i4);
ind_22.push_back(i2);
ind_22.push_back(i6);

process pro_23=process(h1,h3,h5,h2,h4,h6);
ind_23.push_back(i1);
ind_23.push_back(i3);
ind_23.push_back(i5);
ind_23.push_back(i2);
ind_23.push_back(i4);
ind_23.push_back(i6);

process pro_24=process(h1,h2,h3,h5,h4,h6);
ind_24.push_back(i1);
ind_24.push_back(i2);
ind_24.push_back(i3);
ind_24.push_back(i5);
ind_24.push_back(i4);
ind_24.push_back(i6);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

PA->add(pro_1  ,leading_color     ,ind_1,1,1);
PA->add(pro_2  ,leading_color     ,ind_2,1,1);
PA->add(pro_3  ,leading_color     ,ind_3,1,1);
PA->add(pro_4  ,leading_color     ,ind_4,1,1);
PA->add(pro_5  ,leading_color     ,ind_5,1,1);
PA->add(pro_6  ,leading_color     ,ind_6,1,1);
PA->add(pro_7  ,sub_leading_color ,ind_7,1,1);
PA->add(pro_8  ,sub_leading_color ,ind_8,1,1);
PA->add(pro_9  ,sub_leading_color ,ind_9,1,1);
PA->add(pro_10 ,sub_leading_color ,ind_10,1,1);
PA->add(pro_11 ,sub_leading_color ,ind_11,1,1);
PA->add(pro_12 ,sub_leading_color ,ind_12,1,1);

/*tadpole
PA->add(pro_13,nf,ind_13,-1*n_f,n_c);
PA->add(pro_14,nf,ind_14,-1*n_f,n_c);
PA->add(pro_15,nf,ind_15,-1*n_f,n_c);
PA->add(pro_16,nf,ind_16,-1*n_f,n_c);
PA->add(pro_17,nf,ind_17,-1*n_f,n_c);
PA->add(pro_18,nf,ind_18,-1*n_f,n_c);
*/
/* bubble on external
PA->add(pro_19,nf,ind_19,-1*n_f,n_c);
PA->add(pro_20,nf,ind_20,-1*n_f,n_c);
*/

/* zero sum contribution
PA->add(pro_21,nf,ind_21,-1*n_f,n_c);
PA->add(pro_22,nf,ind_22,-1*n_f,n_c);
PA->add(pro_23,nf,ind_23,-1*n_f,n_c);
PA->add(pro_24,nf,ind_24,-1*n_f,n_c);
*/

return PA;

};

partial_amplitude_cached* A_loop_2q_3g_1y_6_4(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

//permutations from eqn. 4.7 hep-ph/9409393
//permutations generated with mathematica
vector<int> ind_1, ind_2, ind_3, ind_4, ind_5, ind_6, ind_7, ind_8, ind_9, ind_10, ind_11, ind_12, ind_13, ind_14, ind_15, ind_16, ind_17, ind_18, ind_19, ind_20, ind_21, ind_22, ind_23, ind_24;


int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);


ph_type h1=pro.p(1);
ph_type h2=pro.p(2);
ph_type h3=pro.p(3);
ph_type h4=pro.p(4);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);


//permutations from eqn. 4.7 hep-ph/9409393
//permutations generated with mathematica

process pro_1=process(h1,h2,h4,h3,h5,h6);
ind_1.push_back(i1);
ind_1.push_back(i2);
ind_1.push_back(i4);
ind_1.push_back(i3);
ind_1.push_back(i5);
ind_1.push_back(i6);

process pro_2=process(h1,h3,h2,h4,h5,h6);
ind_2.push_back(i1);
ind_2.push_back(i3);
ind_2.push_back(i2);
ind_2.push_back(i4);
ind_2.push_back(i5);
ind_2.push_back(i6);

process pro_3=process(h1,h4,h3,h2,h5,h6);
ind_3.push_back(i1);
ind_3.push_back(i4);
ind_3.push_back(i3);
ind_3.push_back(i2);
ind_3.push_back(i5);
ind_3.push_back(i6);

process pro_4=process(h1,h4,h3,h5,h2,h6);
ind_4.push_back(i1);
ind_4.push_back(i4);
ind_4.push_back(i3);
ind_4.push_back(i5);
ind_4.push_back(i2);
ind_4.push_back(i6);

process pro_5=process(h1,h3,h5,h2,h4,h6);
ind_5.push_back(i1);
ind_5.push_back(i3);
ind_5.push_back(i5);
ind_5.push_back(i2);
ind_5.push_back(i4);
ind_5.push_back(i6);

process pro_6=process(h1,h5,h2,h4,h3,h6);
ind_6.push_back(i1);
ind_6.push_back(i5);
ind_6.push_back(i2);
ind_6.push_back(i4);
ind_6.push_back(i3);
ind_6.push_back(i6);

process pro_7=process(h1,h2,h4,h5,h3,h6);
ind_7.push_back(i1);
ind_7.push_back(i2);
ind_7.push_back(i4);
ind_7.push_back(i5);
ind_7.push_back(i3);
ind_7.push_back(i6);

process pro_8=process(h1,h4,h5,h3,h2,h6);
ind_8.push_back(i1);
ind_8.push_back(i4);
ind_8.push_back(i5);
ind_8.push_back(i3);
ind_8.push_back(i2);
ind_8.push_back(i6);

process pro_9=process(h1,h5,h3,h2,h4,h6);
ind_9.push_back(i1);
ind_9.push_back(i5);
ind_9.push_back(i3);
ind_9.push_back(i2);
ind_9.push_back(i4);
ind_9.push_back(i6);

process pro_10=process(h1,h3,h2,h5,h4,h6);
ind_10.push_back(i1);
ind_10.push_back(i3);
ind_10.push_back(i2);
ind_10.push_back(i5);
ind_10.push_back(i4);
ind_10.push_back(i6);

process pro_11=process(h1,h2,h5,h4,h3,h6);
ind_11.push_back(i1);
ind_11.push_back(i2);
ind_11.push_back(i5);
ind_11.push_back(i4);
ind_11.push_back(i3);
ind_11.push_back(i6);

process pro_12=process(h1,h5,h4,h3,h2,h6);
ind_12.push_back(i1);
ind_12.push_back(i5);
ind_12.push_back(i4);
ind_12.push_back(i3);
ind_12.push_back(i2);
ind_12.push_back(i6);

/* tp
process pro_13=process(h1,h5,h3,h4,h2,h6);
ind_13.push_back(i1);
ind_13.push_back(i5);
ind_13.push_back(i3);
ind_13.push_back(i4);
ind_13.push_back(i2);
ind_13.push_back(i6);

process pro_14=process(h1,h5,h4,h2,h3,h6);
ind_14.push_back(i1);
ind_14.push_back(i5);
ind_14.push_back(i4);
ind_14.push_back(i2);
ind_14.push_back(i3);
ind_14.push_back(i6);

process pro_15=process(h1,h5,h2,h3,h4,h6);
ind_15.push_back(i1);
ind_15.push_back(i5);
ind_15.push_back(i2);
ind_15.push_back(i3);
ind_15.push_back(i4);
ind_15.push_back(i6);
*/
/*boe
process pro_16=process(h1,h2,h5,h3,h4,h6);
ind_16.push_back(i1);
ind_16.push_back(i2);
ind_16.push_back(i5);
ind_16.push_back(i3);
ind_16.push_back(i4);
ind_16.push_back(i6);
*/

process pro_17=process(h1,h4,h2,h5,h3,h6);
ind_17.push_back(i1);
ind_17.push_back(i4);
ind_17.push_back(i2);
ind_17.push_back(i5);
ind_17.push_back(i3);
ind_17.push_back(i6);

process pro_18=process(h1,h3,h4,h2,h5,h6);
ind_18.push_back(i1);
ind_18.push_back(i3);
ind_18.push_back(i4);
ind_18.push_back(i2);
ind_18.push_back(i5);
ind_18.push_back(i6);
/* boe
process pro_19=process(h1,h3,h5,h4,h2,h6);
ind_19.push_back(i1);
ind_19.push_back(i3);
ind_19.push_back(i5);
ind_19.push_back(i4);
ind_19.push_back(i2);
ind_19.push_back(i6);
*/
process pro_20=process(h1,h2,h3,h5,h4,h6);
ind_20.push_back(i1);
ind_20.push_back(i2);
ind_20.push_back(i3);
ind_20.push_back(i5);
ind_20.push_back(i4);
ind_20.push_back(i6);

process pro_21=process(h1,h4,h2,h3,h5,h6);
ind_21.push_back(i1);
ind_21.push_back(i4);
ind_21.push_back(i2);
ind_21.push_back(i3);
ind_21.push_back(i5);
ind_21.push_back(i6);
/* boe
process pro_22=process(h1,h4,h5,h2,h3,h6);
ind_22.push_back(i1);
ind_22.push_back(i4);
ind_22.push_back(i5);
ind_22.push_back(i2);
ind_22.push_back(i3);
ind_22.push_back(i6);
*/
process pro_23=process(h1,h3,h4,h5,h2,h6);
ind_23.push_back(i1);
ind_23.push_back(i3);
ind_23.push_back(i4);
ind_23.push_back(i5);
ind_23.push_back(i2);
ind_23.push_back(i6);

process pro_24=process(h1,h2,h3,h4,h5,h6);
ind_24.push_back(i1);
ind_24.push_back(i2);
ind_24.push_back(i3);
ind_24.push_back(i4);
ind_24.push_back(i5);
ind_24.push_back(i6);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

PA->add(pro_1  ,leading_color,ind_1,-1,1);
PA->add(pro_2  ,leading_color,ind_2,-1,1);
PA->add(pro_3  ,leading_color,ind_3,-1,1);
PA->add(pro_4  ,sub_leading_color,ind_4,-1,1);
PA->add(pro_5  ,sub_leading_color,ind_5,-1,1);
PA->add(pro_6  ,sub_leading_color,ind_6,-1,1);
PA->add(pro_7  ,sub_leading_color,ind_7,-1,1);
PA->add(pro_8  ,sub_leading_color,ind_8,-1,1);
PA->add(pro_9  ,sub_leading_color,ind_9,-1,1);
PA->add(pro_10 ,sub_leading_color,ind_10,-1,1);
PA->add(pro_11 ,sub_leading_color,ind_11,-1,1);
PA->add(pro_12 ,sub_leading_color,ind_12,-1,1);

// tadpoles of bubbles on external legs
//PA->add(pro_13 ,nf,ind_13,-n_f,n_c);
//PA->add(pro_14 ,nf,ind_14,-n_f,n_c);
//PA->add(pro_15 ,nf,ind_15,-n_f,n_c);
//PA->add(pro_16 ,nf,ind_16,-n_f,n_c);
///////////////////////////////////

PA->add(pro_17 ,nf,ind_17,-n_f,n_c);
PA->add(pro_18 ,nf,ind_18,-n_f,n_c);
///////////////////////////////////
// tadpoles of bubbles on external legs
//PA->add(pro_19 ,nf,ind_19,-n_f,n_c);
///////////////////////////////////
PA->add(pro_20 ,nf,ind_20,-n_f,n_c);
PA->add(pro_21 ,nf,ind_21,-n_f,n_c);
///////////////////////////////////
// tadpoles of bubbles on external legs
//PA->add(pro_22 ,nf,ind_22,-n_f,n_c);
///////////////////////////////////
PA->add(pro_23 ,nf,ind_23,-n_f,n_c);
PA->add(pro_24 ,nf,ind_24,-n_f,n_c);
///////////////////////////////////

return PA;

};





Squared_ME* A_loop_2q_3g_1y_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, int color, int tree_color, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);


multi_precision_fraction n2m1(n_c*n_c-1);
multi_precision_fraction n2m2(n_c*n_c-2);
multi_precision_fraction n2p1(n_c*n_c+1);
multi_precision_fraction _over_nc(1,n_c);
multi_precision_fraction _8(8);
multi_precision_fraction _4(4);
multi_precision_fraction _2(2);
multi_precision_fraction ext=_4*n2m1*_over_nc;
multi_precision_fraction ext_tree=_2*n2m1*_over_nc*_over_nc;



for(int l1 = 1; l1 < 4;++l1){
	for(int l2 = 1; l2 < 4; ++l2){
		for(int l3 = 1; l3 < 4; ++l3){
			if((l1!=l2)&&(l1!=l3)&&(l2!=l3)){


int i1=ind.at(0);
int i2=ind.at(l1);
int i3=ind.at(l2);
int i4=ind.at(l3);
int i5=ind.at(4);
int i6=ind.at(5);


ph_type h1=pro.p(1);
ph_type h2=pro.p(l1+1);
ph_type h3=pro.p(l2+1);
ph_type h4=pro.p(l3+1);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);



process pro_12345=process(h1,h2,h3,h4,h5,h6);
vector<int> ind_12345;
ind_12345.push_back(i1);
ind_12345.push_back(i2);
ind_12345.push_back(i3);
ind_12345.push_back(i4);
ind_12345.push_back(i5);
ind_12345.push_back(i6);

process pro_13245=process(h1,h3,h2,h4,h5,h6);
vector<int> ind_13245;
ind_13245.push_back(i1);
ind_13245.push_back(i3);
ind_13245.push_back(i2);
ind_13245.push_back(i4);
ind_13245.push_back(i5);
ind_13245.push_back(i6);

process pro_12435=process(h1,h2,h4,h3,h5,h6);
vector<int> ind_12435;
ind_12435.push_back(i1);
ind_12435.push_back(i2);
ind_12435.push_back(i4);
ind_12435.push_back(i3);
ind_12435.push_back(i5);
ind_12435.push_back(i6);

process pro_14325=process(h1,h4,h3,h2,h5,h6);
vector<int> ind_14325;
ind_14325.push_back(i1);
ind_14325.push_back(i4);
ind_14325.push_back(i3);
ind_14325.push_back(i2);
ind_14325.push_back(i5);
ind_14325.push_back(i6);

process pro_13425=process(h1,h3,h4,h2,h5,h6);
vector<int> ind_13425;
ind_13425.push_back(i1);
ind_13425.push_back(i3);
ind_13425.push_back(i4);
ind_13425.push_back(i2);
ind_13425.push_back(i5);
ind_13425.push_back(i6);

process pro_14235=process(h1,h4,h2,h3,h5,h6);
vector<int> ind_14235;
ind_14235.push_back(i1);
ind_14235.push_back(i4);
ind_14235.push_back(i2);
ind_14235.push_back(i3);
ind_14235.push_back(i5);
ind_14235.push_back(i6);

//----in addition needed for A_7;3
// position 4 in call to A_7;3 is the isolated leg
// example: A_loop_2q_3g_1y_6_3(1,2,3,4,5)~A_7;3(1,4,5,2,3)

process pro_12543=process(h1,h4,h3,h2,h5,h6);
vector<int> ind_12543;
ind_12543.push_back(i1);
ind_12543.push_back(i4);
ind_12543.push_back(i3);
ind_12543.push_back(i2);
ind_12543.push_back(i5);
ind_12543.push_back(i6);

process pro_13542=process(h1,h4,h2,h3,h5,h6);
vector<int> ind_13542;
ind_13542.push_back(i1);
ind_13542.push_back(i4);
ind_13542.push_back(i2);
ind_13542.push_back(i3);
ind_13542.push_back(i5);
ind_13542.push_back(i6);

process pro_14532=process(h1,h3,h2,h4,h5,h6);
vector<int> ind_14532;
ind_14532.push_back(i1);
ind_14532.push_back(i3);
ind_14532.push_back(i2);
ind_14532.push_back(i4);
ind_14532.push_back(i5);
ind_14532.push_back(i6);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i6,_ph_type);
//-------------------------------------------

size_t T_12345=SM->add(new CTree_with_prefactor(pro_12345,ind_12345  ,_prop_hel_fn));
SM->add_tree(cached_cross_term_md(T_12345,T_12345, ext_tree*n2m1*n2m1));

size_t L1_12345=SM->add(A_loop_2q_3g_1y_6_1(pro_12345  ,ind_12345  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));
SM->add_loop(cached_cross_term_md(L1_12345,T_12345, ext*n2m1*n2m1));



size_t T_13245;
size_t T_12435;
size_t T_14325;
size_t T_13425;
size_t T_14235;


if(tree_color==0){ 
T_13245=SM->add(new CTree_with_prefactor(pro_13245,ind_13245  ,_prop_hel_fn));
T_12435=SM->add(new CTree_with_prefactor(pro_12435,ind_12435  ,_prop_hel_fn));
T_14325=SM->add(new CTree_with_prefactor(pro_14325,ind_14325  ,_prop_hel_fn));
T_13425=SM->add(new CTree_with_prefactor(pro_13425,ind_13425  ,_prop_hel_fn));
T_14235=SM->add(new CTree_with_prefactor(pro_14235,ind_14235  ,_prop_hel_fn));
}

if((color==0)||(color==2)){
//use full color form of A71:
size_t L1_13245=SM->add(A_loop_2q_3g_1y_6_1(pro_13245  ,ind_13245  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,0,lo_or_nlo));
size_t L1_12435=SM->add(A_loop_2q_3g_1y_6_1(pro_12435  ,ind_12435  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,0,lo_or_nlo));
size_t L1_14325=SM->add(A_loop_2q_3g_1y_6_1(pro_14325  ,ind_14325  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,0,lo_or_nlo));
size_t L1_13425=SM->add(A_loop_2q_3g_1y_6_1(pro_13425  ,ind_13425  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,0,lo_or_nlo));
size_t L1_14235=SM->add(A_loop_2q_3g_1y_6_1(pro_14235  ,ind_14235  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,0,lo_or_nlo));

size_t L3_12543=SM->add(A_loop_2q_3g_1y_6_3(pro_12543  ,ind_12543  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
size_t L3_14532=SM->add(A_loop_2q_3g_1y_6_3(pro_14532  ,ind_14532  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
size_t L3_13542=SM->add(A_loop_2q_3g_1y_6_3(pro_13542  ,ind_13542  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));

size_t L4_12345=SM->add(A_loop_2q_3g_1y_6_4(pro_12345,ind_12345,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
size_t L4_14325=SM->add(A_loop_2q_3g_1y_6_4(pro_14325,ind_14325,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));

SM->add_loop(cached_cross_term_md(L1_13245,T_12345,-ext*n2m1     ));
SM->add_loop(cached_cross_term_md(L1_12435,T_12345,-ext*n2m1     ));
SM->add_loop(cached_cross_term_md(L1_14325,T_12345, ext*n2p1     ));
SM->add_loop(cached_cross_term_md(L1_13425,T_12345, ext          ));
SM->add_loop(cached_cross_term_md(L1_14235,T_12345, ext          ));

SM->add_loop(cached_cross_term_md(L3_12543,T_12345,ext*n2m1      ));
SM->add_loop(cached_cross_term_md(L3_14532,T_12345,ext*n2m1      ));
SM->add_loop(cached_cross_term_md(L3_13542,T_12345,-ext          ));

SM->add_loop(cached_cross_term_md(L4_12345,T_12345,ext*n2m2       ));
SM->add_loop(cached_cross_term_md(L4_14325,T_12345,-2*ext        ));

}

if(tree_color==0){
SM->add_tree(cached_cross_term_md(T_13245,T_12345,-ext_tree*n2m1  ));
SM->add_tree(cached_cross_term_md(T_12435,T_12345,-ext_tree*n2m1  ));
SM->add_tree(cached_cross_term_md(T_14325,T_12345, ext_tree*n2p1  ));
SM->add_tree(cached_cross_term_md(T_13425,T_12345, ext_tree       ));
SM->add_tree(cached_cross_term_md(T_14235,T_12345, ext_tree       ));
}

// AXIAL, Vectorial

};
};
};
};

return SM;
}



Virtual_SME* vsme_2q3g1y(const std::vector<int>& indext,int ns,int nf,int nc, bool up_down_quark, int photonZW, int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);


	ind_b.push_back(indext[4]);
	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[5]);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.


clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;

	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,p,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	//VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,m,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,m,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,m,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,p,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,p,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,p,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,m,qbm,ym),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	//VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

	//VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,p,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,m,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,m,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,m,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,p,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,p,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,p,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,m,qbm,yp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	
	//VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,m,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,m,p,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,m,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_3g_1y_M2(process(qp,p,p,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));


after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}


BH_Ampl_2q3g1y::BH_Ampl_2q3g1y(bool up_down_quark, int photonZW, int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
		BH_Ampl_impl(
				vsme_2q3g1y(mom_assignment,
					0,
					settings::BH_interface_settings::s_nf,
					settings::BH_interface_settings::s_nc,
					up_down_quark,
					photonZW*settings::BH_interface_settings::s_photon_only,
					color,
					tree_color,
					lo_or_nlo),
				bhi,  // parent
				6,    // NbrExtParticles
				3,	  // NbrPowersOfAlphaS
				1,    // NbrPowersOfAlphaQED
				4,   // GeVdim
				1.  // factor
		),
		momenta_assignment(mom_assignment), d_color(color)
{}


double BH_Ampl_2q3g1y::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(d_color){
		//Full color
		case 0 : return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
		//Leading color
		case 1: return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* 33/6. /*beta_0*/ *get_double_pole()/(-12.) /* = lc born tree/full tree */  );
		//Full-Leading color
		case 2: return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* ( 23/6. -  33/6. /*beta_0*/ *(35/36. - get_double_pole()/(-12.)  /* = lc born tree */ ))  );
		}
	case 2:	return get_double_pole();
	}

};



}
