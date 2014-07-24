/*
 * matrix elements for 2q2Q1g2l
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

// various switches: photonZW:  0 for photon
// 				1 for photon+Z
// 				2 for Z->neutrinos only
// 				3 for W only
//#define _PHOTON_ONLY 1     // 0 photon only; 1 for all other replaced by a setting
#define LC_includes_nf 1 //0 corresponde to W+3jet-paper setup. 1 is new setup for Z+3 and W+4 jets.

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {


//expected indices:
// q g Gb G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_2l_7_1(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B3
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);

	process pro_12345=process(h1,h2,h3,h4,h5,h6,h7);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);
	ind_12345.push_back(i6);
	ind_12345.push_back(i7);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i7,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(3,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_12345,leading_color,ind_12345,1,1);
#if LC_includes_nf==1
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
#endif
}

if(color==0){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4+r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,leading_color,ind_12345,1,1);
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}

if(color==2){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r2+r3),-1);
// primitive amplitue
#if LC_includes_nf==0
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
#endif
}

if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(5-1);
	 i2=ind.at(1-1);
	 i3=ind.at(3-1);
	 i4=ind.at(4-1);
	 i5=ind.at(2-1);
	 i6=ind.at(6-1);
	 i7=ind.at(7-1);

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(3);
	 h4=pro.p(4);
	 h5=pro.p(2);
	 h6=pro.p(6);
	 h7=pro.p(7);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);


//a-type ~ sub_leading_color contributions:

	process pro_25341=process(h2,h5,h3,h4,h1,h6,h7);
	vector<int> ind_25341;
	ind_25341.push_back(i2);
	ind_25341.push_back(i5);
	ind_25341.push_back(i3);
	ind_25341.push_back(i4);
	ind_25341.push_back(i1);
	ind_25341.push_back(i6);
	ind_25341.push_back(i7);

	process pro_23415=process(h2,h3,h4,h1,h5,h6,h7);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);
	ind_23415.push_back(i6);
	ind_23415.push_back(i7);

	process pro_24315=process(h2,hp4,hp3,h1,h5,h6,h7);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);
	ind_24315.push_back(i6);
	ind_24315.push_back(i7);

	process pro_23541=process(h2,h3,h5,h4,h1,h6,h7);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);
	ind_23541.push_back(i6);
	ind_23541.push_back(i7);

	process pro_24531=process(h2,hp4,h5,hp3,h1,h6,h7);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);
	ind_24531.push_back(i6);
	ind_24531.push_back(i7);


	process pro_23451=process(h2,h3,h4,h5,h1,h6,h7);
	vector<int> ind_23451;
	ind_23451.push_back(i2);
	ind_23451.push_back(i3);
	ind_23451.push_back(i4);
	ind_23451.push_back(i5);
	ind_23451.push_back(i1);
	ind_23451.push_back(i6);
	ind_23451.push_back(i7);


	process pro_24351=process(h2,hp4,hp3,h5,h1,h6,h7);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);
	ind_24351.push_back(i6);
	ind_24351.push_back(i7);


	process pro_25431=process(h2,h5,hp4,hp3,h1,h6,h7);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);
	ind_25431.push_back(i6);
	ind_25431.push_back(i7);

	//leading_color: PA->add(pro_25341,sub_leading_color,ind_25341,1,1);
	PA->add(pro_25341,leading_color,ind_25341,-1,n_c*n_c);
	PA->add(pro_23415,sub_leading_color,ind_23415,1,n_c*n_c);
	PA->add(pro_24315,sub_leading_color,ind_24315,1,n_c*n_c);
	/*EGKMZ*/ PA->add(pro_23541,sub_leading_color,ind_23541,1,n_c*n_c);
	//EGKMZ PA->add(pro_24531,sub_leading_color,ind_24531,1,n_c*n_c);
	/* corr compared to EGKMZ */ PA->add(pro_24531,sub_leading_color,ind_24531,-1,n_c*n_c);
	PA->add(pro_23451,leading_color,ind_23451,1,n_c*n_c);
	PA->add(pro_24351,leading_color,ind_24351,1,n_c*n_c);
	PA->add(pro_25431,leading_color,ind_25431,-1,n_c*n_c);



//b-type ~ sub_leading_color contributions:

	process pro_21435=process(h2,h1,h4,h3,h5,h6,h7);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);
	ind_21435.push_back(i6);
	ind_21435.push_back(i7);

	PA->add(pro_21435,slc_q,ind_21435,1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	PA->add(pro_21435,slc_G,ind_21435,1,n_c*n_c);

}

	return PA;
}


//expected indices:
// q Gb G g qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_2l_7_3(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B1
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);

	process pro_12345=process(h1,h2,h3,h4,h5,h6,h7);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);
	ind_12345.push_back(i6);
	ind_12345.push_back(i7);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i7,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(3,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_12345,leading_color,ind_12345,1,1);
#if LC_includes_nf==1
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
#endif
}

if(color==0){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4+r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,leading_color,ind_12345,1,1);
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}

if(color==2){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r2+r3),-1);
// primitive amplitue
#if LC_includes_nf==0
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
#endif
}



if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(5-1);
	 i2=ind.at(1-1);
	 i3=ind.at(2-1);
	 i4=ind.at(3-1);
	 i5=ind.at(4-1);
	 i6=ind.at(6-1);
	 i7=ind.at(7-1);

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(3);
	 h5=pro.p(4);
	 h6=pro.p(6);
	 h7=pro.p(7);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:

	process pro_23451=process(h2,h3,h4,h5,h1,h6,h7);
	vector<int> ind_23451;
	ind_23451.push_back(i2);
	ind_23451.push_back(i3);
	ind_23451.push_back(i4);
	ind_23451.push_back(i5);
	ind_23451.push_back(i1);
	ind_23451.push_back(i6);
	ind_23451.push_back(i7);

	process pro_23415=process(h2,h3,h4,h1,h5,h6,h7);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);
	ind_23415.push_back(i6);
	ind_23415.push_back(i7);

	process pro_24315=process(h2,hp4,hp3,h1,h5,h6,h7);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);
	ind_24315.push_back(i6);
	ind_24315.push_back(i7);

	process pro_25341=process(h2,h5,h3,h4,h1,h6,h7);
	vector<int> ind_25341;
	ind_25341.push_back(i2);
	ind_25341.push_back(i5);
	ind_25341.push_back(i3);
	ind_25341.push_back(i4);
	ind_25341.push_back(i1);
	ind_25341.push_back(i6);
	ind_25341.push_back(i7);

	process pro_25431=process(h2,h5,hp4,hp3,h1,h6,h7);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);
	ind_25431.push_back(i6);
	ind_25431.push_back(i7);


	process pro_23541=process(h2,h3,h5,h4,h1,h6,h7);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);
	ind_23541.push_back(i6);
	ind_23541.push_back(i7);


	process pro_24531=process(h2,hp4,h5,hp3,h1,h6,h7);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);
	ind_24531.push_back(i6);
	ind_24531.push_back(i7);


	process pro_24351=process(h2,hp4,hp3,h5,h1,h6,h7);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);
	ind_24351.push_back(i6);
	ind_24351.push_back(i7);

	//leading_color: PA->add(pro_23451,sub_leading_color,ind_23451,1,1);
	PA->add(pro_23451,leading_color,ind_23451,-1,n_c*n_c);
	PA->add(pro_23415,sub_leading_color,ind_23415,1,n_c*n_c);
	PA->add(pro_24315,sub_leading_color,ind_24315,1,n_c*n_c);
	PA->add(pro_25341,leading_color,ind_25341,1,n_c*n_c); 
	PA->add(pro_25431,leading_color,ind_25431,1,n_c*n_c);
	PA->add(pro_23541,sub_leading_color,ind_23541,1,n_c*n_c);
	// EGKMZ PA->add(pro_24531,sub_leading_color,ind_24531,1,n_c*n_c);
	/*corr sign*/ PA->add(pro_24531,sub_leading_color,ind_24531,-1,n_c*n_c);
	PA->add(pro_24351,leading_color,ind_24351,-1,n_c*n_c);



//b-type ~ sub_leading_color contributions:

	process pro_21543=process(h2,h1,h5,h4,h3,h6,h7);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);
	ind_21543.push_back(i6);
	ind_21543.push_back(i7);

	PA->add(pro_21543,slc_q,ind_21543,1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	PA->add(pro_21543,slc_G,ind_21543,1,n_c*n_c);

}

	return PA;
}


//expected indices:
// q Gb g G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_2l_7_2(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B4
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);

	process pro_12345=process(h1,h2,h3,h4,h5,h6,h7);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);
	ind_12345.push_back(i6);
	ind_12345.push_back(i7);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i7,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(3,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!

if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*r1,-1);
// primitive amplitue
}
if(color==0){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4+r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}
if(color==2){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}


if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(5-1);
	 i2=ind.at(1-1);
	 i3=ind.at(2-1);
	 i4=ind.at(4-1);
	 i5=ind.at(3-1);
	 i6=ind.at(6-1);
	 i7=ind.at(7-1);

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(4);
	 h5=pro.p(3);
	 h6=pro.p(6);
	 h7=pro.p(7);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:

	process pro_24315=process(h2,hp4,hp3,h1,h5,h6,h7);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);
	ind_24315.push_back(i6);
	ind_24315.push_back(i7);

	process pro_25431=process(h2,h5,hp4,hp3,h1,h6,h7);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);
	ind_25431.push_back(i6);
	ind_25431.push_back(i7);

	process pro_24351=process(h2,hp4,hp3,h5,h1,h6,h7);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);
	ind_24351.push_back(i6);
	ind_24351.push_back(i7);

	process pro_24531=process(h2,hp4,h5,hp3,h1,h6,h7);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);
	ind_24531.push_back(i6);
	ind_24531.push_back(i7);

	process pro_23541=process(h2,h3,h5,h4,h1,h6,h7);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);
	ind_23541.push_back(i6);
	ind_23541.push_back(i7);

	PA->add(pro_24315,sub_leading_color,ind_24315,-1,1);
	PA->add(pro_25431,leading_color,ind_25431,-1,1);
	PA->add(pro_24351,leading_color,ind_24351,-1,1);
	PA->add(pro_23541,sub_leading_color,ind_23541,-1,n_c*n_c);
	// EGKMZ sign: PA->add(pro_24531,sub_leading_color,ind_24531,-1,n_c*n_c);
	/* corr sign:*/ PA->add(pro_24531,sub_leading_color,ind_24531,1,n_c*n_c);


//b-type ~ sub_leading_color contributions:

	process pro_21435=process(h2,h1,h4,h3,h5,h6,h7);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);
	ind_21435.push_back(i6);
	ind_21435.push_back(i7);

	process pro_21543=process(h2,h1,h5,h4,h3,h6,h7);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);
	ind_21543.push_back(i6);
	ind_21543.push_back(i7);

	process pro_25143=process(h2,h5,h1,h4,h3,h6,h7);
	vector<int> ind_25143;
	ind_25143.push_back(i2);
	ind_25143.push_back(i5);
	ind_25143.push_back(i1);
	ind_25143.push_back(i4);
	ind_25143.push_back(i3);
	ind_25143.push_back(i6);
	ind_25143.push_back(i7);

	process pro_21453=process(h2,h1,h4,h5,h3,h6,h7);
	vector<int> ind_21453;
	ind_21453.push_back(i2);
	ind_21453.push_back(i1);
	ind_21453.push_back(i4);
	ind_21453.push_back(i5);
	ind_21453.push_back(i3);
	ind_21453.push_back(i6);
	ind_21453.push_back(i7);

	PA->add(pro_21543,slc_q,ind_21543,-1,1);
	PA->add(pro_21435,slc_q,ind_21435,-1,1);
	PA->add(pro_25143,slc_q,ind_25143,-1,1);
	PA->add(pro_21453,slc_q,ind_21453,-1,1);
	PA->add(pro_21453,slc_q,ind_21453,1,n_c*n_c);
//c-type ~ sub_leading_color contributions:
	PA->add(pro_21543,slc_G,ind_21543,-1,n_c*n_c);
	PA->add(pro_21435,slc_G,ind_21435,-1,n_c*n_c);
	PA->add(pro_25143,slc_G,ind_25143,-1,n_c*n_c);
}

	return PA;
}


//expected indices:
// q Gb G qb g l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_2l_7_4(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type,int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B2
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);

	process pro_12345=process(h1,h2,h3,h4,h5,h6,h7);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);
	ind_12345.push_back(i6);
	ind_12345.push_back(i7);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i6,i7,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(3,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!

if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*r1,-1);
// primitive amplitue
}
if(color==0){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r4+r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}
if(color==2){
// SCHEME-SHIFT
        PA->add_subtraction(pro_12345,ind_12345,r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro_12345,ind_12345,nm4_over_2*(r2+r3),-1);
// primitive amplitue
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}

if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(4-1);
	 i2=ind.at(1-1);
	 i3=ind.at(2-1);
	 i4=ind.at(3-1);
	 i5=ind.at(5-1);
	 i6=ind.at(6-1);
	 i7=ind.at(7-1);

	 h1=pro.p(4);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(3);
	 h5=pro.p(5);
	 h6=pro.p(6);
	 h7=pro.p(7);


// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:


	process pro_25431=process(h2,h5,hp4,hp3,h1,h6,h7);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);
	ind_25431.push_back(i6);
	ind_25431.push_back(i7);


	process pro_24531=process(h2,hp4,h5,hp3,h1,h6,h7);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);
	ind_24531.push_back(i6);
	ind_24531.push_back(i7);


	process pro_24351=process(h2,hp4,hp3,h5,h1,h6,h7);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);
	ind_24351.push_back(i6);
	ind_24351.push_back(i7);


	process pro_23415=process(h2,h3,h4,h1,h5,h6,h7);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);
	ind_23415.push_back(i6);
	ind_23415.push_back(i7);

	process pro_24315=process(h2,hp4,hp3,h1,h5,h6,h7);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);
	ind_24315.push_back(i6);
	ind_24315.push_back(i7);

	PA->add(pro_25431,leading_color,ind_25431,1,1);
	PA->add(pro_24531,sub_leading_color,ind_24531,1,1); // relative sign compared ot EGKMZ
	PA->add(pro_24351,leading_color,ind_24351,1,1);
	PA->add(pro_23415,sub_leading_color,ind_23415,-1,n_c*n_c);
	PA->add(pro_24315,sub_leading_color,ind_24315,-1,n_c*n_c);


//b-type ~ sub_leading_color contributions:


	process pro_21543=process(h2,h1,h5,h4,h3,h6,h7);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);
	ind_21543.push_back(i6);
	ind_21543.push_back(i7);


	process pro_21453=process(h2,h1,h4,h5,h3,h6,h7);
	vector<int> ind_21453;
	ind_21453.push_back(i2);
	ind_21453.push_back(i1);
	ind_21453.push_back(i4);
	ind_21453.push_back(i5);
	ind_21453.push_back(i3);
	ind_21453.push_back(i6);
	ind_21453.push_back(i7);

	process pro_21435=process(h2,h1,h4,h3,h5,h6,h7);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);
	ind_21435.push_back(i6);
	ind_21435.push_back(i7);


	PA->add(pro_21543,slc_q,ind_21543,-1,n_c*n_c);
	PA->add(pro_21453,slc_q,ind_21453,-1,n_c*n_c);
	PA->add(pro_21435,slc_q,ind_21435,-1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	process pro_25143=process(h2,h5,h1,h4,h3,h6,h7);
	vector<int> ind_25143;
	ind_25143.push_back(i2);
	ind_25143.push_back(i5);
	ind_25143.push_back(i1);
	ind_25143.push_back(i4);
	ind_25143.push_back(i3);
	ind_25143.push_back(i6);
	ind_25143.push_back(i7);

	PA->add(pro_21543,slc_G,ind_21543,-1,1);
	PA->add(pro_21453,slc_G,ind_21453,-1,1);
	PA->add(pro_21435,slc_G,ind_21435,-1,1);
	PA->add(pro_25143,slc_G,ind_25143,-1,1);
	PA->add(pro_25143,slc_G,ind_25143,1,n_c*n_c);
}
	return PA;
}



Squared_ME* A_loop_2q_2Q_1g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int photonZW, const ph_type e_ph_type,int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);
int i7=ind.at(6);


// assume input with identical flavor quarks
// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors

//quarks labels i1-i4
//extra gluon called i5
//leptons i6-i7


ph_type h1=pro.p(1);
ph_type h2=particle_ID(gluino,pro.p(2).helicity(),1,true);
ph_type h3=particle_ID(gluino,pro.p(3).helicity(),1);
ph_type h4=pro.p(4);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);
ph_type h7=pro.p(7);

ph_type hp1=pro.p(1);
ph_type hp2=pro.p(2);
ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),1);
ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),1,true);
ph_type hp5=pro.p(5);
ph_type hp6=pro.p(6);
ph_type hp7=pro.p(7);

ph_type hf1=particle_ID(gluino,pro.p(1).helicity(),1);
ph_type hf2=pro.p(2);
ph_type hf3=pro.p(3);
ph_type hf4=particle_ID(gluino,pro.p(4).helicity(),1,true);

ph_type hpf1=particle_ID(gluino,pro.p(1).helicity(),1);
ph_type hpf2=particle_ID(gluino,pro.p(2).helicity(),1,true);
ph_type hpf3=pro.p(3);
ph_type hpf4=pro.p(4);

//-------------------------------------------------------
// tree labels

process pro_15234=process(h1,h5,h2,h3,h4,h6,h7);
vector<int> ind_15234;
ind_15234.push_back(i1);
ind_15234.push_back(i5);
ind_15234.push_back(i2);
ind_15234.push_back(i3);
ind_15234.push_back(i4);
ind_15234.push_back(i6);
ind_15234.push_back(i7);

process pro_12534=process(h1,h2,h5,h3,h4,h6,h7);
vector<int> ind_12534;
ind_12534.push_back(i1);
ind_12534.push_back(i2);
ind_12534.push_back(i5);
ind_12534.push_back(i3);
ind_12534.push_back(i4);
ind_12534.push_back(i6);
ind_12534.push_back(i7);

process pro_12354=process(h1,h2,h3,h5,h4,h6,h7);
vector<int> ind_12354;
ind_12354.push_back(i1);
ind_12354.push_back(i2);
ind_12354.push_back(i3);
ind_12354.push_back(i5);
ind_12354.push_back(i4);
ind_12354.push_back(i6);
ind_12354.push_back(i7);

process pro_12345=process(h1,h2,h3,h4,h5,h6,h7);
vector<int> ind_12345;
ind_12345.push_back(i1);
ind_12345.push_back(i2);
ind_12345.push_back(i3);
ind_12345.push_back(i4);
ind_12345.push_back(i5);
ind_12345.push_back(i6);
ind_12345.push_back(i7);

process pro_35412=process(hf3,h5,hf4,hf1,hf2,h6,h7);
vector<int> ind_35412;
ind_35412.push_back(i3);
ind_35412.push_back(i5);
ind_35412.push_back(i4);
ind_35412.push_back(i1);
ind_35412.push_back(i2);
ind_35412.push_back(i6);
ind_35412.push_back(i7);

process pro_34512=process(hf3,hf4,h5,hf1,hf2,h6,h7);
vector<int> ind_34512;
ind_34512.push_back(i3);
ind_34512.push_back(i4);
ind_34512.push_back(i5);
ind_34512.push_back(i1);
ind_34512.push_back(i2);
ind_34512.push_back(i6);
ind_34512.push_back(i7);

process pro_34152=process(hf3,hf4,hf1,h5,hf2,h6,h7);
vector<int> ind_34152;
ind_34152.push_back(i3);
ind_34152.push_back(i4);
ind_34152.push_back(i1);
ind_34152.push_back(i5);
ind_34152.push_back(i2);
ind_34152.push_back(i6);
ind_34152.push_back(i7);

process pro_34125=process(hf3,hf4,hf1,hf2,h5,h6,h7);
vector<int> ind_34125;
ind_34125.push_back(i3);
ind_34125.push_back(i4);
ind_34125.push_back(i1);
ind_34125.push_back(i2);
ind_34125.push_back(i5);
ind_34125.push_back(i6);
ind_34125.push_back(i7);

//---------------- second quark line arrangement for case of identical quarks

process pro_35214=process(hpf3,hp5,hpf2,hpf1,hpf4,hp6,hp7);
vector<int> ind_35214;
ind_35214.push_back(i3);
ind_35214.push_back(i5);
ind_35214.push_back(i2);
ind_35214.push_back(i1);
ind_35214.push_back(i4);
ind_35214.push_back(i6);
ind_35214.push_back(i7);

process pro_32514=process(hpf3,hpf2,hp5,hpf1,hpf4,hp6,hp7);
vector<int> ind_32514;
ind_32514.push_back(i3);
ind_32514.push_back(i2);
ind_32514.push_back(i5);
ind_32514.push_back(i1);
ind_32514.push_back(i4);
ind_32514.push_back(i6);
ind_32514.push_back(i7);

process pro_32154=process(hpf3,hpf2,hpf1,hp5,hpf4,hp6,hp7);
vector<int> ind_32154;
ind_32154.push_back(i3);
ind_32154.push_back(i2);
ind_32154.push_back(i1);
ind_32154.push_back(i5);
ind_32154.push_back(i4);
ind_32154.push_back(i6);
ind_32154.push_back(i7);

process pro_32145=process(hpf3,hpf2,hpf1,hpf4,hp5,hp6,hp7);
vector<int> ind_32145;
ind_32145.push_back(i3);
ind_32145.push_back(i2);
ind_32145.push_back(i1);
ind_32145.push_back(i4);
ind_32145.push_back(i5);
ind_32145.push_back(i6);
ind_32145.push_back(i7);

process pro_15432=process(hp1,hp5,hp4,hp3,hp2,hp6,hp7);
vector<int> ind_15432;
ind_15432.push_back(i1);
ind_15432.push_back(i5);
ind_15432.push_back(i4);
ind_15432.push_back(i3);
ind_15432.push_back(i2);
ind_15432.push_back(i6);
ind_15432.push_back(i7);

process pro_14532=process(hp1,hp4,hp5,hp3,hp2,hp6,hp7);
vector<int> ind_14532;
ind_14532.push_back(i1);
ind_14532.push_back(i4);
ind_14532.push_back(i5);
ind_14532.push_back(i3);
ind_14532.push_back(i2);
ind_14532.push_back(i6);
ind_14532.push_back(i7);

process pro_14352=process(hp1,hp4,hp3,hp5,hp2,hp6,hp7);
vector<int> ind_14352;
ind_14352.push_back(i1);
ind_14352.push_back(i4);
ind_14352.push_back(i3);
ind_14352.push_back(i5);
ind_14352.push_back(i2);
ind_14352.push_back(i6);
ind_14352.push_back(i7);

process pro_14325=process(hp1,hp4,hp3,hp2,hp5,hp6,hp7);
vector<int> ind_14325;
ind_14325.push_back(i1);
ind_14325.push_back(i4);
ind_14325.push_back(i3);
ind_14325.push_back(i2);
ind_14325.push_back(i5);
ind_14325.push_back(i6);
ind_14325.push_back(i7);




vector<ph_type> _ph_type_1;
_ph_type_1.push_back(h1);
_ph_type_1.push_back(e_ph_type);

vector<ph_type> _ph_type_3;
_ph_type_3.push_back(h3);
_ph_type_3.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn_up_1(0,photonZW,0,i6,i7,_ph_type_1);
	prop_hel_fn _prop_hel_fn_up_3(0,photonZW,0,i6,i7,_ph_type_3);
	prop_hel_fn _prop_hel_fn_down_1(1,photonZW,0,i6,i7,_ph_type_1);
	prop_hel_fn _prop_hel_fn_down_3(1,photonZW,0,i6,i7,_ph_type_3);
//-------------------------------------------



//-------------------------------------------


	size_t T_15234_up;
	size_t T_12354_up;
	size_t T_34152_up;
	size_t T_35412_up;

//if(tree_color!=1){
	size_t T_12534_up;
	size_t T_12345_up;
	size_t T_34512_up;
	size_t T_34125_up;
//}

	size_t L1_15234_up;
	size_t L3_12354_up;
	size_t L1_35412_up;
	size_t L3_34152_up;

	size_t L1_color_15234_up;
	size_t L3_color_12354_up;
	size_t L1_color_35412_up;
	size_t L3_color_34152_up;

//if(color!=1){
	size_t L2_12534_up;
	size_t L4_12345_up;
	size_t L2_34512_up;
	size_t L4_34125_up;
//}


	size_t T_15432_up;
	size_t T_14352_up;
	size_t T_35214_up;
	size_t T_32154_up;

//if(tree_color!=1){
	size_t T_14532_up;
	size_t T_14325_up;
	size_t T_32514_up;
	size_t T_32145_up;
//}

	size_t L1_15432_up;
	size_t L3_14352_up;
	size_t L1_35214_up;
	size_t L3_32154_up;

	size_t L1_color_15432_up;
	size_t L3_color_14352_up;
	size_t L1_color_35214_up;
	size_t L3_color_32154_up;

//if(color!=1){
	size_t L2_14532_up;
	size_t L4_14325_up;
	size_t L2_32514_up;
	size_t L4_32145_up;
//}

	size_t T_15234_down;
	size_t T_12354_down;
	size_t T_34152_down;
	size_t T_35412_down;

//if(tree_color!=1){
	size_t T_12534_down;
	size_t T_12345_down;
	size_t T_34512_down;
	size_t T_34125_down;
//}


	size_t L1_15234_down;
	size_t L3_12354_down;
	size_t L1_35412_down;
	size_t L3_34152_down;

	size_t L1_color_15234_down;
	size_t L3_color_12354_down;
	size_t L1_color_35412_down;
	size_t L3_color_34152_down;

//if(color!=1){
	size_t L2_12534_down;
	size_t L4_12345_down;
	size_t L2_34512_down;
	size_t L4_34125_down;
//}


	size_t T_15432_down;
	size_t T_14352_down;
	size_t T_35214_down;
	size_t T_32154_down;

//if(tree_color!=1){
	size_t T_14532_down;
	size_t T_14325_down;
	size_t T_32514_down;
	size_t T_32145_down;
//}



	size_t L1_15432_down;
	size_t L3_14352_down;
	size_t L1_35214_down;
	size_t L3_32154_down;

	size_t L1_color_15432_down;
	size_t L3_color_14352_down;
	size_t L1_color_35214_down;
	size_t L3_color_32154_down;

//if(color!=1){
	size_t L2_14532_down;
	size_t L4_14325_down;
	size_t L2_32514_down;
	size_t L4_32145_down;
//}

	if(h1.helicity()+h4.helicity()==0){
		T_15234_up=SM->add(new CTree_with_prefactor(pro_15234,ind_15234,_prop_hel_fn_up_1));
		T_12354_up=SM->add(new CTree_with_prefactor(pro_12354,ind_12354,_prop_hel_fn_up_1));
		T_35412_up=SM->add(new CTree_with_prefactor(pro_35412,ind_35412,_prop_hel_fn_up_3));
		T_34152_up=SM->add(new CTree_with_prefactor(pro_34152,ind_34152,_prop_hel_fn_up_3));

if(tree_color!=1){
		T_12534_up=SM->add(new CTree_with_prefactor(pro_12534,ind_12534,_prop_hel_fn_up_1));
		T_12345_up=SM->add(new CTree_with_prefactor(pro_12345,ind_12345,_prop_hel_fn_up_1));
		T_34512_up=SM->add(new CTree_with_prefactor(pro_34512,ind_34512,_prop_hel_fn_up_3));
		T_34125_up=SM->add(new CTree_with_prefactor(pro_34125,ind_34125,_prop_hel_fn_up_3));
}

if(color!=1){
		L1_15234_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15234,ind_15234,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L3_12354_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_12354,ind_12354,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L1_35412_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35412,ind_35412,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L3_34152_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_34152,ind_34152,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
}
		L1_color_15234_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15234,ind_15234,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_12354_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_12354,ind_12354,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L1_color_35412_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35412,ind_35412,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_34152_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_34152,ind_34152,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));

if(color!=1){
		L2_12534_up=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_12534,ind_12534,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L4_12345_up=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_12345,ind_12345,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L2_34512_up=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_34512,ind_34512,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L4_34125_up=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_34125,ind_34125,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
}

	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_15432_up=SM->add(new CTree_with_prefactor(pro_15432,ind_15432,_prop_hel_fn_up_1));
		T_14352_up=SM->add(new CTree_with_prefactor(pro_14352,ind_14352,_prop_hel_fn_up_1));
		T_35214_up=SM->add(new CTree_with_prefactor(pro_35214,ind_35214,_prop_hel_fn_up_3));
		T_32154_up=SM->add(new CTree_with_prefactor(pro_32154,ind_32154,_prop_hel_fn_up_3));

if(tree_color!=1){
		T_14532_up=SM->add(new CTree_with_prefactor(pro_14532,ind_14532,_prop_hel_fn_up_1));
		T_14325_up=SM->add(new CTree_with_prefactor(pro_14325,ind_14325,_prop_hel_fn_up_1));
		T_32514_up=SM->add(new CTree_with_prefactor(pro_32514,ind_32514,_prop_hel_fn_up_3));
		T_32145_up=SM->add(new CTree_with_prefactor(pro_32145,ind_32145,_prop_hel_fn_up_3));
}



if(color!=1){
		L1_15432_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15432,ind_15432,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L3_14352_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_14352,ind_14352,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L1_35214_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35214,ind_35214,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L3_32154_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_32154,ind_32154,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
}
		L1_color_15432_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15432,ind_15432,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_14352_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_14352,ind_14352,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L1_color_35214_up=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35214,ind_35214,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_32154_up=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_32154,ind_32154,n_s,n_f,n_c,0,photonZW,e_ph_type,color, lo_or_nlo));


if(color!=1){
		L2_14532_up=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_14532,ind_14532,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L4_14325_up=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_14325,ind_14325,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L2_32514_up=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_32514,ind_32514,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
		L4_32145_up=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_32145,ind_32145,n_s,n_f,n_c,0,photonZW,e_ph_type,0, lo_or_nlo));
}

	};


	if(h1.helicity()+h4.helicity()==0){
		T_15234_down=SM->add(new CTree_with_prefactor(pro_15234,ind_15234,_prop_hel_fn_down_1));
		T_12354_down=SM->add(new CTree_with_prefactor(pro_12354,ind_12354,_prop_hel_fn_down_1));
		T_35412_down=SM->add(new CTree_with_prefactor(pro_35412,ind_35412,_prop_hel_fn_down_3));
		T_34152_down=SM->add(new CTree_with_prefactor(pro_34152,ind_34152,_prop_hel_fn_down_3));

if(tree_color!=1){
		T_12534_down=SM->add(new CTree_with_prefactor(pro_12534,ind_12534,_prop_hel_fn_down_1));
		T_12345_down=SM->add(new CTree_with_prefactor(pro_12345,ind_12345,_prop_hel_fn_down_1));
		T_34512_down=SM->add(new CTree_with_prefactor(pro_34512,ind_34512,_prop_hel_fn_down_3));
		T_34125_down=SM->add(new CTree_with_prefactor(pro_34125,ind_34125,_prop_hel_fn_down_3));
}

if(color!=1){
		L1_15234_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15234,ind_15234,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L3_12354_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_12354,ind_12354,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L1_35412_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35412,ind_35412,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L3_34152_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_34152,ind_34152,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
}
		L1_color_15234_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15234,ind_15234,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_12354_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_12354,ind_12354,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L1_color_35412_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35412,ind_35412,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_34152_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_34152,ind_34152,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
if(color!=1){
		L2_12534_down=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_12534,ind_12534,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L4_12345_down=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_12345,ind_12345,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L2_34512_down=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_34512,ind_34512,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L4_34125_down=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_34125,ind_34125,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
}

	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_15432_down=SM->add(new CTree_with_prefactor(pro_15432,ind_15432,_prop_hel_fn_down_1));
		T_14352_down=SM->add(new CTree_with_prefactor(pro_14352,ind_14352,_prop_hel_fn_down_1));
		T_35214_down=SM->add(new CTree_with_prefactor(pro_35214,ind_35214,_prop_hel_fn_down_3));
		T_32154_down=SM->add(new CTree_with_prefactor(pro_32154,ind_32154,_prop_hel_fn_down_3));
if(tree_color!=1){
		T_14532_down=SM->add(new CTree_with_prefactor(pro_14532,ind_14532,_prop_hel_fn_down_1));
		T_14325_down=SM->add(new CTree_with_prefactor(pro_14325,ind_14325,_prop_hel_fn_down_1));
		T_32514_down=SM->add(new CTree_with_prefactor(pro_32514,ind_32514,_prop_hel_fn_down_3));
		T_32145_down=SM->add(new CTree_with_prefactor(pro_32145,ind_32145,_prop_hel_fn_down_3));
}

if(color!=1){
		L1_15432_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15432,ind_15432,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L3_14352_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_14352,ind_14352,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L1_35214_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35214,ind_35214,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L3_32154_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_32154,ind_32154,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
}
		L1_color_15432_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_15432,ind_15432,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_14352_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_14352,ind_14352,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L1_color_35214_down=SM->add(A_loop_2q_2Q_1g_2l_7_1(pro_35214,ind_35214,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
		L3_color_32154_down=SM->add(A_loop_2q_2Q_1g_2l_7_3(pro_32154,ind_32154,n_s,n_f,n_c,1,photonZW,e_ph_type,color, lo_or_nlo));
if(color!=1){
		L2_14532_down=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_14532,ind_14532,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L4_14325_down=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_14325,ind_14325,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L2_32514_down=SM->add(A_loop_2q_2Q_1g_2l_7_2(pro_32514,ind_32514,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
		L4_32145_down=SM->add(A_loop_2q_2Q_1g_2l_7_4(pro_32145,ind_32145,n_s,n_f,n_c,1,photonZW,e_ph_type,0, lo_or_nlo));
}
	};


//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction nf(n_f,1);
multi_precision_fraction n2m1(n_c*n_c-1,1);
multi_precision_fraction _over_nc(1,n_c);

multi_precision_fraction ext(2*4*n_c*n_c*(n_c*n_c-1),1);
multi_precision_fraction ext_tree(4*n_c*(n_c*n_c-1),1);


multi_precision_fraction Qu(2,3);
multi_precision_fraction Qd(-1,3);


//--------------------------------------------------------
// mulitplicities of processes

multi_precision_fraction uu_multiplicity(1,1);//2/2
multi_precision_fraction dd_multiplicity(1,1);//3/2
multi_precision_fraction uup_multiplicity(1,1);//2
multi_precision_fraction ddp_multiplicity(1,1);//6
multi_precision_fraction ud_multiplicity(1,1);//6
multi_precision_fraction du_multiplicity(1,1);//6


// W cases
if(photonZW==3){
	//--------------------------------------------------------
	// W case: case4q==0 -> identical quarks
	// Notice that we use _up to built W amps!
	//--------------------------------------------------------
	if(case4q==0){
		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=-uu*_over_nc;
		multi_precision_fraction uu_EE=-uu_E*_over_nc;
		multi_precision_fraction uu_EEE=-uu_EE*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
		};
		if(hp3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_35214_up,T_35214_up,uu));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_32154_up,uu));
		};

if(tree_color!=1){

		if((h1.helicity()==-1&&h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,uu_EE));
		};

		if(hp3.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35214_up,T_35214_up,uu));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_32514_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35214_up,T_32154_up,0));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_32145_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_up,T_35214_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_32514_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_32154_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32514_up,T_32145_up,0));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_32154_up,T_35214_up,0));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_32514_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_32154_up,T_32154_up,uu));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_32145_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32145_up,T_35214_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_32145_up,T_32514_up,0));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_32154_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_32145_up,uu_EE));
		};

		if(h3.helicity()==-1&&h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35214_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_12354_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_35214_up,T_12345_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_12534_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_32514_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_12345_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_32154_up,T_15234_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_32154_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_12345_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_32145_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12534_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12345_up,uu_EEE));

		};
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_32514_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_32154_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_15234_up,T_32145_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_32514_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_12534_up,T_32154_up,0));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_32145_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_up,T_35214_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_32514_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_32145_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_12345_up,T_35214_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32514_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32145_up,uu_EEE));

		};
}
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;
		multi_precision_fraction uu_EE_virtual=-uu_E_virtual*_over_nc;
		multi_precision_fraction uu_EEE_virtual=-uu_EE_virtual*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_color_15234_up,T_15234_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L3_color_12354_up,T_12354_up,uu_virtual));
		};

		if(h3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_color_35214_up,T_35214_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L3_color_32154_up,T_32154_up,uu_virtual));
		};
if(color!=1){
		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_up,T_15234_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,uu_E_virtual));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,uu_EE_virtual));
		};
		if(h3.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35214_up,T_35214_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_32514_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35214_up,T_32154_up,0));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_32145_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_35214_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_32514_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_32154_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_32514_up,T_32145_up,0));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_32154_up,T_35214_up,0));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_32514_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_32154_up,T_32154_up,uu));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_32145_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_32145_up,T_35214_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L4_32145_up,T_32514_up,0));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_32154_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_32145_up,uu_EE_virtual));
		};
		if(h3.helicity()==-1&&h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_12354_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35214_up,T_12345_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_12534_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_32514_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_12345_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_15234_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_32154_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_12345_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_32145_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12534_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12345_up,uu_EEE_virtual));
		};
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_32514_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_32154_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_up,T_32145_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_32514_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_12534_up,T_32154_up,0));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_32145_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_35214_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_32514_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_32145_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_12345_up,T_35214_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32514_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32145_up,uu_EEE_virtual));
		};
}
	}
	//--------------------------------------------------------
	// W case: case4q==1 -> identical anti-quarks
	// Notice that we use _up to built W amps!
	//--------------------------------------------------------
	if(case4q==1){

#if 0
_PRINT("-------------------------");
_PRINT("identical anti-quarks");
_PRINT("-------------------------");
#endif
		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=-uu*_over_nc;
		multi_precision_fraction uu_EE=-uu_E*_over_nc;
		multi_precision_fraction uu_EEE=-uu_EE*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
		};
		if(hp1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_15432_up,T_15432_up,uu));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_14352_up,uu));
		};

if(tree_color!=1){
		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,uu_EE));
		};

		if(hp1.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15432_up,T_15432_up,uu));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_14532_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15432_up,T_14352_up,0));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_14325_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532_up,T_15432_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_14532_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_14352_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532_up,T_14325_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_14352_up,T_15432_up,0));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_14532_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_14352_up,T_14352_up,uu));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_14325_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_up,T_15432_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325_up,T_14532_up,0));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_14352_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_14325_up,uu_EE));
		};
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432_up,T_15234_up,uu_E));
// 3			SM->add_tree(cached_cross_term_md(T_15432_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_12345_up,uu_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_14532_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12534_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12345_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_12354_up,uu_E));
// 1			SM->add_tree(cached_cross_term_md(T_14352_up,T_12345_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_12534_up,uu_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_14325_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_12345_up,uu_EEE));
		};
		if(hp1.helicity()==-1&&(h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_up,T_15432_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_15234_up,T_14532_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_14325_up,uu_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_12534_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14532_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14325_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_14532_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_14352_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_12354_up,T_14325_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_14532_up,uu_EEE));
//			SM->add_tree(cached_cross_term_md(T_12345_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_14325_up,uu_EEE));
		};
}
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;
		multi_precision_fraction uu_EE_virtual=-uu_E_virtual*_over_nc;
		multi_precision_fraction uu_EEE_virtual=-uu_EE_virtual*_over_nc;

		if(hp1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_color_15234_up,T_15234_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L3_color_12354_up,T_12354_up,uu_virtual));
		};

		if(hp1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_color_15432_up,T_15432_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L3_color_14352_up,T_14352_up,uu_virtual));
		};

if(color!=1){
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,uu_E_virtual));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,uu_EE_virtual));
		};
		if(hp1.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15432_up,T_15432_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_14532_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_15432_up,T_14352_up,0));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_14325_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_15432_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_14532_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_14352_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532_up,T_14325_up,uu_E_virtual));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_14352_up,T_15432_up,0));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_14532_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_14352_up,T_14352_up,uu_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_14325_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_up,T_15432_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325_up,T_14532_up,0));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_14352_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_14325_up,uu_EE_virtual));
		};
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_15432_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_up,T_14532_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_14325_up,uu_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_12534_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14532_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14325_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_14532_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_14352_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_12354_up,T_14325_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_14532_up,uu_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_12345_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_14325_up,uu_EEE_virtual));
		};
		if(hp1.helicity()==-1&&(h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_15234_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_15432_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_12345_up,uu_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_14532_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12534_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12345_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_12354_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_14352_up,T_12345_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_12534_up,uu_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_14325_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_12345_up,uu_EEE_virtual));
		};

}
	}

	//--------------------------------------------------------
	// W case: case4q==2 -> no identical quarks
	// Notice that we use _up to built W amps!
	if(case4q==2){
#if 0
_PRINT("-------------------------");
_PRINT("no identical quarks");
_PRINT("-------------------------");
#endif
			multi_precision_fraction ud=ud_multiplicity*ext_tree;
			multi_precision_fraction ud_EE=ud_multiplicity*ext_tree*_over_nc*_over_nc;
			multi_precision_fraction ud_virtual=ud_multiplicity*ext;
			multi_precision_fraction ud_EE_virtual=ud_multiplicity*ext*_over_nc*_over_nc;


		if(hp1.helicity()==-1&&h1.helicity()+h4.helicity()==0){

			SM->add_tree(cached_cross_term_md(T_15234_up  ,T_15234_up  ,ud));
			SM->add_tree(cached_cross_term_md(T_12354_up  ,T_12354_up  ,ud));

			SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_15234_up  ,ud_virtual));
			SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_12354_up  ,ud_virtual));
		};
if(tree_color!=1){
		
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,ud));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,ud));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,ud_EE));
		}
}
if(color!=1){
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,ud_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,ud_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,ud_E_virtual));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,ud_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,ud_EE_virtual));
		};
}
	
}
	return SM;
}


////////////////////////////////
////////////////////////////////
// Z/gamma cases
////////////////////////////////
////////////////////////////////
//--------------------------------------------------------
// (u,d)
{
if(case4q==2 || case4q==-1){
		multi_precision_fraction ud=ud_multiplicity*ext_tree;
		multi_precision_fraction ud_E=-ud*_over_nc;
		multi_precision_fraction ud_EE=-ud_E*_over_nc;
		multi_precision_fraction ud_EEE=-ud_EE*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_up  ,T_15234_up  ,ud));
		SM->add_tree(cached_cross_term_md(T_12354_up  ,T_12354_up  ,ud));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,ud));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,ud));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_up  ,ud));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_up  ,ud));
		
		SM->add_tree(cached_cross_term_md(T_15234_up  ,T_34152_down,ud));
		SM->add_tree(cached_cross_term_md(T_12354_up  ,T_35412_down,ud));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,ud_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,uu));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34512_down,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34125_down,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_down,T_35412_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34512_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34152_down,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_34125_down,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34512_down,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,uu));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34125_down,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_down,T_35412_down,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34152_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34125_down,ud_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12534_up,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12345_up,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_down,T_15234_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12534_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12354_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_12345_up,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12534_up,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12345_up,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_down,T_15234_up,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12354_up,ud_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12345_up,ud_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_up,T_35412_down,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34512_down,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34125_down,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_up,T_35412_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34512_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34152_down,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_34125_down,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34512_down,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_up,T_34152_down,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34125_down,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_up,T_35412_down,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34152_down,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34125_down,ud_EE));
		
		};
}

		multi_precision_fraction ud_virtual=ud_multiplicity*ext;
		multi_precision_fraction ud_EE_virtual=ud_virtual*_over_nc*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_15234_up  ,ud_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_12354_up  ,ud_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_34152_down,ud_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_35412_down,ud_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_15234_up  ,ud_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_12354_up  ,ud_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_34152_down,ud_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_35412_down,ud_virtual));
		};

if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_up,T_15234_up,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,ud_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,ud_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,ud_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_down,T_35412_down,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34512_down,ud_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34125_down,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_35412_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34512_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34152_down,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_34125_down,ud_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34512_down,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_down,T_34152_down,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34125_down,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_down,T_35412_down,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34152_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34125_down,ud_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_down,T_15234_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12534_up,ud_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12345_up,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_15234_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12534_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12354_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_12345_up,ud_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12534_up,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_down,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12345_up,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_down,T_15234_up,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12354_up,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12345_up,ud_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_up,T_35412_down,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34512_down,ud_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34125_down,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_35412_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34512_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34152_down,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_34125_down,ud_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34512_down,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_up,T_34152_down,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34125_down,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_up,T_35412_down,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34152_down,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34125_down,ud_EE_virtual));
		
		};
}


	};

//--------------------------------------------------------
// (d,u)
if(case4q==3 || case4q==-1){

		multi_precision_fraction du=du_multiplicity*ext_tree;
		multi_precision_fraction du_E=-du*_over_nc;
		multi_precision_fraction du_EE=-du_E*_over_nc;
		multi_precision_fraction du_EEE=-du_EE*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_down  ,T_15234_down  ,du));
		SM->add_tree(cached_cross_term_md(T_12354_down  ,T_12354_down  ,du));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,du));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,du));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_down  ,du));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_down  ,du));
		
		SM->add_tree(cached_cross_term_md(T_15234_down  ,T_34152_up,du));
		SM->add_tree(cached_cross_term_md(T_12354_down  ,T_35412_up,du));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_down,T_15234_down,uu));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12534_down,du_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12345_down,du_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_down,T_15234_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12534_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12354_down,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_12345_down,du_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12534_down,du_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_down,T_12354_down,uu));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12345_down,du_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_down,T_15234_down,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12354_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12345_down,du_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34512_up,du_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34125_up,du_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_up,T_35412_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34512_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34152_up,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_34125_up,du_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34512_up,du_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34125_up,du_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_up,T_35412_up,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34152_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34125_up,du_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_down,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12534_down,du_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12345_down,du_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_up,T_15234_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12534_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12354_down,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_12345_down,du_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12534_down,du_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_down,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12345_down,du_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_up,T_15234_down,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12354_down,du_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12345_down,du_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_down,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34512_up,du_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34125_up,du_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_down,T_35412_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34512_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34152_up,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_34125_up,du_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34512_up,du_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_down,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34125_up,du_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_down,T_35412_up,du_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34152_up,du_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34125_up,du_EE));
		
		};
}

		multi_precision_fraction du_virtual=du_multiplicity*ext;
		multi_precision_fraction du_EE_virtual=du_virtual*_over_nc*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_15234_down  ,du_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_12354_down  ,du_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_34152_up,du_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_35412_up,du_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_15234_down  ,du_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_12354_down  ,du_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_34152_up,du_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_35412_up,du_virtual));
		};
if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_down,T_15234_down,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234_down,T_12534_down,du_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_12345_down,du_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_15234_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12534_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12354_down,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_12345_down,du_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12534_down,du_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_down,T_12354_down,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12345_down,du_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_down,T_15234_down,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12354_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12345_down,du_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_up,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34512_up,du_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34125_up,du_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_35412_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34512_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34152_up,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_34125_up,du_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34512_up,du_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_up,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34125_up,du_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_up,T_35412_up,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34152_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34125_up,du_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_up,T_15234_down,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12534_down,du_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12345_down,du_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_15234_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12534_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12354_down,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_12345_down,du_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12534_down,du_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_up,T_12354_down,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12345_down,du_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_up,T_15234_down,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12354_down,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12345_down,du_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_down,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34512_up,du_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34125_up,du_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_35412_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34512_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34152_up,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_34125_up,du_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34512_up,du_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_down,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34125_up,du_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_down,T_35412_up,du_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34152_up,du_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34125_up,du_EE_virtual));
		
		};
}

	};

//--------------------------------------------------------
// (u,u')
if(case4q==1 || case4q==-1){

		multi_precision_fraction uup=uup_multiplicity*ext_tree;
		multi_precision_fraction uup_E=-uup*_over_nc;
		multi_precision_fraction uup_EE=-uup_E*_over_nc;
		multi_precision_fraction uup_EEE=-uup_EE*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_up  ,T_15234_up  ,uup));
		SM->add_tree(cached_cross_term_md(T_12354_up  ,T_12354_up  ,uup));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,uup));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,uup));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_up  ,uup));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_up  ,uup));
		
		SM->add_tree(cached_cross_term_md(T_15234_up  ,T_34152_up,uup));
		SM->add_tree(cached_cross_term_md(T_12354_up  ,T_35412_up,uup));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,uup_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,uup_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,uup_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,uup_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,uup_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,uup_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34512_up,uup_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34125_up,uup_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_up,T_35412_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34512_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34152_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_34125_up,uup_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34512_up,uup_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34125_up,uup_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_up,T_35412_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34152_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34125_up,uup_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12534_up,uup_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12345_up,uup_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_up,T_15234_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12534_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12354_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_12345_up,uup_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12534_up,uup_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12345_up,uup_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_up,T_15234_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12354_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12345_up,uup_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_up,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34512_up,uup_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34125_up,uup_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_up,T_35412_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34512_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34152_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_34125_up,uup_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34512_up,uup_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_up,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34125_up,uup_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_up,T_35412_up,uup_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34152_up,uup_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34125_up,uup_EE));
		
		};
}

		multi_precision_fraction uup_virtual=uup_multiplicity*ext;
		multi_precision_fraction uup_EE_virtual=uup_virtual*_over_nc*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_15234_up  ,uup_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_12354_up  ,uup_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_34152_up,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_35412_up,uup_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_15234_up  ,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_12354_up  ,uup_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_34152_up,uup_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_35412_up,uup_virtual));
		};
if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_up,T_15234_up,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,uup_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,uup_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,uup_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,uup_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,uup_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,uup_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_up,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34512_up,uup_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34125_up,uup_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_35412_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34512_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34152_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_34125_up,uup_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34512_up,uup_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_up,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34125_up,uup_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_up,T_35412_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34152_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34125_up,uup_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_up,T_15234_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12534_up,uup_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12345_up,uup_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_15234_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12534_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12354_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_12345_up,uup_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12534_up,uup_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12345_up,uup_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_up,T_15234_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12354_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12345_up,uup_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_up,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34512_up,uup_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34125_up,uup_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_35412_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34512_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34152_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_34125_up,uup_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34512_up,uup_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_up,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34125_up,uup_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_up,T_35412_up,uup_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34152_up,uup_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34125_up,uup_EE_virtual));
		
		};
}
}

//--------------------------------------------------------
// (d,d')
if(case4q==4 || case4q==-1){

		multi_precision_fraction ddp=ddp_multiplicity*ext_tree;
		multi_precision_fraction ddp_E=-ddp*_over_nc;
		multi_precision_fraction ddp_EE=-ddp_E*_over_nc;
		multi_precision_fraction ddp_EEE=-ddp_EE*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_down  ,T_15234_down  ,ddp));
		SM->add_tree(cached_cross_term_md(T_12354_down  ,T_12354_down  ,ddp));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,ddp));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,ddp));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_down  ,ddp));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_down  ,ddp));
		
		SM->add_tree(cached_cross_term_md(T_15234_down  ,T_34152_down,ddp));
		SM->add_tree(cached_cross_term_md(T_12354_down  ,T_35412_down,ddp));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_down,T_15234_down,uu));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12534_down,ddp_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12345_down,ddp_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_down,T_15234_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12534_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12354_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_12345_down,ddp_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12534_down,ddp_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_down,T_12354_down,uu));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12345_down,ddp_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_down,T_15234_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12354_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12345_down,ddp_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,uu));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34512_down,ddp_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34125_down,ddp_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_down,T_35412_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34512_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34152_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_34125_down,ddp_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34512_down,ddp_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,uu));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34125_down,ddp_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_down,T_35412_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34152_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34125_down,ddp_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_down,uu));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12534_down,ddp_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12345_down,ddp_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_down,T_15234_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12534_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12354_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_12345_down,ddp_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12534_down,ddp_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_down,uu));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12345_down,ddp_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_down,T_15234_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12354_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12345_down,ddp_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_down,T_35412_down,uu));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34512_down,ddp_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34125_down,ddp_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_down,T_35412_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34512_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34152_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_34125_down,ddp_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34512_down,ddp_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_down,T_34152_down,uu));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34125_down,ddp_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_down,T_35412_down,ddp_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34152_down,ddp_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34125_down,ddp_EE));
		
		};
}

		multi_precision_fraction ddp_virtual=ddp_multiplicity*ext;
		multi_precision_fraction ddp_EE_virtual=ddp_virtual*_over_nc*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_15234_down  ,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_12354_down  ,ddp_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_34152_down,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_35412_down,ddp_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_15234_down  ,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_12354_down  ,ddp_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_34152_down,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_35412_down,ddp_virtual));
		};

if(color!=1){
		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_down,T_15234_down,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234_down,T_12534_down,ddp_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_12345_down,ddp_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_15234_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12534_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12354_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_12345_down,ddp_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12534_down,ddp_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_down,T_12354_down,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12345_down,ddp_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_down,T_15234_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12354_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12345_down,ddp_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_down,T_35412_down,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34512_down,ddp_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34125_down,ddp_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_35412_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34512_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34152_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_34125_down,ddp_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34512_down,ddp_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_down,T_34152_down,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34125_down,ddp_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_down,T_35412_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34152_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34125_down,ddp_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_down,T_15234_down,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12534_down,ddp_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12345_down,ddp_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_15234_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12534_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12354_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_12345_down,ddp_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12534_down,ddp_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_down,T_12354_down,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12345_down,ddp_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_down,T_15234_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12354_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12345_down,ddp_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_down,T_35412_down,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34512_down,ddp_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34125_down,ddp_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_35412_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34512_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34152_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_34125_down,ddp_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34512_down,ddp_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_down,T_34152_down,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34125_down,ddp_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_down,T_35412_down,ddp_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34152_down,ddp_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34125_down,ddp_EE_virtual));
		
		};
}

}

//--------------------------------------------------------
// (u,u)
//--------------------------------------------------------
if(case4q==0 || case4q==-1){

		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=-uu*_over_nc;
		multi_precision_fraction uu_EE=-uu_E*_over_nc;
		multi_precision_fraction uu_EEE=-uu_EE*_over_nc;
		
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;
		multi_precision_fraction uu_EE_virtual=-uu_E_virtual*_over_nc;
		multi_precision_fraction uu_EEE_virtual=-uu_EE_virtual*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
		SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,uu));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_up,uu));
		SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_15234_up,T_34152_up,uu));
		SM->add_tree(cached_cross_term_md(T_12354_up,T_35412_up,uu));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12534_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_12345_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_15234_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12534_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_12354_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12534_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_12345_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15234_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12354_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_12345_up,uu_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_up,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34512_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_34125_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_up,T_35412_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34512_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_34152_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_34125_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34512_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_up,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_34125_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_up,T_35412_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34152_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_34125_up,uu_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_up,T_15234_up,uu));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12534_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_12345_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_up,T_15234_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12534_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_12354_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12534_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_up,T_12354_up,uu));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_12345_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_up,T_15234_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12354_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_12345_up,uu_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_up,T_35412_up,uu));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34512_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_34125_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_up,T_35412_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34512_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_34152_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_up,T_34125_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34512_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_up,T_34152_up,uu));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_34125_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_up,T_35412_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34152_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_34125_up,uu_EE));
		
		};
}


	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_15234_up  ,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_12354_up  ,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_34152_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_35412_up,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_up,T_15234_up  ,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_up,T_12354_up  ,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_up  ,T_34152_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_up  ,T_35412_up,uu_virtual));
		};

if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_up,T_15234_up,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234_up,T_12534_up,uu_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_12345_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_15234_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12534_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_12354_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12534_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_12345_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15234_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12354_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_12345_up,uu_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_up,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34512_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_34125_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_35412_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34512_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_34152_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_34125_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34512_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_up,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_34125_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_up,T_35412_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34152_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_34125_up,uu_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_up,T_15234_up,uu));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12534_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_12345_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_15234_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12534_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_12354_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_up,T_12345_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12534_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_up,T_12354_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_12345_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_up,T_15234_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12354_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_12345_up,uu_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_up,T_35412_up,uu));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34512_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_34125_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_35412_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34512_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_34152_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_up,T_34125_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34512_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_up,T_34152_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_34125_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_up,T_35412_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34152_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_34125_up,uu_EE_virtual));
		
		};
}

	if(h1.helicity()+h2.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15432_up,T_15432_up,uu));
		SM->add_tree(cached_cross_term_md(T_14352_up,T_14352_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_32154_up,T_32154_up,uu));
		SM->add_tree(cached_cross_term_md(T_35214_up,T_35214_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_32154_up,T_15432_up,uu));
		SM->add_tree(cached_cross_term_md(T_35214_up,T_14352_up,uu));
		
		SM->add_tree(cached_cross_term_md(T_15432_up,T_32154_up,uu));
		SM->add_tree(cached_cross_term_md(T_14352_up,T_35214_up,uu));
		}	

if(tree_color!=1){

		if(h1.helicity()+h2.helicity()==0){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15432_up,T_15432_up,uu));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_14532_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15432_up,T_14352_up,0));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_14325_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532_up,T_15432_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_14532_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_14352_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532_up,T_14325_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_14352_up,T_15432_up,0));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_14532_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_14352_up,T_14352_up,uu));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_14325_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_up,T_15432_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325_up,T_14532_up,0));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_14352_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_14325_up,uu_EE));
	
		//block 4
		//relative to block 1: 1432->3214	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35214_up,T_35214_up,uu));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_32514_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35214_up,T_32154_up,0));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_32145_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_up,T_35214_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_32514_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_32154_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32514_up,T_32145_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_32154_up,T_35214_up,0));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_32514_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_32154_up,T_32154_up,uu));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_32145_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32145_up,T_35214_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32145_up,T_32514_up,0));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_32154_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_32145_up,uu_EE));
		
		//block 3
		//relative to block 1: 1432_L->3214_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_32154_up,T_15432_up,uu));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_14532_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_32154_up,T_14352_up,0));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_14325_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_32145_up,T_15432_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_14532_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_14352_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32145_up,T_14325_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35214_up,T_15432_up,0));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_14532_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35214_up,T_14352_up,uu));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_14325_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32514_up,T_15432_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32514_up,T_14532_up,0));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_14352_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_14325_up,uu_EE));

		
		//block 2
		//relative to block 3: 1432->3214
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_14352_up,T_35214_up,uu));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_32514_up,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_14352_up,T_32154_up,0));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_32145_up,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14325_up,T_35214_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_32514_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_32154_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325_up,T_32145_up,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15432_up,T_35214_up,0));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_32514_up,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15432_up,T_32154_up,uu));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_32145_up,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14532_up,T_35214_up,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532_up,T_32514_up,0));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_32154_up,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_32145_up,uu_EE));
		
		};
}


	if(h1.helicity()+h2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15432_up,T_15432_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_14352_up,T_14352_up,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_32154_up,T_32154_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35214_up,T_35214_up,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_32154_up,T_15432_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35214_up,T_14352_up,uu_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15432_up,T_32154_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_14352_up,T_35214_up,uu_virtual));
		};

if(color!=1){

		
		if(h1.helicity()+h2.helicity()==0){
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15432_up,T_15432_up,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15432_up,T_14532_up,uu_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15432_up,T_14352_up,0));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_14325_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_15432_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_14532_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_14352_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532_up,T_14325_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_14352_up,T_15432_up,0));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_14532_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_14352_up,T_14352_up,uu));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_14325_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_up,T_15432_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325_up,T_14532_up,0));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_14352_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_14325_up,uu_EE_virtual));
		
		//block 4
		//relative to block 1: 1432->3214	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35214_up,T_35214_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_32514_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35214_up,T_32154_up,0));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_32145_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_35214_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_32514_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_32154_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_32514_up,T_32145_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_32154_up,T_35214_up,0));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_32514_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_32154_up,T_32154_up,uu));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_32145_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_32145_up,T_35214_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_32145_up,T_32514_up,0));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_32154_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_32145_up,uu_EE_virtual));
		
		//block 3
		//relative to block 1: 1432_L->3214_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_32154_up,T_15432_up,uu));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_14532_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_32154_up,T_14352_up,0));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_14325_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_15432_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_14532_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_14352_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_32145_up,T_14325_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35214_up,T_15432_up,0));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_14532_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35214_up,T_14352_up,uu));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_14325_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_32514_up,T_15432_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_32514_up,T_14532_up,0));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_14352_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_14325_up,uu_EE_virtual));

		
		//block 2
		//relative to block 3: 1432->3214
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_14352_up,T_35214_up,uu));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_32514_up,uu_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_14352_up,T_32154_up,0));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_32145_up,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_35214_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_32514_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_32154_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325_up,T_32145_up,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15432_up,T_35214_up,0));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_32514_up,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15432_up,T_32154_up,uu));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_32145_up,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_14532_up,T_35214_up,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532_up,T_32514_up,0));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_32154_up,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_32145_up,uu_EE_virtual));
		
		};
}

////////////////////////////////
///identical quark exchange term
/// block A
////////////////////////////////

if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35214_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_12354_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_35214_up,T_12345_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_12534_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_32514_up,T_12354_up,0));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_12345_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_32154_up,T_15234_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_32154_up,T_12534_up,0));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_12345_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_32145_up,T_15234_up,0));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12534_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_12345_up,uu_EEE));

		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_32514_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_32154_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_15234_up,T_32145_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_32514_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_12534_up,T_32154_up,0));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_32145_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_up,T_35214_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_12354_up,T_32514_up,0));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_32145_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_12345_up,T_35214_up,0));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32514_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_32145_up,uu_EEE));

		};
}

if(color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_12354_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35214_up,T_12345_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_12534_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_32514_up,T_12354_up,0));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_12345_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_15234_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_32154_up,T_12534_up,0));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_12345_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_32145_up,T_15234_up,0));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12534_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_12345_up,uu_EEE_virtual));
		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_32514_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_32154_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_up,T_32145_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_32514_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_12534_up,T_32154_up,0));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_32145_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_35214_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_up,T_32514_up,0));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_32145_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_12345_up,T_35214_up,0));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32514_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_32145_up,uu_EEE_virtual));
		};
}

////////////////////////////
//// 2<->4 relative to block A
////////////////////////////
/// to be done


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35412_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_14532_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_14352_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_35412_up,T_14325_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_14532_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_34512_up,T_14352_up,0));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_14325_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_34152_up,T_15432_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_34152_up,T_14532_up,0));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_14325_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_34125_up,T_15432_up,0));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_14532_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_14325_up,uu_EEE));

		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432_up,T_35412_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_34512_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_34152_up,uu_E));
			//SM->add_tree(cached_cross_term_md(T_15432_up,T_34125_up,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532_up,T_35412_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_34512_up,uu_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_14532_up,T_34152_up,0));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_34125_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352_up,T_35412_up,uu_E));
			//0: SM->add_tree(cached_cross_term_md(T_14352_up,T_34512_up,0));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_34152_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_34125_up,uu_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_14325_up,T_35412_up,0));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_34512_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_34152_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_34125_up,uu_EEE));

		};
}

if(color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_14532_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_14352_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35412_up,T_14325_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_14532_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_34512_up,T_14352_up,0));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_14325_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_15432_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_up,T_14532_up,0));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_14325_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_34125_up,T_15432_up,0));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_14532_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_14325_up,uu_EEE_virtual));
		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_35412_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_34512_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_34152_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15432_up,T_34125_up,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_35412_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_34512_up,uu_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_14532_up,T_34152_up,0));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_34125_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_35412_up,uu_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_14352_up,T_34512_up,0));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_34152_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_34125_up,uu_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_14325_up,T_35412_up,0));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_34512_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_34152_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_34125_up,uu_EEE_virtual));
		};
}

////////////////////////////////
///identical anti-quark exchange term 
///block B
////////////////////////////////


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432_up,T_15234_up,uu_E));
// 3			SM->add_tree(cached_cross_term_md(T_15432_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432_up,T_12345_up,uu_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_14532_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12534_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532_up,T_12345_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_12534_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352_up,T_12354_up,uu_E));
// 1			SM->add_tree(cached_cross_term_md(T_14352_up,T_12345_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_up,T_15234_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_12534_up,uu_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_14325_up,T_12354_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325_up,T_12345_up,uu_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_up,T_15432_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_15234_up,T_14532_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234_up,T_14325_up,uu_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_12534_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14532_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534_up,T_14325_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_14532_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354_up,T_14352_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_12354_up,T_14325_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_up,T_15432_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_14532_up,uu_EEE));
//			SM->add_tree(cached_cross_term_md(T_12345_up,T_14352_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345_up,T_14325_up,uu_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_15432_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_up,T_14532_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_up,T_14325_up,uu_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_12534_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14532_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_up,T_14325_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_14532_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_up,T_14352_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_12354_up,T_14325_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_up,T_15432_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_14532_up,uu_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_12345_up,T_14352_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_up,T_14325_up,uu_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_15234_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_15432_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_up,T_12345_up,uu_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_14532_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12534_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_up,T_12345_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_12534_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_up,T_12354_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_14352_up,T_12345_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_up,T_15234_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_12534_up,uu_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_14325_up,T_12354_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_up,T_12345_up,uu_EEE_virtual));
		};

}
////////////////////////////
//// 1<->3 relative to block B
////////////////////////////
/// to be done


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35412_up,T_35214_up,uu_E));
// 1			SM->add_tree(cached_cross_term_md(T_35412_up,T_32514_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35412_up,T_32145_up,uu_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_34512_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_32514_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34512_up,T_32145_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_34152_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_32514_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34152_up,T_32154_up,uu_E));
// 3			SM->add_tree(cached_cross_term_md(T_34152_up,T_32145_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_up,T_35214_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_32514_up,uu_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_34125_up,T_32154_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_34125_up,T_32145_up,uu_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35214_up,T_35412_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_35214_up,T_34512_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_34152_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_35214_up,T_34125_up,uu_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_32514_up,T_35412_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_34512_up,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_34152_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32514_up,T_34125_up,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_32154_up,T_35412_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_34512_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32154_up,T_34152_up,uu_E));
//			SM->add_tree(cached_cross_term_md(T_32154_up,T_34125_up,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32145_up,T_35412_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_34512_up,uu_EEE));
//			SM->add_tree(cached_cross_term_md(T_32145_up,T_34152_up,uu_E));
			SM->add_tree(cached_cross_term_md(T_32145_up,T_34125_up,uu_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_35412_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35214_up,T_34512_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_34152_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_up,T_34125_up,uu_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_32514_up,T_35412_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_34512_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_34152_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_up,T_34125_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_35412_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_34512_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_up,T_34152_up,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_32154_up,T_34125_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_32145_up,T_35412_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_34512_up,uu_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_32145_up,T_34152_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_up,T_34125_up,uu_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_35214_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_35412_up,T_32514_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_up,T_32145_up,uu_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_34512_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_32514_up,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_up,T_32145_up,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_32514_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_up,T_32154_up,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_34152_up,T_32145_up,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_up,T_35214_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_32514_up,uu_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_34125_up,T_32154_up,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_up,T_32145_up,uu_EEE_virtual));
		};

}


}


// (d,d)
//--------------------------------------------------------
if(case4q==5 || case4q==-1){


		multi_precision_fraction dd=dd_multiplicity*ext_tree;
		multi_precision_fraction dd_E=-dd*_over_nc;
		multi_precision_fraction dd_EE=-dd_E*_over_nc;
		multi_precision_fraction dd_EEE=-dd_EE*_over_nc;
		
		multi_precision_fraction dd_virtual=dd_multiplicity*ext;
		multi_precision_fraction dd_E_virtual=-dd_virtual*_over_nc;
		multi_precision_fraction dd_EE_virtual=-dd_E_virtual*_over_nc;
		multi_precision_fraction dd_EEE_virtual=-dd_EE_virtual*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234_down,T_15234_down,dd));
		SM->add_tree(cached_cross_term_md(T_12354_down,T_12354_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,dd));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_down,dd));
		SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_15234_down,T_34152_down,dd));
		SM->add_tree(cached_cross_term_md(T_12354_down,T_35412_down,dd));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234_down,T_15234_down,dd));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12534_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_12345_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_down,T_15234_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12534_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_12354_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_12345_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12534_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354_down,T_12354_down,dd));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_12345_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_down,T_15234_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12354_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_12345_down,dd_EE));
	
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35412_down,T_35412_down,dd));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34512_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_34125_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_down,T_35412_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34512_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_34152_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_34125_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34512_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_34152_down,T_34152_down,dd));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_34125_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_down,T_35412_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34152_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_34125_down,dd_EE));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_34152_down,T_15234_down,dd));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12534_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_12345_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_34125_down,T_15234_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12534_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_12354_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34125_down,T_12345_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35412_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12534_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35412_down,T_12354_down,dd));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_12345_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34512_down,T_15234_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_34512_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12354_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_12345_down,dd_EE));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_12354_down,T_35412_down,dd));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34512_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_34125_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12345_down,T_35412_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34512_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_34152_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345_down,T_34125_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15234_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34512_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15234_down,T_34152_down,dd));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_34125_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12534_down,T_35412_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34152_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_34125_down,dd_EE));
		
		};
}


	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_15234_down  ,dd_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_12354_down  ,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_34152_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_35412_down,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_34152_down,T_15234_down  ,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35412_down,T_12354_down  ,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15234_down  ,T_34152_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354_down  ,T_35412_down,dd_virtual));
		};

if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234_down,T_15234_down,dd));
		  	SM->add_loop(cached_cross_term_md(L1_15234_down,T_12534_down,dd_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_12345_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_15234_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12534_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_12354_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_12345_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12534_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354_down,T_12354_down,dd));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_12345_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_down,T_15234_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12354_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_12345_down,dd_EE_virtual));
		
		//block 2
		//relative to block 1: 1234->3412	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35412_down,T_35412_down,dd));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34512_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_34125_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_35412_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34512_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_34152_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_34125_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34512_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_34152_down,T_34152_down,dd));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_34125_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_down,T_35412_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34152_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_34125_down,dd_EE_virtual));
		
		//block 3
		//relative to block 1: 1234_L->3412_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_34152_down,T_15234_down,dd));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12534_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_12345_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_15234_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12534_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_12354_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_34125_down,T_12345_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35412_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12534_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35412_down,T_12354_down,dd));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_12345_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_34512_down,T_15234_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_34512_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12354_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_12345_down,dd_EE_virtual));

		
		//block 4
		//relative to block 3: 1234->3412
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_12354_down,T_35412_down,dd));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34512_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_34125_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_35412_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34512_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_34152_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345_down,T_34125_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15234_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34512_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15234_down,T_34152_down,dd));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_34125_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_12534_down,T_35412_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34152_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_34125_down,dd_EE_virtual));
		
		};
}

	if(h1.helicity()+h2.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15432_down,T_15432_down,dd));
		SM->add_tree(cached_cross_term_md(T_14352_down,T_14352_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_32154_down,T_32154_down,dd));
		SM->add_tree(cached_cross_term_md(T_35214_down,T_35214_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_32154_down,T_15432_down,dd));
		SM->add_tree(cached_cross_term_md(T_35214_down,T_14352_down,dd));
		
		SM->add_tree(cached_cross_term_md(T_15432_down,T_32154_down,dd));
		SM->add_tree(cached_cross_term_md(T_14352_down,T_35214_down,dd));
		}	

if(tree_color!=1){

		if(h1.helicity()+h2.helicity()==0){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15432_down,T_15432_down,dd));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_14532_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15432_down,T_14352_down,0));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_14325_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532_down,T_15432_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_14532_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_14352_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532_down,T_14325_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_14352_down,T_15432_down,0));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_14532_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_14352_down,T_14352_down,dd));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_14325_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_down,T_15432_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325_down,T_14532_down,0));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_14352_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_14325_down,dd_EE));
	
		//block 4
		//relative to block 1: 1432->3214	
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_35214_down,T_35214_down,dd));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_32514_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_35214_down,T_32154_down,0));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_32145_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_down,T_35214_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_32514_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_32154_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32514_down,T_32145_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_32154_down,T_35214_down,0));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_32514_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_32154_down,T_32154_down,dd));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_32145_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32145_down,T_35214_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32145_down,T_32514_down,0));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_32154_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_32145_down,dd_EE));
		
		//block 3
		//relative to block 1: 1432_L->3214_L care with gluon
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_32154_down,T_15432_down,dd));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_14532_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_32154_down,T_14352_down,0));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_14325_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_32145_down,T_15432_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_14532_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_14352_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32145_down,T_14325_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_35214_down,T_15432_down,0));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_14532_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_35214_down,T_14352_down,dd));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_14325_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32514_down,T_15432_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_32514_down,T_14532_down,0));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_14352_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_14325_down,dd_EE));

		
		//block 2
		//relative to block 3: 1432->3214
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_14352_down,T_35214_down,dd));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_32514_down,dd_EE));
			//0: SM->add_tree(cached_cross_term_md(T_14352_down,T_32154_down,0));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_32145_down,dd_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14325_down,T_35214_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_32514_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_32154_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325_down,T_32145_down,dd_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_15432_down,T_35214_down,0));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_32514_down,dd_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_15432_down,T_32154_down,dd));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_32145_down,dd_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14532_down,T_35214_down,dd_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532_down,T_32514_down,0));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_32154_down,dd_EE));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_32145_down,dd_EE));
		
		};
}


	if(h1.helicity()+h2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15432_down,T_15432_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_14352_down,T_14352_down,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_32154_down,T_32154_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35214_down,T_35214_down,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L3_color_32154_down,T_15432_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_color_35214_down,T_14352_down,dd_virtual));
		
		SM->add_loop(cached_cross_term_md(L1_color_15432_down,T_32154_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_14352_down,T_35214_down,dd_virtual));
		};

if(color!=1){

		
		if(h1.helicity()+h2.helicity()==0){
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15432_down,T_15432_down,dd));
		  	SM->add_loop(cached_cross_term_md(L1_15432_down,T_14532_down,dd_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15432_down,T_14352_down,0));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_14325_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_15432_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_14532_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_14352_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532_down,T_14325_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_14352_down,T_15432_down,0));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_14532_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_14352_down,T_14352_down,dd));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_14325_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_down,T_15432_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325_down,T_14532_down,0));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_14352_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_14325_down,dd_EE_virtual));
		
		//block 4
		//relative to block 1: 1432->3214	
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_35214_down,T_35214_down,dd));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_32514_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L1_35214_down,T_32154_down,0));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_32145_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_35214_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_32514_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_32154_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_32514_down,T_32145_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_32154_down,T_35214_down,0));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_32514_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_32154_down,T_32154_down,dd));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_32145_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_32145_down,T_35214_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_32145_down,T_32514_down,0));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_32154_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_32145_down,dd_EE_virtual));
		
		//block 3
		//relative to block 1: 1432_L->3214_L care with gluon
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_32154_down,T_15432_down,dd));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_14532_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_32154_down,T_14352_down,0));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_14325_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_15432_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_14532_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_14352_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_32145_down,T_14325_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_35214_down,T_15432_down,0));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_14532_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_35214_down,T_14352_down,dd));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_14325_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_32514_down,T_15432_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_32514_down,T_14532_down,0));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_14352_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_14325_down,dd_EE_virtual));

		
		//block 2
		//relative to block 3: 1432->3214
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L3_14352_down,T_35214_down,dd));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_32514_down,dd_EE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_14352_down,T_32154_down,0));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_32145_down,dd_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_35214_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_32514_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_32154_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325_down,T_32145_down,dd_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L1_15432_down,T_35214_down,0));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_32514_down,dd_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L1_15432_down,T_32154_down,dd));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_32145_down,dd_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L2_14532_down,T_35214_down,dd_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532_down,T_32514_down,0));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_32154_down,dd_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_32145_down,dd_EE_virtual));
		
		};
}

////////////////////////////////
///identical quark exchange term
/// block A
////////////////////////////////

if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35214_down,T_15234_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_12534_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_12354_down,dd_E));
			//SM->add_tree(cached_cross_term_md(T_35214_down,T_12345_down,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_32514_down,T_15234_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_12534_down,dd_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_32514_down,T_12354_down,0));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_12345_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_32154_down,T_15234_down,dd_E));
			//0: SM->add_tree(cached_cross_term_md(T_32154_down,T_12534_down,0));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_12354_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_12345_down,dd_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_32145_down,T_15234_down,0));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_12534_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_12354_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_12345_down,dd_EEE));

		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_down,T_35214_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_32514_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_32154_down,dd_E));
			//SM->add_tree(cached_cross_term_md(T_15234_down,T_32145_down,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534_down,T_35214_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_32514_down,dd_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_12534_down,T_32154_down,0));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_32145_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_down,T_35214_down,dd_E));
			//0: SM->add_tree(cached_cross_term_md(T_12354_down,T_32514_down,0));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_32154_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_32145_down,dd_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_12345_down,T_35214_down,0));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_32514_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_32154_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_32145_down,dd_EEE));

		};
}

if(color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_15234_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_12534_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_12354_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35214_down,T_12345_down,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_15234_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_12534_down,dd_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_32514_down,T_12354_down,0));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_12345_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_15234_down,dd_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_32154_down,T_12534_down,0));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_12354_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_12345_down,dd_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_32145_down,T_15234_down,0));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_12534_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_12354_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_12345_down,dd_EEE_virtual));
		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_35214_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_32514_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_32154_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_down,T_32145_down,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_35214_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_32514_down,dd_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_12534_down,T_32154_down,0));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_32145_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_35214_down,dd_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_12354_down,T_32514_down,0));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_32154_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_32145_down,dd_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_12345_down,T_35214_down,0));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_32514_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_32154_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_32145_down,dd_EEE_virtual));
		};

}

////////////////////////////
//// 2<->4 relative to block A
////////////////////////////
/// to be done


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35412_down,T_15432_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_14532_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_14352_down,dd_E));
			//SM->add_tree(cached_cross_term_md(T_35412_down,T_14325_down,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_34512_down,T_15432_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_14532_down,dd_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_34512_down,T_14352_down,0));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_14325_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_34152_down,T_15432_down,dd_E));
			//0: SM->add_tree(cached_cross_term_md(T_34152_down,T_14532_down,0));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_14352_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_14325_down,dd_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_34125_down,T_15432_down,0));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_14532_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_14352_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_14325_down,dd_EEE));

		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432_down,T_35412_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_34512_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_34152_down,dd_E));
			//SM->add_tree(cached_cross_term_md(T_15432_down,T_34125_down,0));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532_down,T_35412_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_34512_down,dd_EEE));
			//0: SM->add_tree(cached_cross_term_md(T_14532_down,T_34152_down,0));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_34125_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352_down,T_35412_down,dd_E));
			//0: SM->add_tree(cached_cross_term_md(T_14352_down,T_34512_down,0));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_34152_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_34125_down,dd_E));
		//fourth row
		 	//0: SM->add_tree(cached_cross_term_md(T_14325_down,T_35412_down,0));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_34512_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_34152_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_34125_down,dd_EEE));

		};
}

if(color!=1){
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_15432_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_14532_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_14352_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35412_down,T_14325_down,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_15432_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_14532_down,dd_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_34512_down,T_14352_down,0));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_14325_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_15432_down,dd_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_34152_down,T_14532_down,0));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_14352_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_14325_down,dd_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_34125_down,T_15432_down,0));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_14532_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_14352_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_14325_down,dd_EEE_virtual));
		};
		if((h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_35412_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_34512_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_34152_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15432_down,T_34125_down,0));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_35412_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_34512_down,dd_EEE_virtual));
			//0: SM->add_loop(cached_cross_term_md(L2_14532_down,T_34152_down,0));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_34125_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_35412_down,dd_E_virtual));
			//0: SM->add_loop(cached_cross_term_md(L3_14352_down,T_34512_down,0));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_34152_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_34125_down,dd_E_virtual));
		//fourth row
		 	//0: SM->add_loop(cached_cross_term_md(L4_14325_down,T_35412_down,0));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_34512_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_34152_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_34125_down,dd_EEE_virtual));
		};
}

////////////////////////////////
///identical anti-quark exchange term 
///block B
////////////////////////////////


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432_down,T_15234_down,dd_E));
// 3			SM->add_tree(cached_cross_term_md(T_15432_down,T_12534_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_12354_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15432_down,T_12345_down,dd_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_14532_down,T_15234_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_12534_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_12354_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14532_down,T_12345_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352_down,T_15234_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_12534_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14352_down,T_12354_down,dd_E));
// 1			SM->add_tree(cached_cross_term_md(T_14352_down,T_12345_down,dd_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325_down,T_15234_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_12534_down,dd_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_14325_down,T_12354_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_14325_down,T_12345_down,dd_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234_down,T_15432_down,dd_E));
//			SM->add_tree(cached_cross_term_md(T_15234_down,T_14532_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_14352_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_15234_down,T_14325_down,dd_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_12534_down,T_15432_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_14532_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_14352_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12534_down,T_14325_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354_down,T_15432_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_14532_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12354_down,T_14352_down,dd_E));
//			SM->add_tree(cached_cross_term_md(T_12354_down,T_14325_down,dd_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345_down,T_15432_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_14532_down,dd_EEE));
//			SM->add_tree(cached_cross_term_md(T_12345_down,T_14352_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_12345_down,T_14325_down,dd_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_15432_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234_down,T_14532_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_14352_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234_down,T_14325_down,dd_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_12534_down,T_15432_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_14532_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_14352_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534_down,T_14325_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_15432_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_14532_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354_down,T_14352_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_12354_down,T_14325_down,dd_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345_down,T_15432_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_14532_down,dd_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_12345_down,T_14352_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345_down,T_14325_down,dd_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_15234_down,dd_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_15432_down,T_12534_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_12354_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432_down,T_12345_down,dd_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_14532_down,T_15234_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_12534_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_12354_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532_down,T_12345_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_15234_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_12534_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352_down,T_12354_down,dd_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_14352_down,T_12345_down,dd_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325_down,T_15234_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_12534_down,dd_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_14325_down,T_12354_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325_down,T_12345_down,dd_EEE_virtual));
		};


}

////////////////////////////
//// 1<->3 relative to block B
////////////////////////////
/// to be done


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35412_down,T_35214_down,dd_E));
// 1			SM->add_tree(cached_cross_term_md(T_35412_down,T_32514_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_32154_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35412_down,T_32145_down,dd_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_34512_down,T_35214_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_32514_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_32154_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34512_down,T_32145_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_34152_down,T_35214_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_32514_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34152_down,T_32154_down,dd_E));
// 3			SM->add_tree(cached_cross_term_md(T_34152_down,T_32145_down,dd_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_34125_down,T_35214_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_32514_down,dd_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_34125_down,T_32154_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_34125_down,T_32145_down,dd_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_35214_down,T_35412_down,dd_E));
//			SM->add_tree(cached_cross_term_md(T_35214_down,T_34512_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_34152_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_35214_down,T_34125_down,dd_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_32514_down,T_35412_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_34512_down,dd_EEE));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_34152_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32514_down,T_34125_down,dd_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_32154_down,T_35412_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_34512_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32154_down,T_34152_down,dd_E));
//			SM->add_tree(cached_cross_term_md(T_32154_down,T_34125_down,dd_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_32145_down,T_35412_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_34512_down,dd_EEE));
//			SM->add_tree(cached_cross_term_md(T_32145_down,T_34152_down,dd_E));
			SM->add_tree(cached_cross_term_md(T_32145_down,T_34125_down,dd_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_35412_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_35214_down,T_34512_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_34152_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35214_down,T_34125_down,dd_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_32514_down,T_35412_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_34512_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_34152_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_32514_down,T_34125_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_35412_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_34512_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_32154_down,T_34152_down,dd_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_32154_down,T_34125_down,dd_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_32145_down,T_35412_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_34512_down,dd_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_32145_down,T_34152_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_32145_down,T_34125_down,dd_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_35214_down,dd_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_35412_down,T_32514_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_32154_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_35412_down,T_32145_down,dd_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_34512_down,T_35214_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_32514_down,dd_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_32154_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_34512_down,T_32145_down,dd_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_35214_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_32514_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_34152_down,T_32154_down,dd_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_34152_down,T_32145_down,dd_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_34125_down,T_35214_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_32514_down,dd_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_34125_down,T_32154_down,dd_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_34125_down,T_32145_down,dd_EEE_virtual));
		};
}


}


}

return SM;
}

/*
Virtual_SME* vsme_2q2Q1g2l(int ns,int nf,int nc,int photonZW, int case4q,int color, int tree_color){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);
	ind.push_back(6);
	ind.push_back(7);

	vector<int> ind76;
	ind76.push_back(1);
	ind76.push_back(2);
	ind76.push_back(3);
	ind76.push_back(4);
	ind76.push_back(5);
	ind76.push_back(7);
	ind76.push_back(6);


clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
}	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
}
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
}	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color));
}

after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}
*/


Virtual_SME* vsme_2q2Q1g2l(std::vector<int> indext,int ns,int nf,int nc,int photonZW, int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);
	ind.push_back(indext[6]);

	vector<int> ind76;
	ind76.push_back(indext[0]);
	ind76.push_back(indext[1]);
	ind76.push_back(indext[2]);
	ind76.push_back(indext[3]);
	ind76.push_back(indext[4]);
	ind76.push_back(indext[6]);
	ind76.push_back(indext[5]);

clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,m,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,m,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,p,lm,lbp),ind,ns,nf,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbp,qm,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qp,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qp,qbm,qm,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbm,qp,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qm,qbp,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_2l_M2(process(qm,qbp,qp,qbm,p,lm,lbp),ind76,ns,nf,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}
after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}



BH_Ampl_2q2Q1g2l::BH_Ampl_2q2Q1g2l(int photonZW,int case4q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q1g2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				settings::BH_interface_settings::s_photon_only*photonZW,
				case4q,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			7,    // NbrExtParticles
			3,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			6,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment), d_color(color)
{}

double BH_Ampl_2q2Q1g2l::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(d_color){
		//Full color
		case 0 : return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
		//Leading color
		case 1: return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* 33/6. /*beta_0*/ *get_double_pole()/(-9.) /* = lc born tree */  );
		//Full-Leading color
		case 2: return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* ( 23/6. -  33/6. /*beta_0*/ *(25/3./9. - get_double_pole()/(-9.)  /* = lc born tree */ ))  );
		}
	case 2:	return get_double_pole();
	}

};




}
