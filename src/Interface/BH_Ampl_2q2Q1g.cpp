/*
 * matrix elements for 2q2Q1g
 *
 *  Created on: Apr 29, 2009
 *
 */


#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include <ctime>
#include "BH_interface_impl.h"

#define _VERBOSE 0

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {


//expected indices:
// q g Gb G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_5_1(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B3
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);

	process pro_12345=process(h1,h2,h3,h4,h5);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);

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
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}

if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(5-1);
	 i2=ind.at(1-1);
	 i3=ind.at(3-1);
	 i4=ind.at(4-1);
	 i5=ind.at(2-1);

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(3);
	 h4=pro.p(4);
	 h5=pro.p(2);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);


//a-type ~ sub_leading_color contributions:

	process pro_25341=process(h2,h5,h3,h4,h1);
	vector<int> ind_25341;
	ind_25341.push_back(i2);
	ind_25341.push_back(i5);
	ind_25341.push_back(i3);
	ind_25341.push_back(i4);
	ind_25341.push_back(i1);

	process pro_23415=process(h2,h3,h4,h1,h5);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);

	process pro_24315=process(h2,hp4,hp3,h1,h5);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);

	process pro_23541=process(h2,h3,h5,h4,h1);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);

	process pro_24531=process(h2,hp4,h5,hp3,h1);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);


	process pro_23451=process(h2,h3,h4,h5,h1);
	vector<int> ind_23451;
	ind_23451.push_back(i2);
	ind_23451.push_back(i3);
	ind_23451.push_back(i4);
	ind_23451.push_back(i5);
	ind_23451.push_back(i1);


	process pro_24351=process(h2,hp4,hp3,h5,h1);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);


	process pro_25431=process(h2,h5,hp4,hp3,h1);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);

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

	process pro_21435=process(h2,h1,h4,h3,h5);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);

	PA->add(pro_21435,slc_q,ind_21435,1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	PA->add(pro_21435,slc_G,ind_21435,1,n_c*n_c);

}

	return PA;
}


//expected indices:
// q Gb G g qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_5_3(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B1
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);

	process pro_12345=process(h1,h2,h3,h4,h5);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);



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
	PA->add(pro_12345,nf,ind_12345,n_f,n_c);
}



if((color==0)||(color==2)){
// For convenience we line up here with the EGKMZ-paper
// Apart from lining up particle labels we differ in putting the quark into the first position for the primitive amplitudes.
	 i1=ind.at(5-1);
	 i2=ind.at(1-1);
	 i3=ind.at(2-1);
	 i4=ind.at(3-1);
	 i5=ind.at(4-1);

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(3);
	 h5=pro.p(4);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:

	process pro_23451=process(h2,h3,h4,h5,h1);
	vector<int> ind_23451;
	ind_23451.push_back(i2);
	ind_23451.push_back(i3);
	ind_23451.push_back(i4);
	ind_23451.push_back(i5);
	ind_23451.push_back(i1);

	process pro_23415=process(h2,h3,h4,h1,h5);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);

	process pro_24315=process(h2,hp4,hp3,h1,h5);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);

	process pro_25341=process(h2,h5,h3,h4,h1);
	vector<int> ind_25341;
	ind_25341.push_back(i2);
	ind_25341.push_back(i5);
	ind_25341.push_back(i3);
	ind_25341.push_back(i4);
	ind_25341.push_back(i1);

	process pro_25431=process(h2,h5,hp4,hp3,h1);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);


	process pro_23541=process(h2,h3,h5,h4,h1);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);


	process pro_24531=process(h2,hp4,h5,hp3,h1);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);


	process pro_24351=process(h2,hp4,hp3,h5,h1);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);

	//leading_color: PA->add(pro_23451,sub_leading_color,ind_25341,1,1);
	PA->add(pro_23451,sub_leading_color,ind_23451,-1,n_c*n_c);
	PA->add(pro_23415,sub_leading_color,ind_23415,1,n_c*n_c);
	PA->add(pro_24315,sub_leading_color,ind_24315,1,n_c*n_c);
	PA->add(pro_25341,sub_leading_color,ind_25341,1,n_c*n_c);
	PA->add(pro_25431,sub_leading_color,ind_25431,1,n_c*n_c);
	PA->add(pro_23541,sub_leading_color,ind_23541,1,n_c*n_c);
	// EGKMZ PA->add(pro_24531,sub_leading_color,ind_24531,1,n_c*n_c);
	/*corr sign*/ PA->add(pro_24531,sub_leading_color,ind_24531,-1,n_c*n_c);
	PA->add(pro_24351,sub_leading_color,ind_24351,-1,n_c*n_c);



//b-type ~ sub_leading_color contributions:

	process pro_21543=process(h2,h1,h5,h4,h3);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);

	PA->add(pro_21543,slc_q,ind_21543,1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	PA->add(pro_21543,slc_G,ind_21543,1,n_c*n_c);

}

	return PA;
}


//expected indices:
// q Gb g G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_1g_5_2(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B4
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);

	process pro_12345=process(h1,h2,h3,h4,h5);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);



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

	 h1=pro.p(5);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(4);
	 h5=pro.p(3);

// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:

	process pro_24315=process(h2,hp4,hp3,h1,h5);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);

	process pro_25431=process(h2,h5,hp4,hp3,h1);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);

	process pro_24351=process(h2,hp4,hp3,h5,h1);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);

	process pro_24531=process(h2,hp4,h5,hp3,h1);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);
	
	process pro_23541=process(h2,h3,h5,h4,h1);
	vector<int> ind_23541;
	ind_23541.push_back(i2);
	ind_23541.push_back(i3);
	ind_23541.push_back(i5);
	ind_23541.push_back(i4);
	ind_23541.push_back(i1);

	PA->add(pro_24315,sub_leading_color,ind_24315,-1,1);
	PA->add(pro_25431,leading_color,ind_25431,-1,1);
	PA->add(pro_24351,leading_color,ind_24351,-1,1);
	PA->add(pro_23541,sub_leading_color,ind_23541,-1,n_c*n_c);
	// EGKMZ sign: PA->add(pro_24531,sub_leading_color,ind_24531,-1,n_c*n_c);
	/* corr sign:*/ PA->add(pro_24531,sub_leading_color,ind_24531,1,n_c*n_c);


//b-type ~ sub_leading_color contributions:

	process pro_21435=process(h2,h1,h4,h3,h5);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);

	process pro_21543=process(h2,h1,h5,h4,h3);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);

	process pro_25143=process(h2,h5,h1,h4,h3);
	vector<int> ind_25143;
	ind_25143.push_back(i2);
	ind_25143.push_back(i5);
	ind_25143.push_back(i1);
	ind_25143.push_back(i4);
	ind_25143.push_back(i3);

	process pro_21453=process(h2,h1,h4,h5,h3);
	vector<int> ind_21453;
	ind_21453.push_back(i2);
	ind_21453.push_back(i1);
	ind_21453.push_back(i4);
	ind_21453.push_back(i5);
	ind_21453.push_back(i3);

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
partial_amplitude_cached* A_loop_2q_2Q_1g_5_4(process pro,vector<int> & ind, int n_s, int n_f, int n_c,int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
// equivalent to B2
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);

	process pro_12345=process(h1,h2,h3,h4,h5);
	vector<int> ind_12345;
	ind_12345.push_back(i1);
	ind_12345.push_back(i2);
	ind_12345.push_back(i3);
	ind_12345.push_back(i4);
	ind_12345.push_back(i5);

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

	 h1=pro.p(4);
	 h2=pro.p(1);
	 h3=pro.p(2);
	 h4=pro.p(3);
	 h5=pro.p(5);


// exchange particle and anti particle for reversed order of the gluinos.
	ph_type hp3=particle_ID(gluino,h3.helicity(),1);
	ph_type hp4=particle_ID(gluino,h4.helicity(),1,true);

//a-type ~ sub_leading_color contributions:


	process pro_25431=process(h2,h5,hp4,hp3,h1);
	vector<int> ind_25431;
	ind_25431.push_back(i2);
	ind_25431.push_back(i5);
	ind_25431.push_back(i4);
	ind_25431.push_back(i3);
	ind_25431.push_back(i1);


	process pro_24531=process(h2,hp4,h5,hp3,h1);
	vector<int> ind_24531;
	ind_24531.push_back(i2);
	ind_24531.push_back(i4);
	ind_24531.push_back(i5);
	ind_24531.push_back(i3);
	ind_24531.push_back(i1);


	process pro_24351=process(h2,hp4,hp3,h5,h1);
	vector<int> ind_24351;
	ind_24351.push_back(i2);
	ind_24351.push_back(i4);
	ind_24351.push_back(i3);
	ind_24351.push_back(i5);
	ind_24351.push_back(i1);


	process pro_23415=process(h2,h3,h4,h1,h5);
	vector<int> ind_23415;
	ind_23415.push_back(i2);
	ind_23415.push_back(i3);
	ind_23415.push_back(i4);
	ind_23415.push_back(i1);
	ind_23415.push_back(i5);

	process pro_24315=process(h2,hp4,hp3,h1,h5);
	vector<int> ind_24315;
	ind_24315.push_back(i2);
	ind_24315.push_back(i4);
	ind_24315.push_back(i3);
	ind_24315.push_back(i1);
	ind_24315.push_back(i5);

	PA->add(pro_25431,leading_color,ind_25431,1,1);
	PA->add(pro_24531,sub_leading_color,ind_24531,1,1); // relative sign compared ot EGKMZ
	PA->add(pro_24351,leading_color,ind_24351,1,1);
	PA->add(pro_23415,sub_leading_color,ind_23415,-1,n_c*n_c);
	PA->add(pro_24315,sub_leading_color,ind_24315,-1,n_c*n_c);


//b-type ~ sub_leading_color contributions:


	process pro_21543=process(h2,h1,h5,h4,h3);
	vector<int> ind_21543;
	ind_21543.push_back(i2);
	ind_21543.push_back(i1);
	ind_21543.push_back(i5);
	ind_21543.push_back(i4);
	ind_21543.push_back(i3);


	process pro_21453=process(h2,h1,h4,h5,h3);
	vector<int> ind_21453;
	ind_21453.push_back(i2);
	ind_21453.push_back(i1);
	ind_21453.push_back(i4);
	ind_21453.push_back(i5);
	ind_21453.push_back(i3);

	process pro_21435=process(h2,h1,h4,h3,h5);
	vector<int> ind_21435;
	ind_21435.push_back(i2);
	ind_21435.push_back(i1);
	ind_21435.push_back(i4);
	ind_21435.push_back(i3);
	ind_21435.push_back(i5);


	PA->add(pro_21543,slc_q,ind_21543,-1,n_c*n_c);
	PA->add(pro_21453,slc_q,ind_21453,-1,n_c*n_c);
	PA->add(pro_21435,slc_q,ind_21435,-1,n_c*n_c);

//c-type ~ sub_leading_color contributions:

	process pro_25143=process(h2,h5,h1,h4,h3);
	vector<int> ind_25143;
	ind_25143.push_back(i2);
	ind_25143.push_back(i5);
	ind_25143.push_back(i1);
	ind_25143.push_back(i4);
	ind_25143.push_back(i3);

	PA->add(pro_21543,slc_G,ind_21543,-1,1);
	PA->add(pro_21453,slc_G,ind_21453,-1,1);
	PA->add(pro_21435,slc_G,ind_21435,-1,1);
	PA->add(pro_25143,slc_G,ind_25143,-1,1);
	PA->add(pro_25143,slc_G,ind_25143,1,n_c*n_c);
}
	return PA;
}



Squared_ME* A_loop_2q_2Q_1g_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);


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

ph_type hp1=pro.p(1);
ph_type hp2=pro.p(2);
ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),1);
ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),1,true);
ph_type hp5=pro.p(5);

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

process pro_15234=process(h1,h5,h2,h3,h4);
vector<int> ind_15234;
ind_15234.push_back(i1);
ind_15234.push_back(i5);
ind_15234.push_back(i2);
ind_15234.push_back(i3);
ind_15234.push_back(i4);

process pro_12534=process(h1,h2,h5,h3,h4);
vector<int> ind_12534;
ind_12534.push_back(i1);
ind_12534.push_back(i2);
ind_12534.push_back(i5);
ind_12534.push_back(i3);
ind_12534.push_back(i4);

process pro_12354=process(h1,h2,h3,h5,h4);
vector<int> ind_12354;
ind_12354.push_back(i1);
ind_12354.push_back(i2);
ind_12354.push_back(i3);
ind_12354.push_back(i5);
ind_12354.push_back(i4);

process pro_12345=process(h1,h2,h3,h4,h5);
vector<int> ind_12345;
ind_12345.push_back(i1);
ind_12345.push_back(i2);
ind_12345.push_back(i3);
ind_12345.push_back(i4);
ind_12345.push_back(i5);

process pro_35412=process(hf3,h5,hf4,hf1,hf2);
vector<int> ind_35412;
ind_35412.push_back(i3);
ind_35412.push_back(i5);
ind_35412.push_back(i4);
ind_35412.push_back(i1);
ind_35412.push_back(i2);

process pro_34512=process(hf3,hf4,h5,hf1,hf2);
vector<int> ind_34512;
ind_34512.push_back(i3);
ind_34512.push_back(i4);
ind_34512.push_back(i5);
ind_34512.push_back(i1);
ind_34512.push_back(i2);

process pro_34152=process(hf3,hf4,hf1,h5,hf2);
vector<int> ind_34152;
ind_34152.push_back(i3);
ind_34152.push_back(i4);
ind_34152.push_back(i1);
ind_34152.push_back(i5);
ind_34152.push_back(i2);

process pro_34125=process(hf3,hf4,hf1,hf2,h5);
vector<int> ind_34125;
ind_34125.push_back(i3);
ind_34125.push_back(i4);
ind_34125.push_back(i1);
ind_34125.push_back(i2);
ind_34125.push_back(i5);

//---------------- second quark line arrangement for case of identical quarks

process pro_35214=process(hpf3,hp5,hpf2,hpf1,hpf4);
vector<int> ind_35214;
ind_35214.push_back(i3);
ind_35214.push_back(i5);
ind_35214.push_back(i2);
ind_35214.push_back(i1);
ind_35214.push_back(i4);

process pro_32514=process(hpf3,hpf2,hp5,hpf1,hpf4);
vector<int> ind_32514;
ind_32514.push_back(i3);
ind_32514.push_back(i2);
ind_32514.push_back(i5);
ind_32514.push_back(i1);
ind_32514.push_back(i4);

process pro_32154=process(hpf3,hpf2,hpf1,hp5,hpf4);
vector<int> ind_32154;
ind_32154.push_back(i3);
ind_32154.push_back(i2);
ind_32154.push_back(i1);
ind_32154.push_back(i5);
ind_32154.push_back(i4);

process pro_32145=process(hpf3,hpf2,hpf1,hpf4,hp5);
vector<int> ind_32145;
ind_32145.push_back(i3);
ind_32145.push_back(i2);
ind_32145.push_back(i1);
ind_32145.push_back(i4);
ind_32145.push_back(i5);

process pro_15432=process(hp1,hp5,hp4,hp3,hp2);
vector<int> ind_15432;
ind_15432.push_back(i1);
ind_15432.push_back(i5);
ind_15432.push_back(i4);
ind_15432.push_back(i3);
ind_15432.push_back(i2);

process pro_14532=process(hp1,hp4,hp5,hp3,hp2);
vector<int> ind_14532;
ind_14532.push_back(i1);
ind_14532.push_back(i4);
ind_14532.push_back(i5);
ind_14532.push_back(i3);
ind_14532.push_back(i2);

process pro_14352=process(hp1,hp4,hp3,hp5,hp2);
vector<int> ind_14352;
ind_14352.push_back(i1);
ind_14352.push_back(i4);
ind_14352.push_back(i3);
ind_14352.push_back(i5);
ind_14352.push_back(i2);

process pro_14325=process(hp1,hp4,hp3,hp2,hp5);
vector<int> ind_14325;
ind_14325.push_back(i1);
ind_14325.push_back(i4);
ind_14325.push_back(i3);
ind_14325.push_back(i2);
ind_14325.push_back(i5);


//-------------------------------------------


	size_t T_15234;
	size_t T_12354;
	size_t T_34152;
	size_t T_35412;

//if(tree_color!=1){
	size_t T_12534;
	size_t T_12345;
	size_t T_34512;
	size_t T_34125;
//}

	size_t L1_15234;
	size_t L3_12354;
	size_t L1_35412;
	size_t L3_34152;

	size_t L1_color_15234;
	size_t L3_color_12354;
	size_t L1_color_35412;
	size_t L3_color_34152;

//if(color!=1){
	size_t L2_12534;
	size_t L4_12345;
	size_t L2_34512;
	size_t L4_34125;
//}


	size_t T_15432;
	size_t T_14352;
	size_t T_35214;
	size_t T_32154;

//if(tree_color!=1){
	size_t T_14532;
	size_t T_14325;
	size_t T_32514;
	size_t T_32145;
//}

	size_t L1_15432;
	size_t L3_14352;
	size_t L1_35214;
	size_t L3_32154;

	size_t L1_color_15432;
	size_t L3_color_14352;
	size_t L1_color_35214;
	size_t L3_color_32154;

//if(color!=1){
	size_t L2_14532;
	size_t L4_14325;
	size_t L2_32514;
	size_t L4_32145;
//}


	if(h1.helicity()+h4.helicity()==0){
		T_15234=SM->add(new CTree_with_prefactor(pro_15234,ind_15234));
		T_12354=SM->add(new CTree_with_prefactor(pro_12354,ind_12354));
		T_35412=SM->add(new CTree_with_prefactor(pro_35412,ind_35412));
		T_34152=SM->add(new CTree_with_prefactor(pro_34152,ind_34152));

if(tree_color!=1){
		T_12534=SM->add(new CTree_with_prefactor(pro_12534,ind_12534));
		T_12345=SM->add(new CTree_with_prefactor(pro_12345,ind_12345));
		T_34512=SM->add(new CTree_with_prefactor(pro_34512,ind_34512));
		T_34125=SM->add(new CTree_with_prefactor(pro_34125,ind_34125));
}

if(color!=1){
		L1_15234=SM->add(A_loop_2q_2Q_1g_5_1(pro_15234,ind_15234,n_s,n_f,n_c,0,lo_or_nlo));
		L3_12354=SM->add(A_loop_2q_2Q_1g_5_3(pro_12354,ind_12354,n_s,n_f,n_c,0,lo_or_nlo));
		L1_35412=SM->add(A_loop_2q_2Q_1g_5_1(pro_35412,ind_35412,n_s,n_f,n_c,0,lo_or_nlo));
		L3_34152=SM->add(A_loop_2q_2Q_1g_5_3(pro_34152,ind_34152,n_s,n_f,n_c,0,lo_or_nlo));
}
		L1_color_15234=SM->add(A_loop_2q_2Q_1g_5_1(pro_15234,ind_15234,n_s,n_f,n_c,color,lo_or_nlo));
		L3_color_12354=SM->add(A_loop_2q_2Q_1g_5_3(pro_12354,ind_12354,n_s,n_f,n_c,color,lo_or_nlo));
		L1_color_35412=SM->add(A_loop_2q_2Q_1g_5_1(pro_35412,ind_35412,n_s,n_f,n_c,color,lo_or_nlo));
		L3_color_34152=SM->add(A_loop_2q_2Q_1g_5_3(pro_34152,ind_34152,n_s,n_f,n_c,color,lo_or_nlo));

if(color!=1){
		L2_12534=SM->add(A_loop_2q_2Q_1g_5_2(pro_12534,ind_12534,n_s,n_f,n_c,0,lo_or_nlo));
		L4_12345=SM->add(A_loop_2q_2Q_1g_5_4(pro_12345,ind_12345,n_s,n_f,n_c,0,lo_or_nlo));
		L2_34512=SM->add(A_loop_2q_2Q_1g_5_2(pro_34512,ind_34512,n_s,n_f,n_c,0,lo_or_nlo));
		L4_34125=SM->add(A_loop_2q_2Q_1g_5_4(pro_34125,ind_34125,n_s,n_f,n_c,0,lo_or_nlo));
}

	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_15432=SM->add(new CTree_with_prefactor(pro_15432,ind_15432));
		T_14352=SM->add(new CTree_with_prefactor(pro_14352,ind_14352));
		T_35214=SM->add(new CTree_with_prefactor(pro_35214,ind_35214));
		T_32154=SM->add(new CTree_with_prefactor(pro_32154,ind_32154));

if(tree_color!=1){
		T_14532=SM->add(new CTree_with_prefactor(pro_14532,ind_14532));
		T_14325=SM->add(new CTree_with_prefactor(pro_14325,ind_14325));
		T_32514=SM->add(new CTree_with_prefactor(pro_32514,ind_32514));
		T_32145=SM->add(new CTree_with_prefactor(pro_32145,ind_32145));
}



if(color!=1){
		L1_15432=SM->add(A_loop_2q_2Q_1g_5_1(pro_15432,ind_15432,n_s,n_f,n_c,0,lo_or_nlo));
		L3_14352=SM->add(A_loop_2q_2Q_1g_5_3(pro_14352,ind_14352,n_s,n_f,n_c,0,lo_or_nlo));
		L1_35214=SM->add(A_loop_2q_2Q_1g_5_1(pro_35214,ind_35214,n_s,n_f,n_c,0,lo_or_nlo));
		L3_32154=SM->add(A_loop_2q_2Q_1g_5_3(pro_32154,ind_32154,n_s,n_f,n_c,0,lo_or_nlo));
}
		L1_color_15432=SM->add(A_loop_2q_2Q_1g_5_1(pro_15432,ind_15432,n_s,n_f,n_c,color,lo_or_nlo));
		L3_color_14352=SM->add(A_loop_2q_2Q_1g_5_3(pro_14352,ind_14352,n_s,n_f,n_c,color,lo_or_nlo));
		L1_color_35214=SM->add(A_loop_2q_2Q_1g_5_1(pro_35214,ind_35214,n_s,n_f,n_c,color,lo_or_nlo));
		L3_color_32154=SM->add(A_loop_2q_2Q_1g_5_3(pro_32154,ind_32154,n_s,n_f,n_c,color,lo_or_nlo));


if(color!=1){
		L2_14532=SM->add(A_loop_2q_2Q_1g_5_2(pro_14532,ind_14532,n_s,n_f,n_c,0,lo_or_nlo));
		L4_14325=SM->add(A_loop_2q_2Q_1g_5_4(pro_14325,ind_14325,n_s,n_f,n_c,0,lo_or_nlo));
		L2_32514=SM->add(A_loop_2q_2Q_1g_5_2(pro_32514,ind_32514,n_s,n_f,n_c,0,lo_or_nlo));
		L4_32145=SM->add(A_loop_2q_2Q_1g_5_4(pro_32145,ind_32145,n_s,n_f,n_c,0,lo_or_nlo));
}

	};


//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction nf(n_f,1);
multi_precision_fraction n2m1(n_c*n_c-1,1);
multi_precision_fraction _over_nc(1,n_c);

multi_precision_fraction ext(2*n_c*n_c*(n_c*n_c-1),1);
multi_precision_fraction ext_tree(n_c*(n_c*n_c-1),1);



multi_precision_fraction Qu(2,3);
multi_precision_fraction Qd(-1,3);


//--------------------------------------------------------
// mulitplicities of processes

multi_precision_fraction uu_multiplicity(1,1);//2/2
multi_precision_fraction ud_multiplicity(1,1);//6



////////////////////////////////
////////////////////////////////
// Z/gamma cases
////////////////////////////////
////////////////////////////////
//--------------------------------------------------------
// (u,d), (u,u'), (d,d'), (d,u)
_PRINT(case4q);

{
if(case4q==1 || case4q==2 || case4q==3 || case4q==4 || case4q==-1){
cout<< "new case ud"<<endl;
		multi_precision_fraction ud=ud_multiplicity*ext_tree;
		multi_precision_fraction ud_E=-ud*_over_nc;
		multi_precision_fraction ud_EE=-ud_E*_over_nc;
		multi_precision_fraction ud_EEE=-ud_EE*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234 ,T_15234 ,ud));
		SM->add_tree(cached_cross_term_md(T_12354 ,T_12354 ,ud));
		
		}	

if(tree_color!=1){
_PRINT("subleading tree");

		if((h1.helicity()+h4.helicity()==0)){
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234,T_15234,uu));
			SM->add_tree(cached_cross_term_md(T_15234,T_12534,ud_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234,T_12354,0));
			SM->add_tree(cached_cross_term_md(T_15234,T_12345,ud_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534,T_15234,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534,T_12534,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12534,T_12354,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534,T_12345,ud_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354,T_15234,0));
			SM->add_tree(cached_cross_term_md(T_12354,T_12534,ud_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354,T_12354,uu));
			SM->add_tree(cached_cross_term_md(T_12354,T_12345,ud_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345,T_15234,ud_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345,T_12534,0));
			SM->add_tree(cached_cross_term_md(T_12345,T_12354,ud_EE));
			SM->add_tree(cached_cross_term_md(T_12345,T_12345,ud_EE));
		};
}

		multi_precision_fraction ud_virtual=ud_multiplicity*ext;
		multi_precision_fraction ud_EE_virtual=ud_virtual*_over_nc*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234 ,T_15234 ,ud_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354 ,T_12354 ,ud_virtual));
		};

if(color!=1){
_PRINT("subleading loop");

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234,T_15234,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234,T_12534,ud_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234,T_12354,0));
			SM->add_loop(cached_cross_term_md(L1_15234,T_12345,ud_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534,T_15234,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_12534,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_12354,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534,T_12345,ud_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354,T_15234,0));
			SM->add_loop(cached_cross_term_md(L3_12354,T_12534,ud_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354,T_12354,uu));
			SM->add_loop(cached_cross_term_md(L3_12354,T_12345,ud_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345,T_15234,ud_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345,T_12534,0));
			SM->add_loop(cached_cross_term_md(L4_12345,T_12354,ud_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345,T_12345,ud_EE_virtual));
		};
}


	};


//--------------------------------------------------------
// (u,u)
//--------------------------------------------------------
if(case4q==0 || case4q==5|| case4q==-1){
cout<< "new case uu"<<endl;

		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=-uu*_over_nc;
		multi_precision_fraction uu_EE=-uu_E*_over_nc;
		multi_precision_fraction uu_EEE=-uu_EE*_over_nc;
		
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;
		multi_precision_fraction uu_EE_virtual=-uu_E_virtual*_over_nc;
		multi_precision_fraction uu_EEE_virtual=-uu_EE_virtual*_over_nc;

	if(h1.helicity()+h4.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15234,T_15234,uu));
		SM->add_tree(cached_cross_term_md(T_12354,T_12354,uu));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15234,T_15234,uu));
			SM->add_tree(cached_cross_term_md(T_15234,T_12534,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15234,T_12354,0));
			SM->add_tree(cached_cross_term_md(T_15234,T_12345,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_12534,T_15234,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534,T_12534,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12534,T_12354,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12534,T_12345,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_12354,T_15234,0));
			SM->add_tree(cached_cross_term_md(T_12354,T_12534,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_12354,T_12354,uu));
			SM->add_tree(cached_cross_term_md(T_12354,T_12345,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345,T_15234,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_12345,T_12534,0));
			SM->add_tree(cached_cross_term_md(T_12345,T_12354,uu_EE));
			SM->add_tree(cached_cross_term_md(T_12345,T_12345,uu_EE));
		};
}


	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15234 ,T_15234 ,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_12354 ,T_12354 ,uu_virtual));
		
		};

if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15234,T_15234,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15234,T_12534,uu_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15234,T_12354,0));
			SM->add_loop(cached_cross_term_md(L1_15234,T_12345,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_12534,T_15234,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_12534,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_12354,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_12534,T_12345,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_12354,T_15234,0));
			SM->add_loop(cached_cross_term_md(L3_12354,T_12534,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_12354,T_12354,uu));
			SM->add_loop(cached_cross_term_md(L3_12354,T_12345,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345,T_15234,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_12345,T_12534,0));
			SM->add_loop(cached_cross_term_md(L4_12345,T_12354,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345,T_12345,uu_EE_virtual));
		
		};
}

	if(h1.helicity()+h2.helicity()==0){

		SM->add_tree(cached_cross_term_md(T_15432,T_15432,uu));
		SM->add_tree(cached_cross_term_md(T_14352,T_14352,uu));
		
		}	

if(tree_color!=1){

		if(h1.helicity()+h2.helicity()==0){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cached_cross_term_md(T_15432,T_15432,uu));
			SM->add_tree(cached_cross_term_md(T_15432,T_14532,uu_EE));
			//0: SM->add_tree(cached_cross_term_md(T_15432,T_14352,0));
			SM->add_tree(cached_cross_term_md(T_15432,T_14325,uu_EE));
		//second row
			SM->add_tree(cached_cross_term_md(T_14532,T_15432,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532,T_14532,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14532,T_14352,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14532,T_14325,uu_E));
		//third row
			//0: SM->add_tree(cached_cross_term_md(T_14352,T_15432,0));
			SM->add_tree(cached_cross_term_md(T_14352,T_14532,uu_EE));
			//leading-color: SM->add_tree(cached_cross_term_md(T_14352,T_14352,uu));
			SM->add_tree(cached_cross_term_md(T_14352,T_14325,uu_EE));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325,T_15432,uu_EE));
			//0 SM->add_tree(cached_cross_term_md(T_14325,T_14532,0));
			SM->add_tree(cached_cross_term_md(T_14325,T_14352,uu_EE));
			SM->add_tree(cached_cross_term_md(T_14325,T_14325,uu_EE));
	
		};
}


	if(h1.helicity()+h2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_color_15432,T_15432,uu_virtual));
		SM->add_loop(cached_cross_term_md(L3_color_14352,T_14352,uu_virtual));
		};

if(color!=1){

		
		if(h1.helicity()+h2.helicity()==0){
		
		//block 1
		//first row
			//leading-color SM->add_loop(cached_cross_term_md(L1_15432,T_15432,uu));
		  	SM->add_loop(cached_cross_term_md(L1_15432,T_14532,uu_EE_virtual));
        		//0: SM->add_loop(cached_cross_term_md(L1_15432,T_14352,0));
			SM->add_loop(cached_cross_term_md(L1_15432,T_14325,uu_EE_virtual));
		//second row
			SM->add_loop(cached_cross_term_md(L2_14532,T_15432,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532,T_14532,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532,T_14352,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L2_14532,T_14325,uu_E));
		//third row
			//0: SM->add_loop(cached_cross_term_md(L3_14352,T_15432,0));
			SM->add_loop(cached_cross_term_md(L3_14352,T_14532,uu_EE_virtual));
			//leading-color: SM->add_loop(cached_cross_term_md(L3_14352,T_14352,uu));
			SM->add_loop(cached_cross_term_md(L3_14352,T_14325,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325,T_15432,uu_EE_virtual));
			//0 SM->add_loop(cached_cross_term_md(L4_14325,T_14532,0));
			SM->add_loop(cached_cross_term_md(L4_14325,T_14352,uu_EE_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325,T_14325,uu_EE_virtual));
		};
}

////////////////////////////////
///identical anti-quark exchange term 
///block B
////////////////////////////////


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15432,T_15234,uu_E));
// 3			SM->add_tree(cached_cross_term_md(T_15432,T_12534,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432,T_12354,uu_E));
			SM->add_tree(cached_cross_term_md(T_15432,T_12345,uu_E));
		//second row
// 2			SM->add_tree(cached_cross_term_md(T_14532,T_15234,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532,T_12534,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_14532,T_12354,uu_E));
			SM->add_tree(cached_cross_term_md(T_14532,T_12345,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_14352,T_15234,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352,T_12534,uu_E));
			SM->add_tree(cached_cross_term_md(T_14352,T_12354,uu_E));
// 1			SM->add_tree(cached_cross_term_md(T_14352,T_12345,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_14325,T_15234,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325,T_12534,uu_EEE));
// 4			SM->add_tree(cached_cross_term_md(T_14325,T_12354,uu_E));
			SM->add_tree(cached_cross_term_md(T_14325,T_12345,uu_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cached_cross_term_md(T_15234,T_15432,uu_E));
//			SM->add_tree(cached_cross_term_md(T_15234,T_14532,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234,T_14352,uu_E));
			SM->add_tree(cached_cross_term_md(T_15234,T_14325,uu_E));
		//second row
//			SM->add_tree(cached_cross_term_md(T_12534,T_15432,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534,T_14532,uu_EEE));
			SM->add_tree(cached_cross_term_md(T_12534,T_14352,uu_E));
			SM->add_tree(cached_cross_term_md(T_12534,T_14325,uu_EEE));
		//third row
			SM->add_tree(cached_cross_term_md(T_12354,T_15432,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354,T_14532,uu_E));
			SM->add_tree(cached_cross_term_md(T_12354,T_14352,uu_E));
//			SM->add_tree(cached_cross_term_md(T_12354,T_14325,uu_E));
		//fourth row
		 	SM->add_tree(cached_cross_term_md(T_12345,T_15432,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345,T_14532,uu_EEE));
//			SM->add_tree(cached_cross_term_md(T_12345,T_14352,uu_E));
			SM->add_tree(cached_cross_term_md(T_12345,T_14325,uu_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15234,T_15432,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L1_15234,T_14532,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234,T_14352,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15234,T_14325,uu_E_virtual));
		//second row
			//SM->add_loop(cached_cross_term_md(L2_12534,T_15432,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_14532,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_14352,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_12534,T_14325,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_12354,T_15432,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354,T_14532,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_12354,T_14352,uu_E_virtual));
			//SM->add_loop(cached_cross_term_md(L3_12354,T_14325,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_12345,T_15432,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345,T_14532,uu_EEE_virtual));
			//SM->add_loop(cached_cross_term_md(L4_12345,T_14352,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_12345,T_14325,uu_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cached_cross_term_md(L1_15432,T_15234,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L1_15432,T_12534,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432,T_12354,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L1_15432,T_12345,uu_E_virtual));
		//second row
		//	SM->add_loop(cached_cross_term_md(L2_14532,T_15234,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532,T_12534,uu_EEE_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532,T_12354,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L2_14532,T_12345,uu_EEE_virtual));
		//third row
			SM->add_loop(cached_cross_term_md(L3_14352,T_15234,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352,T_12534,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L3_14352,T_12354,uu_E_virtual));
		//	SM->add_loop(cached_cross_term_md(L3_14352,T_12345,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cached_cross_term_md(L4_14325,T_15234,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325,T_12534,uu_EEE_virtual));
		//	SM->add_loop(cached_cross_term_md(L4_14325,T_12354,uu_E_virtual));
			SM->add_loop(cached_cross_term_md(L4_14325,T_12345,uu_EEE_virtual));
		};

}


}



}

return SM;
}

/*
Virtual_SME* vsme_2q2Q1g(int ns,int nf,int nc, int case4q,int color, int tree_color){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);


clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbp,qm,qbm,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qp,qbm,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qm,qbp,m),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbm,qp,qbp,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qm,qbp,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qp,qbm,m),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbp,qm,qbm,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qp,qbm,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qm,qbp,p),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbm,qp,qbp,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qm,qbp,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qp,qbm,p),ind,ns,nf,nc,case4q,color,tree_color));

after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}
*/
Virtual_SME* vsme_2q2Q1g(std::vector<int> indext,int ns,int nf,int nc, int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);

clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;

	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbp,qm,qbm,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qp,qbm,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qm,qbp,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbm,qp,qbp,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qm,qbp,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qp,qbm,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbp,qm,qbm,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qp,qbm,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qp,qbm,qm,qbp,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbm,qp,qbp,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qm,qbp,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1g_M2(process(qm,qbp,qp,qbm,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}


BH_Ampl_2q2Q1g::BH_Ampl_2q2Q1g(int case4q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q1g(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				case4q,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			5,    // NbrExtParticles
			3,	  // NbrPowersOfAlphaS
			0,    // NbrPowersOfAlphaQED
			2,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
		{}

}
