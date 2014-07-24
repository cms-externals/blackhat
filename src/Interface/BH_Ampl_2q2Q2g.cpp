/*
 * matrix elements for 2q2Q2g
 *
 *  Created on: Oct Apr, 2008
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
// q g g Gb G qb
partial_amplitude_cached* A_loop_2q_2Q_2g_6_12(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}



//expected indices:
// q Gb g g G qb
partial_amplitude_cached* A_loop_2q_2Q_2g_6_23(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}

//expected indices:
// q Gb G g g qb
partial_amplitude_cached* A_loop_2q_2Q_2g_6_34(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}


//expected indices:
// q Gb G qb g g
partial_amplitude_cached* A_loop_2q_2Q_2g_6_41(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}



//expected indices:
// q g Gb G g qb
partial_amplitude_cached* A_loop_2q_2Q_2g_6_13(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}



//expected indices:
// q Gb g G qb g
partial_amplitude_cached* A_loop_2q_2Q_2g_6_24(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}





Squared_ME* A_loop_2q_2Q_2g_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);


// assume input with identical flavor quarks
// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors

//quarks labels i1-i4
//extra gluons called i5-i6


ph_type h1=pro.p(1);
ph_type h2=particle_ID(gluino,pro.p(2).helicity(),1,true);
ph_type h3=particle_ID(gluino,pro.p(3).helicity(),1);
ph_type h4=pro.p(4);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);



ph_type hp1=pro.p(1);
ph_type hp2=pro.p(2);
ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),1);
ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),1,true);
ph_type hp5=pro.p(5);
ph_type hp6=pro.p(6);

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

process pro_156234=process(h1,h5,h6,h2,h3,h4);
vector<int> ind_156234;
ind_156234.push_back(i1);
ind_156234.push_back(i5);
ind_156234.push_back(i6);
ind_156234.push_back(i2);
ind_156234.push_back(i3);
ind_156234.push_back(i4);

process pro_165234=process(h1,h6,h5,h2,h3,h4);
vector<int> ind_165234;
ind_165234.push_back(i1);
ind_165234.push_back(i6);
ind_165234.push_back(i5);
ind_165234.push_back(i2);
ind_165234.push_back(i3);
ind_165234.push_back(i4);

//-

process pro_125634=process(h1,h2,h5,h6,h3,h4);
vector<int> ind_125634;
ind_125634.push_back(i1);
ind_125634.push_back(i2);
ind_125634.push_back(i5);
ind_125634.push_back(i6);
ind_125634.push_back(i3);
ind_125634.push_back(i4);

process pro_126534=process(h1,h2,h6,h5,h3,h4);
vector<int> ind_126534;
ind_126534.push_back(i1);
ind_126534.push_back(i2);
ind_126534.push_back(i6);
ind_126534.push_back(i5);
ind_126534.push_back(i3);
ind_126534.push_back(i4);

//-

process pro_123564=process(h1,h2,h3,h5,h6,h4);
vector<int> ind_123564;
ind_123564.push_back(i1);
ind_123564.push_back(i2);
ind_123564.push_back(i3);
ind_123564.push_back(i5);
ind_123564.push_back(i6);
ind_123564.push_back(i4);

process pro_123654=process(h1,h2,h3,h6,h5,h4);
vector<int> ind_123654;
ind_123654.push_back(i1);
ind_123654.push_back(i2);
ind_123654.push_back(i3);
ind_123654.push_back(i6);
ind_123654.push_back(i5);
ind_123654.push_back(i4);

//-

process pro_123456=process(h1,h2,h3,h4,h5,h6);
vector<int> ind_123456;
ind_123456.push_back(i1);
ind_123456.push_back(i2);
ind_123456.push_back(i3);
ind_123456.push_back(i4);
ind_123456.push_back(i5);
ind_123456.push_back(i6);

process pro_123465=process(h1,h2,h3,h4,h6,h5);
vector<int> ind_123465;
ind_123465.push_back(i1);
ind_123465.push_back(i2);
ind_123465.push_back(i3);
ind_123465.push_back(i4);
ind_123465.push_back(i6);
ind_123465.push_back(i5);

//-

process pro_152364=process(h1,h5,h2,h3,h6,h4);
vector<int> ind_152364;
ind_152364.push_back(i1);
ind_152364.push_back(i5);
ind_152364.push_back(i2);
ind_152364.push_back(i3);
ind_152364.push_back(i6);
ind_152364.push_back(i4);

process pro_162354=process(h1,h6,h2,h3,h5,h4);
vector<int> ind_162354;
ind_162354.push_back(i1);
ind_162354.push_back(i6);
ind_162354.push_back(i2);
ind_162354.push_back(i3);
ind_162354.push_back(i5);
ind_162354.push_back(i4);

//-

process pro_125346=process(h1,h2,h5,h3,h4,h6);
vector<int> ind_125346;
ind_125346.push_back(i1);
ind_125346.push_back(i2);
ind_125346.push_back(i5);
ind_125346.push_back(i3);
ind_125346.push_back(i4);
ind_125346.push_back(i6);

process pro_126345=process(h1,h2,h6,h3,h4,h5);
vector<int> ind_126345;
ind_126345.push_back(i1);
ind_126345.push_back(i2);
ind_126345.push_back(i6);
ind_126345.push_back(i3);
ind_126345.push_back(i4);
ind_126345.push_back(i5);

//---------------- second quark line arrangement for case of identical quarks
//not checked
process pro_356214=process(hpf3,h5,h6,hpf2,hpf1,hpf4);
vector<int> ind_356214;
ind_356214.push_back(i3);
ind_356214.push_back(i5);
ind_356214.push_back(i6);
ind_356214.push_back(i2);
ind_356214.push_back(i1);
ind_356214.push_back(i4);

process pro_365214=process(hpf3,h6,h5,hpf2,hpf1,hpf4);
vector<int> ind_365214;
ind_365214.push_back(i3);
ind_365214.push_back(i6);
ind_365214.push_back(i5);
ind_365214.push_back(i2);
ind_365214.push_back(i1);
ind_365214.push_back(i4);

//-

process pro_325614=process(hpf3,hpf2,h5,h6,hpf1,hpf4);
vector<int> ind_325614;
ind_325614.push_back(i3);
ind_325614.push_back(i2);
ind_325614.push_back(i5);
ind_325614.push_back(i6);
ind_325614.push_back(i1);
ind_325614.push_back(i4);

process pro_326514=process(hpf3,hpf2,h6,h5,hpf1,hpf4);
vector<int> ind_326514;
ind_326514.push_back(i3);
ind_326514.push_back(i2);
ind_326514.push_back(i6);
ind_326514.push_back(i5);
ind_326514.push_back(i1);
ind_326514.push_back(i4);

//-

process pro_321564=process(hpf3,hpf2,hpf1,h5,h6,hpf4);
vector<int> ind_321564;
ind_321564.push_back(i3);
ind_321564.push_back(i2);
ind_321564.push_back(i1);
ind_321564.push_back(i5);
ind_321564.push_back(i6);
ind_321564.push_back(i4);

process pro_321654=process(hpf3,hpf2,hpf1,h6,h5,hpf4);
vector<int> ind_321654;
ind_321654.push_back(i3);
ind_321654.push_back(i2);
ind_321654.push_back(i1);
ind_321654.push_back(i6);
ind_321654.push_back(i5);
ind_321654.push_back(i4);

//-


process pro_321456=process(hpf3,hpf2,hpf1,hpf4,h5,h6);
vector<int> ind_321456;
ind_321456.push_back(i3);
ind_321456.push_back(i2);
ind_321456.push_back(i1);
ind_321456.push_back(i4);
ind_321456.push_back(i5);
ind_321456.push_back(i6);


process pro_321465=process(hpf3,hpf2,hpf1,hpf4,h6,h5);
vector<int> ind_321465;
ind_321465.push_back(i3);
ind_321465.push_back(i2);
ind_321465.push_back(i1);
ind_321465.push_back(i4);
ind_321465.push_back(i6);
ind_321465.push_back(i5);


//-

process pro_352164=process(h3,h5,h2,h1,h6,h4);
vector<int> ind_352164;
ind_352164.push_back(i3);
ind_352164.push_back(i5);
ind_352164.push_back(i2);
ind_352164.push_back(i1);
ind_352164.push_back(i6);
ind_352164.push_back(i4);

process pro_362154=process(h3,h6,h2,h1,h5,h4);
vector<int> ind_362154;
ind_362154.push_back(i3);
ind_362154.push_back(i6);
ind_362154.push_back(i2);
ind_362154.push_back(i1);
ind_362154.push_back(i5);
ind_362154.push_back(i4);

//-

process pro_325146=process(h3,h2,h5,h1,h4,h6);
vector<int> ind_325146;
ind_325146.push_back(i3);
ind_325146.push_back(i2);
ind_325146.push_back(i5);
ind_325146.push_back(i1);
ind_325146.push_back(i4);
ind_325146.push_back(i6);

process pro_326145=process(h3,h2,h6,h1,h4,h5);
vector<int> ind_326145;
ind_326145.push_back(i3);
ind_326145.push_back(i2);
ind_326145.push_back(i6);
ind_326145.push_back(i1);
ind_326145.push_back(i4);
ind_326145.push_back(i5);

//-------------------------------------------


	size_t T_156234;
	size_t T_165234;
	size_t T_123564;
	size_t T_123654;

	size_t T_125634;
	size_t T_126534;
	size_t T_123456;
	size_t T_123465;

	size_t T_152364;
	size_t T_162354;
	size_t T_125346;
	size_t T_126345;

//-

	size_t L12_156234;
	size_t L12_165234;
	size_t L34_123564;
	size_t L34_123654;

	size_t L12_color_156234;
	size_t L12_color_165234;
	size_t L34_color_123564;
	size_t L34_color_123654;

	size_t L23_125634;
	size_t L23_126534;
	size_t L41_123456;
	size_t L41_123465;

	size_t L13_152364;
	size_t L13_162354;
	size_t L24_125346;
	size_t L24_126345;


	size_t L13_color_152364;
	size_t L13_color_162354;
//------------------

	size_t T_356214;
	size_t T_365214;
	size_t T_321564;
	size_t T_321654;

	size_t T_325614;
	size_t T_326514;
	size_t T_321456;
	size_t T_321465;

	size_t T_352164;
	size_t T_362154;
	size_t T_325146;
	size_t T_326145;

//-

	size_t L12_356214;
	size_t L12_365214;
	size_t L34_321564;
	size_t L34_321654;

	size_t L12_color_356214;
	size_t L12_color_365214;
	size_t L34_color_321564;
	size_t L34_color_321654;

	size_t L23_325614;
	size_t L23_326514;
	size_t L41_321456;
	size_t L41_321465;

	size_t L13_352164;
	size_t L13_362154;
	size_t L24_325146;
	size_t L24_326145;

	size_t L13_color_352164;
	size_t L13_color_362154;
	size_t L24_color_325146;
	size_t L24_color_326145;
//------------------

	if(h1.helicity()+h4.helicity()==0){
		T_156234=SM->add(new TreeHelAmpl(pro_156234),ind_156234);
		T_165234=SM->add(new TreeHelAmpl(pro_165234),ind_165234);
		T_123564=SM->add(new TreeHelAmpl(pro_123564),ind_123564);
		T_123654=SM->add(new TreeHelAmpl(pro_123654),ind_123654);
		
		T_152364=SM->add(new TreeHelAmpl(pro_152364),ind_152364);
		T_162354=SM->add(new TreeHelAmpl(pro_162354),ind_162354);

if(tree_color!=1){
		T_125634=SM->add(new TreeHelAmpl(pro_125634),ind_125634);
		T_126534=SM->add(new TreeHelAmpl(pro_126534),ind_126534);
		T_123456=SM->add(new TreeHelAmpl(pro_123456),ind_123456);
		T_123465=SM->add(new TreeHelAmpl(pro_123465),ind_123465);

		T_125346=SM->add(new TreeHelAmpl(pro_125346),ind_126345);
		T_126345=SM->add(new TreeHelAmpl(pro_126345),ind_125346);


}
		L12_color_156234=SM->add(A_loop_2q_2Q_2g_6_12(pro_156234,ind_156234,n_s,n_f,n_c,color,lo_or_nlo));
		L12_color_165234=SM->add(A_loop_2q_2Q_2g_6_12(pro_165234,ind_165234,n_s,n_f,n_c,color,lo_or_nlo));
		L34_color_123564=SM->add(A_loop_2q_2Q_2g_6_34(pro_123564,ind_123564,n_s,n_f,n_c,color,lo_or_nlo));
		L34_color_123654=SM->add(A_loop_2q_2Q_2g_6_34(pro_123654,ind_123654,n_s,n_f,n_c,color,lo_or_nlo));

if(color!=1){
		L12_156234=SM->add(A_loop_2q_2Q_2g_6_12(pro_156234,ind_156234,n_s,n_f,n_c,0,lo_or_nlo));
		L12_165234=SM->add(A_loop_2q_2Q_2g_6_12(pro_165234,ind_165234,n_s,n_f,n_c,0,lo_or_nlo));
		L34_123564=SM->add(A_loop_2q_2Q_2g_6_34(pro_123564,ind_123564,n_s,n_f,n_c,0,lo_or_nlo));
		L34_123654=SM->add(A_loop_2q_2Q_2g_6_34(pro_123654,ind_123654,n_s,n_f,n_c,0,lo_or_nlo));

		L23_125634=SM->add(A_loop_2q_2Q_2g_6_23(pro_125634,ind_125634,n_s,n_f,n_c,0,lo_or_nlo));
		L23_126534=SM->add(A_loop_2q_2Q_2g_6_23(pro_126534,ind_126534,n_s,n_f,n_c,0,lo_or_nlo));
		L41_123456=SM->add(A_loop_2q_2Q_2g_6_41(pro_123456,ind_123456,n_s,n_f,n_c,0,lo_or_nlo));
		L41_123465=SM->add(A_loop_2q_2Q_2g_6_41(pro_123465,ind_123465,n_s,n_f,n_c,0,lo_or_nlo));
}

		L13_color_152364=SM->add(A_loop_2q_2Q_2g_6_13(pro_152364,ind_152364,n_s,n_f,n_c,color,lo_or_nlo));
		L13_color_162354=SM->add(A_loop_2q_2Q_2g_6_13(pro_162354,ind_162354,n_s,n_f,n_c,color,lo_or_nlo));
if(color!=1){
		L13_152364=SM->add(A_loop_2q_2Q_2g_6_13(pro_152364,ind_152364,n_s,n_f,n_c,0,lo_or_nlo));
		L13_162354=SM->add(A_loop_2q_2Q_2g_6_13(pro_162354,ind_162354,n_s,n_f,n_c,0,lo_or_nlo));
		L24_125346=SM->add(A_loop_2q_2Q_2g_6_24(pro_125346,ind_125346,n_s,n_f,n_c,0,lo_or_nlo));
		L24_126345=SM->add(A_loop_2q_2Q_2g_6_24(pro_126345,ind_126345,n_s,n_f,n_c,0,lo_or_nlo));
}

	};

//-------------------------

	if(hp1.helicity()+hp2.helicity()==0){
		T_356214=SM->add(new TreeHelAmpl(pro_356214),ind_356214);
		T_365214=SM->add(new TreeHelAmpl(pro_365214),ind_365214);
		T_321564=SM->add(new TreeHelAmpl(pro_321564),ind_321564);
		T_321654=SM->add(new TreeHelAmpl(pro_321654),ind_321654);

if(tree_color!=1){
		T_325614=SM->add(new TreeHelAmpl(pro_325614),ind_325614);
		T_326514=SM->add(new TreeHelAmpl(pro_326514),ind_326514);
		T_321456=SM->add(new TreeHelAmpl(pro_321456),ind_321456);
		T_321465=SM->add(new TreeHelAmpl(pro_321465),ind_321465);

		T_352164=SM->add(new TreeHelAmpl(pro_352164),ind_352164);
		T_362154=SM->add(new TreeHelAmpl(pro_362154),ind_362154);
		T_325146=SM->add(new TreeHelAmpl(pro_325146),ind_325146);
		T_326145=SM->add(new TreeHelAmpl(pro_326145),ind_326145);

}



if(color!=1){
		L12_356214=SM->add(A_loop_2q_2Q_2g_6_12(pro_356214,ind_356214,n_s,n_f,n_c,0,lo_or_nlo));
		L12_365214=SM->add(A_loop_2q_2Q_2g_6_12(pro_365214,ind_365214,n_s,n_f,n_c,0,lo_or_nlo));
		L34_321564=SM->add(A_loop_2q_2Q_2g_6_34(pro_321564,ind_321564,n_s,n_f,n_c,0,lo_or_nlo));
		L34_321654=SM->add(A_loop_2q_2Q_2g_6_34(pro_321654,ind_321654,n_s,n_f,n_c,0,lo_or_nlo));
}
		L12_color_356214=SM->add(A_loop_2q_2Q_2g_6_12(pro_356214,ind_356214,n_s,n_f,n_c,color,lo_or_nlo));
		L12_color_365214=SM->add(A_loop_2q_2Q_2g_6_12(pro_365214,ind_365214,n_s,n_f,n_c,color,lo_or_nlo));
		L34_color_321564=SM->add(A_loop_2q_2Q_2g_6_34(pro_321564,ind_321564,n_s,n_f,n_c,color,lo_or_nlo));
		L34_color_321654=SM->add(A_loop_2q_2Q_2g_6_34(pro_321654,ind_321654,n_s,n_f,n_c,color,lo_or_nlo));


if(color!=1){
		L23_325614=SM->add(A_loop_2q_2Q_2g_6_23(pro_325614,ind_325614,n_s,n_f,n_c,0,lo_or_nlo));
		L23_326514=SM->add(A_loop_2q_2Q_2g_6_23(pro_326514,ind_326514,n_s,n_f,n_c,0,lo_or_nlo));
		L41_321456=SM->add(A_loop_2q_2Q_2g_6_41(pro_321456,ind_321456,n_s,n_f,n_c,0,lo_or_nlo));
		L41_321465=SM->add(A_loop_2q_2Q_2g_6_41(pro_321465,ind_321465,n_s,n_f,n_c,0,lo_or_nlo));

		L13_352164=SM->add(A_loop_2q_2Q_2g_6_13(pro_352164,ind_352164,n_s,n_f,n_c,0,lo_or_nlo));
		L13_362154=SM->add(A_loop_2q_2Q_2g_6_13(pro_362154,ind_362154,n_s,n_f,n_c,0,lo_or_nlo));
		L24_325146=SM->add(A_loop_2q_2Q_2g_6_24(pro_325146,ind_325146,n_s,n_f,n_c,0,lo_or_nlo));
		L24_326145=SM->add(A_loop_2q_2Q_2g_6_24(pro_326145,ind_326145,n_s,n_f,n_c,0,lo_or_nlo));

}

		L13_color_352164=SM->add(A_loop_2q_2Q_2g_6_13(pro_352164,ind_352164,n_s,n_f,n_c,color,lo_or_nlo));
		L13_color_362154=SM->add(A_loop_2q_2Q_2g_6_13(pro_362154,ind_362154,n_s,n_f,n_c,color,lo_or_nlo));

	};


//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction nf(n_f,1);
multi_precision_fraction n2m1(n_c*n_c-1,1);
multi_precision_fraction _over_nc(1,n_c);

multi_precision_fraction n2=nc*nc;
multi_precision_fraction n3=nc*n2;
multi_precision_fraction one(1,1);
multi_precision_fraction mone(-1,1);
multi_precision_fraction mn2=mone*n2;

multi_precision_fraction n2_times_n2m1=n2*n2m1;
multi_precision_fraction n3_times_n2m1=n3*n2m1;

/* notice: problem with integers becoming too large in the numerical n_c->infinity limit
?? factor (n2m1/n2) suppressed in born part
?? factor (n2m1/n2) suppressed in born part
*/



//--------------------------------------------------------
// mulitplicities of processes


////////////////////////////////
////////////////////////////////
//--------------------------------------------------------
// (u,d), (u,u'), (d,d'), (d,u)
{
if(case4q==1 || case4q==2 || case4q==3 || case4q==4 || case4q==-1){
cout<< "new case ud"<<endl;


	if(h1.helicity()+h4.helicity()==0){

SM->add_tree(cross_term(T_156234,T_156234,n2_times_n2m1));
SM->add_tree(cross_term(T_165234,T_165234,n2_times_n2m1));

SM->add_tree(cross_term(T_123564,T_123564,n2_times_n2m1));
SM->add_tree(cross_term(T_123654,T_123654,n2_times_n2m1));

SM->add_tree(cross_term(T_152364,T_152364,n2_times_n2m1));
SM->add_tree(cross_term(T_162354,T_162354,n2_times_n2m1));

		
		}	

if(tree_color!=1){
//subleading tree;

		if((h1.helicity()+h4.helicity()==0)){

//  SM->add_tree(cross_term(T_156234,T_156234,n2_times_n2m1));
  SM->add_tree(cross_term(T_156234,T_165234,mn2));
  SM->add_tree(cross_term(T_156234,T_123564,n2));
  SM->add_tree(cross_term(T_156234,T_123654,n2));
  SM->add_tree(cross_term(T_156234,T_125634,n2m1));
  SM->add_tree(cross_term(T_156234,T_126534,mone));
  SM->add_tree(cross_term(T_156234,T_123456,n2m1));
  SM->add_tree(cross_term(T_156234,T_123465,mone));
//SM->add_tree(cross_term(T_156234,T_152364,0));
//SM->add_tree(cross_term(T_156234,T_162354,0));
  SM->add_tree(cross_term(T_156234,T_125346,mone));
  SM->add_tree(cross_term(T_156234,T_126345,n2m1));

  SM->add_tree(cross_term(T_165234,T_156234,mn2));
//  SM->add_tree(cross_term(T_165234,T_165234,n2_times_n2m1));
  SM->add_tree(cross_term(T_165234,T_123564,n2));
  SM->add_tree(cross_term(T_165234,T_123654,n2));
  SM->add_tree(cross_term(T_165234,T_125634,mone));
  SM->add_tree(cross_term(T_165234,T_126534,n2m1));
  SM->add_tree(cross_term(T_165234,T_123456,mone));
  SM->add_tree(cross_term(T_165234,T_123465,n2m1));
//SM->add_tree(cross_term(T_165234,T_152364,0));
//SM->add_tree(cross_term(T_165234,T_162354,0));
  SM->add_tree(cross_term(T_165234,T_125346,n2m1));
  SM->add_tree(cross_term(T_165234,T_126345,mone));

//-

  SM->add_tree(cross_term(T_123564,T_156234,n2));
  SM->add_tree(cross_term(T_123564,T_165234,n2));
//SM->add_tree(cross_term(T_123564,T_123564,n2_times_n2m1));
  SM->add_tree(cross_term(T_123564,T_123654,mn2));
  SM->add_tree(cross_term(T_123564,T_125634,n2m1));
  SM->add_tree(cross_term(T_123564,T_126534,mone));
  SM->add_tree(cross_term(T_123564,T_123456,n2m1));
  SM->add_tree(cross_term(T_123564,T_123465,mone));
//SM->add_tree(cross_term(T_123564,T_152364,0));
//SM->add_tree(cross_term(T_123564,T_162354,0));
  SM->add_tree(cross_term(T_123564,T_125346,n2m1));
  SM->add_tree(cross_term(T_123564,T_126345,mone));

  SM->add_tree(cross_term(T_123654,T_156234,n2));
  SM->add_tree(cross_term(T_123654,T_165234,n2));
  SM->add_tree(cross_term(T_123654,T_123564,mn2));
//SM->add_tree(cross_term(T_123654,T_123654,n2_times_n2m1));
  SM->add_tree(cross_term(T_123654,T_125634,mone));
  SM->add_tree(cross_term(T_123654,T_126534,n2m1));
  SM->add_tree(cross_term(T_123654,T_123456,mone));
  SM->add_tree(cross_term(T_123654,T_123465,n2m1));
//SM->add_tree(cross_term(T_123654,T_152364,0));
//SM->add_tree(cross_term(T_123654,T_162354,0));
  SM->add_tree(cross_term(T_123654,T_125346,mone));
  SM->add_tree(cross_term(T_123654,T_126345,n2m1));

//--

  SM->add_tree(cross_term(T_125634,T_156234,n2m1));
  SM->add_tree(cross_term(T_125634,T_165234,mone));
  SM->add_tree(cross_term(T_125634,T_123564,n2m1));
  SM->add_tree(cross_term(T_125634,T_123654,mone));
  SM->add_tree(cross_term(T_125634,T_125634,n2m1));
  SM->add_tree(cross_term(T_125634,T_126534,mone));
  SM->add_tree(cross_term(T_125634,T_123456,one));
  SM->add_tree(cross_term(T_125634,T_123465,one));
  SM->add_tree(cross_term(T_125634,T_152364,mone));
  SM->add_tree(cross_term(T_125634,T_162354,n2m1));
//SM->add_tree(cross_term(T_125634,T_125346,0));
//SM->add_tree(cross_term(T_125634,T_126345,0));

  SM->add_tree(cross_term(T_126534,T_156234,mone));
  SM->add_tree(cross_term(T_126534,T_165234,n2m1));
  SM->add_tree(cross_term(T_126534,T_123564,mone));
  SM->add_tree(cross_term(T_126534,T_123654,n2m1));
  SM->add_tree(cross_term(T_126534,T_125634,mone));
  SM->add_tree(cross_term(T_126534,T_126534,n2m1));
  SM->add_tree(cross_term(T_126534,T_123456,one));
  SM->add_tree(cross_term(T_126534,T_123465,one));
  SM->add_tree(cross_term(T_126534,T_152364,n2m1));
  SM->add_tree(cross_term(T_126534,T_162354,mone));
//SM->add_tree(cross_term(T_126534,T_125346,0));
//SM->add_tree(cross_term(T_126534,T_126345,0));

//---------------------

  SM->add_tree(cross_term(T_123456,T_156234,n2m1));
  SM->add_tree(cross_term(T_123456,T_165234,mone));
  SM->add_tree(cross_term(T_123456,T_123564,n2m1));
  SM->add_tree(cross_term(T_123456,T_123654,mone));
  SM->add_tree(cross_term(T_123456,T_125634,one));
  SM->add_tree(cross_term(T_123456,T_126534,one));
  SM->add_tree(cross_term(T_123456,T_123456,n2m1));
  SM->add_tree(cross_term(T_123456,T_123465,mone));
  SM->add_tree(cross_term(T_123456,T_152364,n2m1));
  SM->add_tree(cross_term(T_123456,T_162354,mone));
//SM->add_tree(cross_term(T_123456,T_125346,0));
//SM->add_tree(cross_term(T_123456,T_126345,0));

  SM->add_tree(cross_term(T_123465,T_156234,mone));
  SM->add_tree(cross_term(T_123465,T_165234,n2m1));
  SM->add_tree(cross_term(T_123465,T_123564,mone));
  SM->add_tree(cross_term(T_123465,T_123654,n2m1));
  SM->add_tree(cross_term(T_123465,T_125634,one));
  SM->add_tree(cross_term(T_123465,T_126534,one));
  SM->add_tree(cross_term(T_123465,T_123456,mone));
  SM->add_tree(cross_term(T_123465,T_123465,n2m1));
  SM->add_tree(cross_term(T_123465,T_152364,mone));
  SM->add_tree(cross_term(T_123465,T_162354,n2m1));
//SM->add_tree(cross_term(T_123465,T_125346,0));
//SM->add_tree(cross_term(T_123465,T_126345,0));

//--

//SM->add_tree(cross_term(T_152364,T_156234,0));
//SM->add_tree(cross_term(T_152364,T_165234,0));
//SM->add_tree(cross_term(T_152364,T_123564,0));
//SM->add_tree(cross_term(T_152364,T_123654,0));
  SM->add_tree(cross_term(T_152364,T_125634,mone));
  SM->add_tree(cross_term(T_152364,T_126534,n2m1));
  SM->add_tree(cross_term(T_152364,T_123456,n2m1));
  SM->add_tree(cross_term(T_152364,T_123465,mone));
//SM->add_tree(cross_term(T_152364,T_152364,n2_times_n2m1));
  SM->add_tree(cross_term(T_152364,T_162354,n2));
  SM->add_tree(cross_term(T_152364,T_125346,n2m1));
  SM->add_tree(cross_term(T_152364,T_126345,n2m1));

//SM->add_tree(cross_term(T_162354,T_156234,0));
//SM->add_tree(cross_term(T_162354,T_165234,0));
//SM->add_tree(cross_term(T_162354,T_123564,0));
//SM->add_tree(cross_term(T_162354,T_123654,0));
  SM->add_tree(cross_term(T_162354,T_125634,n2m1));
  SM->add_tree(cross_term(T_162354,T_126534,mone));
  SM->add_tree(cross_term(T_162354,T_123456,mone));
  SM->add_tree(cross_term(T_162354,T_123465,n2m1));
  SM->add_tree(cross_term(T_162354,T_152364,n2));
//SM->add_tree(cross_term(T_162354,T_162354,n2_times_n2m1));
  SM->add_tree(cross_term(T_162354,T_125346,n2m1));
  SM->add_tree(cross_term(T_162354,T_126345,n2m1));

//--

  SM->add_tree(cross_term(T_125346,T_156234,mone));
  SM->add_tree(cross_term(T_125346,T_165234,n2m1));
  SM->add_tree(cross_term(T_125346,T_123564,n2m1));
  SM->add_tree(cross_term(T_125346,T_123654,mone));
//SM->add_tree(cross_term(T_125346,T_125634,0));
//SM->add_tree(cross_term(T_125346,T_126534,0));
//SM->add_tree(cross_term(T_125346,T_123456,0));
//SM->add_tree(cross_term(T_125346,T_123465,0));
  SM->add_tree(cross_term(T_125346,T_152364,n2m1));
  SM->add_tree(cross_term(T_125346,T_162354,n2m1));
  SM->add_tree(cross_term(T_125346,T_125346,n2m1));
  SM->add_tree(cross_term(T_125346,T_126345,one));

  SM->add_tree(cross_term(T_126345,T_156234,n2m1));
  SM->add_tree(cross_term(T_126345,T_165234,mone));
  SM->add_tree(cross_term(T_126345,T_123564,mone));
  SM->add_tree(cross_term(T_126345,T_123654,n2m1));
//SM->add_tree(cross_term(T_126345,T_125634,0));
//SM->add_tree(cross_term(T_126345,T_126534,0));
//SM->add_tree(cross_term(T_126345,T_123456,0));
//SM->add_tree(cross_term(T_126345,T_123465,0));
  SM->add_tree(cross_term(T_126345,T_152364,n2m1));
  SM->add_tree(cross_term(T_126345,T_162354,n2m1));
  SM->add_tree(cross_term(T_126345,T_125346,one));
  SM->add_tree(cross_term(T_126345,T_126345,n2m1));


};
}


	if(h1.helicity()+h4.helicity()==0){

  SM->add_loop(cross_term(L12_color_156234,T_156234,n3_times_n2m1));
  SM->add_loop(cross_term(L12_color_165234,T_165234,n3_times_n2m1));

  SM->add_loop(cross_term(L34_color_123564,T_123564,n3_times_n2m1));
  SM->add_loop(cross_term(L34_color_123654,T_123654,n3_times_n2m1));

  SM->add_loop(cross_term(L13_color_152364,T_152364,n3_times_n2m1));
  SM->add_loop(cross_term(L13_color_162354,T_162354,n3_times_n2m1));


		};

if(color!=1){
//subleading loop

		if((h1.helicity()+h4.helicity()==0)){

//  SM->add_loop(cross_term(L12_156234,T_156234,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L12_156234,T_165234,mn2*nc));
  SM->add_loop(cross_term(L12_156234,T_123564,n2*nc));
  SM->add_loop(cross_term(L12_156234,T_123654,n2*nc));
  SM->add_loop(cross_term(L12_156234,T_125634,n2m1*nc));
  SM->add_loop(cross_term(L12_156234,T_126534,mone*nc));
  SM->add_loop(cross_term(L12_156234,T_123456,n2m1*nc));
  SM->add_loop(cross_term(L12_156234,T_123465,mone*nc));
//SM->add_loop(cross_term(L12_156234,T_152364,0*nc));
//SM->add_loop(cross_term(L12_156234,T_162354,0*nc));
  SM->add_loop(cross_term(L12_156234,T_125346,mone*nc));
  SM->add_loop(cross_term(L12_156234,T_126345,n2m1*nc));
  
  SM->add_loop(cross_term(L12_165234,T_156234,mn2*nc));
//  SM->add_loop(cross_term(L12_165234,T_165234,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L12_165234,T_123564,n2*nc));
  SM->add_loop(cross_term(L12_165234,T_123654,n2*nc));
  SM->add_loop(cross_term(L12_165234,T_125634,mone*nc));
  SM->add_loop(cross_term(L12_165234,T_126534,n2m1*nc));
  SM->add_loop(cross_term(L12_165234,T_123456,mone*nc));
  SM->add_loop(cross_term(L12_165234,T_123465,n2m1*nc));
//SM->add_loop(cross_term(L12_165234,T_152364,0*nc));
//SM->add_loop(cross_term(L12_165234,T_162354,0*nc));
  SM->add_loop(cross_term(L12_165234,T_125346,n2m1*nc));
  SM->add_loop(cross_term(L12_165234,T_126345,mone*nc));
  
//-
 
  SM->add_loop(cross_term(L34_123564,T_156234,n2*nc));
  SM->add_loop(cross_term(L34_123564,T_165234,n2*nc));
//  SM->add_loop(cross_term(L34_123564,T_123564,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L34_123564,T_123654,mn2*nc));
  SM->add_loop(cross_term(L34_123564,T_125634,n2m1*nc));
  SM->add_loop(cross_term(L34_123564,T_126534,mone*nc));
  SM->add_loop(cross_term(L34_123564,T_123456,n2m1*nc));
  SM->add_loop(cross_term(L34_123564,T_123465,mone*nc));
//SM->add_loop(cross_term(L34_123564,T_152364,0*nc));
//SM->add_loop(cross_term(L34_123564,T_162354,0*nc));
  SM->add_loop(cross_term(L34_123564,T_125346,n2m1*nc));
  SM->add_loop(cross_term(L34_123564,T_126345,mone*nc));


  SM->add_loop(cross_term(L34_123654,T_156234,n2*nc));
  SM->add_loop(cross_term(L34_123654,T_165234,n2*nc));
  SM->add_loop(cross_term(L34_123654,T_123564,mn2*nc));
// SM->add_loop(cross_term(L34_123654,T_123654,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L34_123654,T_125634,mone*nc));
  SM->add_loop(cross_term(L34_123654,T_126534,n2m1*nc));
  SM->add_loop(cross_term(L34_123654,T_123456,mone*nc));
  SM->add_loop(cross_term(L34_123654,T_123465,n2m1*nc));
//SM->add_loop(cross_term(L34_123654,T_152364,0*nc));
//SM->add_loop(cross_term(L34_123654,T_162354,0*nc));
  SM->add_loop(cross_term(L34_123654,T_125346,mone*nc));
  SM->add_loop(cross_term(L34_123654,T_126345,n2m1*nc));

//-

  SM->add_loop(cross_term(L23_125634,T_156234,n2m1*nc));
  SM->add_loop(cross_term(L23_125634,T_165234,mone*nc));
  SM->add_loop(cross_term(L23_125634,T_123564,n2m1*nc));
  SM->add_loop(cross_term(L23_125634,T_123654,mone*nc));
  SM->add_loop(cross_term(L23_125634,T_125634,n2m1*nc));
  SM->add_loop(cross_term(L23_125634,T_126534,mone*nc));
  SM->add_loop(cross_term(L23_125634,T_123456,one*nc));
  SM->add_loop(cross_term(L23_125634,T_123465,one*nc));
  SM->add_loop(cross_term(L23_125634,T_152364,mone*nc));
  SM->add_loop(cross_term(L23_125634,T_162354,n2m1*nc));
//SM->add_loop(cross_term(L23_125634,T_125346,0*nc));
//SM->add_loop(cross_term(L23_125634,T_126345,0*nc));

  SM->add_loop(cross_term(L23_126534,T_156234,mone*nc));
  SM->add_loop(cross_term(L23_126534,T_165234,n2m1*nc));
  SM->add_loop(cross_term(L23_126534,T_123564,mone*nc));
  SM->add_loop(cross_term(L23_126534,T_123654,n2m1*nc));
  SM->add_loop(cross_term(L23_126534,T_125634,mone*nc));
  SM->add_loop(cross_term(L23_126534,T_126534,n2m1*nc));
  SM->add_loop(cross_term(L23_126534,T_123456,one*nc));
  SM->add_loop(cross_term(L23_126534,T_123465,one*nc));
  SM->add_loop(cross_term(L23_126534,T_152364,n2m1*nc));
  SM->add_loop(cross_term(L23_126534,T_162354,mone*nc));
//SM->add_loop(cross_term(L23_126534,T_125346,0*nc));
//SM->add_loop(cross_term(L23_126534,T_126345,0*nc));

//-

  SM->add_loop(cross_term(L41_123456,T_156234,n2m1*nc));
  SM->add_loop(cross_term(L41_123456,T_165234,mone*nc));
  SM->add_loop(cross_term(L41_123456,T_123564,n2m1*nc));
  SM->add_loop(cross_term(L41_123456,T_123654,mone*nc));
  SM->add_loop(cross_term(L41_123456,T_125634,one*nc));
  SM->add_loop(cross_term(L41_123456,T_126534,one*nc));
  SM->add_loop(cross_term(L41_123456,T_123456,n2m1*nc));
  SM->add_loop(cross_term(L41_123456,T_123465,mone*nc));
  SM->add_loop(cross_term(L41_123456,T_152364,n2m1*nc));
  SM->add_loop(cross_term(L41_123456,T_162354,mone*nc));
//SM->add_loop(cross_term(L41_123456,T_125346,0*nc));
//SM->add_loop(cross_term(L41_123456,T_126345,0*nc));

  SM->add_loop(cross_term(L41_123465,T_156234,mone*nc));
  SM->add_loop(cross_term(L41_123465,T_165234,n2m1*nc));
  SM->add_loop(cross_term(L41_123465,T_123564,mone*nc));
  SM->add_loop(cross_term(L41_123465,T_123654,n2m1*nc));
  SM->add_loop(cross_term(L41_123465,T_125634,one*nc));
  SM->add_loop(cross_term(L41_123465,T_126534,one*nc));
  SM->add_loop(cross_term(L41_123465,T_123456,mone*nc));
  SM->add_loop(cross_term(L41_123465,T_123465,n2m1*nc));
  SM->add_loop(cross_term(L41_123465,T_152364,mone*nc));
  SM->add_loop(cross_term(L41_123465,T_162354,n2m1*nc));
//SM->add_loop(cross_term(L41_123465,T_125346,0*nc));
//SM->add_loop(cross_term(L41_123465,T_126345,0*nc));

//-

//SM->add_loop(cross_term(L13_152364,T_156234,0*nc));
//SM->add_loop(cross_term(L13_152364,T_165234,0*nc));
//SM->add_loop(cross_term(L13_152364,T_123564,0*nc));
//SM->add_loop(cross_term(L13_152364,T_123654,0*nc));
  SM->add_loop(cross_term(L13_152364,T_125634,mone*nc));
  SM->add_loop(cross_term(L13_152364,T_126534,n2m1*nc));
  SM->add_loop(cross_term(L13_152364,T_123456,n2m1*nc));
  SM->add_loop(cross_term(L13_152364,T_123465,mone*nc));
//  SM->add_loop(cross_term(L13_152364,T_152364,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L13_152364,T_162354,n2*nc));
  SM->add_loop(cross_term(L13_152364,T_125346,n2m1*nc));
  SM->add_loop(cross_term(L13_152364,T_126345,n2m1*nc));

//SM->add_loop(cross_term(L13_162354,T_156234,0*nc));
//SM->add_loop(cross_term(L13_162354,T_165234,0*nc));
//SM->add_loop(cross_term(L13_162354,T_123564,0*nc));
//SM->add_loop(cross_term(L13_162354,T_123654,0*nc));
  SM->add_loop(cross_term(L13_162354,T_125634,n2m1*nc));
  SM->add_loop(cross_term(L13_162354,T_126534,mone*nc));
  SM->add_loop(cross_term(L13_162354,T_123456,mone*nc));
  SM->add_loop(cross_term(L13_162354,T_123465,n2m1*nc));
  SM->add_loop(cross_term(L13_162354,T_152364,n2*nc));
//  SM->add_loop(cross_term(L13_162354,T_162354,n2_times_n2m1*nc));
  SM->add_loop(cross_term(L13_162354,T_125346,n2m1*nc));
  SM->add_loop(cross_term(L13_162354,T_126345,n2m1*nc));

//---

  SM->add_loop(cross_term(L24_125346,T_156234,mone*nc));
  SM->add_loop(cross_term(L24_125346,T_165234,n2m1*nc));
  SM->add_loop(cross_term(L24_125346,T_123564,n2m1*nc));
  SM->add_loop(cross_term(L24_125346,T_123654,mone*nc));
//SM->add_loop(cross_term(L24_125346,T_125634,0*nc));
//SM->add_loop(cross_term(L24_125346,T_126534,0*nc));
//SM->add_loop(cross_term(L24_125346,T_123456,0*nc));
//SM->add_loop(cross_term(L24_125346,T_123465,0*nc));
  SM->add_loop(cross_term(L24_125346,T_152364,n2m1*nc));
  SM->add_loop(cross_term(L24_125346,T_162354,n2m1*nc));
  SM->add_loop(cross_term(L24_125346,T_125346,n2m1*nc));
  SM->add_loop(cross_term(L24_125346,T_126345,one*nc));

  SM->add_loop(cross_term(L24_126345,T_156234,n2m1*nc));
  SM->add_loop(cross_term(L24_126345,T_165234,mone*nc));
  SM->add_loop(cross_term(L24_126345,T_123564,mone*nc));
  SM->add_loop(cross_term(L24_126345,T_123654,n2m1*nc));
//SM->add_loop(cross_term(L24_126345,T_125634,0*nc));
//SM->add_loop(cross_term(L24_126345,T_126534,0*nc));
//SM->add_loop(cross_term(L24_126345,T_123456,0*nc));
//SM->add_loop(cross_term(L24_126345,T_123465,0*nc));
  SM->add_loop(cross_term(L24_126345,T_152364,n2m1*nc));
  SM->add_loop(cross_term(L24_126345,T_162354,n2m1*nc));
  SM->add_loop(cross_term(L24_126345,T_125346,one*nc));
  SM->add_loop(cross_term(L24_126345,T_126345,n2m1*nc));

	};
}


	};

#if 0
//CONTINUE HERE
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

		SM->add_tree(cross_term(T_15234,T_15234,uu));
		SM->add_tree(cross_term(T_12354,T_12354,uu));
		}	

if(tree_color!=1){

		if((h1.helicity()+h4.helicity()==0)){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cross_term(T_15234,T_15234,uu));
			SM->add_tree(cross_term(T_15234,T_12534,uu_EE));
			//0: SM->add_tree(cross_term(T_15234,T_12354,0));
			SM->add_tree(cross_term(T_15234,T_12345,uu_EE));
		//second row
			SM->add_tree(cross_term(T_12534,T_15234,uu_EE));
			SM->add_tree(cross_term(T_12534,T_12534,uu_EE));
			SM->add_tree(cross_term(T_12534,T_12354,uu_EE));
			//0 SM->add_tree(cross_term(T_12534,T_12345,uu_E));
		//third row
			//0: SM->add_tree(cross_term(T_12354,T_15234,0));
			SM->add_tree(cross_term(T_12354,T_12534,uu_EE));
			//leading-color: SM->add_tree(cross_term(T_12354,T_12354,uu));
			SM->add_tree(cross_term(T_12354,T_12345,uu_EE));
		//fourth row
		 	SM->add_tree(cross_term(T_12345,T_15234,uu_EE));
			//0 SM->add_tree(cross_term(T_12345,T_12534,0));
			SM->add_tree(cross_term(T_12345,T_12354,uu_EE));
			SM->add_tree(cross_term(T_12345,T_12345,uu_EE));
		};
}


	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cross_term(L1_color_15234 ,T_15234 ,uu_virtual));
		SM->add_loop(cross_term(L3_color_12354 ,T_12354 ,uu_virtual));
		
		};

if(color!=1){

		if((h1.helicity()+h4.helicity()==0)){
		
		
		//block 1
		//first row
			//leading-color SM->add_loop(cross_term(L1_15234,T_15234,uu));
		  	SM->add_loop(cross_term(L1_15234,T_12534,uu_EE_virtual));
        		//0: SM->add_loop(cross_term(L1_15234,T_12354,0));
			SM->add_loop(cross_term(L1_15234,T_12345,uu_EE_virtual));
		//second row
			SM->add_loop(cross_term(L2_12534,T_15234,uu_EE_virtual));
			SM->add_loop(cross_term(L2_12534,T_12534,uu_EE_virtual));
			SM->add_loop(cross_term(L2_12534,T_12354,uu_EE_virtual));
			//0 SM->add_loop(cross_term(L2_12534,T_12345,uu_E));
		//third row
			//0: SM->add_loop(cross_term(L3_12354,T_15234,0));
			SM->add_loop(cross_term(L3_12354,T_12534,uu_EE_virtual));
			//leading-color: SM->add_loop(cross_term(L3_12354,T_12354,uu));
			SM->add_loop(cross_term(L3_12354,T_12345,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cross_term(L4_12345,T_15234,uu_EE_virtual));
			//0 SM->add_loop(cross_term(L4_12345,T_12534,0));
			SM->add_loop(cross_term(L4_12345,T_12354,uu_EE_virtual));
			SM->add_loop(cross_term(L4_12345,T_12345,uu_EE_virtual));
		
		};
}

	if(h1.helicity()+h2.helicity()==0){

		SM->add_tree(cross_term(T_15432,T_15432,uu));
		SM->add_tree(cross_term(T_14352,T_14352,uu));
		
		}	

if(tree_color!=1){

		if(h1.helicity()+h2.helicity()==0){
	
		//block 1
		//first row
			//leading-color SM->add_tree(cross_term(T_15432,T_15432,uu));
			SM->add_tree(cross_term(T_15432,T_14532,uu_EE));
			//0: SM->add_tree(cross_term(T_15432,T_14352,0));
			SM->add_tree(cross_term(T_15432,T_14325,uu_EE));
		//second row
			SM->add_tree(cross_term(T_14532,T_15432,uu_EE));
			SM->add_tree(cross_term(T_14532,T_14532,uu_EE));
			SM->add_tree(cross_term(T_14532,T_14352,uu_EE));
			//0 SM->add_tree(cross_term(T_14532,T_14325,uu_E));
		//third row
			//0: SM->add_tree(cross_term(T_14352,T_15432,0));
			SM->add_tree(cross_term(T_14352,T_14532,uu_EE));
			//leading-color: SM->add_tree(cross_term(T_14352,T_14352,uu));
			SM->add_tree(cross_term(T_14352,T_14325,uu_EE));
		//fourth row
		 	SM->add_tree(cross_term(T_14325,T_15432,uu_EE));
			//0 SM->add_tree(cross_term(T_14325,T_14532,0));
			SM->add_tree(cross_term(T_14325,T_14352,uu_EE));
			SM->add_tree(cross_term(T_14325,T_14325,uu_EE));
	
		};
}


	if(h1.helicity()+h2.helicity()==0){
		SM->add_loop(cross_term(L1_color_15432,T_15432,uu_virtual));
		SM->add_loop(cross_term(L3_color_14352,T_14352,uu_virtual));
		};

if(color!=1){

		
		if(h1.helicity()+h2.helicity()==0){
		
		//block 1
		//first row
			//leading-color SM->add_loop(cross_term(L1_15432,T_15432,uu));
		  	SM->add_loop(cross_term(L1_15432,T_14532,uu_EE_virtual));
        		//0: SM->add_loop(cross_term(L1_15432,T_14352,0));
			SM->add_loop(cross_term(L1_15432,T_14325,uu_EE_virtual));
		//second row
			SM->add_loop(cross_term(L2_14532,T_15432,uu_EE_virtual));
			SM->add_loop(cross_term(L2_14532,T_14532,uu_EE_virtual));
			SM->add_loop(cross_term(L2_14532,T_14352,uu_EE_virtual));
			//0 SM->add_loop(cross_term(L2_14532,T_14325,uu_E));
		//third row
			//0: SM->add_loop(cross_term(L3_14352,T_15432,0));
			SM->add_loop(cross_term(L3_14352,T_14532,uu_EE_virtual));
			//leading-color: SM->add_loop(cross_term(L3_14352,T_14352,uu));
			SM->add_loop(cross_term(L3_14352,T_14325,uu_EE_virtual));
		//fourth row
		 	SM->add_loop(cross_term(L4_14325,T_15432,uu_EE_virtual));
			//0 SM->add_loop(cross_term(L4_14325,T_14532,0));
			SM->add_loop(cross_term(L4_14325,T_14352,uu_EE_virtual));
			SM->add_loop(cross_term(L4_14325,T_14325,uu_EE_virtual));
		};
}

////////////////////////////////
///identical anti-quark exchange term 
///block B
////////////////////////////////


if(tree_color!=1){
		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cross_term(T_15432,T_15234,uu_E));
// 3			SM->add_tree(cross_term(T_15432,T_12534,uu_E));
			SM->add_tree(cross_term(T_15432,T_12354,uu_E));
			SM->add_tree(cross_term(T_15432,T_12345,uu_E));
		//second row
// 2			SM->add_tree(cross_term(T_14532,T_15234,uu_E));
			SM->add_tree(cross_term(T_14532,T_12534,uu_EEE));
			SM->add_tree(cross_term(T_14532,T_12354,uu_E));
			SM->add_tree(cross_term(T_14532,T_12345,uu_EEE));
		//third row
			SM->add_tree(cross_term(T_14352,T_15234,uu_E));
			SM->add_tree(cross_term(T_14352,T_12534,uu_E));
			SM->add_tree(cross_term(T_14352,T_12354,uu_E));
// 1			SM->add_tree(cross_term(T_14352,T_12345,uu_E));
		//fourth row
		 	SM->add_tree(cross_term(T_14325,T_15234,uu_E));
			SM->add_tree(cross_term(T_14325,T_12534,uu_EEE));
// 4			SM->add_tree(cross_term(T_14325,T_12354,uu_E));
			SM->add_tree(cross_term(T_14325,T_12345,uu_EEE));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_tree(cross_term(T_15234,T_15432,uu_E));
//			SM->add_tree(cross_term(T_15234,T_14532,uu_E));
			SM->add_tree(cross_term(T_15234,T_14352,uu_E));
			SM->add_tree(cross_term(T_15234,T_14325,uu_E));
		//second row
//			SM->add_tree(cross_term(T_12534,T_15432,uu_E));
			SM->add_tree(cross_term(T_12534,T_14532,uu_EEE));
			SM->add_tree(cross_term(T_12534,T_14352,uu_E));
			SM->add_tree(cross_term(T_12534,T_14325,uu_EEE));
		//third row
			SM->add_tree(cross_term(T_12354,T_15432,uu_E));
			SM->add_tree(cross_term(T_12354,T_14532,uu_E));
			SM->add_tree(cross_term(T_12354,T_14352,uu_E));
//			SM->add_tree(cross_term(T_12354,T_14325,uu_E));
		//fourth row
		 	SM->add_tree(cross_term(T_12345,T_15432,uu_E));
			SM->add_tree(cross_term(T_12345,T_14532,uu_EEE));
//			SM->add_tree(cross_term(T_12345,T_14352,uu_E));
			SM->add_tree(cross_term(T_12345,T_14325,uu_EEE));
		};
}

if(color!=1){

		if((h1.helicity()+h4.helicity()==0) && (h3.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cross_term(L1_15234,T_15432,uu_E_virtual));
			//SM->add_loop(cross_term(L1_15234,T_14532,uu_E_virtual));
			SM->add_loop(cross_term(L1_15234,T_14352,uu_E_virtual));
			SM->add_loop(cross_term(L1_15234,T_14325,uu_E_virtual));
		//second row
			//SM->add_loop(cross_term(L2_12534,T_15432,uu_E_virtual));
			SM->add_loop(cross_term(L2_12534,T_14532,uu_EEE_virtual));
			SM->add_loop(cross_term(L2_12534,T_14352,uu_E_virtual));
			SM->add_loop(cross_term(L2_12534,T_14325,uu_EEE_virtual));
		//third row
			SM->add_loop(cross_term(L3_12354,T_15432,uu_E_virtual));
			SM->add_loop(cross_term(L3_12354,T_14532,uu_E_virtual));
			SM->add_loop(cross_term(L3_12354,T_14352,uu_E_virtual));
			//SM->add_loop(cross_term(L3_12354,T_14325,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cross_term(L4_12345,T_15432,uu_E_virtual));
			SM->add_loop(cross_term(L4_12345,T_14532,uu_EEE_virtual));
			//SM->add_loop(cross_term(L4_12345,T_14352,uu_E_virtual));
			SM->add_loop(cross_term(L4_12345,T_14325,uu_EEE_virtual));
		};
		if((h3.helicity()+h4.helicity()==0) && (h1.helicity()+h4.helicity()==0)){
		//first row
			SM->add_loop(cross_term(L1_15432,T_15234,uu_E_virtual));
		//	SM->add_loop(cross_term(L1_15432,T_12534,uu_E_virtual));
			SM->add_loop(cross_term(L1_15432,T_12354,uu_E_virtual));
			SM->add_loop(cross_term(L1_15432,T_12345,uu_E_virtual));
		//second row
		//	SM->add_loop(cross_term(L2_14532,T_15234,uu_E_virtual));
			SM->add_loop(cross_term(L2_14532,T_12534,uu_EEE_virtual));
			SM->add_loop(cross_term(L2_14532,T_12354,uu_E_virtual));
			SM->add_loop(cross_term(L2_14532,T_12345,uu_EEE_virtual));
		//third row
			SM->add_loop(cross_term(L3_14352,T_15234,uu_E_virtual));
			SM->add_loop(cross_term(L3_14352,T_12534,uu_E_virtual));
			SM->add_loop(cross_term(L3_14352,T_12354,uu_E_virtual));
		//	SM->add_loop(cross_term(L3_14352,T_12345,uu_E_virtual));
		//fourth row
		 	SM->add_loop(cross_term(L4_14325,T_15234,uu_E_virtual));
			SM->add_loop(cross_term(L4_14325,T_12534,uu_EEE_virtual));
		//	SM->add_loop(cross_term(L4_14325,T_12354,uu_E_virtual));
			SM->add_loop(cross_term(L4_14325,T_12345,uu_EEE_virtual));
		};

}

}

#endif


}

return SM;
}

/*
Virtual_SME* vsme_2q2Q2g(int ns,int nf,int nc, int case4q,int color, int tree_color){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);
	ind.push_back(6);


clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color));
;
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color));


	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color));

after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}
*/

Virtual_SME* vsme_2q2Q2g(std::vector<int> indext,int ns,int nf,int nc,int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);

clock_t before, after;
before=clock();
//cout<< "start construction:"<< endl;

	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,m,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,m,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,p,m),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbp,qm,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qp,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qp,qbm,qm,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbm,qp,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qm,qbp,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2g_M2(process(qm,qbp,qp,qbm,p,p),ind,ns,nf,nc,case4q,color,tree_color,lo_or_nlo));

after=clock();
//cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;

return VSM;
}




BH_Ampl_2q2Q2g::BH_Ampl_2q2Q2g(int case4q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q2g(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				case4q,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			6,    // NbrExtParticles
			4,	  // NbrPowersOfAlphaS
			0,    // NbrPowersOfAlphaQED
			4,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
{}






}
