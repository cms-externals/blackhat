/*
 * matrix elements for 2q2Q2g2l
 *
 *  Created on: Aug 20, 2009
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
// #define _PHOTON_ONLY 1     // 0 photon only; 1 for all other replaced by a setting

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {


//expected indices:
// q g g Gb G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_12(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;
}


//expected indices:
// q g Gb G g qb l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_12_34(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;
}


//expected indices:
// q Gb G g g qb l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_34(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;
}


//expected indices:
// q Gb g g G qb l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_23(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	//subleading contributions because of overall 1/Nc for this partial amplitude
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;
}



//expected indices:
// q Gb g G qb g l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_23_41(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	//subleading contributions because of overall 1/Nc for this partial amplitude
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;
}


//expected indices:
// q Gb G qb g g l lb
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_41(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	//subleading contributions because of overall 1/Nc for this partial amplitude
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;

}


//expected indices:
// q Gb G qb ; g g ; l lb
// gluons are contracted in color trace without fundamentals
// fundamentals appear in Kronecker deltas: {1,2} & {3,4}
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_1(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	//subleading contributions because of overall 1/Nc for this partial amplitude
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;

}


//expected indices:
// q Gb G qb ; g g ; l lb
// gluons are contracted in color trace without fundamentals
// fundamentals appear in Kronecker deltas: {3,2} & {1,4}
partial_amplitude_cached* A_loop_2q_2Q_2g_2l_8_2(process pro,vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);
	int i7=ind.at(6);
	int i8=ind.at(7);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
	ph_type h7=pro.p(7);
	ph_type h8=pro.p(8);

	process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);
	ind_123456.push_back(i7);
	ind_123456.push_back(i8);

vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm4_over_2*(r1+r2+r3),-1);
// primitive amplitue
	//subleading contributions because of overall 1/Nc for this partial amplitude
	PA->add(pro_123456,leading_color,ind_123456,1,1);
	PA->add(pro_123456,nf,ind_123456,n_f,n_c);
}

	return PA;

}




Squared_ME* A_loop_2q_2Q_2g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int photonZW, const ph_type e_ph_type,int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);
int i7=ind.at(6);
int i8=ind.at(7);

// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors

//quarks labels i1-i4
//extra gluon called i5-i6
//leptons i7-i8

int flavor(1);

ph_type h1=pro.p(1);
ph_type h2=particle_ID(gluino,pro.p(2).helicity(),flavor,true);
ph_type h3=particle_ID(gluino,pro.p(3).helicity(),flavor);
ph_type h4=pro.p(4);
ph_type h5=pro.p(5);
ph_type h6=pro.p(6);
ph_type h7=pro.p(7);
ph_type h8=pro.p(8);

ph_type hp1=pro.p(1);
ph_type hp2=pro.p(2);
ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),flavor);
ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),flavor,true);
ph_type hp5=pro.p(5);
ph_type hp6=pro.p(6);
ph_type hp7=pro.p(7);
ph_type hp8=pro.p(8);

ph_type hf1=particle_ID(gluino,pro.p(1).helicity(),flavor);
ph_type hf2=pro.p(2);
ph_type hf3=pro.p(3);
ph_type hf4=particle_ID(gluino,pro.p(4).helicity(),flavor,true);

ph_type hpf1=particle_ID(gluino,pro.p(1).helicity(),flavor);
ph_type hpf2=particle_ID(gluino,pro.p(2).helicity(),flavor,true);
ph_type hpf3=pro.p(3);
ph_type hpf4=pro.p(4);

//-------------------------------------------------------
// tree labels

process pro_156234=process(h1,h5,h6,h2,h3,h4,h7,h8);
vector<int> ind_156234;
ind_156234.push_back(i1);
ind_156234.push_back(i5);
ind_156234.push_back(i6);
ind_156234.push_back(i2);
ind_156234.push_back(i3);
ind_156234.push_back(i4);
ind_156234.push_back(i7);
ind_156234.push_back(i8);

process pro_165234=process(h1,h6,h5,h2,h3,h4,h7,h8);
vector<int> ind_165234;
ind_165234.push_back(i1);
ind_165234.push_back(i6);
ind_165234.push_back(i5);
ind_165234.push_back(i2);
ind_165234.push_back(i3);
ind_165234.push_back(i4);
ind_165234.push_back(i7);
ind_165234.push_back(i8);

process pro_152364=process(h1,h5,h2,h3,h6,h4,h7,h8);
vector<int> ind_152364;
ind_152364.push_back(i1);
ind_152364.push_back(i5);
ind_152364.push_back(i2);
ind_152364.push_back(i3);
ind_152364.push_back(i6);
ind_152364.push_back(i4);
ind_152364.push_back(i7);
ind_152364.push_back(i8);

process pro_162354=process(h1,h6,h2,h3,h5,h4,h7,h8);
vector<int> ind_162354;
ind_162354.push_back(i1);
ind_162354.push_back(i6);
ind_162354.push_back(i2);
ind_162354.push_back(i3);
ind_162354.push_back(i5);
ind_162354.push_back(i4);
ind_162354.push_back(i7);
ind_162354.push_back(i8);

process pro_125634=process(h1,h2,h5,h6,h3,h4,h7,h8);
vector<int> ind_125634;
ind_125634.push_back(i1);
ind_125634.push_back(i2);
ind_125634.push_back(i5);
ind_125634.push_back(i6);
ind_125634.push_back(i3);
ind_125634.push_back(i4);
ind_125634.push_back(i7);
ind_125634.push_back(i8);

process pro_126534=process(h1,h2,h6,h5,h3,h4,h7,h8);
vector<int> ind_126534;
ind_126534.push_back(i1);
ind_126534.push_back(i2);
ind_126534.push_back(i6);
ind_126534.push_back(i5);
ind_126534.push_back(i3);
ind_126534.push_back(i4);
ind_126534.push_back(i7);
ind_126534.push_back(i8);

process pro_125346=process(h1,h2,h5,h3,h4,h6,h7,h8);
vector<int> ind_125346;
ind_125346.push_back(i1);
ind_125346.push_back(i2);
ind_125346.push_back(i5);
ind_125346.push_back(i3);
ind_125346.push_back(i4);
ind_125346.push_back(i6);
ind_125346.push_back(i7);
ind_125346.push_back(i8);

process pro_126345=process(h1,h2,h6,h3,h4,h5,h7,h8);
vector<int> ind_126345;
ind_126345.push_back(i1);
ind_126345.push_back(i2);
ind_126345.push_back(i6);
ind_126345.push_back(i3);
ind_126345.push_back(i4);
ind_126345.push_back(i5);
ind_126345.push_back(i7);
ind_126345.push_back(i8);

process pro_123564=process(h1,h2,h3,h5,h6,h4,h7,h8);
vector<int> ind_123564;
ind_123564.push_back(i1);
ind_123564.push_back(i2);
ind_123564.push_back(i3);
ind_123564.push_back(i5);
ind_123564.push_back(i6);
ind_123564.push_back(i4);
ind_123564.push_back(i7);
ind_123564.push_back(i8);

process pro_123654=process(h1,h2,h3,h6,h5,h4,h7,h8);
vector<int> ind_123654;
ind_123654.push_back(i1);
ind_123654.push_back(i2);
ind_123654.push_back(i3);
ind_123654.push_back(i6);
ind_123654.push_back(i5);
ind_123654.push_back(i4);
ind_123654.push_back(i7);
ind_123654.push_back(i8);

process pro_123456=process(h1,h2,h3,h4,h5,h6,h7,h8);
vector<int> ind_123456;
ind_123456.push_back(i1);
ind_123456.push_back(i2);
ind_123456.push_back(i3);
ind_123456.push_back(i4);
ind_123456.push_back(i5);
ind_123456.push_back(i6);
ind_123456.push_back(i7);
ind_123456.push_back(i8);

process pro_123465=process(h1,h2,h3,h4,h6,h5,h7,h8);
vector<int> ind_123465;
ind_123465.push_back(i1);
ind_123465.push_back(i2);
ind_123465.push_back(i3);
ind_123465.push_back(i4);
ind_123465.push_back(i6);
ind_123465.push_back(i5);
ind_123465.push_back(i7);
ind_123465.push_back(i8);

// 1<->3 & 2<->4   

process pro_356412=process(hf3,h5,h6,hf4,hf1,hf2,h7,h8);
vector<int> ind_356412;
ind_356412.push_back(i3);
ind_356412.push_back(i5);
ind_356412.push_back(i6);
ind_356412.push_back(i4);
ind_356412.push_back(i1);
ind_356412.push_back(i2);
ind_356412.push_back(i7);
ind_356412.push_back(i8);

process pro_365412=process(hf3,h6,h5,hf4,hf1,hf2,h7,h8);
vector<int> ind_365412;
ind_365412.push_back(i3);
ind_365412.push_back(i6);
ind_365412.push_back(i5);
ind_365412.push_back(i4);
ind_365412.push_back(i1);
ind_365412.push_back(i2);
ind_365412.push_back(i7);
ind_365412.push_back(i8);

process pro_354162=process(hf3,h5,hf4,hf1,h6,hf2,h7,h8);
vector<int> ind_354162;
ind_354162.push_back(i3);
ind_354162.push_back(i5);
ind_354162.push_back(i4);
ind_354162.push_back(i1);
ind_354162.push_back(i6);
ind_354162.push_back(i2);
ind_354162.push_back(i7);
ind_354162.push_back(i8);

process pro_364152=process(hf3,h6,hf4,hf1,h5,hf2,h7,h8);
vector<int> ind_364152;
ind_364152.push_back(i3);
ind_364152.push_back(i6);
ind_364152.push_back(i4);
ind_364152.push_back(i1);
ind_364152.push_back(i5);
ind_364152.push_back(i2);
ind_364152.push_back(i7);
ind_364152.push_back(i8);

process pro_345612=process(hf3,hf4,h5,h6,hf1,hf2,h7,h8);
vector<int> ind_345612;
ind_345612.push_back(i3);
ind_345612.push_back(i4);
ind_345612.push_back(i5);
ind_345612.push_back(i6);
ind_345612.push_back(i1);
ind_345612.push_back(i2);
ind_345612.push_back(i7);
ind_345612.push_back(i8);

process pro_346512=process(hf3,hf4,h6,h5,hf1,hf2,h7,h8);
vector<int> ind_346512;
ind_346512.push_back(i3);
ind_346512.push_back(i4);
ind_346512.push_back(i6);
ind_346512.push_back(i5);
ind_346512.push_back(i1);
ind_346512.push_back(i2);
ind_346512.push_back(i7);
ind_346512.push_back(i8);

process pro_345126=process(hf3,hf4,h5,hf1,hf2,h6,h7,h8);
vector<int> ind_345126;
ind_345126.push_back(i3);
ind_345126.push_back(i4);
ind_345126.push_back(i5);
ind_345126.push_back(i1);
ind_345126.push_back(i2);
ind_345126.push_back(i6);
ind_345126.push_back(i7);
ind_345126.push_back(i8);

process pro_346125=process(hf3,hf4,h6,hf1,hf2,h5,h7,h8);
vector<int> ind_346125;
ind_346125.push_back(i3);
ind_346125.push_back(i4);
ind_346125.push_back(i6);
ind_346125.push_back(i1);
ind_346125.push_back(i2);
ind_346125.push_back(i5);
ind_346125.push_back(i7);
ind_346125.push_back(i8);

process pro_341562=process(hf3,hf4,hf1,h5,h6,hf2,h7,h8);
vector<int> ind_341562;
ind_341562.push_back(i3);
ind_341562.push_back(i4);
ind_341562.push_back(i1);
ind_341562.push_back(i5);
ind_341562.push_back(i6);
ind_341562.push_back(i2);
ind_341562.push_back(i7);
ind_341562.push_back(i8);

process pro_341652=process(hf3,hf4,hf1,h6,h5,hf2,h7,h8);
vector<int> ind_341652;
ind_341652.push_back(i3);
ind_341652.push_back(i4);
ind_341652.push_back(i1);
ind_341652.push_back(i6);
ind_341652.push_back(i5);
ind_341652.push_back(i2);
ind_341652.push_back(i7);
ind_341652.push_back(i8);

process pro_341256=process(hf3,hf4,hf1,hf2,h5,h6,h7,h8);
vector<int> ind_341256;
ind_341256.push_back(i3);
ind_341256.push_back(i4);
ind_341256.push_back(i1);
ind_341256.push_back(i2);
ind_341256.push_back(i5);
ind_341256.push_back(i6);
ind_341256.push_back(i7);
ind_341256.push_back(i8);

process pro_341265=process(hf3,hf4,hf1,hf2,h6,h5,h7,h8);
vector<int> ind_341265;
ind_341265.push_back(i3);
ind_341265.push_back(i4);
ind_341265.push_back(i1);
ind_341265.push_back(i2);
ind_341265.push_back(i6);
ind_341265.push_back(i5);
ind_341265.push_back(i7);
ind_341265.push_back(i8);

//---------------- second quark line arrangement for case of identical quarks
// 1<->3

process pro_356214=process(hpf3,h5,h6,hpf2,hpf1,hpf4,h7,h8);
vector<int> ind_356214;
ind_356214.push_back(i3);
ind_356214.push_back(i5);
ind_356214.push_back(i6);
ind_356214.push_back(i2);
ind_356214.push_back(i1);
ind_356214.push_back(i4);
ind_356214.push_back(i7);
ind_356214.push_back(i8);

process pro_365214=process(hpf3,h6,h5,hpf2,hpf1,hpf4,h7,h8);
vector<int> ind_365214;
ind_365214.push_back(i3);
ind_365214.push_back(i6);
ind_365214.push_back(i5);
ind_365214.push_back(i2);
ind_365214.push_back(i1);
ind_365214.push_back(i4);
ind_365214.push_back(i7);
ind_365214.push_back(i8);

process pro_352164=process(hpf3,h5,hpf2,hpf1,h6,hpf4,h7,h8);
vector<int> ind_352164;
ind_352164.push_back(i3);
ind_352164.push_back(i5);
ind_352164.push_back(i2);
ind_352164.push_back(i1);
ind_352164.push_back(i6);
ind_352164.push_back(i4);
ind_352164.push_back(i7);
ind_352164.push_back(i8);

process pro_362154=process(hpf3,h6,hpf2,hpf1,h5,hpf4,h7,h8);
vector<int> ind_362154;
ind_362154.push_back(i3);
ind_362154.push_back(i6);
ind_362154.push_back(i2);
ind_362154.push_back(i1);
ind_362154.push_back(i5);
ind_362154.push_back(i4);
ind_362154.push_back(i7);
ind_362154.push_back(i8);

process pro_325614=process(hpf3,hpf2,h5,h6,hpf1,hpf4,h7,h8);
vector<int> ind_325614;
ind_325614.push_back(i3);
ind_325614.push_back(i2);
ind_325614.push_back(i5);
ind_325614.push_back(i6);
ind_325614.push_back(i1);
ind_325614.push_back(i4);
ind_325614.push_back(i7);
ind_325614.push_back(i8);

process pro_326514=process(hpf3,hpf2,h6,h5,hpf1,hpf4,h7,h8);
vector<int> ind_326514;
ind_326514.push_back(i3);
ind_326514.push_back(i2);
ind_326514.push_back(i6);
ind_326514.push_back(i5);
ind_326514.push_back(i1);
ind_326514.push_back(i4);
ind_326514.push_back(i7);
ind_326514.push_back(i8);

process pro_325146=process(hpf3,hpf2,h5,hpf1,hpf4,h6,h7,h8);
vector<int> ind_325146;
ind_325146.push_back(i3);
ind_325146.push_back(i2);
ind_325146.push_back(i5);
ind_325146.push_back(i1);
ind_325146.push_back(i4);
ind_325146.push_back(i6);
ind_325146.push_back(i7);
ind_325146.push_back(i8);

process pro_326145=process(hpf3,hpf2,h6,hpf1,hpf4,h5,h7,h8);
vector<int> ind_326145;
ind_326145.push_back(i3);
ind_326145.push_back(i2);
ind_326145.push_back(i6);
ind_326145.push_back(i1);
ind_326145.push_back(i4);
ind_326145.push_back(i5);
ind_326145.push_back(i7);
ind_326145.push_back(i8);

process pro_321564=process(hpf3,hpf2,hpf1,h5,h6,hpf4,h7,h8);
vector<int> ind_321564;
ind_321564.push_back(i3);
ind_321564.push_back(i2);
ind_321564.push_back(i1);
ind_321564.push_back(i5);
ind_321564.push_back(i6);
ind_321564.push_back(i4);
ind_321564.push_back(i7);
ind_321564.push_back(i8);

process pro_321654=process(hpf3,hpf2,hpf1,h6,h5,hpf4,h7,h8);
vector<int> ind_321654;
ind_321654.push_back(i3);
ind_321654.push_back(i2);
ind_321654.push_back(i1);
ind_321654.push_back(i6);
ind_321654.push_back(i5);
ind_321654.push_back(i4);
ind_321654.push_back(i7);
ind_321654.push_back(i8);

process pro_321456=process(hpf3,hpf2,hpf1,hpf4,h5,h6,h7,h8);
vector<int> ind_321456;
ind_321456.push_back(i3);
ind_321456.push_back(i2);
ind_321456.push_back(i1);
ind_321456.push_back(i4);
ind_321456.push_back(i5);
ind_321456.push_back(i6);
ind_321456.push_back(i7);
ind_321456.push_back(i8);

process pro_321465=process(hpf3,hpf2,hpf1,hpf4,h6,h5,h7,h8);
vector<int> ind_321465;
ind_321465.push_back(i3);
ind_321465.push_back(i2);
ind_321465.push_back(i1);
ind_321465.push_back(i4);
ind_321465.push_back(i6);
ind_321465.push_back(i5);
ind_321465.push_back(i7);
ind_321465.push_back(i8);

// 2<->4
process pro_156432=process(hp1,h5,h6,hp4,hp3,hp2,h7,h8);
vector<int> ind_156432;
ind_156432.push_back(i1);
ind_156432.push_back(i5);
ind_156432.push_back(i6);
ind_156432.push_back(i4);
ind_156432.push_back(i3);
ind_156432.push_back(i2);
ind_156432.push_back(i7);
ind_156432.push_back(i8);

process pro_165432=process(hp1,h6,h5,hp4,hp3,hp2,h7,h8);
vector<int> ind_165432;
ind_165432.push_back(i1);
ind_165432.push_back(i6);
ind_165432.push_back(i5);
ind_165432.push_back(i4);
ind_165432.push_back(i3);
ind_165432.push_back(i2);
ind_165432.push_back(i7);
ind_165432.push_back(i8);

process pro_154362=process(hp1,h5,hp4,hp3,h6,hp2,h7,h8);
vector<int> ind_154362;
ind_154362.push_back(i1);
ind_154362.push_back(i5);
ind_154362.push_back(i4);
ind_154362.push_back(i3);
ind_154362.push_back(i6);
ind_154362.push_back(i2);
ind_154362.push_back(i7);
ind_154362.push_back(i8);

process pro_164352=process(hp1,h6,hp4,hp3,h5,hp2,h7,h8);
vector<int> ind_164352;
ind_164352.push_back(i1);
ind_164352.push_back(i6);
ind_164352.push_back(i4);
ind_164352.push_back(i3);
ind_164352.push_back(i5);
ind_164352.push_back(i2);
ind_164352.push_back(i7);
ind_164352.push_back(i8);

process pro_145632=process(hp1,hp4,h5,h6,hp3,hp2,h7,h8);
vector<int> ind_145632;
ind_145632.push_back(i1);
ind_145632.push_back(i4);
ind_145632.push_back(i5);
ind_145632.push_back(i6);
ind_145632.push_back(i3);
ind_145632.push_back(i2);
ind_145632.push_back(i7);
ind_145632.push_back(i8);

process pro_146532=process(hp1,hp4,h6,h5,hp3,hp2,h7,h8);
vector<int> ind_146532;
ind_146532.push_back(i1);
ind_146532.push_back(i4);
ind_146532.push_back(i6);
ind_146532.push_back(i5);
ind_146532.push_back(i3);
ind_146532.push_back(i2);
ind_146532.push_back(i7);
ind_146532.push_back(i8);

process pro_145326=process(hp1,hp4,h5,hp3,hp2,h6,h7,h8);
vector<int> ind_145326;
ind_145326.push_back(i1);
ind_145326.push_back(i4);
ind_145326.push_back(i5);
ind_145326.push_back(i3);
ind_145326.push_back(i2);
ind_145326.push_back(i6);
ind_145326.push_back(i7);
ind_145326.push_back(i8);

process pro_146325=process(hp1,hp4,h6,hp3,hp2,h5,h7,h8);
vector<int> ind_146325;
ind_146325.push_back(i1);
ind_146325.push_back(i4);
ind_146325.push_back(i6);
ind_146325.push_back(i3);
ind_146325.push_back(i2);
ind_146325.push_back(i5);
ind_146325.push_back(i7);
ind_146325.push_back(i8);

process pro_143562=process(hp1,hp4,hp3,h5,h6,hp2,h7,h8);
vector<int> ind_143562;
ind_143562.push_back(i1);
ind_143562.push_back(i4);
ind_143562.push_back(i3);
ind_143562.push_back(i5);
ind_143562.push_back(i6);
ind_143562.push_back(i2);
ind_143562.push_back(i7);
ind_143562.push_back(i8);

process pro_143652=process(hp1,hp4,hp3,h6,h5,hp2,h7,h8);
vector<int> ind_143652;
ind_143652.push_back(i1);
ind_143652.push_back(i4);
ind_143652.push_back(i3);
ind_143652.push_back(i6);
ind_143652.push_back(i5);
ind_143652.push_back(i2);
ind_143652.push_back(i7);
ind_143652.push_back(i8);

process pro_143256=process(hp1,hp4,hp3,hp2,h5,h6,h7,h8);
vector<int> ind_143256;
ind_143256.push_back(i1);
ind_143256.push_back(i4);
ind_143256.push_back(i3);
ind_143256.push_back(i2);
ind_143256.push_back(i5);
ind_143256.push_back(i6);
ind_143256.push_back(i7);
ind_143256.push_back(i8);

process pro_143265=process(hp1,hp4,hp3,hp2,h6,h5,h7,h8);
vector<int> ind_143265;
ind_143265.push_back(i1);
ind_143265.push_back(i4);
ind_143265.push_back(i3);
ind_143265.push_back(i2);
ind_143265.push_back(i6);
ind_143265.push_back(i5);
ind_143265.push_back(i7);
ind_143265.push_back(i8);


vector<ph_type> _ph_type_1;
_ph_type_1.push_back(h1);
_ph_type_1.push_back(e_ph_type);

vector<ph_type> _ph_type_3;
_ph_type_3.push_back(h3);
_ph_type_3.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn_up_1(0,photonZW,0,i7,i8,_ph_type_1);
	prop_hel_fn _prop_hel_fn_up_3(0,photonZW,0,i7,i8,_ph_type_3);
	prop_hel_fn _prop_hel_fn_down_1(1,photonZW,0,i7,i8,_ph_type_1);
	prop_hel_fn _prop_hel_fn_down_3(1,photonZW,0,i7,i8,_ph_type_3);
//-------------------------------------------


size_t T_156234_down; size_t T_156234_up; size_t L_156234_down; size_t L_color_156234_down; size_t L_color_156234_up; size_t L_156234_up;
size_t T_165234_down; size_t T_165234_up; size_t L_165234_down; size_t L_color_165234_down; size_t L_color_165234_up; size_t L_165234_up;
size_t T_152364_down; size_t T_152364_up; size_t L_152364_down; size_t L_color_152364_down; size_t L_color_152364_up; size_t L_152364_up;
size_t T_162354_down; size_t T_162354_up; size_t L_162354_down; size_t L_color_162354_down; size_t L_color_162354_up; size_t L_162354_up;
size_t T_125634_down; size_t T_125634_up; size_t L_125634_down; size_t L_color_125634_down; size_t L_color_125634_up; size_t L_125634_up;
size_t T_126534_down; size_t T_126534_up; size_t L_126534_down; size_t L_color_126534_down; size_t L_color_126534_up; size_t L_126534_up;
size_t T_125346_down; size_t T_125346_up; size_t L_125346_down; size_t L_color_125346_down; size_t L_color_125346_up; size_t L_125346_up;
size_t T_126345_down; size_t T_126345_up; size_t L_126345_down; size_t L_color_126345_down; size_t L_color_126345_up; size_t L_126345_up;
size_t T_123564_down; size_t T_123564_up; size_t L_123564_down; size_t L_color_123564_down; size_t L_color_123564_up; size_t L_123564_up;
size_t T_123654_down; size_t T_123654_up; size_t L_123654_down; size_t L_color_123654_down; size_t L_color_123654_up; size_t L_123654_up;
size_t T_123456_down; size_t T_123456_up; size_t L_123456_down; size_t L_color_123456_down; size_t L_color_123456_up; size_t L_123456_up;
size_t T_123465_down; size_t T_123465_up; size_t L_123465_down; size_t L_color_123465_down; size_t L_color_123465_up; size_t L_123465_up;

// 1<->3 & 2<->4   

size_t T_356412_down; size_t T_356412_up; size_t L_356412_down; size_t L_color_356412_down; size_t L_color_356412_up; size_t L_356412_up;
size_t T_365412_down; size_t T_365412_up; size_t L_365412_down; size_t L_color_365412_down; size_t L_color_365412_up; size_t L_365412_up;
size_t T_354162_down; size_t T_354162_up; size_t L_354162_down; size_t L_color_354162_down; size_t L_color_354162_up; size_t L_354162_up;
size_t T_364152_down; size_t T_364152_up; size_t L_364152_down; size_t L_color_364152_down; size_t L_color_364152_up; size_t L_364152_up;
size_t T_345612_down; size_t T_345612_up; size_t L_345612_down; size_t L_color_345612_down; size_t L_color_345612_up; size_t L_345612_up;
size_t T_346512_down; size_t T_346512_up; size_t L_346512_down; size_t L_color_346512_down; size_t L_color_346512_up; size_t L_346512_up;
size_t T_345126_down; size_t T_345126_up; size_t L_345126_down; size_t L_color_345126_down; size_t L_color_345126_up; size_t L_345126_up;
size_t T_346125_down; size_t T_346125_up; size_t L_346125_down; size_t L_color_346125_down; size_t L_color_346125_up; size_t L_346125_up;
size_t T_341562_down; size_t T_341562_up; size_t L_341562_down; size_t L_color_341562_down; size_t L_color_341562_up; size_t L_341562_up;
size_t T_341652_down; size_t T_341652_up; size_t L_341652_down; size_t L_color_341652_down; size_t L_color_341652_up; size_t L_341652_up;
size_t T_341256_down; size_t T_341256_up; size_t L_341256_down; size_t L_color_341256_down; size_t L_color_341256_up; size_t L_341256_up;
size_t T_341265_down; size_t T_341265_up; size_t L_341265_down; size_t L_color_341265_down; size_t L_color_341265_up; size_t L_341265_up;

//---------------- second quark line arrangement for case of identical quarks
// 1<->3

size_t T_356214_down; size_t T_356214_up; size_t L_356214_down; size_t L_color_356214_down; size_t L_color_356214_up; size_t L_356214_up;
size_t T_365214_down; size_t T_365214_up; size_t L_365214_down; size_t L_color_365214_down; size_t L_color_365214_up; size_t L_365214_up;
size_t T_352164_down; size_t T_352164_up; size_t L_352164_down; size_t L_color_352164_down; size_t L_color_352164_up; size_t L_352164_up;
size_t T_362154_down; size_t T_362154_up; size_t L_362154_down; size_t L_color_362154_down; size_t L_color_362154_up; size_t L_362154_up;
size_t T_325614_down; size_t T_325614_up; size_t L_325614_down; size_t L_color_325614_down; size_t L_color_325614_up; size_t L_325614_up;
size_t T_326514_down; size_t T_326514_up; size_t L_326514_down; size_t L_color_326514_down; size_t L_color_326514_up; size_t L_326514_up;
size_t T_325146_down; size_t T_325146_up; size_t L_325146_down; size_t L_color_325146_down; size_t L_color_325146_up; size_t L_325146_up;
size_t T_326145_down; size_t T_326145_up; size_t L_326145_down; size_t L_color_326145_down; size_t L_color_326145_up; size_t L_326145_up;
size_t T_321564_down; size_t T_321564_up; size_t L_321564_down; size_t L_color_321564_down; size_t L_color_321564_up; size_t L_321564_up;
size_t T_321654_down; size_t T_321654_up; size_t L_321654_down; size_t L_color_321654_down; size_t L_color_321654_up; size_t L_321654_up;
size_t T_321456_down; size_t T_321456_up; size_t L_321456_down; size_t L_color_321456_down; size_t L_color_321456_up; size_t L_321456_up;
size_t T_321465_down; size_t T_321465_up; size_t L_321465_down; size_t L_color_321465_down; size_t L_color_321465_up; size_t L_321465_up;

// 2<->4
size_t T_156432_down; size_t T_156432_up; size_t L_156432_down; size_t L_color_156432_down; size_t L_color_156432_up; size_t L_156432_up;
size_t T_165432_down; size_t T_165432_up; size_t L_165432_down; size_t L_color_165432_down; size_t L_color_165432_up; size_t L_165432_up;
size_t T_154362_down; size_t T_154362_up; size_t L_154362_down; size_t L_color_154362_down; size_t L_color_154362_up; size_t L_154362_up;
size_t T_164352_down; size_t T_164352_up; size_t L_164352_down; size_t L_color_164352_down; size_t L_color_164352_up; size_t L_164352_up;
size_t T_145632_down; size_t T_145632_up; size_t L_145632_down; size_t L_color_145632_down; size_t L_color_145632_up; size_t L_145632_up;
size_t T_146532_down; size_t T_146532_up; size_t L_146532_down; size_t L_color_146532_down; size_t L_color_146532_up; size_t L_146532_up;
size_t T_145326_down; size_t T_145326_up; size_t L_145326_down; size_t L_color_145326_down; size_t L_color_145326_up; size_t L_145326_up;
size_t T_146325_down; size_t T_146325_up; size_t L_146325_down; size_t L_color_146325_down; size_t L_color_146325_up; size_t L_146325_up;
size_t T_143562_down; size_t T_143562_up; size_t L_143562_down; size_t L_color_143562_down; size_t L_color_143562_up; size_t L_143562_up;
size_t T_143652_down; size_t T_143652_up; size_t L_143652_down; size_t L_color_143652_down; size_t L_color_143652_up; size_t L_143652_up;
size_t T_143256_down; size_t T_143256_up; size_t L_143256_down; size_t L_color_143256_down; size_t L_color_143256_up; size_t L_143256_up;
size_t T_143265_down; size_t T_143265_up; size_t L_143265_down; size_t L_color_143265_down; size_t L_color_143265_up; size_t L_143265_up;

//-------------------------------------------

	if(h1.helicity()+h4.helicity()==0){
	
		T_156234_up=SM->add(new CTree_with_prefactor(pro_156234,ind_156234,_prop_hel_fn_up_1));
		T_152364_up=SM->add(new CTree_with_prefactor(pro_152364,ind_152364,_prop_hel_fn_up_1));
		T_123564_up=SM->add(new CTree_with_prefactor(pro_123564,ind_123564,_prop_hel_fn_up_1));
		T_165234_up=SM->add(new CTree_with_prefactor(pro_165234,ind_165234,_prop_hel_fn_up_1));
		T_162354_up=SM->add(new CTree_with_prefactor(pro_162354,ind_162354,_prop_hel_fn_up_1));
		T_123654_up=SM->add(new CTree_with_prefactor(pro_123654,ind_123654,_prop_hel_fn_up_1));


		T_356412_up=SM->add(new CTree_with_prefactor(pro_356412,ind_356412,_prop_hel_fn_up_3));
		T_354162_up=SM->add(new CTree_with_prefactor(pro_354162,ind_354162,_prop_hel_fn_up_3));
		T_341562_up=SM->add(new CTree_with_prefactor(pro_341562,ind_341562,_prop_hel_fn_up_3));
		T_365412_up=SM->add(new CTree_with_prefactor(pro_365412,ind_365412,_prop_hel_fn_up_3));
		T_364152_up=SM->add(new CTree_with_prefactor(pro_364152,ind_364152,_prop_hel_fn_up_3));
		T_341652_up=SM->add(new CTree_with_prefactor(pro_341652,ind_341652,_prop_hel_fn_up_3));

if(tree_color!=1){
		T_125634_up=SM->add(new CTree_with_prefactor(pro_125634,ind_125634,_prop_hel_fn_up_1));
		T_125346_up=SM->add(new CTree_with_prefactor(pro_125346,ind_125346,_prop_hel_fn_up_1));
		T_123456_up=SM->add(new CTree_with_prefactor(pro_123456,ind_123456,_prop_hel_fn_up_1));
		T_126534_up=SM->add(new CTree_with_prefactor(pro_126534,ind_126534,_prop_hel_fn_up_1));
		T_126345_up=SM->add(new CTree_with_prefactor(pro_126345,ind_126345,_prop_hel_fn_up_1));
		T_123465_up=SM->add(new CTree_with_prefactor(pro_123465,ind_123465,_prop_hel_fn_up_1));

		T_345612_up=SM->add(new CTree_with_prefactor(pro_345612,ind_345612,_prop_hel_fn_up_3));
		T_345126_up=SM->add(new CTree_with_prefactor(pro_345126,ind_345126,_prop_hel_fn_up_3));
		T_341256_up=SM->add(new CTree_with_prefactor(pro_341256,ind_341256,_prop_hel_fn_up_3));
		T_346512_up=SM->add(new CTree_with_prefactor(pro_346512,ind_346512,_prop_hel_fn_up_3));
		T_346125_up=SM->add(new CTree_with_prefactor(pro_346125,ind_346125,_prop_hel_fn_up_3));
		T_341265_up=SM->add(new CTree_with_prefactor(pro_341265,ind_341265,_prop_hel_fn_up_3));
}


if(color!=1){
		L_156234_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156234,ind_156234,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_152364_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_152364,ind_152364,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_123564_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123564,ind_123564,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_165234_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165234,ind_165234,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_162354_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_162354,ind_162354,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_123654_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123654,ind_123654,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));


		L_356412_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356412,ind_356412,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_354162_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_354162,ind_354162,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_341562_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341562,ind_341562,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_365412_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365412,ind_365412,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_364152_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_364152,ind_364152,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_341652_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341652,ind_341652,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
}

		L_color_156234_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156234,ind_156234,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_152364_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_152364,ind_152364,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_123564_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123564,ind_123564,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_165234_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165234,ind_165234,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_162354_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_162354,ind_162354,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_123654_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123654,ind_123654,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));

		L_color_356412_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356412,ind_356412,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_354162_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_354162,ind_354162,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_341562_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341562,ind_341562,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_365412_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365412,ind_365412,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_364152_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_364152,ind_364152,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_341652_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341652,ind_341652,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));

if(color!=1){
		L_125634_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_125634,ind_125634,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_125346_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_125346,ind_125346,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_123456_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_123456,ind_123456,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_126534_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_126534,ind_126534,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_126345_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_126345,ind_126345,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_123465_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_123465,ind_123465,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
	
		L_345612_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_345612,ind_345612,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_345126_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_345126,ind_345126,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_341256_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_341256,ind_341256,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_346512_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_346512,ind_346512,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_346125_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_346125,ind_346125,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_341265_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_341265,ind_341265,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));

}
	};



/////////////////////////////////////////////
//2<->4
	if(hp1.helicity()+hp2.helicity()==0){
	
		T_156432_up=SM->add(new CTree_with_prefactor(pro_156432,ind_156432,_prop_hel_fn_up_1));
		T_154362_up=SM->add(new CTree_with_prefactor(pro_154362,ind_154362,_prop_hel_fn_up_1));
		T_143562_up=SM->add(new CTree_with_prefactor(pro_143562,ind_143562,_prop_hel_fn_up_1));
		T_165432_up=SM->add(new CTree_with_prefactor(pro_165432,ind_165432,_prop_hel_fn_up_1));
		T_164352_up=SM->add(new CTree_with_prefactor(pro_164352,ind_164352,_prop_hel_fn_up_1));
		T_143652_up=SM->add(new CTree_with_prefactor(pro_143652,ind_143652,_prop_hel_fn_up_1));


		T_356214_up=SM->add(new CTree_with_prefactor(pro_356214,ind_356214,_prop_hel_fn_up_3));
		T_352164_up=SM->add(new CTree_with_prefactor(pro_352164,ind_352164,_prop_hel_fn_up_3));
		T_321564_up=SM->add(new CTree_with_prefactor(pro_321564,ind_321564,_prop_hel_fn_up_3));
		T_365214_up=SM->add(new CTree_with_prefactor(pro_365214,ind_365214,_prop_hel_fn_up_3));
		T_362154_up=SM->add(new CTree_with_prefactor(pro_362154,ind_362154,_prop_hel_fn_up_3));
		T_321654_up=SM->add(new CTree_with_prefactor(pro_321654,ind_321654,_prop_hel_fn_up_3));

if(tree_color!=1){
		T_145632_up=SM->add(new CTree_with_prefactor(pro_145632,ind_145632,_prop_hel_fn_up_1));
		T_145326_up=SM->add(new CTree_with_prefactor(pro_145326,ind_145326,_prop_hel_fn_up_1));
		T_143256_up=SM->add(new CTree_with_prefactor(pro_143256,ind_143256,_prop_hel_fn_up_1));
		T_146532_up=SM->add(new CTree_with_prefactor(pro_146532,ind_146532,_prop_hel_fn_up_1));
		T_146325_up=SM->add(new CTree_with_prefactor(pro_146325,ind_146325,_prop_hel_fn_up_1));
		T_143265_up=SM->add(new CTree_with_prefactor(pro_143265,ind_143265,_prop_hel_fn_up_1));

		T_325614_up=SM->add(new CTree_with_prefactor(pro_325614,ind_325614,_prop_hel_fn_up_3));
		T_325146_up=SM->add(new CTree_with_prefactor(pro_325146,ind_325146,_prop_hel_fn_up_3));
		T_321456_up=SM->add(new CTree_with_prefactor(pro_321456,ind_321456,_prop_hel_fn_up_3));
		T_326514_up=SM->add(new CTree_with_prefactor(pro_326514,ind_326514,_prop_hel_fn_up_3));
		T_326145_up=SM->add(new CTree_with_prefactor(pro_326145,ind_326145,_prop_hel_fn_up_3));
		T_321465_up=SM->add(new CTree_with_prefactor(pro_321465,ind_321465,_prop_hel_fn_up_3));
}


if(color!=1){
		L_156432_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156432,ind_156432,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_154362_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_154362,ind_154362,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_143562_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143562,ind_143562,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_165432_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165432,ind_165432,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_164352_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_164352,ind_164352,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_143652_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143652,ind_143652,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));


		L_356214_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356214,ind_356214,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_352164_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_352164,ind_352164,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_321564_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321564,ind_321564,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_365214_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365214,ind_365214,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_362154_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_362154,ind_362154,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_321654_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321654,ind_321654,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
}
		L_color_156432_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156432,ind_156432,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_154362_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_154362,ind_154362,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_143562_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143562,ind_143562,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_165432_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165432,ind_165432,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_164352_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_164352,ind_164352,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_143652_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143652,ind_143652,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));

		L_color_356214_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356214,ind_356214,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_352164_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_352164,ind_352164,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_321564_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321564,ind_321564,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_365214_up=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365214,ind_365214,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_362154_up=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_362154,ind_362154,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_321654_up=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321654,ind_321654,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));

if(color!=1){
		L_145632_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_145632,ind_145632,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_145326_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_145326,ind_145326,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_143256_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_143256,ind_143256,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_146532_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_146532,ind_146532,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_146325_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_146325,ind_146325,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_143265_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_143265,ind_143265,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
	
		L_325614_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_325614,ind_325614,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_325146_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_325146,ind_325146,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_321456_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_321456,ind_321456,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_326514_up=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_326514,ind_326514,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_326145_up=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_326145,ind_326145,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));
		L_321465_up=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_321465,ind_321465,n_s,n_f,n_c,0,photonZW,e_ph_type,0,lo_or_nlo));

}
	};

/////////////////////////////////////////////
//up<->down

// W cases constructed using up
if(photonZW!=3){

	if(h1.helicity()+h4.helicity()==0){
	
		T_156234_down=SM->add(new CTree_with_prefactor(pro_156234,ind_156234,_prop_hel_fn_down_1));
		T_152364_down=SM->add(new CTree_with_prefactor(pro_152364,ind_152364,_prop_hel_fn_down_1));
		T_123564_down=SM->add(new CTree_with_prefactor(pro_123564,ind_123564,_prop_hel_fn_down_1));
		T_165234_down=SM->add(new CTree_with_prefactor(pro_165234,ind_165234,_prop_hel_fn_down_1));
		T_162354_down=SM->add(new CTree_with_prefactor(pro_162354,ind_162354,_prop_hel_fn_down_1));
		T_123654_down=SM->add(new CTree_with_prefactor(pro_123654,ind_123654,_prop_hel_fn_down_1));


		T_356412_down=SM->add(new CTree_with_prefactor(pro_356412,ind_356412,_prop_hel_fn_down_3));
		T_354162_down=SM->add(new CTree_with_prefactor(pro_354162,ind_354162,_prop_hel_fn_down_3));
		T_341562_down=SM->add(new CTree_with_prefactor(pro_341562,ind_341562,_prop_hel_fn_down_3));
		T_365412_down=SM->add(new CTree_with_prefactor(pro_365412,ind_365412,_prop_hel_fn_down_3));
		T_364152_down=SM->add(new CTree_with_prefactor(pro_364152,ind_364152,_prop_hel_fn_down_3));
		T_341652_down=SM->add(new CTree_with_prefactor(pro_341652,ind_341652,_prop_hel_fn_down_3));

if(tree_color!=1){
		T_125634_down=SM->add(new CTree_with_prefactor(pro_125634,ind_125634,_prop_hel_fn_down_1));
		T_125346_down=SM->add(new CTree_with_prefactor(pro_125346,ind_125346,_prop_hel_fn_down_1));
		T_123456_down=SM->add(new CTree_with_prefactor(pro_123456,ind_123456,_prop_hel_fn_down_1));
		T_126534_down=SM->add(new CTree_with_prefactor(pro_126534,ind_126534,_prop_hel_fn_down_1));
		T_126345_down=SM->add(new CTree_with_prefactor(pro_126345,ind_126345,_prop_hel_fn_down_1));
		T_123465_down=SM->add(new CTree_with_prefactor(pro_123465,ind_123465,_prop_hel_fn_down_1));

		T_345612_down=SM->add(new CTree_with_prefactor(pro_345612,ind_345612,_prop_hel_fn_down_3));
		T_345126_down=SM->add(new CTree_with_prefactor(pro_345126,ind_345126,_prop_hel_fn_down_3));
		T_341256_down=SM->add(new CTree_with_prefactor(pro_341256,ind_341256,_prop_hel_fn_down_3));
		T_346512_down=SM->add(new CTree_with_prefactor(pro_346512,ind_346512,_prop_hel_fn_down_3));
		T_346125_down=SM->add(new CTree_with_prefactor(pro_346125,ind_346125,_prop_hel_fn_down_3));
		T_341265_down=SM->add(new CTree_with_prefactor(pro_341265,ind_341265,_prop_hel_fn_down_3));
}


if(color!=1){
		L_156234_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156234,ind_156234,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_152364_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_152364,ind_152364,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_123564_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123564,ind_123564,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_165234_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165234,ind_165234,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_162354_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_162354,ind_162354,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_123654_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123654,ind_123654,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));


		L_356412_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356412,ind_356412,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_354162_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_354162,ind_354162,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_341562_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341562,ind_341562,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_365412_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365412,ind_365412,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_364152_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_364152,ind_364152,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_341652_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341652,ind_341652,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
}

		L_color_156234_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156234,ind_156234,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_152364_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_152364,ind_152364,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_123564_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123564,ind_123564,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_165234_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165234,ind_165234,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_162354_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_162354,ind_162354,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_123654_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_123654,ind_123654,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		
		L_color_356412_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356412,ind_356412,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_354162_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_354162,ind_354162,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_341562_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341562,ind_341562,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_365412_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365412,ind_365412,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_364152_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_364152,ind_364152,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_341652_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_341652,ind_341652,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));

		if(color!=1){
		L_125634_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_125634,ind_125634,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_125346_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_125346,ind_125346,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_123456_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_123456,ind_123456,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_126534_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_126534,ind_126534,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_126345_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_126345,ind_126345,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_123465_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_123465,ind_123465,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
	
		L_345612_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_345612,ind_345612,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_345126_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_345126,ind_345126,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_341256_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_341256,ind_341256,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_346512_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_346512,ind_346512,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_346125_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_346125,ind_346125,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_341265_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_341265,ind_341265,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));

}
	};



/////////////////////////////////////////////
//2<->4
	if(hp1.helicity()+hp2.helicity()==0){
	
		T_156432_down=SM->add(new CTree_with_prefactor(pro_156432,ind_156432,_prop_hel_fn_down_1));
		T_154362_down=SM->add(new CTree_with_prefactor(pro_154362,ind_154362,_prop_hel_fn_down_1));
		T_143562_down=SM->add(new CTree_with_prefactor(pro_143562,ind_143562,_prop_hel_fn_down_1));
		T_165432_down=SM->add(new CTree_with_prefactor(pro_165432,ind_165432,_prop_hel_fn_down_1));
		T_164352_down=SM->add(new CTree_with_prefactor(pro_164352,ind_164352,_prop_hel_fn_down_1));
		T_143652_down=SM->add(new CTree_with_prefactor(pro_143652,ind_143652,_prop_hel_fn_down_1));


		T_356214_down=SM->add(new CTree_with_prefactor(pro_356214,ind_356214,_prop_hel_fn_down_3));
		T_352164_down=SM->add(new CTree_with_prefactor(pro_352164,ind_352164,_prop_hel_fn_down_3));
		T_321564_down=SM->add(new CTree_with_prefactor(pro_321564,ind_321564,_prop_hel_fn_down_3));
		T_365214_down=SM->add(new CTree_with_prefactor(pro_365214,ind_365214,_prop_hel_fn_down_3));
		T_362154_down=SM->add(new CTree_with_prefactor(pro_362154,ind_362154,_prop_hel_fn_down_3));
		T_321654_down=SM->add(new CTree_with_prefactor(pro_321654,ind_321654,_prop_hel_fn_down_3));

if(tree_color!=1){
		T_145632_down=SM->add(new CTree_with_prefactor(pro_145632,ind_145632,_prop_hel_fn_down_1));
		T_145326_down=SM->add(new CTree_with_prefactor(pro_145326,ind_145326,_prop_hel_fn_down_1));
		T_143256_down=SM->add(new CTree_with_prefactor(pro_143256,ind_143256,_prop_hel_fn_down_1));
		T_146532_down=SM->add(new CTree_with_prefactor(pro_146532,ind_146532,_prop_hel_fn_down_1));
		T_146325_down=SM->add(new CTree_with_prefactor(pro_146325,ind_146325,_prop_hel_fn_down_1));
		T_143265_down=SM->add(new CTree_with_prefactor(pro_143265,ind_143265,_prop_hel_fn_down_1));

		T_325614_down=SM->add(new CTree_with_prefactor(pro_325614,ind_325614,_prop_hel_fn_down_3));
		T_325146_down=SM->add(new CTree_with_prefactor(pro_325146,ind_325146,_prop_hel_fn_down_3));
		T_321456_down=SM->add(new CTree_with_prefactor(pro_321456,ind_321456,_prop_hel_fn_down_3));
		T_326514_down=SM->add(new CTree_with_prefactor(pro_326514,ind_326514,_prop_hel_fn_down_3));
		T_326145_down=SM->add(new CTree_with_prefactor(pro_326145,ind_326145,_prop_hel_fn_down_3));
		T_321465_down=SM->add(new CTree_with_prefactor(pro_321465,ind_321465,_prop_hel_fn_down_3));
}


if(color!=1){
		L_156432_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156432,ind_156432,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_154362_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_154362,ind_154362,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_143562_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143562,ind_143562,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_165432_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165432,ind_165432,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_164352_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_164352,ind_164352,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_143652_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143652,ind_143652,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));


		L_356214_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356214,ind_356214,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_352164_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_352164,ind_352164,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_321564_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321564,ind_321564,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_365214_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365214,ind_365214,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_362154_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_362154,ind_362154,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_321654_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321654,ind_321654,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
}
		
		L_color_156432_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_156432,ind_156432,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_154362_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_154362,ind_154362,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_143562_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143562,ind_143562,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_165432_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_165432,ind_165432,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_164352_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_164352,ind_164352,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_143652_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_143652,ind_143652,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));

		L_color_356214_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_356214,ind_356214,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_352164_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_352164,ind_352164,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_321564_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321564,ind_321564,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_365214_down=SM->add(A_loop_2q_2Q_2g_2l_8_12(   pro_365214,ind_365214,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_362154_down=SM->add(A_loop_2q_2Q_2g_2l_8_12_34(pro_362154,ind_362154,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L_color_321654_down=SM->add(A_loop_2q_2Q_2g_2l_8_34(   pro_321654,ind_321654,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));

if(color!=1){
		L_145632_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_145632,ind_145632,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_145326_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_145326,ind_145326,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_143256_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_143256,ind_143256,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_146532_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_146532,ind_146532,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_146325_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_146325,ind_146325,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_143265_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_143265,ind_143265,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
	
		L_325614_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_325614,ind_325614,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_325146_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_325146,ind_325146,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_321456_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_321456,ind_321456,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_326514_down=SM->add(A_loop_2q_2Q_2g_2l_8_23(   pro_326514,ind_326514,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_326145_down=SM->add(A_loop_2q_2Q_2g_2l_8_23_41(pro_326145,ind_326145,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));
		L_321465_down=SM->add(A_loop_2q_2Q_2g_2l_8_41(   pro_321465,ind_321465,n_s,n_f,n_c,1,photonZW,e_ph_type,0,lo_or_nlo));

}
	};

}

//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction two(2,1);
multi_precision_fraction _8(8,1);
multi_precision_fraction _4(4,1);
multi_precision_fraction _2(2,1);
multi_precision_fraction mone(-1,1);
//multi_precision_fraction nf(n_f,1);
multi_precision_fraction n2m1(n_c*n_c-1,1);
multi_precision_fraction n2m1_2=n2m1*n2m1;
multi_precision_fraction _over_nc(1,n_c);
multi_precision_fraction _over_nc2(1,n_c*n_c);
multi_precision_fraction _over_nc3(1,n_c*n_c*n_c);

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
		multi_precision_fraction c(1,1);
#if 0
_PRINT("-------------------------");
_PRINT("identical quarks");
_PRINT("-------------------------");
#endif


		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,n2m1_2*_4));
		};
		if(hp3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
SM->add_tree(cached_cross_term_md(T_362154_up,T_362154_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_352164_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_356214_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_321564_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_365214_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_321654_up,n2m1_2*_4));
		};

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
SM->add_loop(cached_cross_term_md(L_color_162354_up,T_162354_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_152364_up,T_152364_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_156234_up,T_156234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123564_up,T_123564_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,n2m1_2*nc*_8));
		};
		if(hp3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
SM->add_loop(cached_cross_term_md(L_color_362154_up,T_362154_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_352164_up,T_352164_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_356214_up,T_356214_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_321564_up,T_321564_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_365214_up,T_365214_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_321654_up,T_321654_up,n2m1_2*nc*_8));
		};



if(tree_color!=1){

		if((h1.helicity()==-1&&h1.helicity()+h4.helicity()==0)){
//SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_152364_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_125346_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_126345_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_162354_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123564_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_165234_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123654_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_162354_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123465_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123456_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_152364_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_126534_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_125634_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_156234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_165234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123654_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_156234_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123564_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123654_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_152364_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123465_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_162354_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_126534_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_125634_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_156234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123564_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_165234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,_4*n2m1_2));

	};

		if((h3.helicity()==-1&&h3.helicity()+h4.helicity()==0)){
SM->add_tree(cached_cross_term_md(T_326145_up,T_326145_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_352164_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_362154_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_325146_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_321564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_356214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_321654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_365214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_326145_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_352164_up,T_352164_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_362154_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_325146_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_326514_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_321465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_325614_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_321456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_326145_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_352164_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_362154_up,T_362154_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_325146_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_326514_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_321465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_325614_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_321456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_326145_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_352164_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_362154_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_325146_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_321564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_356214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_321654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_365214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_352164_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_362154_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_326514_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_321564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_356214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_321465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_325614_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_321654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_365214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_321456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_326145_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_325146_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_326514_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_321564_up,T_321564_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_356214_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_321465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_325614_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_321654_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_365214_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_321456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_326145_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_325146_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_326514_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_321564_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_356214_up,T_356214_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_321465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_325614_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_321654_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_365214_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_321456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_352164_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_362154_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_326514_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_321564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_356214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_321465_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_325614_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_321654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_365214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_321456_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_352164_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_362154_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_326514_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_321564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_356214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_321465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_325614_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_321654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_365214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_321456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_326145_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_325146_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_326514_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_321564_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_356214_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_321465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_325614_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_321654_up,T_321654_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_365214_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_321456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_326145_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_325146_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_326514_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_321564_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_356214_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_321465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_325614_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_321654_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_365214_up,T_365214_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_321456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_352164_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_362154_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_326514_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_321564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_356214_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_321465_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_325614_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_321654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_365214_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_321456_up,(_4*n2m1_2)*_over_nc2));


	};



		if(h3.helicity()==-1&&h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
SM->add_tree(cached_cross_term_md(T_326145_up,T_162354_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_125346_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_126345_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_152364_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_123465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_126534_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_123456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326145_up,T_125634_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_162354_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_125346_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_126345_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_152364_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_156234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_123564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_165234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_352164_up,T_123654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_162354_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_125346_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_126345_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_152364_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_156234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_123564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_165234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_362154_up,T_123654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_162354_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_125346_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_126345_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_152364_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_123465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_126534_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_123456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_325146_up,T_125634_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_125346_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_126345_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_156234_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_123465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_126534_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_123564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_165234_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_123456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_125634_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_326514_up,T_123654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_162354_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_152364_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_156234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_123465_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_126534_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_123564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_165234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_123456_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_125634_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321564_up,T_123654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_162354_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_152364_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_156234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_123465_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_126534_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_123564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_165234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_123456_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_125634_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_356214_up,T_123654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_125346_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_126345_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_156234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_123465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_126534_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_123564_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_165234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_123456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_125634_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321465_up,T_123654_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_125346_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_126345_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_156234_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_123465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_126534_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_123564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_165234_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_123456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_125634_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_325614_up,T_123654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_162354_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_152364_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_156234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_123465_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_126534_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_123564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_165234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_123456_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_125634_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321654_up,T_123654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_162354_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_152364_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_156234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_123465_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_126534_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_123564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_165234_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_123456_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_125634_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_123654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_125346_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_126345_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_156234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_123465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_126534_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_123564_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_165234_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_123456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_125634_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_321456_up,T_123654_up,c*_4*n2m1_2*_over_nc));
		};
		
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
SM->add_tree(cached_cross_term_md(T_162354_up,T_326145_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_352164_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_362154_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_325146_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_321564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_356214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_321654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_365214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_326145_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_352164_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_362154_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_325146_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_326514_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_321465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_325614_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_321456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_326145_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_352164_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_362154_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_325146_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_326514_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_321465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_325614_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_321456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_326145_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_352164_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_362154_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_325146_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_321564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_356214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_321654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_365214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_352164_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_362154_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_326514_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_321564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_356214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_321465_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_325614_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_321654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_365214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_321456_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_326145_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_325146_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_326514_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_321564_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_356214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_321465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_325614_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_321654_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_365214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_321456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_326145_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_325146_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_326514_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_321564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_356214_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_321465_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_325614_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_321654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_365214_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_321456_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_352164_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_362154_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_326514_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_321564_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_356214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_321465_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_325614_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_321654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_365214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_321456_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_352164_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_362154_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_326514_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_321564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_356214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_321465_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_325614_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_321654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_365214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_321456_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_326145_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_325146_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_326514_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_321564_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_356214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_321465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_325614_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_321654_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_365214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_321456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_326145_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_325146_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_326514_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_321564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_356214_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_321465_up,c*_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_325614_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_321654_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_365214_up,c*_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_321456_up,mone*(c*_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_352164_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_362154_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_326514_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_321564_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_356214_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_321465_up,mone*(c*_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_325614_up,c*_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_321654_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_365214_up,mone*(c*_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_321456_up,c*_4*n2m1_2*_over_nc));

	
		};
}

if(color!=1){


		if((h1.helicity()==-1&&h1.helicity()+h4.helicity()==0)){
//SM->add_loop(cached_cross_term_md(L_162354_up,T_162354_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125346_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126345_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_152364_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126534_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125634_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_162354_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_125346_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_126345_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_152364_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_156234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_165234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_162354_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_125346_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_126345_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_152364_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_156234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_165234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_162354_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125346_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126345_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_152364_up,T_152364_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126534_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125634_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_125346_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126345_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_156234_up,T_156234_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126534_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123564_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_165234_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_125634_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123654_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_162354_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_152364_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_156234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123465_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_126534_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_165234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123456_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_125634_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_162354_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_152364_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_156234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_126534_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_165234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_125634_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_125346_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126345_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_156234_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126534_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_123564_up,T_123564_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_165234_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_125634_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123654_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_125346_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126345_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_156234_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126534_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123564_up,_8*nc*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_165234_up,T_165234_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_125634_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123654_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_162354_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_152364_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_156234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123465_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_126534_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_165234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123456_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_125634_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_162354_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_152364_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_156234_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_126534_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_165234_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_125634_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_125346_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126345_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_156234_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126534_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123564_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_165234_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_125634_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_123654_up,T_123654_up,_8*nc*n2m1_2));

	};

		if((h3.helicity()==-1&&h3.helicity()+h4.helicity()==0)){
SM->add_loop(cached_cross_term_md(L_326145_up,T_326145_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_352164_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_362154_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_325146_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_321564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_356214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_321654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_365214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_326145_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_352164_up,T_352164_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_362154_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_325146_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_326514_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_321465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_325614_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_321456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_326145_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_352164_up,_8*nc*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_362154_up,T_362154_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_325146_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_326514_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_321465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_325614_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_321456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_326145_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_352164_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_362154_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_325146_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_321564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_356214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_321654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_365214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_352164_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_362154_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_326514_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_321564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_356214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_321465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_325614_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_321654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_365214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_321456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_326145_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_325146_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_326514_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_321564_up,T_321564_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_356214_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_321465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_325614_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_321654_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_365214_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_321456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_326145_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_325146_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_326514_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_321564_up,_8*nc*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_356214_up,T_356214_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_321465_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_325614_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_321654_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_365214_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_321456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_352164_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_362154_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_326514_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_321564_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_356214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_321465_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_325614_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_321654_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_365214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_321456_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_352164_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_362154_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_326514_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_321564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_356214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_321465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_325614_up,(_8*nc*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_321654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_365214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_321456_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_326145_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_325146_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_326514_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_321564_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_356214_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_321465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_325614_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_321654_up,T_321654_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_365214_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_321456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_326145_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_325146_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_326514_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_321564_up,_8*nc*n2m1)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_356214_up,mone*(_8*nc*n2m1))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_321465_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_325614_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_321654_up,_8*nc*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_365214_up,T_365214_up,_8*nc*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_321456_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_352164_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_362154_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_326514_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_321564_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_356214_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_321465_up,mone*((_8*nc*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_325614_up,(_8*nc*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_321654_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_365214_up,mone*((_8*nc*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_321456_up,(_8*nc*n2m1_2)*_over_nc2));


	};



		if(h3.helicity()==-1&&h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
SM->add_loop(cached_cross_term_md(L_326145_up,T_162354_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_125346_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_126345_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_152364_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_123465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_126534_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_123456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326145_up,T_125634_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_162354_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_125346_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_126345_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_152364_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_156234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_123564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_165234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_352164_up,T_123654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_162354_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_125346_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_126345_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_152364_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_156234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_123564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_165234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_362154_up,T_123654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_162354_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_125346_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_126345_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_152364_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_123465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_126534_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_123456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_325146_up,T_125634_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_125346_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_126345_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_156234_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_123465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_126534_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_123564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_165234_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_123456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_125634_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_326514_up,T_123654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_162354_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_152364_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_156234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_123465_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_126534_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_123564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_165234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_123456_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_125634_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321564_up,T_123654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_162354_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_152364_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_156234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_123465_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_126534_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_123564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_165234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_123456_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_125634_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_356214_up,T_123654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_125346_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_126345_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_156234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_123465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_126534_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_123564_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_165234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_123456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_125634_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321465_up,T_123654_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_125346_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_126345_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_156234_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_123465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_126534_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_123564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_165234_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_123456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_125634_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_325614_up,T_123654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_162354_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_152364_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_156234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_123465_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_126534_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_123564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_165234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_123456_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_125634_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321654_up,T_123654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_162354_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_152364_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_156234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_123465_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_126534_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_123564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_165234_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_123456_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_125634_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_365214_up,T_123654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_125346_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_126345_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_156234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_123465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_126534_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_123564_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_165234_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_123456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_125634_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_321456_up,T_123654_up,c*_8*nc*n2m1_2*_over_nc));
		};
		
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h1.helicity()+h2.helicity()==0)){
SM->add_loop(cached_cross_term_md(L_162354_up,T_326145_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_352164_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_362154_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_325146_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_321564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_356214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_321654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_365214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_326145_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_352164_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_362154_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_325146_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_326514_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_321465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_325614_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_321456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_326145_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_352164_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_362154_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_325146_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_326514_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_321465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_325614_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_321456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_326145_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_352164_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_362154_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_325146_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_321564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_356214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_321654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_365214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_352164_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_362154_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_326514_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_321564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_356214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_321465_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_325614_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_321654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_365214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_321456_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_326145_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_325146_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_326514_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_321564_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_356214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_321465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_325614_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_321654_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_365214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_321456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_326145_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_325146_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_326514_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_321564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_356214_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_321465_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_325614_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_321654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_365214_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_321456_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_352164_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_362154_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_326514_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_321564_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_356214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_321465_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_325614_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_321654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_365214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_321456_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_352164_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_362154_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_326514_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_321564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_356214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_321465_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_325614_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_321654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_365214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_321456_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_326145_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_325146_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_326514_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_321564_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_356214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_321465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_325614_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_321654_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_365214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_321456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_326145_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_325146_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_326514_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_321564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_356214_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_321465_up,c*_8*nc*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_325614_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_321654_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_365214_up,c*_8*nc*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_321456_up,mone*(c*_8*nc*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_352164_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_362154_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_326514_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_321564_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_356214_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_321465_up,mone*(c*_8*nc*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_325614_up,c*_8*nc*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_321654_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_365214_up,mone*(c*_8*nc*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_321456_up,c*_8*nc*n2m1_2*_over_nc));


}

}}
	//--------------------------------------------------------
	// W case: case4q==1 -> identical anti-quarks
	// Notice that we use _up to built W amps!
	//--------------------------------------------------------

if(case4q==1){
		multi_precision_fraction c(1,1);


#if 0
_PRINT("-------------------------");
_PRINT("identical anti-quarks");
_PRINT("-------------------------");
#endif

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,n2m1_2*_4));
		};
		if(hp1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
SM->add_tree(cached_cross_term_md(T_164352_up,T_164352_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_154362_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_156432_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_143562_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_165432_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_143652_up,n2m1_2*_4));
		};

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
SM->add_loop(cached_cross_term_md(L_color_162354_up,T_162354_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_152364_up,T_152364_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_156234_up,T_156234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123564_up,T_123564_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,n2m1_2*nc*_8));
		};
		if(hp1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
SM->add_loop(cached_cross_term_md(L_color_164352_up,T_164352_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_154362_up,T_154362_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_156432_up,T_156432_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_143562_up,T_143562_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_165432_up,T_165432_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_143652_up,T_143652_up,n2m1_2*nc*_8));
		};



if(tree_color!=1){

		if((h1.helicity()==-1&&h1.helicity()+h4.helicity()==0)){
//SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_152364_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_125346_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_126345_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_162354_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123564_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_165234_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123654_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_162354_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123465_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123456_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_152364_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_156234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_126534_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123564_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_165234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_125634_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123654_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_156234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123465_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126534_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_165234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123654_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_125346_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126345_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_156234_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123564_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123654_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_162354_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_152364_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123465_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_125634_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_162354_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_152364_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_156234_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_126534_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123564_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_165234_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123456_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_125634_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123654_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_125346_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126345_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_156234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123465_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126534_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123564_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_165234_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123456_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_125634_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,_4*n2m1_2));

		};

		if(hp1.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){

SM->add_tree(cached_cross_term_md(T_145326_up,T_145326_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_164352_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_154362_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_146325_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_156432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_143562_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_165432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_143652_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_145326_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_164352_up,T_164352_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_154362_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_146325_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_143265_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_146532_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_143256_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_145632_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_145326_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_164352_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_154362_up,T_154362_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_146325_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_143265_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_146532_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_143256_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_145632_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_145326_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_164352_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_154362_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_146325_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_156432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_143562_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_165432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_143652_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_164352_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_154362_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_143265_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_156432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_143562_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_146532_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_143256_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_165432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_143652_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_145632_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_145326_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_146325_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_143265_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_156432_up,T_156432_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_143562_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_146532_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_143256_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_165432_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_143652_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_145632_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_145326_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_146325_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_143265_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_156432_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_143562_up,T_143562_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_146532_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_143256_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_165432_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_143652_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_145632_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_164352_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_154362_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_143265_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_156432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_143562_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_146532_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_143256_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_165432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_143652_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_145632_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_164352_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_154362_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143265_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_156432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143562_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_146532_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143256_up,(_4*n2m1_2)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_165432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143652_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_145632_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_145326_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_146325_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_143265_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_156432_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_143562_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_146532_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_143256_up,mone*((_4*n2m1_2)*_over_nc2))); 
//SM->add_tree(cached_cross_term_md(T_165432_up,T_165432_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_143652_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_145632_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_145326_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_146325_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_143265_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_156432_up,_4*n2m1)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_143562_up,mone*(_4*n2m1))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_146532_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_143256_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_165432_up,_4*n2m1)); 
//SM->add_tree(cached_cross_term_md(T_143652_up,T_143652_up,_4*n2m1_2)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_145632_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_164352_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_154362_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_143265_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_156432_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_143562_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_146532_up,mone*((_4*n2m1)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_143256_up,(_4*n2m1)*_over_nc2)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_165432_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_143652_up,mone*((_4*n2m1_2)*_over_nc2))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_145632_up,(_4*n2m1_2)*_over_nc2));
		
		};

		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h3.helicity()+h4.helicity()==0)){
SM->add_tree(cached_cross_term_md(T_162354_up,T_145326_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_164352_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_154362_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_146325_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_156432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_143562_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_165432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_143652_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_145326_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_164352_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_154362_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_146325_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_143265_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_146532_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_143256_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_145632_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_145326_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_164352_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_154362_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_146325_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_143265_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_146532_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_143256_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_145632_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_145326_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_164352_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_154362_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_146325_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_156432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_143562_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_165432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_143652_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_164352_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_154362_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_143265_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_156432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_143562_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_146532_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_143256_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_165432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_143652_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_145632_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_145326_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_146325_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_143265_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_156432_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_143562_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_146532_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_143256_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_165432_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_143652_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_145632_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_145326_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_146325_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_143265_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_156432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_143562_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_146532_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_143256_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_165432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_143652_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_145632_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_164352_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_154362_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_143265_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_156432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_143562_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_146532_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_143256_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_165432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_143652_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_145632_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_164352_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_154362_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_143265_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_156432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_143562_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_146532_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_143256_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_143652_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_145632_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_145326_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_146325_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_143265_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_156432_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_143562_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_146532_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_143256_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_165432_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_143652_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_145632_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_145326_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_146325_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_143265_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_156432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_143562_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_146532_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_143256_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_165432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_143652_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_145632_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_164352_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_154362_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_143265_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_156432_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_143562_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_146532_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_143256_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_165432_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_143652_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_145632_up,_4*n2m1_2*_over_nc));

		};
		if(h1.helicity()==-1&&(h3.helicity()+h4.helicity()==0)&&(h1.helicity()+h4.helicity()==0)){

SM->add_tree(cached_cross_term_md(T_145326_up,T_162354_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_125346_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_126345_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_152364_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_123465_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_126534_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_123456_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145326_up,T_125634_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_162354_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_125346_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_126345_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_152364_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_156234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_123564_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_165234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_164352_up,T_123654_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_162354_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_125346_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_126345_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_152364_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_156234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_123564_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_165234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_154362_up,T_123654_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_162354_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_125346_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_126345_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_152364_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_123465_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_126534_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_123456_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_146325_up,T_125634_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_125346_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_126345_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_156234_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_123465_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_126534_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_123564_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_165234_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_123456_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_125634_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143265_up,T_123654_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_162354_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_152364_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_156234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_123465_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_126534_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_123564_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_165234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_123456_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_125634_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_156432_up,T_123654_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_162354_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_152364_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_156234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_123465_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_126534_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_123564_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_165234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_123456_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_125634_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143562_up,T_123654_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_125346_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_126345_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_156234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_123465_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_126534_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_123564_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_165234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_123456_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_125634_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_146532_up,T_123654_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_125346_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_126345_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_156234_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_123465_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_126534_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_123564_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_165234_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_123456_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_125634_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_123654_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_162354_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_152364_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_156234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_123465_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_126534_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_123564_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_165234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_123456_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_125634_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_123654_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_162354_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_152364_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_156234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_123465_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_126534_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_123564_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_165234_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_123456_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_125634_up,_4*n2m1_2*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_143652_up,T_123654_up,mone*(_4*n2m1_2*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_125346_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_126345_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_156234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_123465_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_126534_up,_4*n2m1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_123564_up,mone*(_4*n2m1*_over_nc))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_165234_up,_4*n2m1*_over_nc)); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_123456_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_125634_up,mone*(_4*n2m1_2*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145632_up,T_123654_up,_4*n2m1_2*_over_nc));

		};

}

if(color!=1){
		if((h1.helicity()==-1&&h1.helicity()+h4.helicity()==0)){
//SM->add_loop(cached_cross_term_md(L_162354_up,T_162354_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125346_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126345_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_152364_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123465_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126534_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123456_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125634_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_162354_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_125346_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_126345_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_152364_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_156234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123564_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_165234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123654_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_162354_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_125346_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_126345_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_152364_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_156234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123564_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_165234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123654_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_162354_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125346_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126345_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_152364_up,T_152364_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123465_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126534_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123456_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125634_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_125346_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126345_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_156234_up,T_156234_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123465_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126534_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123564_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_165234_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123456_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_125634_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123654_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_162354_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_152364_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_156234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123465_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_126534_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123564_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_165234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123456_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_125634_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123654_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_162354_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_152364_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_156234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123465_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_126534_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123564_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_165234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123456_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_125634_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123654_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_125346_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126345_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_156234_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123465_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126534_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_123564_up,T_123564_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_165234_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123456_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_125634_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123654_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_125346_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126345_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_156234_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123465_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126534_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123564_up,nc*_8*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_165234_up,T_165234_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123456_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_125634_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123654_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_162354_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_152364_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_156234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123465_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_126534_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123564_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_165234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123456_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_125634_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123654_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_162354_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_152364_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_156234_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123465_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_126534_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123564_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_165234_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123456_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_125634_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123654_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_125346_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126345_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_156234_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123465_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126534_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123564_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_165234_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123456_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_125634_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_123654_up,T_123654_up,nc*_8*n2m1_2));

		};

		if(hp1.helicity()==-1&&(hp1.helicity()+hp2.helicity()==0)){

SM->add_loop(cached_cross_term_md(L_145326_up,T_145326_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_164352_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_154362_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_146325_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_156432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_143562_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_165432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_143652_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_145326_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_164352_up,T_164352_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_154362_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_146325_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_143265_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_146532_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_143256_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_145632_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_145326_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_164352_up,nc*_8*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_154362_up,T_154362_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_146325_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_143265_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_146532_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_143256_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_145632_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_145326_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_164352_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_154362_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_146325_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_156432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_143562_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_165432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_143652_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_164352_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_154362_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_143265_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_156432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_143562_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_146532_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_143256_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_165432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_143652_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_145632_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_145326_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_146325_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_143265_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_156432_up,T_156432_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_143562_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_146532_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_143256_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_165432_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_143652_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_145632_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_145326_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_146325_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_143265_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_156432_up,nc*_8*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_143562_up,T_143562_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_146532_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_143256_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_165432_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_143652_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_145632_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_164352_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_154362_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_143265_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_156432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_143562_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_146532_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_143256_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_165432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_143652_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_145632_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_164352_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_154362_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_143265_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_156432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_143562_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_146532_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_143256_up,(nc*_8*n2m1_2)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_165432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_143652_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_145632_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_145326_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_146325_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_143265_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_156432_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_143562_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_146532_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_143256_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
//SM->add_loop(cached_cross_term_md(L_165432_up,T_165432_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_143652_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_145632_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_145326_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_146325_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_143265_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_156432_up,nc*_8*n2m1)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_143562_up,mone*(nc*_8*n2m1))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_146532_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_143256_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_165432_up,nc*_8*n2m1)); 
//SM->add_loop(cached_cross_term_md(L_143652_up,T_143652_up,nc*_8*n2m1_2)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_145632_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_164352_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_154362_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_143265_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_156432_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_143562_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_146532_up,mone*((nc*_8*n2m1)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_143256_up,(nc*_8*n2m1)*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_165432_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_143652_up,mone*((nc*_8*n2m1_2)*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_145632_up,(nc*_8*n2m1_2)*_over_nc2));
		
		};

		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(h3.helicity()+h4.helicity()==0)){
SM->add_loop(cached_cross_term_md(L_162354_up,T_145326_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_164352_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_154362_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_146325_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_156432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_143562_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_165432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_143652_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_145326_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_164352_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_154362_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_146325_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_143265_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_146532_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_143256_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_145632_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_145326_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_164352_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_154362_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_146325_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_143265_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_146532_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_143256_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_145632_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_145326_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_164352_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_154362_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_146325_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_156432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_143562_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_165432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_143652_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_164352_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_154362_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_143265_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_156432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_143562_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_146532_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_143256_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_165432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_143652_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_145632_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_145326_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_146325_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_143265_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_156432_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_143562_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_146532_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_143256_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_165432_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_143652_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_145632_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_145326_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_146325_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_143265_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_156432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_143562_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_146532_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_143256_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_165432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_143652_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_145632_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_164352_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_154362_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_143265_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_156432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_143562_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_146532_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_143256_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_165432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_143652_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_145632_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_164352_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_154362_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_143265_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_156432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_143562_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_146532_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_143256_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_165432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_143652_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_145632_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_145326_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_146325_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_143265_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_156432_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_143562_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_146532_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_143256_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_165432_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_143652_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_145632_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_145326_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_146325_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_143265_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_156432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_143562_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_146532_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_143256_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_165432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_143652_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_145632_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_164352_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_154362_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_143265_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_156432_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_143562_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_146532_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_143256_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_165432_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_143652_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_145632_up,nc*_8*n2m1_2*_over_nc));

		};
		if(h1.helicity()==-1&&(h3.helicity()+h4.helicity()==0)&&(h1.helicity()+h4.helicity()==0)){

SM->add_loop(cached_cross_term_md(L_145326_up,T_162354_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_125346_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_126345_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_152364_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_123465_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_126534_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_123456_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145326_up,T_125634_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_162354_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_125346_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_126345_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_152364_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_156234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_123564_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_165234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_164352_up,T_123654_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_162354_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_125346_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_126345_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_152364_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_156234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_123564_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_165234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_154362_up,T_123654_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_162354_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_125346_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_126345_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_152364_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_123465_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_126534_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_123456_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_146325_up,T_125634_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_125346_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_126345_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_156234_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_123465_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_126534_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_123564_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_165234_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_123456_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_125634_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143265_up,T_123654_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_162354_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_152364_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_156234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_123465_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_126534_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_123564_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_165234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_123456_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_125634_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_156432_up,T_123654_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_162354_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_152364_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_156234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_123465_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_126534_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_123564_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_165234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_123456_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_125634_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143562_up,T_123654_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_125346_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_126345_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_156234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_123465_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_126534_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_123564_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_165234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_123456_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_125634_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_146532_up,T_123654_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_125346_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_126345_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_156234_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_123465_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_126534_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_123564_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_165234_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_123456_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_125634_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_143256_up,T_123654_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_162354_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_152364_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_156234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_123465_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_126534_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_123564_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_165234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_123456_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_125634_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_165432_up,T_123654_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_162354_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_152364_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_156234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_123465_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_126534_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_123564_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_165234_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_123456_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_125634_up,nc*_8*n2m1_2*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_143652_up,T_123654_up,mone*(nc*_8*n2m1_2*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_125346_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_126345_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_156234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_123465_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_126534_up,nc*_8*n2m1*_over_nc3)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_123564_up,mone*(nc*_8*n2m1*_over_nc))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_165234_up,nc*_8*n2m1*_over_nc)); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_123456_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_125634_up,mone*(nc*_8*n2m1_2*_over_nc3))); 
SM->add_loop(cached_cross_term_md(L_145632_up,T_123654_up,nc*_8*n2m1_2*_over_nc));

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

		if(hp1.helicity()==-1&&h1.helicity()+h4.helicity()==0){

SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,n2m1*n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,n2m1*n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,n2m1*n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,n2m1*n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,n2m1*n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,n2m1*n2m1*_4));

SM->add_loop(cached_cross_term_md(L_color_162354_up,T_162354_up,n2m1*n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_152364_up,T_152364_up,n2m1*n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_156234_up,T_156234_up,n2m1*n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123564_up,T_123564_up,n2m1*n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,n2m1*n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,n2m1*n2m1*nc*_8));
		
		};
if(tree_color!=1){
		
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){

//SM->add_tree(cached_cross_term_md(T_162354_up,T_162354_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125346_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126345_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_152364_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123465_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_126534_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_123456_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_162354_up,T_125634_up,((n2m1)*_over_nc2)*_4)); 

SM->add_tree(cached_cross_term_md(T_125346_up,T_162354_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_125346_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_126345_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_152364_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_156234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123564_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_165234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125346_up,T_123654_up,((n2m1)*_over_nc2)*_4)); 

SM->add_tree(cached_cross_term_md(T_126345_up,T_162354_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_125346_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_126345_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_152364_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_156234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123564_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_165234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_126345_up,T_123654_up,(mone*n2m1_2)*_over_nc2*_4)); 

SM->add_tree(cached_cross_term_md(T_152364_up,T_162354_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125346_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126345_up,(mone*n2m1_2)*_over_nc2*_4)); 
//SM->add_tree(cached_cross_term_md(T_152364_up,T_152364_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123465_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_126534_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_123456_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_152364_up,T_125634_up,(mone*n2m1_2)*_over_nc2*_4)); 

SM->add_tree(cached_cross_term_md(T_156234_up,T_125346_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126345_up,(mone*n2m1_2)*_over_nc2*_4)); 
//SM->add_tree(cached_cross_term_md(T_156234_up,T_156234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123465_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_126534_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123564_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_165234_up,-n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123456_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_125634_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_156234_up,T_123654_up,n2m1*_4)); 

SM->add_tree(cached_cross_term_md(T_123465_up,T_162354_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_152364_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_156234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123465_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_126534_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123564_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_165234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123456_up,mone*((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_125634_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123465_up,T_123654_up,((n2m1)*_over_nc2)*_4)); 

SM->add_tree(cached_cross_term_md(T_126534_up,T_162354_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_152364_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_156234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123465_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_126534_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123564_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_165234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123456_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_125634_up,mone*((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_126534_up,T_123654_up,((n2m1)*_over_nc2)*_4)); 

SM->add_tree(cached_cross_term_md(T_123564_up,T_125346_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126345_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_156234_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123465_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_126534_up,(mone*n2m1_2)*_over_nc2*_4)); 
//SM->add_tree(cached_cross_term_md(T_123564_up,T_123564_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_165234_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123456_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_125634_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123564_up,T_123654_up,-n2m1*_4)); 

SM->add_tree(cached_cross_term_md(T_165234_up,T_125346_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126345_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_156234_up,-n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123465_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_126534_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123564_up,n2m1*_4)); 
//SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,n2m1_2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123456_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_125634_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_123654_up,n2m1*_4)); 

SM->add_tree(cached_cross_term_md(T_123456_up,T_162354_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_152364_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_156234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123465_up,mone*((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_126534_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123564_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_165234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_125634_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123456_up,T_123654_up,(mone*n2m1_2)*_over_nc2*_4)); 

SM->add_tree(cached_cross_term_md(T_125634_up,T_162354_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_152364_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_156234_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123465_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_126534_up,mone*((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123564_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_165234_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123456_up,(n2m1)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_125634_up,(n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_125634_up,T_123654_up,(mone*n2m1_2)*_over_nc2*_4)); 

SM->add_tree(cached_cross_term_md(T_123654_up,T_125346_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126345_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_156234_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123465_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_126534_up,((n2m1)*_over_nc2)*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123564_up,-n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_165234_up,n2m1*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_123456_up,(mone*n2m1_2)*_over_nc2*_4)); 
SM->add_tree(cached_cross_term_md(T_123654_up,T_125634_up,(mone*n2m1_2)*_over_nc2*_4)); 
//SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,n2m1_2*_4));

		}
}
if(color!=1){
		if(hp1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)){


//SM->add_loop(cached_cross_term_md(L_162354_up,T_162354_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125346_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126345_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_152364_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123465_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_126534_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_123456_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_162354_up,T_125634_up,((n2m1)*_over_nc2)*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_125346_up,T_162354_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_125346_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_126345_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_152364_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_156234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123564_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_165234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125346_up,T_123654_up,((n2m1)*_over_nc2)*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_126345_up,T_162354_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_125346_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_126345_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_152364_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_156234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123564_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_165234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126345_up,T_123654_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_152364_up,T_162354_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125346_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126345_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
//SM->add_loop(cached_cross_term_md(L_152364_up,T_152364_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123465_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_126534_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_123456_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_152364_up,T_125634_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_156234_up,T_125346_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126345_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
//SM->add_loop(cached_cross_term_md(L_156234_up,T_156234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123465_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_126534_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123564_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_165234_up,-n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123456_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_125634_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_156234_up,T_123654_up,n2m1*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_123465_up,T_162354_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_152364_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_156234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123465_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_126534_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123564_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_165234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123456_up,mone*((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_125634_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123465_up,T_123654_up,((n2m1)*_over_nc2)*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_126534_up,T_162354_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_152364_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_156234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123465_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_126534_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123564_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_165234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123456_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_125634_up,mone*((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_126534_up,T_123654_up,((n2m1)*_over_nc2)*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_123564_up,T_125346_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126345_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_156234_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123465_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_126534_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
//SM->add_loop(cached_cross_term_md(L_123564_up,T_123564_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_165234_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123456_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_125634_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123564_up,T_123654_up,-n2m1*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_165234_up,T_125346_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126345_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_156234_up,-n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123465_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_126534_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123564_up,n2m1*nc*_8)); 
//SM->add_loop(cached_cross_term_md(L_165234_up,T_165234_up,n2m1_2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123456_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_125634_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_165234_up,T_123654_up,n2m1*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_123456_up,T_162354_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_152364_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_156234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123465_up,mone*((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_126534_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123564_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_165234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123456_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_125634_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123456_up,T_123654_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_125634_up,T_162354_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_152364_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_156234_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123465_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_126534_up,mone*((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123564_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_165234_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123456_up,(n2m1)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_125634_up,(n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_125634_up,T_123654_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 

SM->add_loop(cached_cross_term_md(L_123654_up,T_125346_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126345_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_156234_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123465_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_126534_up,((n2m1)*_over_nc2)*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123564_up,-n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_165234_up,n2m1*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_123456_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
SM->add_loop(cached_cross_term_md(L_123654_up,T_125634_up,(mone*n2m1_2)*_over_nc2*nc*_8)); 
//SM->add_loop(cached_cross_term_md(L_123654_up,T_123654_up,n2m1_2*nc*_8));


// for subleading color we need to add color structures that do not appear for the born pieces

			};
}
}
return SM;
	}

///////////////
//Z/gamma cases
///////////////

return SM;
}


Virtual_SME* vsme_2q2Q2g2l(std::vector<int> indext,int ns,int n_f,int nc,int photonZW, int case4q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);
	ind.push_back(indext[6]);
	ind.push_back(indext[7]);

	vector<int> ind87;
	ind87.push_back(indext[0]);
	ind87.push_back(indext[1]);
	ind87.push_back(indext[2]);
	ind87.push_back(indext[3]);
	ind87.push_back(indext[4]);
	ind87.push_back(indext[5]);
	ind87.push_back(indext[7]);
	ind87.push_back(indext[6]);

clock_t before, after;
before=clock();

		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,m,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,m,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,p,m,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,p,m,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}

		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,m,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,m,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}		
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
}		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,p,p,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case4q,color,tree_color,lo_or_nlo));

if(photonZW!=3){VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbp,qm,qbm,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qp,qbm,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qp,qbm,qm,qbp,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbm,qp,qbp,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qm,qbp,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_2q_2Q_2g_2l_M2(process(qm,qbp,qp,qbm,p,p,lm,lbp),ind87,ns,n_f,nc,photonZW,lp,case4q,color,tree_color,lo_or_nlo));
}

after=clock();
std::cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<std::endl;

return VSM;
}



BH_Ampl_2q2Q2g2l::BH_Ampl_2q2Q2g2l(int photonZW,int case4q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q2g2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				settings::BH_interface_settings::s_photon_only*photonZW,
				case4q,
				color,
				tree_color,
				lo_or_nlo),
			bhi,	// parent
			8,	// NbrExtParticles
			4,	// NbrPowersOfAlphaS
			2,	// NbrPowersOfAlphaQED
			8,	// GeVdim
			1.	// factor
		),
		momenta_assignment(mom_assignment), d_color(color)
{}


double BH_Ampl_2q2Q2g2l::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(d_color){
		//Full color
		case 0 : return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
		//Leading color
        //notice full color renormalization
		//case 1: return get_single_pole()+(3 /*NbrPowersOfAlphaS*/* 33/6. /*beta_0*/ *get_double_pole()/(-12.) /* = lc born tree */  );
		case 1: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/ *get_double_pole()/(-12.) /* = lc born tree */  );
		//Full-Leading color
		case 2: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* ( 23/6. -  33/6. /*beta_0*/ *(25/3./9. - get_double_pole()/(-12.)  /* = lc born tree */ ))  );
		}
	case 2:	return get_double_pole();
	}

};




}
