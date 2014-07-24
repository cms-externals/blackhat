/*
 * matrix elements for 6q2l
 *
 *  Created on: Sept 6, 2009
 *
 *  only leading color born and NLO terms are available in the moment
 *  to do: (1) full color trees
 *  	(2) correction of numerical instabilities
 *  	(3) correction of certain 6-quark trees which seem off in the moment
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
//#define _PHOTON_ONLY 1     // 0 photon only; 1 for all other    replaced by a setting

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

vector<int> index(vector<int> & indin, int num){
	  int index(num);
	  int n=int(log10(num))+1;
	  vector<int> ind;
	  for(int i=0;i<n;i++) {
	       int j=index/int(std::pow(10.,n-1-i));
	       index-=j*int(std::pow(10.,n-1-i));
	       ind.push_back(indin[j-1]);}
	  return ind;
}


process qQQGGq(process & pro, int num){
	ph_type h1,h2,h3,h4,h5,h6,h7,h8;
        int n=pro.n();
        int index(num);
        vector<int> ind;
        for(int i=0;i<n;i++) {
                int j=index/int(std::pow(10.,n-1-i));
                index-=j*int(std::pow(10.,n-1-i));
                ind.push_back(j);}
	int flavor(1);
	h1=particle_ID(quark,pro.p(ind[1-1]).helicity(),flavor);
	h2=particle_ID(gluino,pro.p(ind[2-1]).helicity(),flavor,true);
	h3=particle_ID(gluino,pro.p(ind[3-1]).helicity(),flavor);
	h4=particle_ID(gluino,pro.p(ind[4-1]).helicity(),flavor+1,true);
	h5=particle_ID(gluino,pro.p(ind[5-1]).helicity(),flavor+1);
	h6=particle_ID(quark,pro.p(ind[6-1]).helicity(),flavor,true);
	h7=pro.p(ind[7-1]);
	h8=pro.p(ind[8-1]);

	return process(h1,h2,h3,h4,h5,h6,h7,h8);	
}


bool helmatch_qQQGGq(process & pro, int num){
	int h1,h2,h3,h4,h5,h6;
        int n=pro.n();
        int index(num);
        vector<int> ind;
        for(int i=0;i<n;i++) {
                int j=index/int(std::pow(10.,n-1-i));
                index-=j*int(std::pow(10.,n-1-i));
                ind.push_back(j);}
	int flavor(1);
	h1=pro.p(ind[1-1]).helicity();
	h2=pro.p(ind[2-1]).helicity();
	h3=pro.p(ind[3-1]).helicity();
	h4=pro.p(ind[4-1]).helicity();
	h5=pro.p(ind[5-1]).helicity();
	h6=pro.p(ind[6-1]).helicity();
		
	return (h1==-1&&h6==1&&h2==-h3&&h4==-h5);
}


//expected indices:
// q Qb Q Gb G qb l lb
// Q and G are neutral distinct quark flavors
partial_amplitude_cached* A_loop_6q_2l_qQQGGq(const process & pro,const vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
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





Squared_ME* A_loop_6q_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int photonZW, const ph_type e_ph_type,int case6q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

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

//quarks labels i1-i6
//leptons i7-i8
//q G1b G1 G2b G2 gb lm lbp

int flavor(1);

ph_type h1=pro.p(1);
ph_type h2=particle_ID(gluino,pro.p(2).helicity(),flavor,true);
ph_type h3=particle_ID(gluino,pro.p(3).helicity(),flavor);
ph_type h4=particle_ID(gluino,pro.p(4).helicity(),flavor+1,true);
ph_type h5=particle_ID(gluino,pro.p(5).helicity(),flavor+1);
ph_type h6=pro.p(6);
ph_type h7=pro.p(7);
ph_type h8=pro.p(8);


//-------------------------------------------------------
// tree labels

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

process pro_145236=process(h1,h4,h5,h2,h3,h6,h7,h8);
vector<int> ind_145236;
ind_145236.push_back(i1);
ind_145236.push_back(i4);
ind_145236.push_back(i5);
ind_145236.push_back(i2);
ind_145236.push_back(i3);
ind_145236.push_back(i6);
ind_145236.push_back(i7);
ind_145236.push_back(i8);


vector<ph_type> _ph_type_1;
_ph_type_1.push_back(h1);
_ph_type_1.push_back(e_ph_type);

vector<ph_type> _ph_type_3;
_ph_type_3.push_back(h3);
_ph_type_3.push_back(e_ph_type);

vector<ph_type> _ph_type_5;
_ph_type_5.push_back(h5);
_ph_type_5.push_back(e_ph_type);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn_up_1(0,photonZW,0,i7,i8,_ph_type_1);
	prop_hel_fn _prop_hel_fn_up_3(0,photonZW,0,i7,i8,_ph_type_3);
	prop_hel_fn _prop_hel_fn_up_5(0,photonZW,0,i7,i8,_ph_type_5);
	prop_hel_fn _prop_hel_fn_down_1(1,photonZW,0,i7,i8,_ph_type_1);
	prop_hel_fn _prop_hel_fn_down_3(1,photonZW,0,i7,i8,_ph_type_3);
	prop_hel_fn _prop_hel_fn_down_5(1,photonZW,0,i7,i8,_ph_type_5);
//-------------------------------------------


size_t T_123456_up;size_t L_color_123456_up;
size_t T_123654_up;size_t L_color_123654_up;
size_t T_125436_up;size_t L_color_125436_up;
size_t T_125634_up;size_t L_color_125634_up;
size_t T_143256_up;size_t L_color_143256_up;
size_t T_143652_up;size_t L_color_143652_up;
size_t T_145236_up;size_t L_color_145236_up;
size_t T_145632_up;size_t L_color_145632_up;
size_t T_163254_up;size_t L_color_163254_up;
size_t T_163452_up;size_t L_color_163452_up;
size_t T_165234_up;size_t L_color_165234_up;
size_t T_165432_up;size_t L_color_165432_up;

size_t T_321456_up;size_t L_color_321456_up;
size_t T_325416_up;size_t L_color_325416_up;
size_t T_341256_up;size_t L_color_341256_up;
size_t T_345216_up;size_t L_color_345216_up;
size_t T_521436_up;size_t L_color_521436_up;
size_t T_523416_up;size_t L_color_523416_up;
size_t T_541236_up;size_t L_color_541236_up;
size_t T_543216_up;size_t L_color_543216_up;

size_t T_321654_up;size_t L_color_321654_up;
size_t T_365214_up;size_t L_color_365214_up;


//-------------------------------------------

/*	
		T_123456_up=SM->add(new Tree_with_prefactor(pro_123456,_prop_hel_fn_up_1),ind_123456);
		T_145236_up=SM->add(new Tree_with_prefactor(pro_145236,_prop_hel_fn_up_1),ind_145236);
*/		
if(photonZW==3){
switch(case6q){
	case 34:{	
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));

L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
	};
	} break;
	case 11: {
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));
L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,12365478)){
T_123654_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12365478),index(ind,12365478),_prop_hel_fn_up_1));
T_165234_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16523478),index(ind,16523478),_prop_hel_fn_up_1));
L_color_123654_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12365478),index(ind,12365478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_165234_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16523478),index(ind,16523478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,12543678)){
T_125436_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12543678),index(ind,12543678),_prop_hel_fn_up_1));
T_143256_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14325678),index(ind,14325678),_prop_hel_fn_up_1));
L_color_125436_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12543678),index(ind,12543678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_143256_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14325678),index(ind,14325678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,12563478)){
T_125634_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12563478),index(ind,12563478),_prop_hel_fn_up_1));
T_163254_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16325478),index(ind,16325478),_prop_hel_fn_up_1));
L_color_125634_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12563478),index(ind,12563478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_163254_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16325478),index(ind,16325478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,14365278)){
T_143652_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14365278),index(ind,14365278),_prop_hel_fn_up_1));
T_165432_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16543278),index(ind,16543278),_prop_hel_fn_up_1));
L_color_143652_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14365278),index(ind,14365278),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_165432_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16543278),index(ind,16543278),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,14563278)){
T_145632_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14563278),index(ind,14563278),_prop_hel_fn_up_1));
T_163452_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16345278),index(ind,16345278),_prop_hel_fn_up_1));
L_color_145632_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14563278),index(ind,14563278),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_163452_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16345278),index(ind,16345278),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};	} break;
	case 22:{
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));
L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,12543678)){
T_125436_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12543678),index(ind,12543678),_prop_hel_fn_up_1));
T_143256_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14325678),index(ind,14325678),_prop_hel_fn_up_1));
L_color_125436_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12543678),index(ind,12543678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_143256_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14325678),index(ind,14325678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,32145678)){
T_321456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,32145678),index(ind,32145678),_prop_hel_fn_up_3));
T_345216_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,34521678),index(ind,34521678),_prop_hel_fn_up_3));
L_color_321456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,32145678),index(ind,32145678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_345216_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,34521678),index(ind,34521678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,32541678)){
T_325416_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,32541678),index(ind,32541678),_prop_hel_fn_up_3));
T_341256_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,34125678),index(ind,34125678),_prop_hel_fn_up_3));
L_color_325416_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,32541678),index(ind,32541678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_341256_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,34125678),index(ind,34125678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,52143678)){
T_521436_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,52143678),index(ind,52143678),_prop_hel_fn_up_5));
T_543216_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,54321678),index(ind,54321678),_prop_hel_fn_up_5));
L_color_521436_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,52143678),index(ind,52143678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_543216_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,54321678),index(ind,54321678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
if(helmatch_qQQGGq(pro,52341678)){
T_523416_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,52341678),index(ind,52341678),_prop_hel_fn_up_5));
T_541236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,54123678),index(ind,54123678),_prop_hel_fn_up_5));
L_color_523416_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,52341678),index(ind,52341678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_541236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,54123678),index(ind,54123678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
	} break;

	case 21: {
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));
L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};

if(helmatch_qQQGGq(pro,12365478)){
T_123654_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12365478),index(ind,12365478),_prop_hel_fn_up_1));
T_165234_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16523478),index(ind,16523478),_prop_hel_fn_up_1));
L_color_123654_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12365478),index(ind,12365478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_165234_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16523478),index(ind,16523478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};

if(helmatch_qQQGGq(pro,32145678)){
T_321456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,32145678),index(ind,32145678),_prop_hel_fn_up_3));
T_345216_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,34521678),index(ind,34521678),_prop_hel_fn_up_3));
L_color_321456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,32145678),index(ind,32145678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_345216_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,34521678),index(ind,34521678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};

if(helmatch_qQQGGq(pro,32165478)){
T_321654_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,32165478),index(ind,32165478),_prop_hel_fn_up_3));
T_365214_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,36521478),index(ind,36521478),_prop_hel_fn_up_3));
L_color_321654_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,32165478),index(ind,32165478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_365214_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,36521478),index(ind,36521478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
};
	} break;


	case 23: {
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));
L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}

if(helmatch_qQQGGq(pro,32145678)){
T_321456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,32145678),index(ind,32145678),_prop_hel_fn_up_3));
T_345216_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,34521678),index(ind,34521678),_prop_hel_fn_up_3));
L_color_321456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,32145678),index(ind,32145678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_345216_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,34521678),index(ind,34521678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}		 
		 
		 } break;
	case 31: {
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));

L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}
		 
if(helmatch_qQQGGq(pro,12365478)){
T_123654_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12365478),index(ind,12365478),_prop_hel_fn_up_1));
T_165234_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,16523478),index(ind,16523478),_prop_hel_fn_up_1));

L_color_123654_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12365478),index(ind,12365478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_165234_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,16523478),index(ind,16523478),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}
		 } break;
	case 33: {
if(helmatch_qQQGGq(pro,12345678)){
T_123456_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12345678),index(ind,12345678),_prop_hel_fn_up_1));
T_145236_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14523678),index(ind,14523678),_prop_hel_fn_up_1));

L_color_123456_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12345678),index(ind,12345678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_145236_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14523678),index(ind,14523678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}
if(helmatch_qQQGGq(pro,12543678)){
T_125436_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,12543678),index(ind,12543678),_prop_hel_fn_up_1));
T_143256_up=SM->add(new CTree_with_prefactor(qQQGGq(pro,14325678),index(ind,14325678),_prop_hel_fn_up_1));

L_color_125436_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,12543678),index(ind,12543678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
L_color_143256_up=SM->add(A_loop_6q_2l_qQQGGq(qQQGGq(pro,14325678),index(ind,14325678),n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
}
	} break;
};
}

//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction nc3=nc*nc*nc;
multi_precision_fraction nc4=nc3*nc;
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

// W cases
if(photonZW==3){
	
	//--------------------------------------------------------
	//see BH_interface.cpp file case 260 for values of case6q and their meaning
	//--------------------------------------------------------

switch(case6q){
	case 34:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3));

SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4));
}
	} break;
	case 23:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4));
}
if(helmatch_qQQGGq(pro,32145678)){
SM->add_tree(cached_cross_term_md(T_321456_up,T_321456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_345216_up,T_345216_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_321456_up,T_321456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_345216_up,T_345216_up,_8*nc4)); 
}
		} break;
	case 31:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,12365478)){
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,_8*nc4));
};
		} break;
	case 33:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,12543678)){
SM->add_tree(cached_cross_term_md(T_125436_up,T_125436_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143256_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_125436_up,T_125436_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_143256_up,T_143256_up,_8*nc4)); 
};
		} break;
	case 21:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4));
};

if(helmatch_qQQGGq(pro,12365478)){
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,32145678)){
SM->add_tree(cached_cross_term_md(T_321456_up,T_321456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_345216_up,T_345216_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_321456_up,T_321456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_345216_up,T_345216_up,_8*nc4));
};

if(helmatch_qQQGGq(pro,36521478)){
SM->add_tree(cached_cross_term_md(T_365214_up,T_365214_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_321654_up,T_321654_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_365214_up,T_365214_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_321654_up,T_321654_up,_8*nc4));
};		
if(helmatch_qQQGGq(pro,14523678)&&helmatch_qQQGGq(pro,36521478)){
SM->add_tree(cached_cross_term_md(T_145236_up,T_365214_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_365214_up,T_145236_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_365214_up,T_145236_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_365214_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,16523478)&&helmatch_qQQGGq(pro,34521678)){
SM->add_tree(cached_cross_term_md(T_165234_up,T_345216_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_345216_up,T_165234_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_345216_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_345216_up,T_165234_up,_8*nc4)); 
};
		} break;
	case 11:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,12365478)){
SM->add_tree(cached_cross_term_md(T_123654_up,T_123654_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_165234_up,T_165234_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_123654_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_165234_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,12543678)){
SM->add_tree(cached_cross_term_md(T_125436_up,T_125436_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_143256_up,T_143256_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_125436_up,T_125436_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_143256_up,T_143256_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,12563478)){
SM->add_tree(cached_cross_term_md(T_125634_up,T_125634_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_163254_up,T_163254_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_125634_up,T_125634_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_163254_up,T_163254_up,_8*nc4)); 
};
if(helmatch_qQQGGq(pro,14365278)){
SM->add_tree(cached_cross_term_md(T_143652_up,T_143652_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_165432_up,T_165432_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_143652_up,T_143652_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_165432_up,T_165432_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,14563278)){
SM->add_tree(cached_cross_term_md(T_145632_up,T_145632_up,_4*nc3)); 
SM->add_tree(cached_cross_term_md(T_163452_up,T_163452_up,_4*nc3)); 
SM->add_loop(cached_cross_term_md(L_color_145632_up,T_145632_up,_8*nc4)); 
SM->add_loop(cached_cross_term_md(L_color_163452_up,T_163452_up,_8*nc4)); 
};
//mixed leading_color terms
if(helmatch_qQQGGq(pro,12345678)&&helmatch_qQQGGq(pro,12563478)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_125634_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_125634_up,T_123456_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_125634_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_125634_up,T_123456_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,12365478)&&helmatch_qQQGGq(pro,12543678)){
SM->add_tree(cached_cross_term_md(T_123654_up,T_125436_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_125436_up,T_123654_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123654_up,T_125436_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_125436_up,T_123654_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,14325678)&&helmatch_qQQGGq(pro,14563278)){
SM->add_tree(cached_cross_term_md(T_143256_up,T_145632_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_145632_up,T_143256_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_143256_up,T_145632_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_145632_up,T_143256_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,14365278)&& helmatch_qQQGGq(pro,14523678)){
SM->add_tree(cached_cross_term_md(T_143652_up,T_145236_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_145236_up,T_143652_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_143652_up,T_145236_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_143652_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,16325478)&&helmatch_qQQGGq(pro,16543278)){
SM->add_tree(cached_cross_term_md(T_163254_up,T_165432_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_165432_up,T_163254_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_163254_up,T_165432_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_165432_up,T_163254_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,16345278)&&helmatch_qQQGGq(pro,16523478)){
SM->add_tree(cached_cross_term_md(T_163452_up,T_165234_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_165234_up,T_163452_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_163452_up,T_165234_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_165234_up,T_163452_up,_8*nc4));
};

		} break;
	case 22:{
if(helmatch_qQQGGq(pro,12345678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_123456_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_145236_up,T_145236_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_123456_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_145236_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,12543678)){
SM->add_tree(cached_cross_term_md(T_125436_up,T_125436_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_143256_up,T_143256_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_125436_up,T_125436_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_143256_up,T_143256_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,32145678)){
SM->add_tree(cached_cross_term_md(T_321456_up,T_321456_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_345216_up,T_345216_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_321456_up,T_321456_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_345216_up,T_345216_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,32541678)){
SM->add_tree(cached_cross_term_md(T_325416_up,T_325416_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_341256_up,T_341256_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_325416_up,T_325416_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_341256_up,T_341256_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,52143678)){
SM->add_tree(cached_cross_term_md(T_521436_up,T_521436_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_543216_up,T_543216_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_521436_up,T_521436_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_543216_up,T_543216_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,52341678)){
SM->add_tree(cached_cross_term_md(T_541236_up,T_541236_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_523416_up,T_523416_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_541236_up,T_541236_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_523416_up,T_523416_up,_8*nc4));
};
//-
if(helmatch_qQQGGq(pro,12345678)&&helmatch_qQQGGq(pro,34125678)){
SM->add_tree(cached_cross_term_md(T_123456_up,T_341256_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_341256_up,T_123456_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_123456_up,T_341256_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_341256_up,T_123456_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,14523678)&&helmatch_qQQGGq(pro,52143678)){
SM->add_tree(cached_cross_term_md(T_145236_up,T_521436_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_521436_up,T_145236_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_145236_up,T_521436_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_521436_up,T_145236_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,12543678)&&helmatch_qQQGGq(pro,54123678)){
SM->add_tree(cached_cross_term_md(T_125436_up,T_541236_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_541236_up,T_125436_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_125436_up,T_541236_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_541236_up,T_125436_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,34521678)&&helmatch_qQQGGq(pro,52341678)){
SM->add_tree(cached_cross_term_md(T_345216_up,T_523416_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_523416_up,T_345216_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_345216_up,T_523416_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_523416_up,T_345216_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,14325678)&&helmatch_qQQGGq(pro,32145678)){
SM->add_tree(cached_cross_term_md(T_143256_up,T_321456_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_321456_up,T_143256_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_143256_up,T_321456_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_321456_up,T_143256_up,_8*nc4));
};
if(helmatch_qQQGGq(pro,32541678)&&helmatch_qQQGGq(pro,54321678)){
SM->add_tree(cached_cross_term_md(T_325416_up,T_543216_up,_4*nc3));
SM->add_tree(cached_cross_term_md(T_543216_up,T_325416_up,_4*nc3));
SM->add_loop(cached_cross_term_md(L_color_325416_up,T_543216_up,_8*nc4));
SM->add_loop(cached_cross_term_md(L_color_543216_up,T_325416_up,_8*nc4));
};
	} break;
}

return SM;
	}
}



Virtual_SME* vsme_6q2l(std::vector<int> indext,int ns,int n_f,int nc,int photonZW, int case6q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
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

/* not needed for W-decays
	vector<int> ind87;
	ind87.push_back(indext[0]);
	ind87.push_back(indext[1]);
	ind87.push_back(indext[2]);
	ind87.push_back(indext[3]);
	ind87.push_back(indext[4]);
	ind87.push_back(indext[5]);
	ind87.push_back(indext[7]);
	ind87.push_back(indext[6]);
*/
clock_t before, after;
before=clock();

//notice: We only use helicity configurations that con contribute to W-decays. 
//	  We will have to add more when we consider Zs&photons.

//no W-case	VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qp,qbm,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));

		VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qp,qbp,qm,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qm,qbp,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbm,qp,qbp,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbp,qp,qbm,qm,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbp,qm,qbm,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qp,qbm,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qp,qbm,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qm,qbm,qp,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbm,qp,qbm,qp,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));

		VSM->add(A_loop_6q_2l_M2(process(qp,qbp,qm,qbp,qm,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qp,qbp,qm,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qm,qbp,qp,qbm,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbp,qm,qbm,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qp,qbm,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qm,qbm,qp,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qp,qbm,qm,qbp,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbm,qp,qbp,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));
		VSM->add(A_loop_6q_2l_M2(process(qm,qbm,qm,qbp,qp,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));

		VSM->add(A_loop_6q_2l_M2(process(qm,qbp,qm,qbp,qm,qbp,lm,lbp),ind,ns,n_f,nc,photonZW,lm,case6q,color,tree_color,lo_or_nlo));


after=clock();
std::cout<< "total construction time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<std::endl;

return VSM;
}



BH_Ampl_6q2l::BH_Ampl_6q2l(int photonZW,int case6q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_6q2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				settings::BH_interface_settings::s_photon_only*photonZW,
				case6q,
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


double BH_Ampl_6q2l::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(d_color){
		//Full color
		case 0 : return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
		//Leading color\
        //notice full color renormalization
		//case 1: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 33/6. /*beta_0*/ *get_double_pole()/(-9.) /* = lc born tree */  );
		case 1: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/ *get_double_pole()/(-9.) /* = lc born tree */  );
		//Full-Leading color
		case 2: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* ( 23/6. -  33/6. /*beta_0*/ *(25/3./9. - get_double_pole()/(-9.)  /* = lc born tree */ ))  );
		}
	case 2:	return get_double_pole();
	}

};




}
