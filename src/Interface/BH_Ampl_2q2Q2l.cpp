/*
 * matrix elements for 2q2Q2l
 *
 *  Created on: Oct 21, 2008
 *
 */



#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include "BH_interface_impl.h"


#define _VERBOSE 0

// various switches: photonZW:  0 for photon
// 				1 for photon+Z
// 				2 for Z->neutrinos only
// 				3 for W only
// #define _PHOTON_ONLY 1     // 0 photon only; 1 for all other     replaced by a setting


// include nf_top pieces
#define Include_nf_top_Pieces 0
#define Include_axial_Pieces 0

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;
namespace BH {

partial_amplitude_cached* A_loop_2q_2Q_2l_6_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark,int photonZW, const ph_type e_ph_type,int color,QCDorder lo_or_nlo=nlo){
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



vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);


if(color==1){
//Leading color
// SCHEME-SHIFT
        PA->add_subtraction(pro,ind,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro,ind,r1,-1);
}
if(color==0){
//Full color
// SCHEME-SHIFT
        PA->add_subtraction(pro,ind,r4+r5,0);
// RENORMALIZATION
	PA->add_subtraction(pro,ind,(r1+r2+r3),-1);
}

	switch(h1.helicity()-h2.helicity()){
	case 0:{
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	process pro_1324=process(h1,hp3,hp2,h4,h5,h6);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);
	ind_1324.push_back(i6);

//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_2314=process(h1,h4,h2,h3,h5,h6);
	vector<int> ind_2314;
	ind_2314.push_back(i1);
	ind_2314.push_back(i4);
	ind_2314.push_back(i2);
	ind_2314.push_back(i3);
	ind_2314.push_back(i5);
	ind_2314.push_back(i6);

if(color==1){
	PA->add(pro,leading_color,ind,1,1);
}
if(color==0){
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro,leading_color,ind,-2,n_c*n_c);
	PA->add(pro_1324,leading_color,ind_1324,-2,n_c*n_c);

	//PA->add(pro_2314,sub_leading_color,ind_2314,r2,1);
	PA->add(pro_2314,sub_leading_color,ind_2314,1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
#if Include_nf_top_Pieces
	PA->add(pro,nf_top,ind,1,n_c);
#endif
}
	       };break;
	default: {

	// with this I flip quark<->anti-quark for particle_2 and particle_3
	ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	process pro_1324=process(h1,hp3,hp2,h4,h5,h6);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);
	ind_1324.push_back(i6);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_3214=process(h1,h4,h3,h2,h5,h6);
	vector<int> ind_3214;
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i5);
	ind_3214.push_back(i6);

if(color==1){
	PA->add(pro,leading_color,ind,1,1);
}
if(color==0){
	PA->add(pro,leading_color,ind,1,1);

	PA->add(pro,leading_color,ind,-2,n_c*n_c);
	PA->add(pro_1324,leading_color,ind_1324,-2,n_c*n_c);

	//PA->add(pro_3214,sub_leading_color,ind_3214,r2,1);
	PA->add(pro_3214,sub_leading_color,ind_3214,-1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
#if Include_nf_top_Pieces
	PA->add(pro,nf_top,ind,1,n_c);
#endif
}
		 };break;};

	return PA;
}


partial_amplitude_cached* A_loop_2q_2Q_2l_6_2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const ph_type e_ph_type, int color,QCDorder lo_or_nlo=nlo){
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


vector<ph_type> _ph_type;
_ph_type.push_back(h1);
_ph_type.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(-11,3), r2(2*n_f,3*n_c),r3(n_s,3*n_c);
        multi_precision_fraction r4(-2,3), r5(1,n_c*n_c);

if(color==1){
//Leading color
// SCHEME-SHIFT
        PA->add_subtraction(pro,ind,r4,0); // notice minus sign: compare 2.10 vs 2.11 sign infront of Kron. deltas
// RENORMALIZATION
	PA->add_subtraction(pro,ind,r1,-1); // notice minus sign: compare 2.10 vs 2.11 sign infront of Kron. deltas
}
if(color==0){
//Full color
// SCHEME-SHIFT
        PA->add_subtraction(pro,ind,(r4+r5),0); // notice minus sign: compare 2.10 vs 2.11 sign infront of Kron. deltas
// RENORMALIZATION
	PA->add_subtraction(pro,ind,(r1+r2+r3),-1); // notice minus sign: compare 2.10 vs 2.11 sign infront of Kron. deltas
}



	switch(h1.helicity()-h2.helicity()){
	case 0:{
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	process pro_1324=process(h1,hp3,hp2,h4,h5,h6);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);
	ind_1324.push_back(i6);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_2314=process(h1,h4,h2,h3,h5,h6);
	vector<int> ind_2314;
	ind_2314.push_back(i1);
	ind_2314.push_back(i4);
	ind_2314.push_back(i2);
	ind_2314.push_back(i3);
	ind_2314.push_back(i5);
	ind_2314.push_back(i6);

if(color==1){
	PA->add(pro_1324,leading_color,ind_1324,1,1);
}
if(color==0){
	PA->add(pro_1324,leading_color,ind_1324,1,1);
	PA->add(pro_1324,leading_color,ind_1324,1,n_c*n_c);
	PA->add(pro,leading_color,ind,1,n_c*n_c);

	PA->add(pro_2314,sub_leading_color,ind_2314,-1,n_c*n_c);

	PA->add(pro,nf,ind,-n_f,n_c);
#if Include_nf_top_Pieces
	PA->add(pro,nf_top,ind,-1,n_c);
#endif
}

	       };break;
	default:{

	// with this I flip quark<->anti-quark for particle_2 and particle_3
	ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	process pro_1324=process(h1,hp3,hp2,h4,h5,h6);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);
	ind_1324.push_back(i6);


	process pro_3214=process(h1,h4,h3,h2,h5,h6);
	vector<int> ind_3214;
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i5);
	ind_3214.push_back(i6);

if(color==1){
	PA->add(pro_1324,leading_color,ind_1324,1,1);
}
if(color==0){
	PA->add(pro_1324,leading_color,ind_1324,1,1);
	PA->add(pro_1324,leading_color,ind_1324,1,n_c*n_c);
	PA->add(pro,leading_color,ind,1,n_c*n_c);

	PA->add(pro_3214,sub_leading_color,ind_3214,1,n_c*n_c);

	PA->add(pro,nf,ind,-n_f,n_c);
#if Include_nf_top_Pieces
	PA->add(pro,nf_top,ind,-1,n_c);
#endif
}
		};break;
		};
return PA;
}


partial_amplitude_cached* A_loop_2q_2Q_2l_6_3_AX(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA = new partial_amplitude_cached(lo_or_nlo); ;

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

	switch(h1.helicity()-h2.helicity()){
	case 0:{
		process pro_1423=process(h1,h4,h2,h3,h5,h6);
		vector<int> ind_1423;
		ind_1423.push_back(i1);
		ind_1423.push_back(i4);
		ind_1423.push_back(i2);
		ind_1423.push_back(i3);
		ind_1423.push_back(i5);
		ind_1423.push_back(i6);

		int leading_vect_ax=2;
		//--------------------- propagator correction
			prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i5,i6,_ph_type);
			PA->set_prefactor(_prop_hel_fn);
		//-------------------------------------------

		PA->add(pro_1423,AX,ind_1423,1,1);

	} break;
	default:{
		process pro_1432=process(h1,h4,h3,h2,h5,h6);
		vector<int> ind_1432;
		ind_1432.push_back(i1);
		ind_1432.push_back(i4);
		ind_1432.push_back(i3);
		ind_1432.push_back(i2);
		ind_1432.push_back(i5);
		ind_1432.push_back(i6);

		int leading_vect_ax=2;
		//--------------------- propagator correction
			prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i5,i6,_ph_type);
			PA->set_prefactor(_prop_hel_fn);
		//-------------------------------------------

		PA->add(pro_1432,AX,ind_1432,-1,1);

	} break;
	}

return PA;
}


Squared_ME* A_loop_2q_2Q_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int photonZW, int color, int tree_color,const ph_type e_ph_type, int case4q,QCDorder lo_or_nlo=nlo){

	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);
	int i6=ind.at(5);


// assume input with identical flavor quarks
// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors
	ph_type h1=pro.p(1);
	ph_type h2=particle_ID(gluino,pro.p(2).helicity(),1,true);
	ph_type h3=particle_ID(gluino,pro.p(3).helicity(),1);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
	ph_type h6=pro.p(6);
// flipped case
	ph_type hf1=particle_ID(gluino,pro.p(1).helicity(),1);
	ph_type hf2=pro.p(2);
	ph_type hf3=pro.p(3);
	ph_type hf4=particle_ID(gluino,pro.p(4).helicity(),1,true);


	process pro_1234=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_1234;
	ind_1234.push_back(i1);
	ind_1234.push_back(i2);
	ind_1234.push_back(i3);
	ind_1234.push_back(i4);
	ind_1234.push_back(i5);
	ind_1234.push_back(i6);

	process pro_3412=process(hf3,hf4,hf1,hf2,h5,h6);
	vector<int> ind_3412;
	ind_3412.push_back(i3);
	ind_3412.push_back(i4);
	ind_3412.push_back(i1);
	ind_3412.push_back(i2);
	ind_3412.push_back(i5);
	ind_3412.push_back(i6);

// flavor exchange helicity labels
	ph_type hp1=pro.p(1);
	ph_type hp2=pro.p(2);
	ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),1);
	ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),1,true);
	ph_type hp5=pro.p(5);
	ph_type hp6=pro.p(6);
// flipped case
	ph_type hpf1=particle_ID(gluino,pro.p(1).helicity(),1);
	ph_type hpf2=particle_ID(gluino,pro.p(2).helicity(),1,true);
	ph_type hpf3=pro.p(3);
	ph_type hpf4=pro.p(4);

	process pro_1432=process(hp1,hp4,hp3,hp2,hp5,hp6);
	vector<int> ind_1432;
	ind_1432.push_back(i1);
	ind_1432.push_back(i4);
	ind_1432.push_back(i3);
	ind_1432.push_back(i2);
	ind_1432.push_back(i5);
	ind_1432.push_back(i6);

	process pro_3214=process(hpf3,hpf2,hpf1,hpf4,hp5,hp6);
	vector<int> ind_3214;
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i5);
	ind_3214.push_back(i6);

	Squared_ME* SM = new Squared_ME(lo_or_nlo);


vector<ph_type> _ph_type_1;
_ph_type_1.push_back(h1);
_ph_type_1.push_back(e_ph_type);

vector<ph_type> _ph_type_3;
_ph_type_3.push_back(h3);
_ph_type_3.push_back(e_ph_type);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn_up_1(0,photonZW,0,i5,i6,_ph_type_1);
	prop_hel_fn _prop_hel_fn_up_3(0,photonZW,0,i5,i6,_ph_type_3);
	prop_hel_fn _prop_hel_fn_down_1(1,photonZW,0,i5,i6,_ph_type_1);
	prop_hel_fn _prop_hel_fn_down_3(1,photonZW,0,i5,i6,_ph_type_3);
//-------------------------------------------

//------------------------------------------------------
// constants
	multi_precision_fraction nc(n_c,1);
	multi_precision_fraction nf(n_f,1);
	multi_precision_fraction n2m1(n_c*n_c-1,1);
	multi_precision_fraction _over_nc(1,n_c);

	multi_precision_fraction ext(4*2*nc*n2m1,1);
	multi_precision_fraction ext_tree(4*n2m1,1);
	
	multi_precision_fraction two(2,1);
	multi_precision_fraction Qu(2,3);
	multi_precision_fraction Qd(-1,3);

//------------------------------------------------------
//


	size_t T_1234_up;
	size_t T_3412_up;

	size_t L1_1234_up;
	size_t L1_3412_up;
	size_t L2_1234_up;
	size_t L2_3412_up;

	size_t T_1432_up;
	size_t T_3214_up;

	size_t L1_1432_up;
	size_t L1_3214_up;
	size_t L2_1432_up;
	size_t L2_3214_up;

	size_t T_1234_down;
	size_t T_3412_down;

	size_t L1_1234_down;
	size_t L1_3412_down;
	size_t L2_1234_down;
	size_t L2_3412_down;

	size_t T_1432_down;
	size_t T_3214_down;

	size_t L1_1432_down;
	size_t L1_3214_down;
	size_t L2_1432_down;
	size_t L2_3214_down;




	if(h1.helicity()+h4.helicity()==0){
		T_1234_up=SM->add(new CTree_with_prefactor(pro_1234,ind_1234,_prop_hel_fn_up_1));
		T_3412_up=SM->add(new CTree_with_prefactor(pro_3412,ind_3412,_prop_hel_fn_up_3));

		L1_1234_up=SM->add(A_loop_2q_2Q_2l_6_1(pro_1234,ind_1234,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L1_3412_up=SM->add(A_loop_2q_2Q_2l_6_1(pro_3412,ind_3412,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L2_1234_up=SM->add(A_loop_2q_2Q_2l_6_2(pro_1234,ind_1234,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L2_3412_up=SM->add(A_loop_2q_2Q_2l_6_2(pro_3412,ind_3412,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_1432_up=SM->add(new CTree_with_prefactor(pro_1432,ind_1432,_prop_hel_fn_up_1));
		T_3214_up=SM->add(new CTree_with_prefactor(pro_3214,ind_3214,_prop_hel_fn_up_3));

		L1_1432_up=SM->add(A_loop_2q_2Q_2l_6_1(pro_1432,ind_1432,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L1_3214_up=SM->add(A_loop_2q_2Q_2l_6_1(pro_3214,ind_3214,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L2_1432_up=SM->add(A_loop_2q_2Q_2l_6_2(pro_1432,ind_1432,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
		L2_3214_up=SM->add(A_loop_2q_2Q_2l_6_2(pro_3214,ind_3214,n_s,n_f,n_c,0,photonZW,e_ph_type,color,lo_or_nlo));
	};


	if(h1.helicity()+h4.helicity()==0){
		T_1234_down=SM->add(new CTree_with_prefactor(pro_1234,ind_1234,_prop_hel_fn_down_1));
		T_3412_down=SM->add(new CTree_with_prefactor(pro_3412,ind_3412,_prop_hel_fn_down_3));

		L1_1234_down=SM->add(A_loop_2q_2Q_2l_6_1(pro_1234,ind_1234,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L1_3412_down=SM->add(A_loop_2q_2Q_2l_6_1(pro_3412,ind_3412,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L2_1234_down=SM->add(A_loop_2q_2Q_2l_6_2(pro_1234,ind_1234,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L2_3412_down=SM->add(A_loop_2q_2Q_2l_6_2(pro_3412,ind_3412,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_1432_down=SM->add(new CTree_with_prefactor(pro_1432,ind_1432,_prop_hel_fn_down_1));
		T_3214_down=SM->add(new CTree_with_prefactor(pro_3214,ind_3214,_prop_hel_fn_down_3));

		L1_1432_down=SM->add(A_loop_2q_2Q_2l_6_1(pro_1432,ind_1432,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L1_3214_down=SM->add(A_loop_2q_2Q_2l_6_1(pro_3214,ind_3214,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L2_1432_down=SM->add(A_loop_2q_2Q_2l_6_2(pro_1432,ind_1432,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
		L2_3214_down=SM->add(A_loop_2q_2Q_2l_6_2(pro_3214,ind_3214,n_s,n_f,n_c,1,photonZW,e_ph_type,color,lo_or_nlo));
	};
//--------------------------------------------------------
// mulitplicities of processes

	multi_precision_fraction uu_multiplicity(1,1); //v
	multi_precision_fraction dd_multiplicity(1,1); //v
	multi_precision_fraction uup_multiplicity(1,1); //v
	multi_precision_fraction ddp_multiplicity(1,1); //v
	multi_precision_fraction ud_multiplicity(1,1); //v
	multi_precision_fraction du_multiplicity(1,1); //v


// W cases
if(photonZW==3){
	//--------------------------------------------------------
	// W case: case4q==0 -> identical quarks
	// Notice that we use _up to built W amps!
	//--------------------------------------------------------
	if(case4q==0){

		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=uu*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_1234_up,T_1234_up,uu));
		};
		if(h3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_3214_up,T_3214_up,uu));
		};
if(tree_color==0){
//Full color
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
			SM->add_tree(cached_cross_term_md(T_1234_up,T_3214_up,uu_E));
		};
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
			SM->add_tree(cached_cross_term_md(T_3214_up,T_1234_up,uu_E));
		};
}
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_1234_up,T_1234_up,uu_virtual));
		};

		if(h3.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_3214_up,T_3214_up,uu_virtual));
		};

if(color==0){
//Full color
		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
			SM->add_loop(cached_cross_term_md(L2_1234_up,T_3214_up,uu_E_virtual));
		};

		if(h1.helicity()==-1&&h3.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
			SM->add_loop(cached_cross_term_md(L2_3214_up,T_1234_up,uu_E_virtual));
		}
	}
	}
	//--------------------------------------------------------
	// W case: case4q==1 -> identical anti-quarks
	// Notice that we use _up to built W amps!
	//--------------------------------------------------------
	if(case4q==1){

		multi_precision_fraction uu=uu_multiplicity*ext_tree;
		multi_precision_fraction uu_E=uu*_over_nc;

		if(h1.helicity()==-1&& h1.helicity()+h4.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_1234_up,T_1234_up,uu));
		};
		if(h1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_tree(cached_cross_term_md(T_1432_up,T_1432_up,uu));
		};
if(tree_color==0){
//Full color
		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
			SM->add_tree(cached_cross_term_md(T_1234_up,T_1432_up,uu_E));
		};
		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
			SM->add_tree(cached_cross_term_md(T_1432_up,T_1234_up,uu_E));
		};
	}	
		multi_precision_fraction uu_virtual=uu_multiplicity*ext;
		multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;

		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_1234_up,T_1234_up,uu_virtual));
		};

		if(h1.helicity()==-1&&hp1.helicity()+hp2.helicity()==0){
			SM->add_loop(cached_cross_term_md(L1_1432_up,T_1432_up,uu_virtual));
		};

if(color==0){
//Full color
		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
			SM->add_loop(cached_cross_term_md(L2_1234_up,T_1432_up,uu_E_virtual));
		};

		if(h1.helicity()==-1&&(h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
			SM->add_loop(cached_cross_term_md(L2_1432_up,T_1234_up,uu_E_virtual));
		}
	}
	}

	//--------------------------------------------------------
	// W case: case4q==2 -> no identical quarks
	// Notice that we use _up to built W amps!
	if(case4q==2){
		if(h1.helicity()==-1&&h1.helicity()+h4.helicity()==0){

			multi_precision_fraction ud=ud_multiplicity*ext_tree;
			multi_precision_fraction ud_virtual=ud_multiplicity*ext;

			SM->add_tree(cached_cross_term_md(T_1234_up  ,T_1234_up  ,ext_tree));

			SM->add_loop(cached_cross_term_md(L1_1234_up  ,T_1234_up  ,ext));
		};
	}

	return SM;
}

// Axial amplitudes
	size_t L3_1234;
	size_t L3_1432;
	if(h1.helicity()+h4.helicity()==0){
		L3_1234=SM->add(A_loop_2q_2Q_2l_6_3_AX(pro_1234,ind_1234,n_s,n_f,n_c,0,photonZW,_ph_type_1,lo_or_nlo));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		L3_1432=SM->add(A_loop_2q_2Q_2l_6_3_AX(pro_1432,ind_1432,n_s,n_f,n_c,0,photonZW,_ph_type_1,lo_or_nlo));
	};


// Now Z/gamma cases
//--------------------------------------------------------
// (u,d)
if(case4q==2 || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){

		multi_precision_fraction ud=ud_multiplicity*ext_tree;
		multi_precision_fraction ud_virtual=ud_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234_up  ,T_1234_up  ,ud));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_3412_down,ud));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_1234_up  ,ud));
		SM->add_tree(cached_cross_term_md(T_1234_up  ,T_3412_down,ud));

		SM->add_loop(cached_cross_term_md(L1_1234_up  ,T_1234_up  ,ud_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_3412_down,ud_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_1234_up  ,ud_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_up  ,T_3412_down,ud_virtual));
		// Axial piece contribution

#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_up  ,ud*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_down,ud*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};
}
//--------------------------------------------------------
// (d,u)
if(case4q==3 || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){

		multi_precision_fraction du=du_multiplicity*ext_tree;
		multi_precision_fraction du_virtual=du_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234_down,T_1234_down,du));
		SM->add_tree(cached_cross_term_md(T_3412_up  ,T_3412_up  ,du));
		SM->add_tree(cached_cross_term_md(T_3412_up  ,T_1234_down,du));
		SM->add_tree(cached_cross_term_md(T_1234_down,T_3412_up  ,du));

		SM->add_loop(cached_cross_term_md(L1_1234_down,T_1234_down,du_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up  ,T_3412_up  ,du_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up  ,T_1234_down,du_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_down,T_3412_up  ,du_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_up  ,du*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_down,du*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};
}
//--------------------------------------------------------
// (u,u')
if(case4q==1 || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){
		multi_precision_fraction uup=uup_multiplicity*ext_tree;
		multi_precision_fraction uup_virtual=uup_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234_up,T_1234_up,uup));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_3412_up,uup));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_1234_up,uup));
		SM->add_tree(cached_cross_term_md(T_1234_up,T_3412_up,uup));

		SM->add_loop(cached_cross_term_md(L1_1234_up,T_1234_up,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up,T_3412_up,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up,T_1234_up,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_up,T_3412_up,uup_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_up,uup*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_up,uup*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};
}
//--------------------------------------------------------
// (d,d')
if(case4q==4 || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){
		multi_precision_fraction ddp=ddp_multiplicity*ext_tree;
		multi_precision_fraction ddp_virtual=ddp_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234_down,T_1234_down,ddp));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_3412_down,ddp));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_1234_down,ddp));
		SM->add_tree(cached_cross_term_md(T_1234_down,T_3412_down,ddp));

		SM->add_loop(cached_cross_term_md(L1_1234_down,T_1234_down,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_3412_down,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_1234_down,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_down,T_3412_down,ddp_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_down,ddp*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_down,ddp*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};
}
//--------------------------------------------------------
// (u,u)
//--------------------------------------------------------
if(case4q==0 || case4q==-1){

	multi_precision_fraction uu=uu_multiplicity*ext_tree;
	multi_precision_fraction uu_E=uu*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1234_up,T_1234_up,uu));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_3412_up,uu));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_1234_up,uu));
		SM->add_tree(cached_cross_term_md(T_1234_up,T_3412_up,uu));
	};
	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1432_up,T_1432_up,uu));
		SM->add_tree(cached_cross_term_md(T_3214_up,T_3214_up,uu));
		SM->add_tree(cached_cross_term_md(T_3214_up,T_1432_up,uu));
		SM->add_tree(cached_cross_term_md(T_1432_up,T_3214_up,uu));
	};
if(tree_color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1234_up,T_1432_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_1234_up,T_3214_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_1432_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_3412_up,T_3214_up,uu_E));
	};
	if((h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1432_up,T_1234_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_3214_up,T_1234_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_1432_up,T_3412_up,uu_E));
		SM->add_tree(cached_cross_term_md(T_3214_up,T_3412_up,uu_E));
	};
	}
	multi_precision_fraction uu_virtual=uu_multiplicity*ext;
	multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1234_up,T_1234_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up,T_3412_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_up,T_1234_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_up,T_3412_up,uu_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_up,uu*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_up,uu*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1432_up,T_1432_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214_up,T_3214_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214_up,T_1432_up,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_1432_up,T_3214_up,uu_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1432,T_1432_up,uu*two));
			SM->add_loop(cached_cross_term_md(L3_1432,T_3214_up,uu*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

if(color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1234_up,T_1432_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1234_up,T_3214_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412_up,T_1432_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412_up,T_3214_up,uu_E_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1432_up,uu_E*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3214_up,uu_E*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1432_up,T_1234_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214_up,T_1234_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1432_up,T_3412_up,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214_up,T_3412_up,uu_E_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1432,T_1234_up,uu_E*two));
			SM->add_loop(cached_cross_term_md(L3_1432,T_3412_up,uu_E*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	}
	}
}
// (d,d)
//--------------------------------------------------------
if(case4q==5 || case4q==-1){

	multi_precision_fraction dd=dd_multiplicity*ext_tree;
	multi_precision_fraction dd_E=dd*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1234_down,T_1234_down,dd));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_3412_down,dd));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_1234_down,dd));
		SM->add_tree(cached_cross_term_md(T_1234_down,T_3412_down,dd));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1432_down,T_1432_down,dd));
		SM->add_tree(cached_cross_term_md(T_3214_down,T_3214_down,dd));
		SM->add_tree(cached_cross_term_md(T_3214_down,T_1432_down,dd));
		SM->add_tree(cached_cross_term_md(T_1432_down,T_3214_down,dd));
	};

if(tree_color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1234_down,T_1432_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_1234_down,T_3214_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_1432_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_3412_down,T_3214_down,dd_E));
	};
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1432_down,T_1234_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_3214_down,T_1234_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_1432_down,T_3412_down,dd_E));
		SM->add_tree(cached_cross_term_md(T_3214_down,T_3412_down,dd_E));
	};
	}

	multi_precision_fraction dd_virtual=dd_multiplicity*ext;
	multi_precision_fraction dd_E_virtual=-dd_virtual*_over_nc;


	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1234_down,T_1234_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_3412_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412_down,T_1234_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234_down,T_3412_down,dd_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1234_down,dd*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3412_down,dd*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1432_down,T_1432_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214_down,T_3214_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214_down,T_1432_down,dd_virtual));
		SM->add_loop(cached_cross_term_md(L1_1432_down,T_3214_down,dd_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1432,T_1432_down,dd*two));
			SM->add_loop(cached_cross_term_md(L3_1432,T_3214_down,dd*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

if(color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1234_down,T_1432_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1234_down,T_3214_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412_down,T_1432_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412_down,T_3214_down,dd_E_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1234,T_1432_down,dd_E*two));
			SM->add_loop(cached_cross_term_md(L3_1234,T_3214_down,dd_E*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};

	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1432_down,T_1234_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214_down,T_1234_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1432_down,T_3412_down,dd_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214_down,T_3412_down,dd_E_virtual));
		// Axial piece contribution
#if Include_axial_Pieces
		if(photonZW!=0){
			SM->add_loop(cached_cross_term_md(L3_1432,T_1234_down,dd_E*two));
			SM->add_loop(cached_cross_term_md(L3_1432,T_3412_down,dd_E*two));
		}
#else
//_WARNING("Axial parts disabled");
#endif
	};
	}
}
	return SM;
}

/*
Virtual_SME* vsme_2q2Q2l(int ns,int nf,int nc,int photonZW,int case4q, int color, int tree_color){
	Virtual_SME* VSM=new Virtual_SME();

	vector<int> ind, ind_123465;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);
	ind.push_back(6);

	ind_123465.push_back(1);
	ind_123465.push_back(2);
	ind_123465.push_back(3);
	ind_123465.push_back(4);
	ind_123465.push_back(6);
	ind_123465.push_back(5);


	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbp,qm,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));
if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qp,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));
}
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qm,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));

	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbm,qp,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qm,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qp,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q));

if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbp,qm,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qp,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qm,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));

	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbm,qp,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qm,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qp,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q));
}
return VSM;
}
*/
Virtual_SME* vsme_2q2Q2l(std::vector<int> indext,int ns,int nf,int nc,int photonZW,int case4q, int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME();

	vector<int> ind, ind_123465;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);

	ind_123465.push_back(indext[0]);
	ind_123465.push_back(indext[1]);
	ind_123465.push_back(indext[2]);
	ind_123465.push_back(indext[3]);
	ind_123465.push_back(indext[5]);
	ind_123465.push_back(indext[4]);


	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbp,qm,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qp,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qm,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbm,qp,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qm,qbp,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qp,qbm,lm,lbp),ind,ns,nf,nc,photonZW,color,tree_color,lm,case4q,lo_or_nlo));

if(photonZW!=3){
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbp,qm,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qp,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qp,qbm,qm,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbm,qp,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qm,qbp,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2l_M2(process(qm,qbp,qp,qbm,lm,lbp),ind_123465,ns,nf,nc,photonZW,color,tree_color,lp,case4q,lo_or_nlo));
}
return VSM;
}


BH_Ampl_2q2Q2l::BH_Ampl_2q2Q2l(int photonZW,int case4q, int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q2l(mom_assignment,0,
							settings::BH_interface_settings::s_nf,
							settings::BH_interface_settings::s_nc,
							settings::BH_interface_settings::s_photon_only*photonZW,
							case4q,
							color,
							tree_color,
							lo_or_nlo),
			bhi,  // parent
			6,    // NbrExtParticles
			2,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			4,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
{}

}
