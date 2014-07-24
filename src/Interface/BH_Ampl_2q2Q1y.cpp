/* 
 * matrix elements for 2q2Q1y processes
 *
 *  Created on: Oct 21, 2008
 *
 */

///BROKEN DIST QUARKS

#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include "BH_interface_impl.h"

#define _VERBOSE 0


using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

partial_amplitude_cached* A_loop_2q_2Q_1y_6_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

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
	process pro_1324=process(h1,h3,h2,h4,h5);
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	//ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	//ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	//process pro_1324=process(h1,hp3,hp2,h4,h5);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_2314=process(h1,h4,h2,h3,h5);
	vector<int> ind_2314;
	ind_2314.push_back(i1);
	ind_2314.push_back(i4);
	ind_2314.push_back(i2);
	ind_2314.push_back(i3);
	ind_2314.push_back(i5);

if(color==1){
//Leading color
	PA->add(pro,leading_color,ind,1,1);
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro,leading_color,ind,-2,n_c*n_c);
	PA->add(pro_1324,leading_color,ind_1324,-2,n_c*n_c);

	//PA->add(pro_2314,sub_leading_color,ind_2314,r2,1);
	PA->add(pro_2314,sub_leading_color,ind_2314,1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
}
	       };break;
	default: {

	process pro_1324=process(h1,h3,h2,h4,h5);
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	//ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	//ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	//process pro_1324=process(h1,hp3,hp2,h4,h5);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_3214=process(h1,h4,h3,h2,h5);
	vector<int> ind_3214;
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i5);

if(color==1){
//Leading color
	PA->add(pro,leading_color,ind,1,1);
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);

	PA->add(pro,leading_color,ind,-2,n_c*n_c);
	PA->add(pro_1324,leading_color,ind_1324,-2,n_c*n_c);

	//PA->add(pro_3214,sub_leading_color,ind_3214,r2,1);
	PA->add(pro_3214,sub_leading_color,ind_3214,-1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
}
		 };break;};

	return PA;
}


partial_amplitude_cached* A_loop_2q_2Q_1y_6_2(process pro, vector<int> & ind, int n_s, int n_f, int n_c,int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

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
	process pro_1324=process(h1,h3,h2,h4,h5);
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	//ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	//ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	//process pro_1324=process(h1,hp3,hp2,h4,h5);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_2314=process(h1,h4,h2,h3,h5);
	vector<int> ind_2314;
	ind_2314.push_back(i1);
	ind_2314.push_back(i4);
	ind_2314.push_back(i2);
	ind_2314.push_back(i3);
	ind_2314.push_back(i5);

if(color==1){
//Leading color
	PA->add(pro_1324,leading_color,ind_1324,1,1);
}
if(color==0){
//Full color
	PA->add(pro_1324,leading_color,ind_1324,1,1);
	PA->add(pro_1324,leading_color,ind_1324,1,n_c*n_c);
	PA->add(pro,leading_color,ind,1,n_c*n_c);

	PA->add(pro_2314,sub_leading_color,ind_2314,-1,n_c*n_c);

	PA->add(pro,nf,ind,-n_f,n_c);
}

	       };break;
	default:{

	process pro_1324=process(h1,h3,h2,h4,h5);
	// with this I flip quark<->anti-quark for particle_2 and particle_3
	//ph_type hp2=particle_ID(gluino,h2.helicity(),1);
	//ph_type hp3=particle_ID(gluino,h3.helicity(),1,true);
	//process pro_1324=process(h1,hp3,hp2,h4,h5);
	vector<int> ind_1324;
	ind_1324.push_back(i1);
	ind_1324.push_back(i3);
	ind_1324.push_back(i2);
	ind_1324.push_back(i4);
	ind_1324.push_back(i5);


//NOTICE CHANGE OF NOTATION HERE compared to eqn. 2.15
	process pro_3214=process(h1,h4,h3,h2,h5);
	vector<int> ind_3214;
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i5);


if(color==1){
//Leading color
	PA->add(pro_1324,leading_color,ind_1324,1,1);
}
if(color==0){
//Full color
	PA->add(pro_1324,leading_color,ind_1324,1,1);
	PA->add(pro_1324,leading_color,ind_1324,1,n_c*n_c);
	PA->add(pro,leading_color,ind,1,n_c*n_c);

	PA->add(pro_3214,sub_leading_color,ind_3214,1,n_c*n_c);

	PA->add(pro,nf,ind,-n_f,n_c);
}
		};break;
		};
return PA;
}



Squared_ME* A_loop_2q_2Q_1y_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c,int case4q,bool up_down_quark,int color,int tree_color,QCDorder lo_or_nlo=nlo){

	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);
	int i5=ind.at(4);


// assume input with identical flavor quarks
// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors
	ph_type h1=pro.p(1);
//	ph_type h2=particle_ID(quark,pro.p(2).helicity(),2,true);
//	ph_type h3=particle_ID(quark,pro.p(3).helicity(),2);
	ph_type h2=particle_ID(gluino,pro.p(2).helicity(),1,true);
	ph_type h3=particle_ID(gluino,pro.p(3).helicity(),1);
	ph_type h4=pro.p(4);
	ph_type h5=pro.p(5);
// flipped case
	ph_type hf1=particle_ID(gluino,pro.p(1).helicity(),1);
	ph_type hf2=pro.p(2);
	ph_type hf3=pro.p(3);
	ph_type hf4=particle_ID(gluino,pro.p(4).helicity(),1,true);


	process pro_1234=process(h1,h2,h3,h4,h5);
	vector<int> ind_1234;
	ind_1234.push_back(i1);
	ind_1234.push_back(i2);
	ind_1234.push_back(i3);
	ind_1234.push_back(i4);
	ind_1234.push_back(i5);

//	process pro_3412=process(h3,h4,h1,h2,h5);
	process pro_3412=process(hf3,hf4,hf1,hf2,h5);
	vector<int> ind_3412;
	ind_3412.push_back(i3);
	ind_3412.push_back(i4);
	ind_3412.push_back(i1);
	ind_3412.push_back(i2);
	ind_3412.push_back(i5);

// flavor exchange helicity labels
	ph_type hp1=pro.p(1);
	ph_type hp2=pro.p(2);
//	ph_type hp3=particle_ID(quark,pro.p(3).helicity(),2);
//	ph_type hp4=particle_ID(quark,pro.p(4).helicity(),2,true);
	ph_type hp3=particle_ID(gluino,pro.p(3).helicity(),1);
	ph_type hp4=particle_ID(gluino,pro.p(4).helicity(),1,true);
	ph_type hp5=pro.p(5);
// flipped case
	ph_type hpf1=particle_ID(gluino,pro.p(1).helicity(),1);
	ph_type hpf2=particle_ID(gluino,pro.p(2).helicity(),1,true);
	ph_type hpf3=pro.p(3);
	ph_type hpf4=pro.p(4);

	process pro_1432=process(hpf3,hpf2,hpf1,hpf4,hp5);
	vector<int> ind_1432;
	ind_1432.push_back(i3);
	ind_1432.push_back(i2);
	ind_1432.push_back(i1);
	ind_1432.push_back(i4);
	ind_1432.push_back(i5);

	process pro_3214=process(hp1,hp4,hp3,hp2,hp5);
	vector<int> ind_3214;
	ind_3214.push_back(i1);
	ind_3214.push_back(i4);
	ind_3214.push_back(i3);
	ind_3214.push_back(i2);
	ind_3214.push_back(i5);

	Squared_ME* SM=new Squared_ME(lo_or_nlo);

//------------------------------------------------------
// constants
	multi_precision_fraction nc(n_c,1);
	multi_precision_fraction nf(n_f,1);
	multi_precision_fraction n2m1(n_c*n_c-1,1);
	multi_precision_fraction _over_nc(1,n_c);

	//multi_precision_fraction ext(4*2*nc*n2m1,1);
	multi_precision_fraction ext(4*nc*n2m1,1);
	//multi_precision_fraction ext_tree(4*n2m1,1);
	multi_precision_fraction ext_tree(2*n2m1,1);
	multi_precision_fraction Qu(2,3);
	multi_precision_fraction Qd(-1,3);

//------------------------------------------------------
//


	size_t T_1234;
	size_t T_3412;

	size_t L1_1234;
	size_t L1_3412;
	size_t L2_1234;
	size_t L2_3412;

	size_t T_1432;
	size_t T_3214;

	size_t L1_1432;
	size_t L1_3214;
	size_t L2_1432;
	size_t L2_3214;


	if(h1.helicity()+h4.helicity()==0){
		T_1234=SM->add(new CTree_with_prefactor(pro_1234,ind_1234));
		T_3412=SM->add(new CTree_with_prefactor(pro_3412,ind_3412));

		L1_1234=SM->add(A_loop_2q_2Q_1y_6_1(pro_1234,ind_1234,n_s,n_f,n_c,color,lo_or_nlo));
		L1_3412=SM->add(A_loop_2q_2Q_1y_6_1(pro_3412,ind_3412,n_s,n_f,n_c,color,lo_or_nlo));
		L2_1234=SM->add(A_loop_2q_2Q_1y_6_2(pro_1234,ind_1234,n_s,n_f,n_c,color,lo_or_nlo));
		L2_3412=SM->add(A_loop_2q_2Q_1y_6_2(pro_3412,ind_3412,n_s,n_f,n_c,color,lo_or_nlo));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		T_1432=SM->add(new CTree_with_prefactor(pro_1432,ind_1432));
		T_3214=SM->add(new CTree_with_prefactor(pro_3214,ind_3214));
		L1_1432=SM->add(A_loop_2q_2Q_1y_6_1(pro_1432,ind_1432,n_s,n_f,n_c,color,lo_or_nlo));
		L1_3214=SM->add(A_loop_2q_2Q_1y_6_1(pro_3214,ind_3214,n_s,n_f,n_c,color,lo_or_nlo));
		L2_1432=SM->add(A_loop_2q_2Q_1y_6_2(pro_1432,ind_1432,n_s,n_f,n_c,color,lo_or_nlo));
		L2_3214=SM->add(A_loop_2q_2Q_1y_6_2(pro_3214,ind_3214,n_s,n_f,n_c,color,lo_or_nlo));
	};

//--------------------------------------------------------
// mulitplicities of processes
/*
	multi_precision_fraction uu_multiplicity(2,2); //v
	multi_precision_fraction dd_multiplicity(3,2); //v
	multi_precision_fraction uup_multiplicity(2,1); //v
	multi_precision_fraction ddp_multiplicity(6,1); //v
	multi_precision_fraction ud_multiplicity(6,1); //v
	multi_precision_fraction du_multiplicity(6,1); //v
*/
	multi_precision_fraction uu_multiplicity(1,1); //v
	multi_precision_fraction dd_multiplicity(1,1); //v
	multi_precision_fraction uup_multiplicity(1,1); //v
	multi_precision_fraction ddp_multiplicity(1,1); //v
	multi_precision_fraction ud_multiplicity(1,1); //v
	multi_precision_fraction du_multiplicity(1,1); //v

//--------------------------------------------------------
// (u,d)
if((up_down_quark==0&&case4q==2) || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){

		multi_precision_fraction ud=ud_multiplicity*ext_tree;
		multi_precision_fraction ud_virtual=ud_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234,T_1234,ud*Qu*Qu));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,ud*Qd*Qd));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,ud*Qu*Qd));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,ud*Qu*Qd));

		SM->add_loop(cached_cross_term_md(L1_1234,T_1234,ud_virtual*Qu*Qu));
		SM->add_loop(cached_cross_term_md(L1_3412,T_3412,ud_virtual*Qd*Qd));
		SM->add_loop(cached_cross_term_md(L1_3412,T_1234,ud_virtual*Qu*Qd));
		SM->add_loop(cached_cross_term_md(L1_1234,T_3412,ud_virtual*Qu*Qd));
	};
}
//--------------------------------------------------------
// (d,u)
if((up_down_quark==1&&case4q==2) || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){

		multi_precision_fraction du=du_multiplicity*ext_tree;
		multi_precision_fraction du_virtual=du_multiplicity*ext;

		SM->add_tree(cached_cross_term_md(T_1234,T_1234,du*Qd*Qd));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,du*Qu*Qu));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,du*Qu*Qd));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,du*Qu*Qd));

		SM->add_loop(cached_cross_term_md(L1_1234,T_1234,du_virtual*Qd*Qd));
		SM->add_loop(cached_cross_term_md(L1_3412,T_3412,du_virtual*Qu*Qu));
		SM->add_loop(cached_cross_term_md(L1_3412,T_1234,du_virtual*Qd*Qu));
		SM->add_loop(cached_cross_term_md(L1_1234,T_3412,du_virtual*Qd*Qu));
	};
}
//--------------------------------------------------------
// (u,u')
if((up_down_quark==0 && case4q==1) || case4q==-1){
	if(h1.helicity()+h4.helicity()==0){
		multi_precision_fraction uup=uup_multiplicity*Qu*Qu*ext_tree;
		multi_precision_fraction uup_virtual=uup_multiplicity*Qu*Qu*ext;

		SM->add_tree(cached_cross_term_md(T_1234,T_1234,uup));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,uup));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,uup));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,uup));

		SM->add_loop(cached_cross_term_md(L1_1234,T_1234,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_3412,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_1234,uup_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234,T_3412,uup_virtual));
	};
}
//--------------------------------------------------------
// (d,d')
if((up_down_quark==1 && case4q==1)||case4q==-1){
	if(h1.helicity()+h4.helicity()==0){
		multi_precision_fraction ddp=ddp_multiplicity*Qd*Qd*ext_tree;
		multi_precision_fraction ddp_virtual=ddp_multiplicity*Qd*Qd*ext;

		SM->add_tree(cached_cross_term_md(T_1234,T_1234,ddp));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,ddp));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,ddp));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,ddp));

		SM->add_loop(cached_cross_term_md(L1_1234,T_1234,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_3412,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_1234,ddp_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234,T_3412,ddp_virtual));
	};
}
//--------------------------------------------------------
// (u,u)
//--------------------------------------------------------
if((up_down_quark==0&&case4q==0) || case4q==-1){

	multi_precision_fraction uu=uu_multiplicity*ext_tree*Qu*Qu;
	multi_precision_fraction uu_E=uu*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1234,T_1234,uu));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,uu));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,uu));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,uu));
	};
	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1432,T_1432,uu));
		SM->add_tree(cached_cross_term_md(T_3214,T_3214,uu));
		SM->add_tree(cached_cross_term_md(T_3214,T_1432,uu));
		SM->add_tree(cached_cross_term_md(T_1432,T_3214,uu));
	};

if(tree_color==0){
//Full color

	if((h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1234,T_1432,uu_E));
		SM->add_tree(cached_cross_term_md(T_1234,T_3214,uu_E));
		SM->add_tree(cached_cross_term_md(T_3412,T_1432,uu_E));
		SM->add_tree(cached_cross_term_md(T_3412,T_3214,uu_E));
	};
	if((h1.helicity()+h4.helicity()==0)&&(hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1432,T_1234,uu_E));
		SM->add_tree(cached_cross_term_md(T_3214,T_1234,uu_E));
		SM->add_tree(cached_cross_term_md(T_1432,T_3412,uu_E));
		SM->add_tree(cached_cross_term_md(T_3214,T_3412,uu_E));
	};
}

	multi_precision_fraction uu_virtual=uu_multiplicity*ext*Qu*Qu;
	multi_precision_fraction uu_E_virtual=-uu_virtual*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1234,T_1234,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_3412,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3412,T_1234,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_1234,T_3412,uu_virtual));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_loop(cached_cross_term_md(L1_1432,T_1432,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214,T_3214,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_3214,T_1432,uu_virtual));
		SM->add_loop(cached_cross_term_md(L1_1432,T_3214,uu_virtual));
	};


if(color==0){
//Full color

	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1234,T_1432,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1234,T_3214,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412,T_1432,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3412,T_3214,uu_E_virtual));
	};

	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_loop(cached_cross_term_md(L2_1432,T_1234,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214,T_1234,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_1432,T_3412,uu_E_virtual));
		SM->add_loop(cached_cross_term_md(L2_3214,T_3412,uu_E_virtual));
	}
}

}
// (d,d)
//--------------------------------------------------------
if((up_down_quark==1&&case4q==0)||case4q==-1){

	multi_precision_fraction dd=dd_multiplicity*ext_tree*(Qd*Qd);
	multi_precision_fraction dd_E=dd*_over_nc;

	if(h1.helicity()+h4.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1234,T_1234,dd));
		SM->add_tree(cached_cross_term_md(T_3412,T_3412,dd));
		SM->add_tree(cached_cross_term_md(T_3412,T_1234,dd));
		SM->add_tree(cached_cross_term_md(T_1234,T_3412,dd));
	};

	if(hp1.helicity()+hp2.helicity()==0){
		SM->add_tree(cached_cross_term_md(T_1432,T_1432,dd));
		SM->add_tree(cached_cross_term_md(T_3214,T_3214,dd));
		SM->add_tree(cached_cross_term_md(T_3214,T_1432,dd));
		SM->add_tree(cached_cross_term_md(T_1432,T_3214,dd));
	};

if(tree_color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1234,T_1432,dd_E));
		SM->add_tree(cached_cross_term_md(T_1234,T_3214,dd_E));
		SM->add_tree(cached_cross_term_md(T_3412,T_1432,dd_E));
		SM->add_tree(cached_cross_term_md(T_3412,T_3214,dd_E));
	};
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
		SM->add_tree(cached_cross_term_md(T_1432,T_1234,dd_E));
		SM->add_tree(cached_cross_term_md(T_3214,T_1234,dd_E));
		SM->add_tree(cached_cross_term_md(T_1432,T_3412,dd_E));
		SM->add_tree(cached_cross_term_md(T_3214,T_3412,dd_E));
	};
}

	multi_precision_fraction dd_virtual=dd_multiplicity*ext*(Qd*Qd);
	multi_precision_fraction dd_E_virtual=-dd_virtual*_over_nc;


	if(h1.helicity()+h4.helicity()==0){
	SM->add_loop(cached_cross_term_md(L1_1234,T_1234,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_3412,T_3412,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_3412,T_1234,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_1234,T_3412,dd_virtual));
	};

	if(hp1.helicity()+hp2.helicity()==0){
	SM->add_loop(cached_cross_term_md(L1_1432,T_1432,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_3214,T_3214,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_3214,T_1432,dd_virtual));
	SM->add_loop(cached_cross_term_md(L1_1432,T_3214,dd_virtual));
	};

if(color==0){
//Full color
	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
	SM->add_loop(cached_cross_term_md(L2_1234,T_1432,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_1234,T_3214,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_3412,T_1432,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_3412,T_3214,dd_E_virtual));
	};

	if((h1.helicity()+h4.helicity()==0) && (hp1.helicity()+hp2.helicity()==0)){
	SM->add_loop(cached_cross_term_md(L2_1432,T_1234,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_3214,T_1234,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_1432,T_3412,dd_E_virtual));
	SM->add_loop(cached_cross_term_md(L2_3214,T_3412,dd_E_virtual));
	};
}
}
	return SM;
}

/*
Virtual_SME* vsme_2q2Q1y(int ns,int nf,int nc,int case4q, bool up_down_quark,int color, int tree_color){
	Virtual_SME* VSM=new Virtual_SME();

	vector<int> ind, ind_4321;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);

	ind_4321.push_back(4);
	ind_4321.push_back(3);
	ind_4321.push_back(2);
	ind_4321.push_back(1);
	ind_4321.push_back(5);


	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbp,qm,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qp,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qm,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));

	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbm,qp,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qm,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qp,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color));

return VSM;
}
*/

Virtual_SME* vsme_2q2Q1y(std::vector<int> indext,int ns,int nf,int nc,int case4q, bool up_down_quark,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME();

	vector<int> ind, ind_4321;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);

	ind_4321.push_back(indext[3]);
	ind_4321.push_back(indext[2]);
	ind_4321.push_back(indext[1]);
	ind_4321.push_back(indext[0]);
	ind_4321.push_back(indext[4]);


	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbp,qm,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qp,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qm,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbm,qp,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qm,qbp,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qp,qbm,yp),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	
    VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbp,qm,qbm,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qp,qbm,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qp,qbm,qm,qbp,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));

	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbm,qp,qbp,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qm,qbp,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_1y_M2(process(qm,qbp,qp,qbm,ym),ind,ns,nf,nc,case4q,up_down_quark,color,tree_color,lo_or_nlo));

return VSM;
}


BH_Ampl_2q2Q1y::BH_Ampl_2q2Q1y(bool up_down_quark,int case4q,int color, int tree_color, const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q1y(mom_assignment,
					0,
					settings::BH_interface_settings::s_nf,
					settings::BH_interface_settings::s_nc,
					case4q,
					up_down_quark,
					color,
					tree_color,
					lo_or_nlo),
			bhi,  // parent
			5,    // NbrExtParticles
			2,	  // NbrPowersOfAlphaS
			1,    // NbrPowersOfAlphaQED
			2,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
{}


}
