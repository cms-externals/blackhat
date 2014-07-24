/*

 * matrix elements for lm lbp -> 2q2g
 *
 *  Created on: Oct 21, 2008
 *
 */


#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "BH_interface_impl.h"
#include "cached_OLHA.h"




#define _VERBOSE 0

// various switches: photonZW:  0 for photon
// 				1 for photon+Z
// 				2 for Z->neutrinos only
// 				3 for W only


using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;
namespace BH {

partial_amplitude_cached* A_loop_2q_4g_2l_8_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const vector<ph_type> _ph_type, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_1=process(h1,h6,h5,h4,h3,h2,h7,h8);
	vector<int> ind_1;
	ind_1.push_back(i1);
	ind_1.push_back(i6);
	ind_1.push_back(i5);
	ind_1.push_back(i4);
	ind_1.push_back(i3);
	ind_1.push_back(i2);
	ind_1.push_back(i7);
	ind_1.push_back(i8);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm4_over_2(4,2); // CAREFUL WITH GENERALIZATION TO N-points
	multi_precision_fraction r4(1,2), r5(-1,2*n_c*n_c);
	multi_precision_fraction r6(1,6), r7(-1,6*n_c*n_c);

if(color==1){
//Leading color
	// primitive amplitue
    PA->add(pro,leading_color,ind,1,1);
	if(n_f!=0){
		PA->add(pro,nf,ind,n_f,n_c);
	}
	// Scheme shift
    PA->add_subtraction(pro,ind,r4,0);
	// renormalization
    PA->add_subtraction(pro,ind,(r1+r2+r3)*nm4_over_2,-1);
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro_1,sub_leading_color,ind_1,-1,n_c*n_c);//SIGN
	if(n_f!=0){PA->add(pro,nf,ind,n_f,n_c);}
	// Scheme shift
	PA->add_subtraction(pro,ind,r4+r5,0);
	// renormalization
	PA->add_subtraction(pro,ind,(r1+r2+r3)*nm4_over_2,-1);
}
if(color==2){
//Full minus Leading color
        PA->add(pro_1,sub_leading_color,ind_1,-1,n_c*n_c); //SIGN
        if(n_f!=0){PA->add(pro,nf,ind,n_f,n_c);}
        // Scheme shift
        PA->add_subtraction(pro,ind,r5,0);
        // renormalization
        PA->add_subtraction(pro,ind,(r2+r3)*nm4_over_2,-1);
}


	return PA;
}


Squared_ME* A_loop_2q_4g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, int color, int tree_color, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);


multi_precision_fraction n2m1(n_c*n_c-1);
multi_precision_fraction n2m1_2=n2m1*n2m1;
multi_precision_fraction n2m1_3=n2m1_2*n2m1;
multi_precision_fraction n2m1_4=n2m1_3*n2m1;
multi_precision_fraction n2m2(n_c*n_c-2);
multi_precision_fraction n2p1(n_c*n_c+1);
multi_precision_fraction _over_nc(1,n_c);
multi_precision_fraction _over_nc2=_over_nc*_over_nc;
multi_precision_fraction _over_nc3=_over_nc*_over_nc2;
multi_precision_fraction _8(8);
multi_precision_fraction _4(4);
multi_precision_fraction ext=_8*n2m1*_over_nc;
multi_precision_fraction ext_tree=_4*n2m1*_over_nc*_over_nc;
multi_precision_fraction poly1(-1-n_c*n_c+n_c*n_c*n_c*n_c);
multi_precision_fraction poly2(-1-2*n_c*n_c+n_c*n_c*n_c*n_c);
multi_precision_fraction poly3(1+3*n_c*n_c);


size_t T_123456;
size_t T_123546;
size_t T_124356;
size_t T_124536;
size_t T_125346;
size_t T_125436;
size_t T_132456;
size_t T_132546;
size_t T_134256;
size_t T_134526;
size_t T_135246;
size_t T_135426;
size_t T_142356;
size_t T_142536;
size_t T_143256;
size_t T_143526;
size_t T_145236;
size_t T_145326;
size_t T_152346;
size_t T_152436;
size_t T_153246;
size_t T_153426;
size_t T_154236;
size_t T_154326;


int i1=ind.at(0);
//int i2=ind.at(l1);
//int i3=ind.at(l2);
//int i4=ind.at(l3);
//int i5=ind.at(l4);
int i6=ind.at(5);
int i7=ind.at(6);
int i8=ind.at(7);

ph_type h1=pro.p(1);
//ph_type h2=pro.p(l1+1);
//ph_type h3=pro.p(l2+1);
//ph_type h4=pro.p(l3+1);
//ph_type h5=pro.p(l4+1);
ph_type h6=pro.p(6);
ph_type h7=pro.p(7);
ph_type h8=pro.p(8);
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i7,i8,_ph_type);
//-------------------------------------------

int i2,i3,i4,i5;

for(int l1 = 1; l1 < 5;++l1){
	for(int l2 = 1; l2 < 5; ++l2){
		for(int l3 = 1; l3 < 5; ++l3){
			for(int l4 = 1; l4 < 5; ++l4){
if((l1!=l2)&&(l1!=l3)&&(l2!=l3)&&(l4!=l1&&l4!=l2&&l4!=l3)){

//int i1=ind.at(0);
i2=ind.at(l1);
i3=ind.at(l2);
i4=ind.at(l3);
i5=ind.at(l4);
//int i6=ind.at(5);
//int i7=ind.at(6);
//int i8=ind.at(7);

//ph_type h1=pro.p(1);
ph_type h2=pro.p(l1+1);
ph_type h3=pro.p(l2+1);
ph_type h4=pro.p(l3+1);
ph_type h5=pro.p(l4+1);
//ph_type h6=pro.p(6);
//ph_type h7=pro.p(7);
//ph_type h8=pro.p(8);


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

process pro_123546=process(h1,h2,h3,h5,h4,h6,h7,h8);
vector<int> ind_123546;
ind_123546.push_back(i1);
ind_123546.push_back(i2);
ind_123546.push_back(i3);
ind_123546.push_back(i5);
ind_123546.push_back(i4);
ind_123546.push_back(i6);
ind_123546.push_back(i7);
ind_123546.push_back(i8);

process pro_124356=process(h1,h2,h4,h3,h5,h6,h7,h8);
vector<int> ind_124356;
ind_124356.push_back(i1);
ind_124356.push_back(i2);
ind_124356.push_back(i4);
ind_124356.push_back(i3);
ind_124356.push_back(i5);
ind_124356.push_back(i6);
ind_124356.push_back(i7);
ind_124356.push_back(i8);

process pro_124536=process(h1,h2,h4,h5,h3,h6,h7,h8);
vector<int> ind_124536;
ind_124536.push_back(i1);
ind_124536.push_back(i2);
ind_124536.push_back(i4);
ind_124536.push_back(i5);
ind_124536.push_back(i3);
ind_124536.push_back(i6);
ind_124536.push_back(i7);
ind_124536.push_back(i8);

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

process pro_125436=process(h1,h2,h5,h4,h3,h6,h7,h8);
vector<int> ind_125436;
ind_125436.push_back(i1);
ind_125436.push_back(i2);
ind_125436.push_back(i5);
ind_125436.push_back(i4);
ind_125436.push_back(i3);
ind_125436.push_back(i6);
ind_125436.push_back(i7);
ind_125436.push_back(i8);

process pro_132456=process(h1,h3,h2,h4,h5,h6,h7,h8);
vector<int> ind_132456;
ind_132456.push_back(i1);
ind_132456.push_back(i3);
ind_132456.push_back(i2);
ind_132456.push_back(i4);
ind_132456.push_back(i5);
ind_132456.push_back(i6);
ind_132456.push_back(i7);
ind_132456.push_back(i8);

process pro_132546=process(h1,h3,h2,h5,h4,h6,h7,h8);
vector<int> ind_132546;
ind_132546.push_back(i1);
ind_132546.push_back(i3);
ind_132546.push_back(i2);
ind_132546.push_back(i5);
ind_132546.push_back(i4);
ind_132546.push_back(i6);
ind_132546.push_back(i7);
ind_132546.push_back(i8);

process pro_134256=process(h1,h3,h4,h2,h5,h6,h7,h8);
vector<int> ind_134256;
ind_134256.push_back(i1);
ind_134256.push_back(i3);
ind_134256.push_back(i4);
ind_134256.push_back(i2);
ind_134256.push_back(i5);
ind_134256.push_back(i6);
ind_134256.push_back(i7);
ind_134256.push_back(i8);

process pro_134526=process(h1,h3,h4,h5,h2,h6,h7,h8);
vector<int> ind_134526;
ind_134526.push_back(i1);
ind_134526.push_back(i3);
ind_134526.push_back(i4);
ind_134526.push_back(i5);
ind_134526.push_back(i2);
ind_134526.push_back(i6);
ind_134526.push_back(i7);
ind_134526.push_back(i8);

process pro_135246=process(h1,h3,h5,h2,h4,h6,h7,h8);
vector<int> ind_135246;
ind_135246.push_back(i1);
ind_135246.push_back(i3);
ind_135246.push_back(i5);
ind_135246.push_back(i2);
ind_135246.push_back(i4);
ind_135246.push_back(i6);
ind_135246.push_back(i7);
ind_135246.push_back(i8);

process pro_135426=process(h1,h3,h5,h4,h2,h6,h7,h8);
vector<int> ind_135426;
ind_135426.push_back(i1);
ind_135426.push_back(i3);
ind_135426.push_back(i5);
ind_135426.push_back(i4);
ind_135426.push_back(i2);
ind_135426.push_back(i6);
ind_135426.push_back(i7);
ind_135426.push_back(i8);

process pro_142356=process(h1,h4,h2,h3,h5,h6,h7,h8);
vector<int> ind_142356;
ind_142356.push_back(i1);
ind_142356.push_back(i4);
ind_142356.push_back(i2);
ind_142356.push_back(i3);
ind_142356.push_back(i5);
ind_142356.push_back(i6);
ind_142356.push_back(i7);
ind_142356.push_back(i8);

process pro_142536=process(h1,h4,h2,h5,h3,h6,h7,h8);
vector<int> ind_142536;
ind_142536.push_back(i1);
ind_142536.push_back(i4);
ind_142536.push_back(i2);
ind_142536.push_back(i5);
ind_142536.push_back(i3);
ind_142536.push_back(i6);
ind_142536.push_back(i7);
ind_142536.push_back(i8);

process pro_143256=process(h1,h4,h3,h2,h5,h6,h7,h8);
vector<int> ind_143256;
ind_143256.push_back(i1);
ind_143256.push_back(i4);
ind_143256.push_back(i3);
ind_143256.push_back(i2);
ind_143256.push_back(i5);
ind_143256.push_back(i6);
ind_143256.push_back(i7);
ind_143256.push_back(i8);

process pro_143526=process(h1,h4,h3,h5,h2,h6,h7,h8);
vector<int> ind_143526;
ind_143526.push_back(i1);
ind_143526.push_back(i4);
ind_143526.push_back(i3);
ind_143526.push_back(i5);
ind_143526.push_back(i2);
ind_143526.push_back(i6);
ind_143526.push_back(i7);
ind_143526.push_back(i8);

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

process pro_145326=process(h1,h4,h5,h3,h2,h6,h7,h8);
vector<int> ind_145326;
ind_145326.push_back(i1);
ind_145326.push_back(i4);
ind_145326.push_back(i5);
ind_145326.push_back(i3);
ind_145326.push_back(i2);
ind_145326.push_back(i6);
ind_145326.push_back(i7);
ind_145326.push_back(i8);

process pro_152346=process(h1,h5,h2,h3,h4,h6,h7,h8);
vector<int> ind_152346;
ind_152346.push_back(i1);
ind_152346.push_back(i5);
ind_152346.push_back(i2);
ind_152346.push_back(i3);
ind_152346.push_back(i4);
ind_152346.push_back(i6);
ind_152346.push_back(i7);
ind_152346.push_back(i8);

process pro_152436=process(h1,h5,h2,h4,h3,h6,h7,h8);
vector<int> ind_152436;
ind_152436.push_back(i1);
ind_152436.push_back(i5);
ind_152436.push_back(i2);
ind_152436.push_back(i4);
ind_152436.push_back(i3);
ind_152436.push_back(i6);
ind_152436.push_back(i7);
ind_152436.push_back(i8);

process pro_153246=process(h1,h5,h3,h2,h4,h6,h7,h8);
vector<int> ind_153246;
ind_153246.push_back(i1);
ind_153246.push_back(i5);
ind_153246.push_back(i3);
ind_153246.push_back(i2);
ind_153246.push_back(i4);
ind_153246.push_back(i6);
ind_153246.push_back(i7);
ind_153246.push_back(i8);

process pro_153426=process(h1,h5,h3,h4,h2,h6,h7,h8);
vector<int> ind_153426;
ind_153426.push_back(i1);
ind_153426.push_back(i5);
ind_153426.push_back(i3);
ind_153426.push_back(i4);
ind_153426.push_back(i2);
ind_153426.push_back(i6);
ind_153426.push_back(i7);
ind_153426.push_back(i8);

process pro_154236=process(h1,h5,h4,h2,h3,h6,h7,h8);
vector<int> ind_154236;
ind_154236.push_back(i1);
ind_154236.push_back(i5);
ind_154236.push_back(i4);
ind_154236.push_back(i2);
ind_154236.push_back(i3);
ind_154236.push_back(i6);
ind_154236.push_back(i7);
ind_154236.push_back(i8);

process pro_154326=process(h1,h5,h4,h3,h2,h6,h7,h8);
vector<int> ind_154326;
ind_154326.push_back(i1);
ind_154326.push_back(i5);
ind_154326.push_back(i4);
ind_154326.push_back(i3);
ind_154326.push_back(i2);
ind_154326.push_back(i6);
ind_154326.push_back(i7);
ind_154326.push_back(i8);
//---


if(tree_color==0){ 
//leading color: T_123456=SM->add(new CTree_with_prefactor(pro_123456,ind_123456,_prop_hel_fn), ind_123456);
T_123546=SM->add(new CTree_with_prefactor(pro_123546,ind_123546,_prop_hel_fn), ind_123546);
T_124356=SM->add(new CTree_with_prefactor(pro_124356,ind_124356,_prop_hel_fn), ind_124356);
T_124536=SM->add(new CTree_with_prefactor(pro_124536,ind_124536,_prop_hel_fn), ind_124536);
T_125346=SM->add(new CTree_with_prefactor(pro_125346,ind_125346,_prop_hel_fn), ind_125346);
T_125436=SM->add(new CTree_with_prefactor(pro_125436,ind_125436,_prop_hel_fn), ind_125436);
T_132456=SM->add(new CTree_with_prefactor(pro_132456,ind_132456,_prop_hel_fn), ind_132456);
T_132546=SM->add(new CTree_with_prefactor(pro_132546,ind_132546,_prop_hel_fn), ind_132546);
T_134256=SM->add(new CTree_with_prefactor(pro_134256,ind_134256,_prop_hel_fn), ind_134256);
T_134526=SM->add(new CTree_with_prefactor(pro_134526,ind_134526,_prop_hel_fn), ind_134526);
T_135246=SM->add(new CTree_with_prefactor(pro_135246,ind_135246,_prop_hel_fn), ind_135246);
T_135426=SM->add(new CTree_with_prefactor(pro_135426,ind_135426,_prop_hel_fn), ind_135426);
T_142356=SM->add(new CTree_with_prefactor(pro_142356,ind_142356,_prop_hel_fn), ind_142356);
T_142536=SM->add(new CTree_with_prefactor(pro_142536,ind_142536,_prop_hel_fn), ind_142536);
T_143256=SM->add(new CTree_with_prefactor(pro_143256,ind_143256,_prop_hel_fn), ind_143256);
T_143526=SM->add(new CTree_with_prefactor(pro_143526,ind_143526,_prop_hel_fn), ind_143526);
T_145236=SM->add(new CTree_with_prefactor(pro_145236,ind_145236,_prop_hel_fn), ind_145236);
T_145326=SM->add(new CTree_with_prefactor(pro_145326,ind_145326,_prop_hel_fn), ind_145326);
T_152346=SM->add(new CTree_with_prefactor(pro_152346,ind_152346,_prop_hel_fn), ind_152346);
T_152436=SM->add(new CTree_with_prefactor(pro_152436,ind_152436,_prop_hel_fn), ind_152436);
T_153246=SM->add(new CTree_with_prefactor(pro_153246,ind_153246,_prop_hel_fn), ind_153246);
T_153426=SM->add(new CTree_with_prefactor(pro_153426,ind_153426,_prop_hel_fn), ind_153426);
T_154236=SM->add(new CTree_with_prefactor(pro_154236,ind_154236,_prop_hel_fn), ind_154236);
T_154326=SM->add(new CTree_with_prefactor(pro_154326,ind_154326,_prop_hel_fn), ind_154326);
}

//----------------------


	T_123456=SM->add(new CTree_with_prefactor(pro_123456,ind_123456  ,_prop_hel_fn), ind_123456  );
	SM->add_tree(cached_cross_term_md(T_123456,T_123456,_4*n2m1_4*_over_nc3)); 
	
    size_t L1_123456=SM->add(A_loop_2q_4g_2l_8_1(pro_123456 ,ind_123456  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));
	SM->add_loop(cached_cross_term_md(L1_123456,T_123456,_8*n2m1_4*_over_nc2)); 

//leading color: SM->add_loop(cached_cross_term_md(L1_123456,T_123456,_8*n2m1_4*_over_nc2)); 
//may use somewhat non-standard leading color approxiamtion; switched off here to be consistent with earlier runs
#if 0    
SM->add_loop(cached_cross_term_md(L1_123456,T_123546,-_8*(n2m1_3*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_124356,-_8*(n2m1_3*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_124536,_8*n2m1_2*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_125346,_8*n2m1_2*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_125436,_8*n2m1_2*n2p1*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_132456,-_8*(n2m1_3*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_132546,_8*n2m1_2*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_134256,_8*n2m1_2*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_134526,-_8*(n2m1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_135246,-_8*(n2m1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_135426,-_8*(n2m1*n2p1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_142356,_8*n2m1_2*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_142536,-_8*(n2m1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_143256,_8*n2m1_2*n2p1*_over_nc2)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_143526,-_8*(n2m1*n2p1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_145236,_8*n2m1*poly1*_over_nc2));//(-1-Power(nc,2)+Power(nc,4)); 
SM->add_loop(cached_cross_term_md(L1_123456,T_145326,_8*n2m1*poly2*_over_nc2));//(-1-2*Power(nc,2)+Power(nc,4)) 
SM->add_loop(cached_cross_term_md(L1_123456,T_152346,-_8*(n2m1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_152436,-_8*(n2m1*n2p1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_153246,-_8*(n2m1*n2p1*_over_nc2))); 
SM->add_loop(cached_cross_term_md(L1_123456,T_153426,_8*n2m1*poly2*_over_nc2));//(-1-2*Power(nc,2)+Power(nc,4)) 
SM->add_loop(cached_cross_term_md(L1_123456,T_154236,_8*n2m1*poly2*_over_nc2));//(-1-2*Power(nc,2)+Power(nc,4))
SM->add_loop(cached_cross_term_md(L1_123456,T_154326,-_8*(n2m1*poly3)*_over_nc2));//(1+3*Power(nc,2));
#endif


if(tree_color==0){ 
//leading color: SM->add_tree(cached_cross_term_md(T_123456,T_123456,_4*n2m1_4*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_123546,T_123456,-_4*(n2m1_3*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_124356,T_123456,-_4*(n2m1_3*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_124536,T_123456,_4*n2m1_2*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125346,T_123456,_4*n2m1_2*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_125436,T_123456,_4*n2m1_2*n2p1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_132456,T_123456,-_4*(n2m1_3*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_132546,T_123456,_4*n2m1_2*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_134256,T_123456,_4*n2m1_2*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_134526,T_123456,-_4*(n2m1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_135246,T_123456,-_4*(n2m1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_135426,T_123456,-_4*(n2m1*n2p1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_142356,T_123456,_4*n2m1_2*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_142536,T_123456,-_4*(n2m1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_143256,T_123456,_4*n2m1_2*n2p1*_over_nc3)); 
SM->add_tree(cached_cross_term_md(T_143526,T_123456,-_4*(n2m1*n2p1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_145236,T_123456,_4*n2m1*poly1*_over_nc3));//(-1-Power(nc,2)+Power(nc,4)); 
SM->add_tree(cached_cross_term_md(T_145326,T_123456,_4*n2m1*poly2*_over_nc3));//(-1-2*Power(nc,2)+Power(nc,4)) 
SM->add_tree(cached_cross_term_md(T_152346,T_123456,-_4*(n2m1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_152436,T_123456,-_4*(n2m1*n2p1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_153246,T_123456,-_4*(n2m1*n2p1*_over_nc3))); 
SM->add_tree(cached_cross_term_md(T_153426,T_123456,_4*n2m1*poly2*_over_nc3));//(-1-2*Power(nc,2)+Power(nc,4)) 
SM->add_tree(cached_cross_term_md(T_154236,T_123456,_4*n2m1*poly2*_over_nc3));//(-1-2*Power(nc,2)+Power(nc,4))
SM->add_tree(cached_cross_term_md(T_154326,T_123456,-_4*(n2m1*poly3)*_over_nc3));//(1+3*Power(nc,2));
}

};};};};};
return SM;
}



Virtual_SME* vsme_2q4g2l(std::vector<int> indext,int n_s,int n_f,int n_c, bool up_down_quark, int photonZW,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind87, ind87_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);
	ind.push_back(indext[6]);
	ind.push_back(indext[7]);

	ind87.push_back(indext[0]);
	ind87.push_back(indext[1]);
	ind87.push_back(indext[2]);
	ind87.push_back(indext[3]);
	ind87.push_back(indext[4]);
	ind87.push_back(indext[5]);
	ind87.push_back(indext[7]);
	ind87.push_back(indext[6]);

	ind_b.push_back(indext[5]);
	ind_b.push_back(indext[4]);
	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[6]);
	ind_b.push_back(indext[7]);

	ind87_b.push_back(indext[5]);
	ind87_b.push_back(indext[4]);
	ind87_b.push_back(indext[3]);
	ind87_b.push_back(indext[2]);
	ind87_b.push_back(indext[1]);
	ind87_b.push_back(indext[0]);
	ind87_b.push_back(indext[7]);
	ind87_b.push_back(indext[6]);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type87, _ph_type87_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type87.push_back(qp);
	_ph_type87.push_back(lp);
	_ph_type87_b.push_back(qm);
	_ph_type87_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

if(photonZW!=3){
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
		VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,p,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
}
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,p,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

if(photonZW!=3){
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,p,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));

	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,p,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
}

if(photonZW!=3){ 
    	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,m,qbm,lm,lbp),ind,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
}
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,m,qbm,lm,lbp),ind_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

if(photonZW!=3){
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,m,qbm,lm,lbp),ind87,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87,lo_or_nlo));

	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,p,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,m,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,m,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,m,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,m,p,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,m,p,p,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,m,p,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
	VSM->add( A_loop_2q_4g_2l_M2(process(qp,p,p,m,m,qbm,lm,lbp),ind87_b,n_s,n_f,n_c,up_down_quark,photonZW,color,tree_color,_ph_type87_b,lo_or_nlo));
}

return VSM;
}

BH_Ampl_2q4g2l::BH_Ampl_2q4g2l(bool up_down_quark, int photonZW, int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q4g2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				up_down_quark,
				photonZW*settings::BH_interface_settings::s_photon_only,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			8,    // NbrExtParticles
			4,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			8,   // GeVdim
			1.,  // factor
			0., //additional scheme shift
			0. //additional renormalization shift
			),
		momenta_assignment(mom_assignment), d_color(color)
{}


double BH_Ampl_2q4g2l::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:
		switch(d_color){
		//Full color
		case 0 : return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
		//Leading color
        //notice full color renormalization
		//case 1: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 33/6. /*beta_0*/ *get_double_pole()/(-15.) /* = lc born tree/full tree */  );
		case 1: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/ *get_double_pole()/(-15.) /* = lc born tree/full tree */  );
		//Full-Leading color
		case 2: return get_single_pole()+(4 /*NbrPowersOfAlphaS*/* ( 23/6. -  33/6. /*beta_0*/ *(35/36. - get_double_pole()/(-15.)  /* = lc born tree */ ))  );
		}
	case 2:	return get_double_pole();
	}

};


}
