/*
 * matrix elements for lm lbp -> 2q2g
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
//#define _PHOTON_ONLY 1     // 0 photon only; 1 for all other   replaced by a setting

// include nf_top pieces
#define Include_nf_top_Pieces 0
#define Include_axial_Pieces 0
#define Include_vectorial_Pieces 1


using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;
namespace BH {

partial_amplitude_cached* A_loop_2q_2g_2l_6_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const vector<ph_type> _ph_type, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_1=process(h1,h4,h3,h2,h5,h6);
	vector<int> ind_1;
	ind_1.push_back(i1);
	ind_1.push_back(i4);
	ind_1.push_back(i3);
	ind_1.push_back(i2);
	ind_1.push_back(i5);
	ind_1.push_back(i6);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction r4(1,2), r5(-1,2*n_c*n_c);
	multi_precision_fraction r6(1,6), r7(-1,6*n_c*n_c);

if(color==1){
//Leading color
	PA->add(pro,leading_color,ind,1,1);
	// Scheme shift
	PA->add_subtraction(pro,ind,r4,0);
	// renormalization
	PA->add_subtraction(pro,ind,r1,-1);
}
if(color==0){
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro_1,sub_leading_color,ind_1,-1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
#if Include_nf_top_Pieces
	PA->add(pro,nf_top,ind,1,n_c);
#endif
	// Scheme shift
	PA->add_subtraction(pro,ind,r4+r5,0);
	// renormalization
	PA->add_subtraction(pro,ind,(r1+r2+r3),-1);
}
	return PA;
}

partial_amplitude_cached* A_loop_2q_2g_2l_6_3(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
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


//eqn 2.13 1324
process pro_1=process(h1,h3,h2,h4,h5,h6);
vector<int> ind_1;
ind_1.push_back(i1);
ind_1.push_back(i3);
ind_1.push_back(i2);
ind_1.push_back(i4);
ind_1.push_back(i5);
ind_1.push_back(i6);

//eqn 2.13 1243
process pro_2=process(h1,h2,h4,h3,h5,h6);
vector<int> ind_2;
ind_2.push_back(i1);
ind_2.push_back(i2);
ind_2.push_back(i4);
ind_2.push_back(i3);
ind_2.push_back(i5);
ind_2.push_back(i6);

//eqn 2.13 1342
process pro_3=process(h1,h3,h4,h2,h5,h6);
vector<int> ind_3;
ind_3.push_back(i1);
ind_3.push_back(i3);
ind_3.push_back(i4);
ind_3.push_back(i2);
ind_3.push_back(i5);
ind_3.push_back(i6);

//eqn 2.13 1423
process pro_4=process(h1,h4,h2,h3,h5,h6);
vector<int> ind_4;
ind_4.push_back(i1);
ind_4.push_back(i4);
ind_4.push_back(i2);
ind_4.push_back(i3);
ind_4.push_back(i5);
ind_4.push_back(i6);

//eqn 2.13 1432
process pro_5=process(h1,h4,h3,h2,h5,h6);
vector<int> ind_5;
ind_5.push_back(i1);
ind_5.push_back(i4);
ind_5.push_back(i3);
ind_5.push_back(i2);
ind_5.push_back(i5);
ind_5.push_back(i6);



//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

PA->add(pro  ,leading_color,ind,1,1);
PA->add(pro_1,leading_color,ind_1,1,1);
PA->add(pro_2,sub_leading_color,ind_2,1,1);
PA->add(pro_3,sub_leading_color,ind_3,1,1);
PA->add(pro_4,sub_leading_color,ind_4,1,1);
PA->add(pro_5,sub_leading_color,ind_5,1,1);
return PA;

}



partial_amplitude_cached* A_loop_2q_2g_2l_6_4_V(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
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

int leading_vect_ax=1;

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

#if 1
 PA->add(pro,VECT,ind,1,1);
#else
 _WARNING("VECT cs disabled!");
#endif
return PA;
}

partial_amplitude_cached* A_loop_2q_2g_2l_6_4_AX(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
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


int leading_vect_ax=2;
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


PA->add(pro,AX,ind,1,1);
return PA;
}

partial_amplitude_cached* A_loop_2q_2g_2l_6_5_AX(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA = new partial_amplitude_cached(lo_or_nlo);;

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


int leading_vect_ax=2;
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i5,i6,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------

PA->add(pro,AXSL,ind,1,1);
return PA;
}



Squared_ME* A_loop_2q_2g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, int color, int tree_color,const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){

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

process pro_b=process(h1,h3,h2,h4,h5,h6);
vector<int> ind_b;
ind_b.push_back(i1);
ind_b.push_back(i3);
ind_b.push_back(i2);
ind_b.push_back(i4);
ind_b.push_back(i5);
ind_b.push_back(i6);

process pro_va=process(h1,h4,h2,h3,h5,h6);
vector<int> ind_va;
ind_va.push_back(i1);
ind_va.push_back(i4);
ind_va.push_back(i2);
ind_va.push_back(i3);
ind_va.push_back(i5);
ind_va.push_back(i6);

process pro_vb=process(h1,h4,h3,h2,h5,h6);
vector<int> ind_vb;
ind_vb.push_back(i1);
ind_vb.push_back(i4);
ind_vb.push_back(i3);
ind_vb.push_back(i2);
ind_vb.push_back(i5);
ind_vb.push_back(i6);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i5,i6,_ph_type);
//-------------------------------------------


Squared_ME* SM = new Squared_ME(lo_or_nlo);


multi_precision_fraction ext(4*2*(n_c*n_c-1));
multi_precision_fraction ext_tree(2*2*(n_c*n_c-1),n_c);

size_t T1=SM->add(new CTree_with_prefactor(pro  , ind, _prop_hel_fn));
size_t T2=SM->add(new CTree_with_prefactor(pro_b, ind_b,_prop_hel_fn));

size_t L1=SM->add(A_loop_2q_2g_2l_6_1(pro  ,ind  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));
size_t L2=SM->add(A_loop_2q_2g_2l_6_1(pro_b,ind_b,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));

//Full color
size_t L3=SM->add(A_loop_2q_2g_2l_6_3(pro  ,ind  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
size_t L4=SM->add(A_loop_2q_2g_2l_6_3(pro_b,ind_b,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));


if(color==1){
//Leading color
SM->add_loop(cached_cross_term_md(L1,T1,ext*(n_c*n_c-1),1.));
//SM->add_loop(cached_cross_term_md(L2,T1,-ext));
}
if(color==0){
//Full color
SM->add_loop(cached_cross_term_md(L1,T1,ext*(n_c*n_c-1),1.));
SM->add_loop(cached_cross_term_md(L2,T1,-ext,1.));
SM->add_loop(cached_cross_term_md(L3,T1,ext,1.));
}

if(tree_color==1){
//Leading color
SM->add_tree(cached_cross_term_md(T1,T1,ext_tree*(n_c*n_c-1),1.));
}
if(tree_color==0){
//Full color
SM->add_tree(cached_cross_term_md(T1,T1,ext_tree*(n_c*n_c-1),1.));
SM->add_tree(cached_cross_term_md(T1,T2,-ext_tree,1.));
}

if(color==1){
//Leading color
SM->add_loop(cached_cross_term_md(L2,T2,ext*(n_c*n_c-1),1.));
//SM->add_loop(cached_cross_term_md(L1,T2,-ext));
}
if(color==0){
//Full color
SM->add_loop(cached_cross_term_md(L2,T2,ext*(n_c*n_c-1),1.));
SM->add_loop(cached_cross_term_md(L1,T2,-ext,1.));
SM->add_loop(cached_cross_term_md(L4,T2,ext,1.));
}

if(tree_color==1){
//Leading color
SM->add_tree(cached_cross_term_md(T2,T2,ext_tree*(n_c*n_c-1),1.));
}
if(tree_color==0){
//Full color
SM->add_tree(cached_cross_term_md(T2,T2,ext_tree*(n_c*n_c-1),1.));
SM->add_tree(cached_cross_term_md(T2,T1,-ext_tree,1.));
}

if(color==0){
//Full color
// axial and vectorial parts
// not needed for W case
if(photonZW!=3){
	multi_precision_fraction vect((n_c*n_c-4),n_c);
	multi_precision_fraction ax1((n_c*n_c-2),n_c);
	multi_precision_fraction ax2(-2,n_c);
	multi_precision_fraction ax3(1,n_c);

#if Include_vectorial_Pieces
	size_t LV=     SM->add(A_loop_2q_2g_2l_6_4_V( pro_va,ind_va,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
	size_t LV_b=   SM->add(A_loop_2q_2g_2l_6_4_V( pro_vb,ind_vb,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
#endif

#if Include_axial_Pieces
	size_t LAX4=   SM->add(A_loop_2q_2g_2l_6_4_AX(pro_va,ind_va,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
	size_t LAX4_b= SM->add(A_loop_2q_2g_2l_6_4_AX(pro_vb,ind_vb,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
	size_t LAX5=   SM->add(A_loop_2q_2g_2l_6_5_AX(pro_va,ind_va,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
	size_t LAX5_b= SM->add(A_loop_2q_2g_2l_6_5_AX(pro_vb,ind_vb,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
#endif

#if Include_vectorial_Pieces
	SM->add_loop(cached_cross_term_md(LV,T1  ,ext*vect,1.));
	SM->add_loop(cached_cross_term_md(LV_b,T2,ext*vect,1.));
#else
_WARNING("Vectorial parts disabled");
#endif
#if Include_axial_Pieces
	SM->add_loop(cached_cross_term_md(LAX4,T1  ,ext*ax1,1.));
	SM->add_loop(cached_cross_term_md(LAX4_b,T2  ,ext*ax1,1.));

	SM->add_loop(cached_cross_term_md(LAX4,T2  ,ext*ax2,1.));
	SM->add_loop(cached_cross_term_md(LAX4_b,T1  ,ext*ax2,1.));

	SM->add_loop(cached_cross_term_md(LAX5,T1  ,ext*ax3,1.));
	SM->add_loop(cached_cross_term_md(LAX5_b,T2  ,ext*ax3,1.));
#else
//_WARNING("Axial parts disabled");
	#endif
}
}
return SM;
}

/*
Virtual_SME* vsme_2q2g2l(int ns,int nf,int nc, bool up_down_quark, int photonZW,int color, int tree_color){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind65, ind65_b;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);
	ind.push_back(6);

	ind65.push_back(1);
	ind65.push_back(2);
	ind65.push_back(3);
	ind65.push_back(4);
	ind65.push_back(6);
	ind65.push_back(5);

	ind_b.push_back(4);
	ind_b.push_back(3);
	ind_b.push_back(2);
	ind_b.push_back(1);
	ind_b.push_back(5);
	ind_b.push_back(6);

	ind65_b.push_back(4);
	ind65_b.push_back(3);
	ind65_b.push_back(2);
	ind65_b.push_back(1);
	ind65_b.push_back(6);
	ind65_b.push_back(5);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type65, _ph_type65_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type65.push_back(qp);
	_ph_type65.push_back(lp);
	_ph_type65_b.push_back(qm);
	_ph_type65_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

if(photonZW!=3){
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type));
}
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b));

if(photonZW!=3){
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65));

	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b));
}

return VSM;
}
*/

Virtual_SME* vsme_2q2g2l(std::vector<int> indext,int ns,int nf,int nc, bool up_down_quark, int photonZW, int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind65, ind65_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);

	ind65.push_back(indext[0]);
	ind65.push_back(indext[1]);
	ind65.push_back(indext[2]);
	ind65.push_back(indext[3]);
	ind65.push_back(indext[5]);
	ind65.push_back(indext[4]);

	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[4]);
	ind_b.push_back(indext[5]);

	ind65_b.push_back(indext[3]);
	ind65_b.push_back(indext[2]);
	ind65_b.push_back(indext[1]);
	ind65_b.push_back(indext[0]);
	ind65_b.push_back(indext[5]);
	ind65_b.push_back(indext[4]);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type65, _ph_type65_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type65.push_back(qp);
	_ph_type65.push_back(lp);
	_ph_type65_b.push_back(qm);
	_ph_type65_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

//speedup for W-cases due to left-handed coupling
if(photonZW!=3){
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
}
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

if(photonZW!=3){
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind65,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65,lo_or_nlo));

	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,p,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,m,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,p,m,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b,lo_or_nlo));
	VSM->add( A_loop_2q_2g_2l_M2(process(qp,m,p,qbm,lm,lbp),ind65_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type65_b,lo_or_nlo));
}

return VSM;
}

BH_Ampl_2q2g2l::BH_Ampl_2q2g2l(bool up_down_quark, int photonZW,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2g2l(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				up_down_quark,
				photonZW*settings::BH_interface_settings::s_photon_only,
				color, 
				tree_color,
				lo_or_nlo
				),
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
