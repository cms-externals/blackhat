/*
 * matrix elements for 2q2g1y processes
 *
 *  Created on: Oct 26, 2008
 *
 */



#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include "cached_integral.h"
#include "BH_interface_impl.h"

#define _VERBOSE 0
#define _INCL_VECT 0   // 1 including vectors => contribs have not been coded yet


using namespace std;

using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {



partial_amplitude_cached* A_loop_2q_2g_1y_6_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,int color,QCDorder lo_or_nlo){
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

	process pro_1=process(h1,h4,h3,h2,h5);
	vector<int> ind_1;
	ind_1.push_back(i1);
	ind_1.push_back(i4);
	ind_1.push_back(i3);
	ind_1.push_back(i2);
	ind_1.push_back(i5);

	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction r4(1,2), r5(-1,2*n_c*n_c);

if(color==1){
//Leading color
	PA->add(pro,leading_color,ind,1,1);
	// Scheme shift
	PA->add_subtraction(pro,ind,r4,0);
	// renormalization
	PA->add_subtraction(pro,ind,r1,-1);
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro_1,sub_leading_color,ind_1,-1,n_c*n_c);
	PA->add(pro,nf,ind,n_f,n_c);
	// Scheme shift
	PA->add_subtraction(pro,ind,r4+r5,0);
	// renormalization
	PA->add_subtraction(pro,ind,(r1+r2+r3),-1);
}
	return PA;
}

// partial amplitudes not needed for leading color
partial_amplitude_cached* A_loop_2q_2g_1y_6_3(process pro, vector<int> & ind, int n_s, int n_f, int n_c,QCDorder lo_or_nlo){
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


//eqn 2.13 1324
process pro_1=process(h1,h3,h2,h4,h5);
vector<int> ind_1;
ind_1.push_back(i1);
ind_1.push_back(i3);
ind_1.push_back(i2);
ind_1.push_back(i4);
ind_1.push_back(i5);

//eqn 2.13 1243
process pro_2=process(h1,h2,h4,h3,h5);
vector<int> ind_2;
ind_2.push_back(i1);
ind_2.push_back(i2);
ind_2.push_back(i4);
ind_2.push_back(i3);
ind_2.push_back(i5);

//eqn 2.13 1342
process pro_3=process(h1,h3,h4,h2,h5);
vector<int> ind_3;
ind_3.push_back(i1);
ind_3.push_back(i3);
ind_3.push_back(i4);
ind_3.push_back(i2);
ind_3.push_back(i5);

//eqn 2.13 1423
process pro_4=process(h1,h4,h2,h3,h5);
vector<int> ind_4;
ind_4.push_back(i1);
ind_4.push_back(i4);
ind_4.push_back(i2);
ind_4.push_back(i3);
ind_4.push_back(i5);

//eqn 2.13 1432
process pro_5=process(h1,h4,h3,h2,h5);
vector<int> ind_5;
ind_5.push_back(i1);
ind_5.push_back(i4);
ind_5.push_back(i3);
ind_5.push_back(i2);
ind_5.push_back(i5);

PA->add(pro  ,leading_color,ind,1,1);
PA->add(pro_1,leading_color,ind_1,1,1);
PA->add(pro_2,sub_leading_color,ind_2,1,1);
PA->add(pro_3,sub_leading_color,ind_3,1,1);
PA->add(pro_4,sub_leading_color,ind_4,1,1);
PA->add(pro_5,sub_leading_color,ind_5,1,1);
return PA;

}

partial_amplitude_cached* A_loop_2q_2g_1y_6_4_V(process pro, vector<int> & ind, int n_s, int n_f, int n_c,QCDorder lo_or_nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);

PA->add(pro,VECT,ind,1,1);
return PA;
}




Squared_ME* A_loop_2q_2g_1y_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int color, int tree_color,QCDorder lo_or_nlo=nlo){

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

process pro_b=process(h1,h3,h2,h4,h5);
vector<int> ind_b;
ind_b.push_back(i1);
ind_b.push_back(i3);
ind_b.push_back(i2);
ind_b.push_back(i4);
ind_b.push_back(i5);

process pro_va=process(h1,h4,h2,h3,h5);
vector<int> ind_va;
ind_va.push_back(i1);
ind_va.push_back(i4);
ind_va.push_back(i2);
ind_va.push_back(i3);
ind_va.push_back(i5);

process pro_vb=process(h1,h4,h3,h2,h5);
vector<int> ind_vb;
ind_vb.push_back(i1);
ind_vb.push_back(i4);
ind_vb.push_back(i3);
ind_vb.push_back(i2);
ind_vb.push_back(i5);

Squared_ME* SM=new Squared_ME(lo_or_nlo);


size_t Qnumerator;
if(up_down_quark) {Qnumerator=-1;} else {Qnumerator=2;}

multi_precision_fraction Q(Qnumerator,3);
multi_precision_fraction ext(4*(n_c*n_c-1));
multi_precision_fraction ext_tree(2*(n_c*n_c-1),n_c);


size_t T1=SM->add(new CTree_with_prefactor(pro  ,ind  ));
size_t T2=SM->add(new CTree_with_prefactor(pro_b,ind_b));

size_t L1=SM->add(A_loop_2q_2g_1y_6_1(pro  ,ind  ,n_s,n_f,n_c,color,lo_or_nlo));
size_t L2=SM->add(A_loop_2q_2g_1y_6_1(pro_b,ind_b,n_s,n_f,n_c,color,lo_or_nlo));

size_t L3=SM->add(A_loop_2q_2g_1y_6_3(pro  ,ind  ,n_s,n_f,n_c,lo_or_nlo));
size_t L4=SM->add(A_loop_2q_2g_1y_6_3(pro_b,ind_b,n_s,n_f,n_c,lo_or_nlo));

#if _INCL_VECT
multi_precision_fraction vect((n_c*n_c-4),n_c);

size_t LV=SM->add(A_loop_2q_2g_2l_6_4_V(pro_va,ind_va,n_s,n_f,n_c,lo_or_nlo));
size_t LV_b=SM->add(A_loop_2q_2g_2l_6_4_V(pro_vb,ind_vb,n_s,n_f,n_c,lo_or_nlo));
#endif

if(color==1){
//Leading color
SM->add_loop(cached_cross_term_md(L1,T1,ext*Q*Q*(n_c*n_c-1)));
//SM->add_loop(cached_cross_term_md(L2,T1,-ext*Q*Q));
}
if(color==0){
//Full color
SM->add_loop(cached_cross_term_md(L1,T1,ext*Q*Q*(n_c*n_c-1)));
SM->add_loop(cached_cross_term_md(L2,T1,-ext*Q*Q));
SM->add_loop(cached_cross_term_md(L3,T1,ext*Q*Q));
}

if(tree_color==1){
//Leading color
SM->add_tree(cached_cross_term_md(T1,T1,ext_tree*Q*Q*(n_c*n_c-1)));
}
if(tree_color==0){
//Full color
SM->add_tree(cached_cross_term_md(T1,T1,ext_tree*Q*Q*(n_c*n_c-1)));
SM->add_tree(cached_cross_term_md(T1,T2,-ext_tree*Q*Q));
}

if(color==1){
//Leading color
SM->add_loop(cached_cross_term_md(L2,T2,ext*Q*Q*(n_c*n_c-1)));
//SM->add_loop(cached_cross_term_md(L1,T2,-ext*Q*Q));
}
if(tree_color==0){
//Full color
SM->add_loop(cached_cross_term_md(L2,T2,ext*Q*Q*(n_c*n_c-1)));
SM->add_loop(cached_cross_term_md(L1,T2,-ext*Q*Q));
SM->add_loop(cached_cross_term_md(L4,T2,ext*Q*Q));
}

if(tree_color==1){
//Leading color
SM->add_tree(cached_cross_term_md(T2,T2,ext_tree*Q*Q*(n_c*n_c-1)));
}
if(tree_color==0){
//Leading color
SM->add_tree(cached_cross_term_md(T2,T2,ext_tree*Q*Q*(n_c*n_c-1)));
SM->add_tree(cached_cross_term_md(T2,T1,-ext_tree*Q*Q));
}

#if _INCL_VECT
SM->add_loop(cached_cross_term_md(LV,T1  ,ext*vect*Q*Q));
SM->add_loop(cached_cross_term_md(LV_b,T2,ext*vect*Q*Q));
#endif


return SM;
}

/*
Virtual_SME* vsme_2q2g1y(int ns,int nf,int nc, bool up_down_quark,int color,int tree_color){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);

	ind_b.push_back(4);
	ind_b.push_back(3);
	ind_b.push_back(2);
	ind_b.push_back(1);
	ind_b.push_back(5);
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color));


	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color));

return VSM;
}
*/

Virtual_SME* vsme_2q2g1y(std::vector<int> indext,int ns,int nf,int nc,bool up_down_quark,int color,int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);

	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[4]);

//    VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,yp),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));

//	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,yp),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));

    VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,ym),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
//	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,ym),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,ym),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,ym),ind,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));

	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
//	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,p,m,qbm,ym),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));
	VSM->add( A_loop_2q_2g_1y_M2(process(qp,m,p,qbm,ym),ind_b,ns,nf,nc,up_down_quark,color,tree_color,lo_or_nlo));



return VSM;
}


BH_Ampl_2q2g1y::BH_Ampl_2q2g1y(bool up_down_quark,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2g1y(mom_assignment,0,
					settings::BH_interface_settings::s_nf,
					settings::BH_interface_settings::s_nc,
					up_down_quark,
					settings::BH_interface_settings::s_BH_color_mode,
					settings::BH_interface_settings::s_BH_tree_color_mode,
					lo_or_nlo
					),
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
