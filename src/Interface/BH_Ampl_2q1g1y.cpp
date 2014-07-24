/*
 * matrix elements for 2q1g1y processes
 *
 *  Created on: Dec 12, 2008
 *
 */



#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include "cached_integral.h"
#include "BH_interface_impl.h"

#define _VERBOSE 0


using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {



partial_amplitude_cached* A_loop_2q_1g_1y_5_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);

	process pro_1=process(h1,h3,h2,h4);
	vector<int> ind_1;
	ind_1.push_back(i1);
	ind_1.push_back(i3);
	ind_1.push_back(i2);
	ind_1.push_back(i4);

	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction r4(1,2), r5(-1,2*n_c*n_c);
	multi_precision_fraction r6(1,6), r7(-1,6*n_c*n_c);

if(color==1){
//Leading color
	PA->add(pro,leading_color,ind,1,1);
	// Scheme shift
	PA->add_subtraction(pro,ind,r4,0);
	// renormalization
	PA->add_subtraction(pro,ind,r4*r1,-1);
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro_1,sub_leading_color,ind_1,1,n_c*n_c);  	///CHECK SIGN!!!!
	PA->add(pro,nf,ind,n_f,n_c);
	// Scheme shift
	PA->add_subtraction(pro,ind,r4+r5,0);
	// renormalization
	PA->add_subtraction(pro,ind,r4*(r1+r2+r3),-1);
}
	return PA;
}



Squared_ME* A_loop_2q_1g_1y_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark,int color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM=new Squared_ME(lo_or_nlo);

/*
multi_precision_fraction ext(4*2*(n_c*n_c-1)*n_c);
multi_precision_fraction ext_tree(2*2*(n_c*n_c-1));
*/
size_t Qnumerator;
if(up_down_quark) {Qnumerator=1;} else {Qnumerator=4;}

multi_precision_fraction _over_nc(1,n_c);
multi_precision_fraction _over_2(1,2);
multi_precision_fraction Q2(Qnumerator,9);

//absorbe helicity sum into factor of two here
multi_precision_fraction ext(4*2*n_c*(n_c*n_c-1)*Q2);
multi_precision_fraction ext_tree=ext*_over_nc*_over_2;


size_t T1=SM->add(new CTree_with_prefactor(pro ,ind  ));

size_t L1=SM->add(A_loop_2q_1g_1y_5_1(pro  ,ind  ,n_s,n_f,n_c,color,lo_or_nlo));



SM->add_loop(cached_cross_term_md(L1,T1,ext));

SM->add_tree(cached_cross_term_md(T1,T1,ext_tree));

return SM;
}

/*
Virtual_SME* vsme_2q1g1y(int ns,int nf,int nc, bool up_down_quark,int color){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);

	ind_b.push_back(3);
	ind_b.push_back(2);
	ind_b.push_back(1);
	ind_b.push_back(4);
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,p,qbp,yp),ind,ns,nf,nc,up_down_quark,color));
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,m,qbp,yp),ind,ns,nf,nc,up_down_quark,color));


	VSM->add( A_loop_2q_1g_1y_M2(process(qm,p,qbp,yp),ind_b,ns,nf,nc,up_down_quark,color));
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,m,qbp,yp),ind_b,ns,nf,nc,up_down_quark,color));

return VSM;
}
*/


Virtual_SME* vsme_2q1g1y(std::vector<int> indext,int ns,int nf,int nc,bool up_down_quark,int color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);

	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[3]);
	
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,p,qbp,yp),ind,ns,nf,nc,up_down_quark,color, lo_or_nlo));
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,m,qbp,yp),ind,ns,nf,nc,up_down_quark,color,lo_or_nlo));

	VSM->add( A_loop_2q_1g_1y_M2(process(qm,p,qbp,yp),ind_b,ns,nf,nc,up_down_quark,color,lo_or_nlo));
	VSM->add( A_loop_2q_1g_1y_M2(process(qm,m,qbp,yp),ind_b,ns,nf,nc,up_down_quark,color,lo_or_nlo));

return VSM;
}


BH_Ampl_2q1g1y::BH_Ampl_2q1g1y(bool up_down_quark,int color,const std::vector<int>& mom_assignment, QCDorder lo_or_nlo, BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q1g1y(mom_assignment,0,
						settings::BH_interface_settings::s_nf,
						settings::BH_interface_settings::s_nc,
						up_down_quark,
						color,
						lo_or_nlo
						),
			bhi,  // parent
			4,    // NbrExtParticles
			1,	  // NbrPowersOfAlphaS
			1,    // NbrPowersOfAlphaQED
			0,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
		{}

}
