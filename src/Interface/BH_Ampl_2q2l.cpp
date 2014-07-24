/*
 * matrix elements for lm lbp -> 2q
 *
 *  Created on: Jan 29, 2008
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
// #define _PHOTON_ONLY 1     // 0 photon only; 1 for all other  replaced by a setting

// Leading color vs. subleading color:



using namespace std;

using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

partial_amplitude_cached* A_loop_2q_2l_4_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const vector<ph_type> _ph_type,int color,QCDorder lo_or_nlo=nlo){
	partial_amplitude_cached* PA=new partial_amplitude_cached(lo_or_nlo);
	int i1=ind.at(0);
	int i2=ind.at(1);
	int i3=ind.at(2);
	int i4=ind.at(3);

	ph_type h1=pro.p(1);
	ph_type h2=pro.p(2);
	ph_type h3=pro.p(3);
	ph_type h4=pro.p(4);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i3,i4,_ph_type);
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
	// proportional to (n-4)/2=0 thus absent
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	// Scheme shift
//	PA->add_subtraction(pro,ind,r4+r5,0));
	PA->add_subtraction(pro,ind,r4,0); // SPECIAL here the color the Nc factors in front of A_{4;1}
	// renormalization
	// proportional to (n-4)/2=0 thus absent
}
	return PA;
}


Squared_ME* A_loop_2q_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW,int color,int tree_color,const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);

ph_type h1=pro.p(1);
ph_type h2=pro.p(2);
ph_type h3=pro.p(3);
ph_type h4=pro.p(4);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i3,i4,_ph_type);
//-------------------------------------------


Squared_ME* SM = new Squared_ME(lo_or_nlo);

multi_precision_fraction ext(4*2*(n_c*n_c-1));
multi_precision_fraction ext_tree(2*2*n_c);

size_t T1=SM->add(new CTree_with_prefactor(pro  ,ind,_prop_hel_fn));

size_t L1=SM->add(A_loop_2q_2l_4_1(pro  ,ind  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));

SM->add_loop(cached_cross_term_md(L1,T1,ext));
SM->add_tree(cached_cross_term_md(T1,T1,ext_tree));


return SM;
}

/*
Virtual_SME* vsme_2q2l(int ns,int nf,int nc, bool up_down_quark, int photonZW,int color, int tree_color){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind43, ind43_b;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);

	ind43.push_back(1);
	ind43.push_back(2);
	ind43.push_back(4);
	ind43.push_back(3);

	ind_b.push_back(2);
	ind_b.push_back(1);
	ind_b.push_back(3);
	ind_b.push_back(4);

	ind43_b.push_back(2);
	ind43_b.push_back(1);
	ind43_b.push_back(4);
	ind43_b.push_back(3);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type43, _ph_type43_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type43.push_back(qp);
	_ph_type43.push_back(lp);
	_ph_type43_b.push_back(qm);
	_ph_type43_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind43,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type43));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind43_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type43_b));

return VSM;
}
*/

Virtual_SME* vsme_2q2l(std::vector<int> indext,int ns,int nf,int nc, bool up_down_quark, int photonZW,int color,int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind43, ind43_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);

	ind43.push_back(indext[0]);
	ind43.push_back(indext[1]);
	ind43.push_back(indext[3]);
	ind43.push_back(indext[2]);

	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[3]);

	ind43_b.push_back(indext[1]);
	ind43_b.push_back(indext[0]);
	ind43_b.push_back(indext[3]);
	ind43_b.push_back(indext[2]);

// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type43, _ph_type43_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type43.push_back(qp);
	_ph_type43.push_back(lp);
	_ph_type43_b.push_back(qm);
	_ph_type43_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind43,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type43,lo_or_nlo));

	VSM->add( A_loop_2q_2l_M2(process(qp,qbm,lm,lbp),ind43_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type43_b,lo_or_nlo));


return VSM;
}


BH_Ampl_2q2l::BH_Ampl_2q2l(bool up_down_quark, int photonZW,int color, int tree_color, const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2l(mom_assignment,0,
									settings::BH_interface_settings::s_nf,
									settings::BH_interface_settings::s_nc,
									up_down_quark,
									photonZW*settings::BH_interface_settings::s_photon_only,
									color,
									tree_color,
									lo_or_nlo),
			bhi,  // parent
			4,    // NbrExtParticles
			0,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			0,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
{}



}
