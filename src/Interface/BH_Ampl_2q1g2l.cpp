/*
 * matrix elements for lm lbp -> 2q1g
 *
 *  Created on: Dec 13, 2008
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
// #define _PHOTON_ONLY 0     // 0 photon only; 1 for all other  (replaced by a setting in settings_reader )

#define Include_axial_Pieces 0


using namespace std;

using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

partial_amplitude_cached* A_loop_2q_1g_2l_5_1(process pro, vector<int> & ind, int n_s, int n_f, int n_c,bool up_down_quark,int photonZW, const vector<ph_type> _ph_type,int color,QCDorder lo_or_nlo){
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

	process pro_1=process(h1,h3,h2,h4,h5);
	vector<int> ind_1;
	ind_1.push_back(i1);
	ind_1.push_back(i3);
	ind_1.push_back(i2);
	ind_1.push_back(i4);
	ind_1.push_back(i5);


//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i4,i5,_ph_type);
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
	PA->add_subtraction(pro,ind,r4*r1,-1);  //(n-4)/2=1/2 thus use factor of r4
}
if(color==0){
//Full color
	PA->add(pro,leading_color,ind,1,1);
	PA->add(pro_1,sub_leading_color,ind_1,1,n_c*n_c); ///CHECK SIGN!!!!
	//PA->add(pro,nf,ind,n_f,n_c); not needed since it is a bubble on an external leg
	// Scheme shift
	PA->add_subtraction(pro,ind,r4+r5,0);
	// renormalization
	PA->add_subtraction(pro,ind,r4*(r1+r2+r3),-1);   //(n-4)/2=1/2 thus use factor of r4
}
	return PA;
}

partial_amplitude_cached* A_loop_2q_1g_2l_5_1_AX(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, const vector<ph_type> _ph_type,QCDorder lo_or_nlo){
	partial_amplitude_cached* PA = new partial_amplitude_cached(lo_or_nlo); ;

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


int leading_vect_ax=2;
//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,leading_vect_ax,i4,i5,_ph_type);
	PA->set_prefactor(_prop_hel_fn);
//-------------------------------------------


PA->add(pro,AX,ind,1,1);
return PA;
}


Squared_ME* A_loop_2q_1g_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, bool up_down_quark, int photonZW, int color, int tree_color,const vector<ph_type> _ph_type,QCDorder lo_or_nlo=nlo){

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

process pro_ax=process(h1,h3,h2,h4,h5);
vector<int> ind_ax;
ind_ax.push_back(i1);
ind_ax.push_back(i3);
ind_ax.push_back(i2);
ind_ax.push_back(i4);
ind_ax.push_back(i5);

//--------------------- propagator correction
	prop_hel_fn _prop_hel_fn(up_down_quark,photonZW,0,i4,i5,_ph_type);
//-------------------------------------------


Squared_ME* SM = new Squared_ME(lo_or_nlo);

/*
multi_precision_fraction ext(4*2*n_c*(n_c*n_c-1));
multi_precision_fraction ext_tree(2*2*(n_c*n_c-1));
*/

multi_precision_fraction _over_nc(1,n_c);
multi_precision_fraction _over_2(1,2);
multi_precision_fraction ext(4*2*n_c*(n_c*n_c-1));
multi_precision_fraction ext_tree=ext*_over_nc*_over_2;

size_t T1=SM->add(new CTree_with_prefactor(pro  ,ind,_prop_hel_fn));

size_t L1=SM->add(A_loop_2q_1g_2l_5_1(pro  ,ind  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));



if(!settings::BH_interface_settings::s_use_W_polarization_A567_assembly){
SM->add_loop(cached_cross_term_md(L1,T1,ext));
}

SM->add_tree(cached_cross_term_md(T1,T1,ext_tree));



if(settings::BH_interface_settings::s_use_W_polarization_A567_assembly){

    ph_type h1p=pro.p(1).conjugate();
    ph_type h2p=pro.p(2).conjugate();
    ph_type h3p=pro.p(3).conjugate();
    ph_type h4p=pro.p(4).conjugate();
    ph_type h5p=pro.p(5).conjugate();
    process pro_pflip(h3p,h2p,h1p,h5p,h4p);


    vector<int> ind_p;
    ind_p.push_back(i3);
    ind_p.push_back(i2);
    ind_p.push_back(i1);
    ind_p.push_back(i5);
    ind_p.push_back(i4);



    size_t T1p=SM->add(new CTree_with_prefactor(pro_pflip  ,ind_p,_prop_hel_fn));
    size_t L1p=SM->add(A_loop_2q_1g_2l_5_1(pro_pflip  ,ind_p  ,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,color,lo_or_nlo));

    multi_precision_fraction ext_over2(4*n_c*(n_c*n_c-1));
    multi_precision_fraction m_ext_over2(-4*n_c*(n_c*n_c-1));
    
    SM->add_loop(cached_cross_term_md(L1,T1,ext_over2));
    SM->add_loop(cached_cross_term_md(L1p,T1p,m_ext_over2));
    //SM->add_loop(cached_cross_term_md(L1p,T1p,ext_over2));
}

if(color==0){
//Full color
// axial and vectorial parts
// not needed for W case
if(photonZW!=3){
	multi_precision_fraction ax(4*2*(n_c*n_c-1));

#if Include_axial_Pieces
	size_t LAX= SM->add(A_loop_2q_1g_2l_5_1_AX(pro_ax,ind_ax,n_s,n_f,n_c,up_down_quark,photonZW,_ph_type,lo_or_nlo));
#endif

#if Include_axial_Pieces
	SM->add_loop(cached_cross_term_md(LAX,T1,ax)); // ext_tree used to take into account the Nf'/Nc for Nf'=#t/b-quark type pairs with large mass split
#endif
}
}


return SM;
}

Virtual_SME* vsme_2q1g2l(std::vector<int> indext,int ns,int nf,int nc, bool up_down_quark, int photonZW, int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM=new Virtual_SME() ;
	vector<int> ind,ind_b, ind54, ind54_b;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);

	ind54.push_back(indext[0]);
	ind54.push_back(indext[1]);
	ind54.push_back(indext[2]);
	ind54.push_back(indext[4]);
	ind54.push_back(indext[3]);

	ind_b.push_back(indext[2]);
	ind_b.push_back(indext[1]);
	ind_b.push_back(indext[0]);
	ind_b.push_back(indext[3]);
	ind_b.push_back(indext[4]);

	ind54_b.push_back(indext[2]);
	ind54_b.push_back(indext[1]);
	ind54_b.push_back(indext[0]);
	ind54_b.push_back(indext[4]);
	ind54_b.push_back(indext[3]);


// We need to specify the quark and electron types to make sure that the Zs/Ws couple correctly.
// The primitive amplitudes are charge and parity-conjugation covariant.

	vector<ph_type> _ph_type,_ph_type_b,_ph_type54, _ph_type54_b;
	_ph_type.push_back(qp);
	_ph_type.push_back(lm);
	_ph_type_b.push_back(qm);
	_ph_type_b.push_back(lm);

	_ph_type54.push_back(qp);
	_ph_type54.push_back(lp);
	_ph_type54_b.push_back(qm);
	_ph_type54_b.push_back(lp);


//NOTICE: unpolarized electrons: NEED TO CONSIDER FACTOR OF 1/4 from averaging! NOT IMPLEMENTED.
//NOTICE: for photons only we can use charge and parity to reduce the amount of construction. we do not use this here.

if(photonZW!=3){
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,p,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,m,qbm,lm,lbp),ind,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type,lo_or_nlo));
}
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,p,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,m,qbm,lm,lbp),ind_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type_b,lo_or_nlo));

if(photonZW!=3){
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,p,qbm,lm,lbp),ind54,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type54,lo_or_nlo));
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,m,qbm,lm,lbp),ind54,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type54,lo_or_nlo));

	VSM->add( A_loop_2q_1g_2l_M2(process(qp,p,qbm,lm,lbp),ind54_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type54_b,lo_or_nlo));
	VSM->add( A_loop_2q_1g_2l_M2(process(qp,m,qbm,lm,lbp),ind54_b,ns,nf,nc,up_down_quark,photonZW,color,tree_color,_ph_type54_b,lo_or_nlo));
}

return VSM;
}


BH_Ampl_2q1g2l::BH_Ampl_2q1g2l(bool up_down_quark, int photonZW,int color, int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q1g2l(mom_assignment,0,
						settings::BH_interface_settings::s_nf,
						settings::BH_interface_settings::s_nc,
						up_down_quark,
						photonZW*settings::BH_interface_settings::s_photon_only,
						color,
						tree_color,
						lo_or_nlo),
			bhi,  // parent
			5,    // NbrExtParticles
			1,	  // NbrPowersOfAlphaS
			2,    // NbrPowersOfAlphaQED
			2,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
		{}

}
