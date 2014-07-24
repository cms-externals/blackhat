/*
 * matrix elements for 2q2Q2P
 *
 *  Created on: Apr, 2009
 *
 */


#include "scheme.h"
#include "assembly.h"
#include "BH_Ampl_processes.h"
#include "cached_OLHA.h"
#include <ctime>
#include "BH_interface_impl.h"

#define _VERBOSE 0

using namespace std;
using BH::CachedOLHA::partial_amplitude_cached;

namespace BH {

//flavors assumed distinctin partial amplitudes
//label of partial amplitudes gives exchanges flavors, i.e. exchanged 
//  	anti-quarks to generate subleading color structures

//expected indices:
// q Qb Q Pb P qb  
partial_amplitude_cached* A_loop_2q_2Q_2P_6_00(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}

	return PA;
}


//expected indices:
// q Qb Q qb P Pb  
partial_amplitude_cached* A_loop_2q_2Q_2P_6_46(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
// not used presently!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}
	return PA;
}


//expected indices:
// q qb Q Pb P Qb  
partial_amplitude_cached* A_loop_2q_2Q_2P_6_26(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);
// leading color only!!
// not used presently!!
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}
	return PA;
}


//expected indices:
// q Pb Q Qb P qb  
partial_amplitude_cached* A_loop_2q_2Q_2P_6_24(process pro,vector<int> & ind, int n_s, int n_f, int n_c, int color,QCDorder lo_or_nlo=nlo){
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

	process pro_123456=process(h1,h2,h3,h4,h5,h6);
	vector<int> ind_123456;
	ind_123456.push_back(i1);
	ind_123456.push_back(i2);
	ind_123456.push_back(i3);
	ind_123456.push_back(i4);
	ind_123456.push_back(i5);
	ind_123456.push_back(i6);

//-------------------------------------------


// renormalization & schemeshift
// notice the implicit minus-sign in the subtraction structure
	multi_precision_fraction r1(11,3), r2(-2*n_f,3*n_c),r3(-n_s,3*n_c);
	multi_precision_fraction nm2_over_2(2,1); // CAREFUL WITH GENERALIZATION TO N-points
        multi_precision_fraction r4(2,3), r5(-1,n_c*n_c);

// leading color only!!
// not used presently
if(color==1){
// SCHEME-SHIFT
        PA->add_subtraction(pro_123456,ind_123456,r4,0);
// RENORMALIZATION
	PA->add_subtraction(pro_123456,ind_123456,nm2_over_2*r1,-1);
// primitive amplitue
	PA->add(pro_123456,leading_color,ind_123456,1,1);
}
	return PA;
}




//expected indices:
// qb q qb q qb q 
Squared_ME* A_loop_2q_2Q_2P_2l_M2(process pro, vector<int> & ind, int n_s, int n_f, int n_c, int case6q,int color, int tree_color,QCDorder lo_or_nlo=nlo){

Squared_ME* SM = new Squared_ME(lo_or_nlo);

int i1=ind.at(0);
int i2=ind.at(1);
int i3=ind.at(2);
int i4=ind.at(3);
int i5=ind.at(4);
int i6=ind.at(5);


// assume input with identical flavor quarks
// introduce distinct flavor quarks in order to use primitive amplitudes with distinct flavors

//3 cases:
// 3 distinct flavors
// 2 distinct flavors
// 3 identical flavors
// built from distinct flavor amplitudes 
// h: (12,34,56)
// h12: (14,32,56)
// h13: (16,34,52)
// h23: (12,36,54)

ph_type h_1=particle_ID(quark,pro.p(2).helicity(),1);
ph_type h_2=particle_ID(quark,pro.p(3).helicity(),2,true);
ph_type h_3=particle_ID(quark,pro.p(4).helicity(),2);
ph_type h_4=particle_ID(quark,pro.p(5).helicity(),3,true);
ph_type h_5=particle_ID(quark,pro.p(6).helicity(),3);
ph_type h_6=particle_ID(quark,pro.p(1).helicity(),1,true);


ph_type h12_1=particle_ID(quark,pro.p(2).helicity(),2);
ph_type h12_2=particle_ID(quark,pro.p(3).helicity(),2,true);
ph_type h12_3=particle_ID(quark,pro.p(4).helicity(),1);
ph_type h12_4=particle_ID(quark,pro.p(5).helicity(),3,true);
ph_type h12_5=particle_ID(quark,pro.p(6).helicity(),3);
ph_type h12_6=particle_ID(quark,pro.p(1).helicity(),1,true);


ph_type h13_1=particle_ID(quark,pro.p(2).helicity(),3);
ph_type h13_2=particle_ID(quark,pro.p(3).helicity(),2,true);
ph_type h13_3=particle_ID(quark,pro.p(4).helicity(),2);
ph_type h13_4=particle_ID(quark,pro.p(5).helicity(),3,true);
ph_type h13_5=particle_ID(quark,pro.p(6).helicity(),1);
ph_type h13_6=particle_ID(quark,pro.p(1).helicity(),1,true);


ph_type h23_1=particle_ID(quark,pro.p(2).helicity(),1);
ph_type h23_2=particle_ID(quark,pro.p(3).helicity(),2,true);
ph_type h23_3=particle_ID(quark,pro.p(4).helicity(),3);
ph_type h23_4=particle_ID(quark,pro.p(5).helicity(),3,true);
ph_type h23_5=particle_ID(quark,pro.p(6).helicity(),2);
ph_type h23_6=particle_ID(quark,pro.p(1).helicity(),1,true);


//-------------------------------------------------------
// tree labels

process pro_123456=process(h_1,h_2,h_3,h_4,h_5,h_6);
vector<int> ind_123456;
ind_123456.push_back(i1);
ind_123456.push_back(i2);
ind_123456.push_back(i3);
ind_123456.push_back(i4);
ind_123456.push_back(i5);
ind_123456.push_back(i6);


process pro_123654=process(h_1,h_2,h_3,h_6,h_5,h_4);
vector<int> ind_123654;
ind_123654.push_back(i1);
ind_123654.push_back(i2);
ind_123654.push_back(i3);
ind_123654.push_back(i6);
ind_123654.push_back(i5);
ind_123654.push_back(i4);

process pro_163452=process(h_1,h_6,h_3,h_4,h_5,h_2);
vector<int> ind_163452;
ind_163452.push_back(i1);
ind_163452.push_back(i6);
ind_163452.push_back(i3);
ind_163452.push_back(i4);
ind_163452.push_back(i5);
ind_163452.push_back(i2);

process pro_143256=process(h_1,h_4,h_3,h_2,h_5,h_6);
vector<int> ind_143256;
ind_143256.push_back(i1);
ind_143256.push_back(i4);
ind_143256.push_back(i3);
ind_143256.push_back(i2);
ind_143256.push_back(i5);
ind_143256.push_back(i6);

//-------------------------------------------


	size_t T_123456;
	size_t T_123654;
	size_t T_163452;
	size_t T_143256;

//-


	size_t L00_123456;
	size_t L00_color_123456;

	size_t L46_123654;
	size_t L26_163452;
	size_t L24_143256;

//------------------

	if(h_1.helicity()+h_6.helicity()==0&& h_2.helicity()+h_3.helicity()==0 &&  h_4.helicity()+h_5.helicity()==0 ){

		T_123456=SM->add(new CTree_with_prefactor(pro_123456,ind_123456));

if(tree_color!=1){

		T_123654=SM->add(new CTree_with_prefactor(pro_123654,ind_123654));
		T_163452=SM->add(new CTree_with_prefactor(pro_163452,ind_163452));
		T_143256=SM->add(new CTree_with_prefactor(pro_143256,ind_143256));

}

if(color!=1){
		L00_123456=SM->add(A_loop_2q_2Q_2P_6_00(pro_123456,ind_123456,n_s,n_f,n_c,0,lo_or_nlo));
}
		L00_color_123456=SM->add(A_loop_2q_2Q_2P_6_00(pro_123456,ind_123456,n_s,n_f,n_c,color,lo_or_nlo));

if(color!=1){

		L26_163452=SM->add(A_loop_2q_2Q_2P_6_26(pro_163452,ind_163452,n_s,n_f,n_c,0,lo_or_nlo));
		L46_123654=SM->add(A_loop_2q_2Q_2P_6_46(pro_123654,ind_123654,n_s,n_f,n_c,0,lo_or_nlo));
		L24_143256=SM->add(A_loop_2q_2Q_2P_6_24(pro_143256,ind_143256,n_s,n_f,n_c,0,lo_or_nlo));

}

	};



//------------------------------------------------------
// constants
multi_precision_fraction nc(n_c,1);
multi_precision_fraction nf(n_f,1);
multi_precision_fraction n2m1(n_c*n_c-1,1);
multi_precision_fraction _over_nc(1,n_c);



multi_precision_fraction n2=nc*nc;
multi_precision_fraction n3=nc*n2;
multi_precision_fraction n4=nc*n3;
multi_precision_fraction one(1,1);
multi_precision_fraction mone(-1,1);
multi_precision_fraction mn2=mone*n2;
multi_precision_fraction n2_times_n2m1=n2*n2m1;


/* notice: problem with integers becoming too large in the numerical n_c->infinity limit
multi_precision_fraction ext(4*n_c*n_c*(n_c*n_c-1),1);
multi_precision_fraction ext_tree(4*n_c*(n_c*n_c-1),1);
*/

//--------------------------------------------------------
// multiplicities of processes



////////////////////////////////
////////////////////////////////
//--------------------------------------------------------
// 3 flavor
{
if(case6q==14
||case6q==15
||case6q==16
||case6q==17
||case6q==18
||case6q==19
||case6q==20
||case6q==21
){
cout<< "new case 3 flavor"<<endl;

	if(h_1.helicity()+h_6.helicity()==0&& h_2.helicity()+h_3.helicity()==0 &&  h_4.helicity()+h_5.helicity()==0 ){

  SM->add_tree(cached_cross_term_md(T_123456,T_123456,n3));
		
		}	

if(tree_color!=1){
_PRINT("subleading tree");

	if(h_1.helicity()+h_6.helicity()==0&& h_2.helicity()+h_3.helicity()==0 &&  h_4.helicity()+h_5.helicity()==0 ){

// SM->add_tree(cached_cross_term_md(T_123456,T_123456,n3));
  SM->add_tree(cached_cross_term_md(T_123456,T_123654,mone*nc));
  SM->add_tree(cached_cross_term_md(T_123456,T_163452,mone*nc));
  SM->add_tree(cached_cross_term_md(T_123456,T_143256,mone*nc));

  SM->add_tree(cached_cross_term_md(T_123654,T_123456,mone*nc));
  SM->add_tree(cached_cross_term_md(T_123654,T_123654,nc));
  SM->add_tree(cached_cross_term_md(T_123654,T_163452,_over_nc));
  SM->add_tree(cached_cross_term_md(T_123654,T_143256,_over_nc));

  SM->add_tree(cached_cross_term_md(T_163452,T_123456,mone*nc));
  SM->add_tree(cached_cross_term_md(T_163452,T_123654,_over_nc));
  SM->add_tree(cached_cross_term_md(T_163452,T_163452,nc));
  SM->add_tree(cached_cross_term_md(T_163452,T_143256,_over_nc));

  SM->add_tree(cached_cross_term_md(T_143256,T_123456,mone*nc));
  SM->add_tree(cached_cross_term_md(T_143256,T_123654,_over_nc));
  SM->add_tree(cached_cross_term_md(T_143256,T_163452,_over_nc));
  SM->add_tree(cached_cross_term_md(T_143256,T_143256,nc));

};
}


	if(h_1.helicity()+h_6.helicity()==0&& h_2.helicity()+h_3.helicity()==0 &&  h_4.helicity()+h_5.helicity()==0 ){

SM->add_loop(cached_cross_term_md(L00_color_123456,T_123456,n4));
		};

if(color!=1){
_PRINT("subleading loop");

	if(h_1.helicity()+h_6.helicity()==0&& h_2.helicity()+h_3.helicity()==0 &&  h_4.helicity()+h_5.helicity()==0 ){

//SM->add_loop(cached_cross_term_md(L00_123456,T_123456,n4));
  SM->add_loop(cached_cross_term_md(L00_123456,T_123654,mn2));
  SM->add_loop(cached_cross_term_md(L00_123456,T_163452,mn2));
  SM->add_loop(cached_cross_term_md(L00_123456,T_143256,mn2));
  SM->add_loop(cached_cross_term_md(L46_123654,T_123456,mn2));
  SM->add_loop(cached_cross_term_md(L46_123654,T_123654,n2));
  SM->add_loop(cached_cross_term_md(L46_123654,T_163452,one));
  SM->add_loop(cached_cross_term_md(L46_123654,T_143256,one));
  SM->add_loop(cached_cross_term_md(L26_163452,T_123456,mn2));
  SM->add_loop(cached_cross_term_md(L26_163452,T_123654,one));
  SM->add_loop(cached_cross_term_md(L26_163452,T_163452,n2));
  SM->add_loop(cached_cross_term_md(L26_163452,T_143256,one));
  SM->add_loop(cached_cross_term_md(L24_143256,T_123456,mn2));
  SM->add_loop(cached_cross_term_md(L24_143256,T_123654,one));
  SM->add_loop(cached_cross_term_md(L24_143256,T_163452,one));
  SM->add_loop(cached_cross_term_md(L24_143256,T_143256,n2));

	};
}


	};


}

return SM;
}

/*
Virtual_SME* vsme_2q2Q2P(int ns,int nf,int nc, int case6q,int color, int tree_color){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(1);
	ind.push_back(2);
	ind.push_back(3);
	ind.push_back(4);
	ind.push_back(5);
	ind.push_back(6);


clock_t before, after;
before=clock();
	
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qp,qbm,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qp,qbm,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qm,qbp,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qp,qbm,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color));
	

	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qm,qbp,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qm,qbp,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qp,qbm,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qm,qbp,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color));

after=clock();

return VSM;
}
*/

Virtual_SME* vsme_2q2Q2P(std::vector<int> indext,int ns,int nf,int nc,int case6q,int color, int tree_color,QCDorder lo_or_nlo=nlo){
	Virtual_SME* VSM= new Virtual_SME();

	vector<int> ind;
	ind.push_back(indext[0]);
	ind.push_back(indext[1]);
	ind.push_back(indext[2]);
	ind.push_back(indext[3]);
	ind.push_back(indext[4]);
	ind.push_back(indext[5]);

clock_t before, after;
before=clock();
	
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qp,qbm,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qp,qbm,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qm,qbp,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qp,qbm,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	

	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qm,qbp,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qp,qbm,qm,qbp,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qp,qbm,qm,qbp),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));
	VSM->add(A_loop_2q_2Q_2P_2l_M2(process(qm,qbp,qm,qbp,qp,qbm),ind,ns,nf,nc,case6q,color,tree_color,lo_or_nlo));

after=clock();

return VSM;
}





BH_Ampl_2q2Q2P::BH_Ampl_2q2Q2P(int case6q,int color,int tree_color,const std::vector<int>& mom_assignment,QCDorder lo_or_nlo,BH_interface_impl* bhi):
	BH_Ampl_impl(
			vsme_2q2Q2P(mom_assignment,
				0,
				settings::BH_interface_settings::s_nf,
				settings::BH_interface_settings::s_nc,
				case6q,
				color,
				tree_color,
				lo_or_nlo),
			bhi,  // parent
			6,    // NbrExtParticles
			4,	  // NbrPowersOfAlphaS
			0,    // NbrPowersOfAlphaQED
			4,   // GeVdim
			1.  // factor
		),
		momenta_assignment(mom_assignment)
		{}

}
