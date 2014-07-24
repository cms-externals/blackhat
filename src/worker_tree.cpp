/*
 * worker_tree.cpp
 *
 *  Created on: 10 Jul 2009
 *      Author: daniel
 */

#include "worker_tree.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include <typeinfo>
#include <cassert>
#include <string>
#include "tree_amp.h"
#include "settings.h"
#include "worker_utils.h"
#include "trees_eval/amplitudes_tree_eval.h"
#include <sys/stat.h>
#include "BHpath.h"
#include "BH_debug.h"
#include "from_file.h"
#include <unistd.h>
#include <locale>
#include "BH_typedefs.h"

template <typename A>
inline void Assert(A assertion)
{
    if( !assertion ) throw BH::BHerror("");
}


using BH::worker::write;
using BH::worker::read_process_from_stream;

using namespace std;


namespace BH {

//defined lower
void add_tree_to_lib(const process& PRO,Rec_Tree& RT,bool firstTime=true);


shift_base::shift_base(std::istream& is) {
	std::string label;
	is >> label;
	assert( label == "sh" );
	is >> d_i;
	is >> d_j;
};

//in rational.cpp
template <class T> size_t massless_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr);
template <class T> size_t massive_i_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr);
template <class T> size_t massive_j_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr);
template <class T> size_t massive_ij_shift_ij(momentum_configuration<T>& mc, vector<int>& ind, int i, int j, size_t P, const complex<T>& Psqr);

template <class T> void massless_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j);
template <class T> void massive_i_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j);
template <class T> void massive_j_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j);
template <class T> void massive_ij_shift_ij_ep(const eval_param<T>& ep, int i, int j, int imass, int jmass, Cmom<T>& shifti, Cmom<T>& shiftj, Cmom<T>& Phat, const momentum<complex<T> >& P, const complex<T>& Psqr, const Cmom<T>*& ref_i, const Cmom<T>*& ref_j);


template <class Pair> massive_shift<Pair>::massive_shift(std::istream& is): shift_base(is){
	int massive_shift_type;
	string title;
	is>> title;
	assert(title=="ms");
	is >> massive_shift_type;
	assert(title=="im");
	is >> _imass;
	assert(title=="jm");
	is >> _jmass;

	switch (massive_shift_type){
	case 0: // Both are massless
		shift_ij=&massless_shift_ij<R>;
		shift_ij_HP=&massless_shift_ij<RHP>;
		shift_ij_VHP=&massless_shift_ij<RVHP>;
		shift_ij_ep=&massless_shift_ij_ep<R>;
		shift_ij_ep_HP=&massless_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massless_shift_ij_ep<RVHP>;

		break;
	case 1: // i is massive and j is massless
		shift_ij=&massive_i_shift_ij<R>;
		shift_ij_HP=&massive_i_shift_ij<RHP>;
		shift_ij_VHP=&massive_i_shift_ij<RVHP>;
		shift_ij_ep=&massive_i_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_i_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_i_shift_ij_ep<RVHP>;

		break;
	case 2: // j is massive and i is massless
		shift_ij=&massive_j_shift_ij<R>;
		shift_ij_HP=&massive_j_shift_ij<RHP>;
		shift_ij_VHP=&massive_j_shift_ij<RVHP>;
		shift_ij_ep=&massive_j_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_j_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_j_shift_ij_ep<RVHP>;

		break;
	case 3: // They are both massive
		shift_ij=&massive_ij_shift_ij<R>;
		shift_ij_HP=&massive_ij_shift_ij<RHP>;
		shift_ij_VHP=&massive_ij_shift_ij<RVHP>;
		shift_ij_ep=&massive_ij_shift_ij_ep<R>;
		shift_ij_ep_HP=&massive_ij_shift_ij_ep<RHP>;
		shift_ij_ep_VHP=&massive_ij_shift_ij_ep<RVHP>;

		break;

	}

#ifdef BH_USE_GMP
	switch (massive_shift_type){
	case 0: // Both are massless
		shift_ij_GMP=&massless_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massless_shift_ij_ep<RGMP>;
		break;
	case 1: // i is massive and j is massless
		shift_ij_GMP=&massive_i_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_i_shift_ij_ep<RGMP>;
		break;
	case 2: // j is massive and i is massless
		shift_ij_GMP=&massive_j_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_j_shift_ij_ep<RGMP>;
		break;
	case 3: // They are both massive
		shift_ij_GMP=&massive_ij_shift_ij<RGMP>;
		shift_ij_ep_GMP=&massive_ij_shift_ij_ep<RGMP>;
		break;
	}
#endif


};

template <class Pair> massive_prop_shift<Pair>::massive_prop_shift(std::istream& is): massive_shift<Pair>(is){
	string title;
	is >> title;
	assert( title == "ml");
	is >> _mass_leg;
}

template <class Pair> template <class T> std::complex<T> massless_shift<Pair>::generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa){
	//Get the eval params for the left and right hand side
	Pa.template get_l_eval_param<T>();
	eval_param<T>& epl(Pa.template get_l_eval_param<T>());
	eval_param<T>& epr(Pa.template get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(ep.p(Pa.left_ind(0))->P());
	epl.set(0,ep.p(Pa.left_ind(0)));
	for(size_t psiter=1;psiter<Pa.maxl-1;psiter++){
		size_t ind=Pa.left_ind(psiter);
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}



	// Create the shifted momenta
	momentum<complex<T> > ij_mom;
	Spinor_to_momentum(ep.p(d_j)->Lt(),ep.p(d_i)->L(),ij_mom);
	const complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));

	// Now compute the on-shell propagator momentum
	ij_mom.mult_by(z);
	ij_mom.add_to(Psum_mom);
	ms.PHat.set_to_U(ij_mom);
	ms.mPHat.set_to(ms.PHat.L(),-ms.PHat.Lt());

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&ms.PHat);
	epl.set(Pa.maxl-1,&ms.mPHat);

	for(size_t kr=1;kr<Pa.maxr;kr++){
		size_t ind=Pa.right_ind(kr);
		epr.set(kr,ep.p(ind));


	}

	// Insert the shifted legs
	ms.shifti.set_to(ep.p(d_i)->L(),ep.p(d_i)->Lt()-z*ep.p(d_j)->Lt());
	ms.shiftj.set_to(ep.p(d_j)->L()+z*ep.p(d_i)->L(),ep.p(d_j)->Lt());
	epr.set(Pa.shifted_ind_i,&ms.shifti);  // The shifted legs are massless
	epl.set(Pa.shifted_ind_j,&ms.shiftj);

	return Psqr;
}

template <class Pair> template <class T> std::complex<T> massive_prop_massless_shift<Pair>::generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa){


	//Get the eval params for the left and right hand side
	Pa.template get_l_eval_param<T>();
	eval_param<T>& epl(Pa.template get_l_eval_param<T>());
	eval_param<T>& epr(Pa.template get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(ep.p(Pa.left_ind(0))->P());
	epl.set(0,ep.p(Pa.left_ind(0)));
	for(size_t psiter=1;psiter<Pa.maxl-1;psiter++){
		size_t ind=Pa.left_ind(psiter);
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}



	// Create the shifted momenta
	momentum<complex<T> > ij_mom=PfLLt(ep.p(d_j)->Lt(),ep.p(d_i)->L());
	Spinor_to_momentum(ep.p(d_j)->Lt(),ep.p(d_i)->L(),ij_mom);
	const complex<T> Psqr=Psum_mom.square()-eval_param<T>::mass2(d_mass_leg);
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));


	// Now compute the on-shell propagator momentum
	ms.PHat=Cmom<T>(Psum_mom+z*ij_mom,_mt_massive);
	ms.mPHat=Cmom<T>(-ms.PHat.P(),_mt_massive);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&ms.PHat);
	epl.set(Pa.maxl-1,&ms.mPHat);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r

	for(size_t kr=1;kr<Pa.maxr;kr++){
		size_t ind=Pa.right_ind(kr);
		epr.set(kr,ep.p(ind));


	}


	ms.shifti=Cmom<T>(ep.p(d_i)->L(),ep.p(d_i)->Lt()-z*ep.p(d_j)->Lt());
	ms.shiftj=Cmom<T>(ep.p(d_j)->L()+z*ep.p(d_i)->L(),ep.p(d_j)->Lt());
	epr.set(Pa.shifted_ind_i,&ms.shifti); // The shifted legs are massless
	epl.set(Pa.shifted_ind_j,&ms.shiftj);

	//Pass down reference vector and masses
	epl.set_ref(ep.ref());
	epr.set_ref(ep.ref());

	return Psqr;

}

template <class Pair> template <class T> std::complex<T> massive_prop_shift<Pair>::generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa){
	//Get the eval params for the left and right hand side
	eval_param<T>& epl(Pa.template get_l_eval_param<T>());
	eval_param<T>& epr(Pa.template get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(ep.p(Pa.left_ind(0))->P());
	epl.set(0,ep.p(Pa.left_ind(0)));
	for(size_t psiter=1;psiter<Pa.maxl-1;psiter++){
		size_t ind=Pa.left_ind(psiter);
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}



	// Create the shifted momenta
	const complex<T> Psqr=Psum_mom.square()-eval_param<T>::mass2(_mass_leg);
    const Cmom<T> *ref_i, *ref_j;
	massive_shift<Pair>::get_shifted_ij(ep,ms.shifti,ms.shiftj,ms.PHat,Psum_mom,Psqr,ref_i,ref_j);
	ms.mPHat=-ms.PHat;

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&ms.PHat);
	epl.set(Pa.maxl-1,&ms.mPHat);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r

	for(size_t kr=1;kr<Pa.maxr;kr++){
		size_t ind=Pa.right_ind(kr);
		epr.set(kr,ep.p(ind));
	}


	epr.set(Pa.shifted_ind_i,&ms.shifti);
	epl.set(Pa.shifted_ind_j,&ms.shiftj);

	//Pass down reference vector and masses
	epl.set_ref(ep.ref());
	epr.set_ref(ep.ref());

	return Psqr;

}

template <class Pair> template <class T> std::complex<T> massive_unshifted_shift<Pair>::generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa){
	//Get the eval params for the left and right hand side
	Pa.template get_l_eval_param<T>();
	eval_param<T>& epl(Pa.template get_l_eval_param<T>());
	eval_param<T>& epr(Pa.template get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(ep.p(Pa.left_ind(0))->P());
	epl.set(0,ep.p(Pa.left_ind(0)));
	for(size_t psiter=1;psiter<Pa.maxl-1;psiter++){
		size_t ind=Pa.left_ind(psiter);
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}



	// Create the shifted momenta
	momentum<complex<T> > ij_mom=PfLLt(ep.p(d_j)->Lt(),ep.p(d_i)->L());
	Spinor_to_momentum(ep.p(d_j)->Lt(),ep.p(d_i)->L(),ij_mom);
	const complex<T> Psqr=Psum_mom.square();
	complex<T> z=-Psqr/(T(2)*(Psum_mom*ij_mom));


	// Now compute the on-shell propagator momentum
	ms.PHat=Cmom<T> (Psum_mom+z*ij_mom);
	ms.mPHat=Cmom<T> (ms.PHat.L(),-ms.PHat.Lt());

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&ms.PHat);
	epl.set(Pa.maxl-1,&ms.mPHat);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r

	for(size_t kr=1;kr<Pa.maxr;kr++){
		size_t ind=Pa.right_ind(kr);
		epr.set(kr,ep.p(ind));
	}


	ms.shifti=Cmom<T>(ep.p(d_i)->L(),ep.p(d_i)->Lt()-z*ep.p(d_j)->Lt());
	ms.shiftj=Cmom<T>(ep.p(d_j)->L()+z*ep.p(d_i)->L(),ep.p(d_j)->Lt());
	epr.set(Pa.shifted_ind_i,&ms.shifti); // The shifted legs are massless
	epl.set(Pa.shifted_ind_j,&ms.shiftj);

	//Pass down reference vector and masses
	epl.set_ref(ep.ref());
	epr.set_ref(ep.ref());

	return Psqr;
}

template <class Pair> template <class T> std::complex<T> massive_shift<Pair>::generate_shift(const eval_param<T>& ep,momenta_struct<T>& ms,Pair& Pa){
	//Get the eval params for the left and right hand side
	Pa.template get_l_eval_param<T>();
	eval_param<T>& epl(Pa.template get_l_eval_param<T>());
	eval_param<T>& epr(Pa.template get_r_eval_param<T>());

	// Compute the shifted legs, shiftBA gives an [i<j shift
	momentum<complex<T> > Psum_mom(ep.p(Pa.left_ind(0))->P());
	epl.set(0,ep.p(Pa.left_ind(0)));
	for(size_t psiter=1;psiter<Pa.maxl-1;psiter++){
		size_t ind=Pa.left_ind(psiter);
		Psum_mom.add_to(ep.p(ind)->P());
		epl.set(psiter,ep.p(ind));
	}



	// Create the shifted momenta
	const complex<T> Psqr=Psum_mom.square();
    const Cmom<T> *ref_i, *ref_j;
	get_shifted_ij(ep,ms.shifti,ms.shiftj,ms.PHat,Psum_mom,Psqr,ref_i,ref_j);
	ms.mPHat=Cmom<T>(ms.PHat.L(),-ms.PHat.Lt());


	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r
	epr.set(0,&ms.PHat);
	epl.set(Pa.maxl-1,&ms.mPHat);

	//Set up the inds for the left and right parts of the recursion including the Phat legs
	// Note we have j in l and i in r

	for(size_t kr=1;kr<Pa.maxr;kr++){
		size_t ind=Pa.right_ind(kr);
		epr.set(kr,ep.p(ind));


	}


	Pa.template get_r_eval_param<T>().set(Pa.shifted_ind_i,&ms.shifti);
	Pa.template get_l_eval_param<T>().set(Pa.shifted_ind_j,&ms.shiftj);

	//Pass down reference vector and masses
	Pa.template get_l_eval_param<T>().set_ref(ep.ref());
	Pa.template get_r_eval_param<T>().set_ref(ep.ref());

	return Psqr;

}



#ifndef BH_PUBLIC

void write(Rec_Pair* RP,std::ostream& os){
	part PART=RP->get_part();

	int cs1=PART.corner_size(1);
	int cs2=PART.corner_size(2);

	os << cs1 << " ";
	for (int j=1;j<cs1;j++){  // left part has last index as propagator
		os << PART.corner_ind(1,j) << " ";
	}
	os << cs2 << " ";
	for (int j=2;j<=cs2;j++){ // right part has first index as propagator
		os << PART.corner_ind(2,j) << " ";
	}
	Rec_Tree* RT;
	bool topLevel=false;
	RT=dynamic_cast<Rec_Tree*>(RP->left());
	assert(RT);
	os << "L ";
	write(RT,os,topLevel);
	RT=dynamic_cast<Rec_Tree*>(RP->right());
	assert(RT);
	os << "R ";
	write(RT,os,topLevel);

	// write shift information

	os  << "shi " <<  RP->get_shifted_ind_i() << " " << RP->get_shifted_ind_j() << " " << "sh " << RP->get_i() << " " << RP->get_j() << " " ;

}

void write(Rec_Pair_massive* RPm,std::ostream& os){
	Rec_Pair* RP=RPm;
	write(RP,os);
	os << "ms " << RPm->get_massive_shift();
	os << "im " << RPm->get_imass();
	os << "jm " << RPm->get_jmass();
}

void write(Rec_Pair_massive_prop_massless_shift* RPm,std::ostream& os){
	Rec_Pair* RP=RPm;
	write(RP,os);
	os << RPm->get_mass_leg();
}


void write(Rec_Pair_massive_prop* RPmp,std::ostream& os){
	Rec_Pair_massive* RPm=RPmp;
	write(RPm,os);
	os << "ml " << RPmp->get_mass_leg();
}


void write(Unknown_Rec_Tree* URT,std::ostream& os){
	int nbr_daughters = URT->nbr_daughters() ;
	os << nbr_daughters << "\n";
	for (int i=1;i<=nbr_daughters;i++){
		Rec_BB* dau= URT->get_daughter(i);
		Rec_Pair* RP=dynamic_cast<Rec_Pair*>(dau);
		Rec_Pair_massive* RPM=dynamic_cast<Rec_Pair_massive*>(dau);
		Rec_Pair_massive_prop* RPMP=dynamic_cast<Rec_Pair_massive_prop*>(dau);
		Rec_Pair_massive_prop_massless_shift* RPMPM=dynamic_cast<Rec_Pair_massive_prop_massless_shift*>(dau);
		Rec_Pair_massive_unshifted* RPMU=dynamic_cast<Rec_Pair_massive_unshifted*>(dau);

		if (RPMU){
			os << "Mu ";
			write(RPMU,os);
		} else if (RPMPM){
			os << "Mpm ";
			write(RPMPM,os);
		} else if (RPMP){
			os << "Mp ";
			write(RPMP,os);
		} else if (RPM){
			os << "M ";
			write(RPM,os);
		} else if (RP){
			os << "m ";
			write(RP,os);
		}
		os << "\n";
	}
}

void write(Known_Rec_Tree* KRT,std::ostream& os){
	BH::worker::write(KRT->get_process(),os);
}


void write(Known_Rec_Tree_offset* KRTO,std::ostream& os){
	vector<particle_ID> pp;
	const process& pro = KRTO->get_process();
	int length=pro.n();
	int offset=KRTO->get_offset();
	for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((offset+(i-2))%length+1)) ; }
	process new_pro(pp);

	BH::worker::write(new_pro,os);
	os << offset << " " << length << " ";
}

void write(Known_Rec_Tree_permutation* KRTP,std::ostream& os){
	const process& pro = KRTP->get_process();
	vector<particle_ID> pp;

	const vector<int>& perm_ind=KRTP->get_perm_ind();
	for (int i=0;i<pro.n();i++){
		pp.push_back(pro.p(perm_ind[i])) ;
	}
	process new_pro(pp);

	BH::worker::write(new_pro,os);
	os << perm_ind.size() << " " ;
	copy(perm_ind.begin(),perm_ind.end(),ostream_iterator<int>(os," "));
}


void write(Rec_Tree* RT,std::ostream& os,bool topLevel){
	Unknown_Rec_Tree* URT=dynamic_cast<Unknown_Rec_Tree*>(RT);
	if ( URT){
		if (topLevel) {
			os << "URT " ;
			write (URT,os);
			return;
		}
		else {
			os << "URTH " ;
			stringstream ss;
			ss<< URT->get_process();
			string processString=ss.str();
			string hash=myHash(processString);
	//		make sure the tree is there
			add_tree_to_lib(URT->get_process(),*URT,false);
			os << URT->get_process().n() << " " << hash << " ";
			return;

		}
	}
	Known_Rec_Tree* KRT=dynamic_cast<Known_Rec_Tree*>(RT);
	if ( KRT){
		os << "KRT " ;
		write (KRT,os);
		return;
	}
	Known_Rec_Tree_offset* KRTO=dynamic_cast<Known_Rec_Tree_offset*>(RT);
	if ( KRTO){
		os << "KRTO " ;
		write (KRTO,os);
		return;
	}
	Known_Rec_Tree_permutation* KRTP=dynamic_cast<Known_Rec_Tree_permutation*>(RT);
	if ( KRTP){
		os << "KRTP " ;
		write (KRTP,os);
		return;
	}
	_WARNING2("Unknown type of Rec_Tree: ",typeid(*RT).name());

}

#endif

template <class T> complex<T> worker_tree_known_offset::eval_fn(momentum_configuration<T>& mc, const vector<int>& ind ){
	vector<int> new_ind=ind;
	rotate_copy(ind.begin(),ind.begin()+(d_offset-1),ind.begin()+d_length,new_ind.begin());
	eval_param<T> ep(mc,new_ind);
	return worker_tree_known::eval(ep);

}

template <class T> complex<T> worker_tree_known_offset::eval_fn(const eval_param<T>& ep){
	eval_param<T> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((d_offset+(i-1))%d_length));
	}
	return worker_tree_known::eval(rotated_ep);
}


indices_struct::indices_struct(std::istream& is,bool putMinusOne=false){
	if (putMinusOne==true){
		indices.push_back(-1);
	}
	is >> size;
	for (int j=1;j<size;j++){  // left part has last index as propagator
		int next;
		is  >> next;
		indices.push_back(next);
	}

};


Tree_Pair_base::Tree_Pair_base(std::istream& is):
		d_left_indices(is,false), d_right_indices(is,true),
		d_eps(d_left_indices.size,d_right_indices.size){


	indshiftl.assign(d_left_indices.size,0);
	d_left_indices.indices.push_back(-1);

	d_right_indices.indices.push_back(-1);
	indshiftr.assign(d_right_indices.size,0);
/* from Darren's code, replace with cs1,cs2 if it is ok*/
	maxl=indshiftl.size();
	maxr=indshiftr.size();
	max=maxl+maxr-2;

	std::string label;

	is >> label;

	assert(label == "L");

	d_left=create_worker_tree(is);

	is >> label;
	assert(label == "R");

	d_right=create_worker_tree(is);

	is >> label;
	Assert(label == "shi");


	is >> shifted_ind_i;
	is >> shifted_ind_j;

	assert(shifted_ind_i < maxr);
	assert(shifted_ind_j < maxl);

	indshiftl.push_back(0);
	indshiftl.push_back(0);
	indshiftr.push_back(0);
	indshiftr.push_back(0);



}

worker_tree_unknown::worker_tree_unknown(std::istream& is){
	is >> d_nbr_pairs;
	assert( d_nbr_pairs >= 0 );
	for (int i=0;i<d_nbr_pairs;i++){
		std::string type;
		is >> type ;
		if ( type == "m"){
			massless_pair* newPair=new massless_pair(is);
			d_pairs.push_back(newPair);
		}
		else if ( type == "M"){
			massive_pair* newPair=new massive_pair(is);
			d_pairs.push_back(newPair);
		}
		else if ( type == "Mp"){
			massive_prop_pair* newPair=new massive_prop_pair(is);
			d_pairs.push_back(newPair);
		}
		else if ( type == "Mpm"){
			massive_prop_massless_pair* newPair=new massive_prop_massless_pair(is);
			d_pairs.push_back(newPair);
		}
		else if ( type == "Mu"){
			massive_unshifted_pair* newPair=new massive_unshifted_pair(is);
			d_pairs.push_back(newPair);
		} else {
			throw BHerror("Syntax error in worker data");
		}
	}
}



worker_tree_known::worker_tree_known(const process& pro){
	init(pro);
}

void worker_tree_known::init(const process& pro){
	if (Tree_is_zero(pro)) {

		d_eval_C_ep_ptr=&ZeroF_eval;
		d_eval_CHP_ep_ptr=&ZeroF_eval;
		d_eval_CVHP_ep_ptr=&ZeroF_eval;

#if BH_USE_GMP
		d_eval_CGMP_ep_ptr=&ZeroF_eval;
#endif

	} else {

		d_eval_C_ep_ptr=A_Tree_Ptr_eval<R>(pro);
		d_eval_CHP_ep_ptr=A_Tree_Ptr_eval<RHP>(pro);
		d_eval_CVHP_ep_ptr=A_Tree_Ptr_eval<RVHP>(pro);

#if BH_USE_GMP
		d_eval_CGMP_ep_ptr=A_Tree_Ptr_eval<RGMP>(pro);
#endif
	}

	 _masses=new mass_param_coll(pro);
		assert(d_eval_C_ep_ptr);
		assert(d_eval_CHP_ep_ptr);
		assert(d_eval_CVHP_ep_ptr);


}

// the known tree is initialized with the rotated process
worker_tree_known_offset::worker_tree_known_offset(std::istream& is): worker_tree_known(read_process_from_stream(is)) {

	is >> d_offset;
	is >> d_length;

}
// the known tree is initialized with the rotated process
worker_tree_known_offset::worker_tree_known_offset(const process& pro,int offset): worker_tree_known(pro), d_offset(offset) {
	d_length=pro.n();
}


worker_tree*  new_known_tree(const process& pro);

worker_tree* create_worker_tree(std::istream& is){
	std::string type;
	is >> type;
	if (type == "URT"){
		return new worker_tree_unknown(is);
	} else	if (type == "URTH"){
		int n;
		string hash;
		is >> n;
		is >> hash;
		stringstream ss;
		ss<< get_worker_dir("trees/") << n << "/tree_" << hash << ".dat";
		string path=ss.str();
		ifstream is2;
		is2.open(path.c_str(),ios::in);
		if (!is2.is_open()){throw BHerror("File not found");}
		string title; is2 >> title;
		assert( title == "URT");
		return new worker_tree_unknown(is2);
	} else if ( type == "KRT"){
		process pro;
		read_process_from_stream(pro,is);
		return new_known_tree(pro);
	} else if ( type == "KRTO"){
		return new worker_tree_known_offset(is);
	} else {
		_WARNING3("Unknown type : \"",type,"\" in worker_tree_unknown::worker_tree_unknown.");
	}
	throw BHerror("Syntax error");
}



template <class T> std::complex<T> worker_tree_unknown::eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){
//	std::complex<T> res(0,0);
//	for (int i=0;i<d_nbr_pairs;i++){
//		res+=d_pairs[i]->eval(mc,ind);
//	}
	eval_param<T> ep(mc,ind);
	return eval(ep);
};

template <class T> std::complex<T> worker_tree_unknown::eval_fn(const eval_param<T>& ep){
	static int depth=0;
	std::complex<T> res(0,0);
	depth++;
	int showme=depth; // needed to make it visible in gdb
	for (int i=0;i<d_nbr_pairs;i++){
		res+=d_pairs[i]->eval(ep);
	}
	depth--;
	return res;
};

template <template <class> class ShiftType> template <class T> complex<T> Tree_Pair<ShiftType>::eval_fn(momentum_configuration<T>& mc,const std::vector<int>& ind){
	eval_param<T> ep(mc,ind);
	return eval(ep);
};

template <template <class> class ShiftType> template <class T> complex<T> Tree_Pair<ShiftType>::eval_fn(const eval_param<T>& ep){


	momenta_struct<T> momenta;

	complex<T> Psqr = ShiftType<Tree_Pair_base>::generate_shift(ep,momenta,*this);
	complex<T> A1=eval_left(this->template get_l_eval_param<T>());
	complex<T> A2=eval_right(this->template get_r_eval_param<T>());
	return Amp_safe(complex<T>(0.0,-1.0)*(A1*A2)/Psqr);

};


std::string string_name(const process& pro);  // in Rec_Tree


std::string tree_filename(const process& pro){
	stringstream ss;
	struct stat st;
	std::string dir=get_worker_dir("trees");
	ss<< dir ;
	ss << "/" << pro.n() ;
	if ( access( ss.str().c_str(), 0 ) != 0 ){
		_WARNING3("Data path ",ss.str(),"not present. Please create it. ");
		throw BHerror("Missing path");
	}
	ss << "/tree_";
	stringstream pros;
	pros << pro;
	std::string prostring=pros.str();
	string hash(myHash(prostring));
	ss << hash ;
	//ss << pros.str();
	ss << ".dat";
	return ss.str();
}

struct getMCs {


mom_conf mc4,mc5,mc6,mc7,mc8,mc9;
getMCs(){
Cmom<R> pformomconstruction1,pformomconstruction2,pformomconstruction3,pformomconstruction4,pformomconstruction5,pformomconstruction6,pformomconstruction7,pformomconstruction8,pformomconstruction9,pformomconstruction10,pformomconstruction11,pformomconstruction12;
pformomconstruction1 = Cmom<R>(2.29333281086030751,-0.88305370744680835,-1.96002342239748243,-0.79868624301796119);
pformomconstruction2 = momC(2.22100900105464863,0.83839402809274431,1.97192484347455149,0.584370471629133949);
pformomconstruction3 = momC(-2.20789098262258162,1.64481285677059797,0.82065855183763165,-1.2230669640882361);
pformomconstruction4 = momC(-2.30645082929237452,-1.60015317741653393,-0.83255997291470071,1.43738273547706334);
mc4.insert(pformomconstruction1);
mc4.insert(pformomconstruction2);
mc4.insert(pformomconstruction3);
mc4.insert(pformomconstruction4);
pformomconstruction1 = momC(1.86058205404180962,0.22322341633486098,1.29381396176358134,1.31832557381242467);
pformomconstruction2 = momC(2.22702886899219324,-1.97843090981622047,0.18090821896703733,1.00635030417771722);
pformomconstruction3 = momC(1.77552105445334353,-1.61888957780145734,0.46152052119123634,-0.56450895317284531);
pformomconstruction4 = momC(-4.59937300385497873,4.17419448896838902,-0.96173464232271022,-1.67493249852413707);
pformomconstruction5 = momC(-1.26375897363236766,-0.80009741768557219,-0.97450805959914479,-0.08523442629315951);
mc5.insert(pformomconstruction1);
mc5.insert(pformomconstruction2);
mc5.insert(pformomconstruction3);
mc5.insert(pformomconstruction4);
mc5.insert(pformomconstruction5);

pformomconstruction1 = momC(1.92955137629582831,0.12176375318546036,1.35203710358853989,1.37125408757648828);
pformomconstruction2 = momC(2.3718946178199526,-2.10154141143786404,0.25155580416061273,1.07057342179232001);
pformomconstruction3 = momC(1.86478515042509431,-1.71633845344031231,0.51744205227939982,-0.51367274895690288);
pformomconstruction4 = momC(-4.87600294030279301,3.29796758789333694,-3.20053955542535968,-1.62952785100596814);
pformomconstruction5 = momC(-0.54186109016474053,-0.06227089168380983,0.5327311409203322,-0.07702797269268775);
pformomconstruction6 = momC(-0.74836711407334168,0.46041941548318888,0.54677345447647505,-0.22159893671324953);
mc6.insert(pformomconstruction1);
mc6.insert(pformomconstruction2);
mc6.insert(pformomconstruction3);
mc6.insert(pformomconstruction4);
mc6.insert(pformomconstruction5);
mc6.insert(pformomconstruction6);
pformomconstruction1 = momC(1.37845962162848452,-0.899604391289195437,-0.76146684945939796,-0.71486439609741331);
pformomconstruction2 = momC(1.31550521424004713,0.86154333881518911,0.98581179913844106,0.128343837855031049);
pformomconstruction3 = momC(-0.483429888680200544,0.286893255432157653,-0.277486867061210162,-0.272759520210466919);
pformomconstruction4 = momC(-0.321701443648813055,0.202997164217524277,-0.246144284682511691,-0.041194189924169683);
pformomconstruction5 = momC(-0.471025810729965125,-0.286107020036279202,0.365488906741813979,-0.080162001649501396);
pformomconstruction6 = momC(-0.337726474245189617,0.258662100312875037,-0.217145127094240915,-0.00104069546211045);
pformomconstruction7 = momC(-1.080081218564363314,-0.42438444745227144,0.150942422417105687,0.981676965488630708);
mc7.insert(pformomconstruction1);
mc7.insert(pformomconstruction2);
mc7.insert(pformomconstruction3);
mc7.insert(pformomconstruction4);
mc7.insert(pformomconstruction5);
mc7.insert(pformomconstruction6);
mc7.insert(pformomconstruction7);
pformomconstruction1 = momC(0.677400017074605876,0.250142496747423391,0.218529533647501042,-0.590376453949043703);
pformomconstruction2 = momC(0.890448994420059572,-0.84042287009358795,-0.169878789888167917,-0.240270696992960786);
pformomconstruction3 = momC(0.319059059309799145,0.282169596702112851,-0.132095668503444812,-0.068773078942627912);
pformomconstruction4 = momC(0.591895986238114489,0.075279124768095475,-0.329199014134954976,0.486108959999209894);
pformomconstruction5 = momC(0.695489395606314995,0.128490939776938918,0.58370742287648654,0.355641986096927112);
pformomconstruction6 = momC(0.449781677244689127,-0.067132947945870905,0.07422384062513292,-0.438506038690458123);
pformomconstruction7 = momC(-2.04102697550486404,0.10089801809451397,-1.68939376079989505,1.14085907352903102);
pformomconstruction8 = momC(-1.58304815438871916,0.07057564195037425,1.44410643617734225,-0.6446837510500775);
mc8.insert(pformomconstruction1);
mc8.insert(pformomconstruction2);
mc8.insert(pformomconstruction3);
mc8.insert(pformomconstruction4);
mc8.insert(pformomconstruction5);
mc8.insert(pformomconstruction6);
mc8.insert(pformomconstruction7);
mc8.insert(pformomconstruction8);
pformomconstruction1 = momC(0.724029539784378082,0.232674459973560293,0.24351697489968873,-0.640921877526500715);
pformomconstruction2 = momC(0.925896021004126494,-0.86306259634180863,-0.13749342014590363,-0.305780895288754375);
pformomconstruction3 = momC(0.313397289344500322,0.274286381522237675,-0.120818993354682674,-0.091583910839617861);
pformomconstruction4 = momC(0.545232501029269833,0.061105453223224026,-0.308924051948866364,0.44509609512161666);
pformomconstruction5 = momC(0.689905146431366202,0.111222770507794176,0.608408960360722244,0.305674898128453623);
pformomconstruction6 = momC(0.487395926677445987,-0.078814343395845723,0.090933685373488494,-0.472307265962787639);
pformomconstruction7 = momC(-1.81586645851486555,0.31768258179445032,-1.00852443447653242,1.47625446229554976);
pformomconstruction8 = momC(-1.24170537476965713,0.08802623630309822,1.16631481734921562,-0.41688531550755247);
pformomconstruction9 = momC(-0.62828459098656425,-0.14312094358671035,-0.53341353805713,-0.29954619042040699);
mc9.insert(pformomconstruction1);
mc9.insert(pformomconstruction2);
mc9.insert(pformomconstruction3);
mc9.insert(pformomconstruction4);
mc9.insert(pformomconstruction5);
mc9.insert(pformomconstruction6);
mc9.insert(pformomconstruction7);
mc9.insert(pformomconstruction8);
mc9.insert(pformomconstruction9);
 }
};

momentum_configuration<R>& get_mc(int n){

	static getMCs mcs;
	switch(n){
	case 4: return mcs.mc4;
	case 5: return mcs.mc5;
	case 6: return mcs.mc6;
	case 7: return mcs.mc7;
	case 8: return mcs.mc8;
	case 9: return mcs.mc9;
	}
}


void add_tree_to_lib(const process& PRO,Rec_Tree& RT,bool firstTime){
#ifndef BH_PUBLIC
	string filename=tree_filename(PRO);
		ifstream is;
		is.open(filename.c_str(),ios::in);
		if( is.is_open() ) {
			if (firstTime) {
				_MESSAGE3("File ", filename ," already present. Hash collision!.");
				throw BHerror("Internal error: hash collision in the generation of the tree worker files");
			} else {
				return;
			}
		} else {
			ofstream file;
			file.open(filename.c_str());
			_MESSAGE2("Creating tree file for ",PRO);
			write(&RT,file);
			file.close();
			stringstream ss;
			std::string dir=get_worker_dir("trees");
			ss<< dir ;
			ss << "/" << PRO.n() ;
			ss << "/tree_list.dat";
			ofstream list;
			list.open(ss.str().c_str(),ios::app);
			stringstream pros;
			pros << PRO;
			std::string prostring=pros.str();
			string hash(myHash(prostring));
			list << hash << " " << PRO << std::endl;
			list.close();
			ifstream itest;
			itest.open(filename.c_str(),ios::in);
			string title;
			itest >> title;
			assert(title=="URT");
			worker_tree_unknown wt(itest);

			mom_conf& mc=get_mc(PRO.n());
			int inds[]={1,2,3,4,5,6,7,8,9};
			std::vector<int>  ind(&inds[0],&inds[0] + PRO.n());
			eval_param<R> ep(mc,ind);
			C rwt=wt.eval(ep);
			C rrt=RT.eval(ep);
			assert(rwt==rrt);
		}


#else
		throw BHerror("Not possible in public version!");
#endif
}

worker_tree* new_known_tree(const process& pro)
{
//	static int n=0;
	if (Tree_is_zero(pro)) {
//		_MESSAGE3(++n," new known tree for ",pro);
		return new worker_tree_known(pro);}
	switch(pro.pcode()){

	// cases where all amplitudes are known:
	case odd_nbr_q:	case odd_nbr_q2: case odd_nbr_l: case odd_nbr_l1: case odd_nbr_l2:
		case 3:
#if _TESTING_MODE == 0
		case 4:	case 5: case 6:   // n gluon amplitudes
#endif
		case 21:
//			_MESSAGE3(++n," new known tree for ",pro);
			return new worker_tree_known(pro);
		break;




#if _TESTING_MODE == 1
		case 7: case 8: // n gluon amplitudes
		case 22: case 23: case 24:  // 2 quark (n-2) gluons amplitudes
		case 221: case 222: case 223: case 224:  // 2 quark 2 leptons (n-4) gluons amplitudes
		case 41: case 42: case 43: {// 4 quark (n-4) gluons amplitudes
			return new Unknown_Rec_Tree(pro);
		break;
		}
#endif


//		case 220:
//		{
//		if( A_Tree_Ptr_eval<R>(pro) != 0){
//			_MESSAGE2("new known tree for ",pro);
//
//			return new worker_tree_known(pro);
//
//		}
//		else {
//			vector<int> new_ind;
//			process new_pro=order_2q2e(pro,new_ind);
//			if( A_Tree_Ptr_eval<R>(new_pro) != 0){
//				_MESSAGE2("new known tree per for ",pro);
//				return new Known_Rec_Tree_permutation(pro,new_ind);
//			}
//		}
//		} break;
		default: {
			if( A_Tree_Ptr_eval<R>(pro) != 0){
//				_MESSAGE3(++n," new known tree for ",pro);

				return new worker_tree_known(pro);
			 }
			for (int j=1;j<pro.n();j++){
				vector<particle_ID> pp;
				rotate_copy(pro.begin(),pro.begin()+j,pro.end(),back_inserter(pp));
				//for (int i=1;i<=pro.n();i++){ pp.push_back(pro.p((pro.n()+j+(i-2))%pro.n()+1)) ; }
				process new_pro(pp);
				 if( A_Tree_Ptr_eval<R>(new_pro) != 0){
//						_MESSAGE3(++n,"new known tree offset for ",new_pro);
						return new worker_tree_known_offset(new_pro,j+1);
				 }
			}
		} break;


		}
	return 0;
}

worker_tree* worker_tree_factory::new_tree(const process& PRO){
//	_MESSAGE2("new tree for ",PRO);
#ifndef BH_PUBLIC
	worker_tree* kwt=new_known_tree(fix_flavors(PRO));
#else
	worker_tree* kwt=new_known_tree(PRO);
#endif
	if ( kwt ){
		return kwt;
	}

	string filename=tree_filename(PRO);
		ifstream ifile;
		ifile.open(filename.c_str(),ios::in);
		if( ifile.is_open() ) {
			try {
				worker_tree* WT=create_worker_tree(ifile);
				return WT;
			} catch ( ... ) {
				_WARNING3("\nerror reading ", filename, "in worker_tree_factory::new_tree\n");
				throw BHerror("Syntax error in worker data");
			}
		} else {
#ifndef BH_PUBLIC
			TREE_FACTORY_TYPE TF;
			_MESSAGE3("Process ",PRO," not in the database, creating it now...");
			Rec_Tree* RT = TF.new_tree(PRO);
			add_tree_to_lib(PRO,*RT);
			_MESSAGE("Done.");
			delete RT;
			ifstream ifile;
			ifile.open(filename.c_str(),ios::in);
			if( ifile.is_open() ) {
				worker_tree* WT=create_worker_tree(ifile);
				return WT;
			} else {
				_WARNING("Amplitude not present, even though we just created it.");
				throw BHerror("BH inconsistency");
			}
#else
			throw BHerror("Missing tree in library!");
#endif
		}
}


}

