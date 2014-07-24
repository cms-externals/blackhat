/*
 * cut_part_base.cpp
 *
 *  Created on: 05.07.2012
 *      Author: daniel
 */


#include "cut_part_base.h"

#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include<vector>
#include<complex>

#include "BH_typedefs.h"
#include "amplitudes.h"
#include "BH_A0.h"
#include "rational.h"
#include "mode_dependent_typedefs.h"
#include "process_utils_basic.h"

#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif

#define _VERBOSE 0 // 0 is silent, 1 gives some information

using namespace std;

namespace BH {


Cut_Part_base::~Cut_Part_base(){
	delete d_tree_ptr;
}

Cut_Part_base::Cut_Part_base(const process& pro):
		HelAmpl(pro) ,
#if BH_USE_GMP
		m_MU(R(1),RHP(1),RVHP(1),RGMP(1)),
		d_mu_index(0) , d_mu_index_HP(0) , d_mu_index_VHP(0), d_mu_index_GMP(0)
#else
		m_MU(R(1),RHP(1),RVHP(1)),
		d_mu_index(0) , d_mu_index_HP(0) , d_mu_index_VHP(0)
#endif
		{
	TREE_FACTORY_TYPE TF;
	d_tree_ptr=TF.new_tree(fix_flavors(pro));

	for(int i=1;i<=d_process.n();i++){
		_masses.push_back(d_process.p(i).mass_label());
	}
}

template <class T> complex<T> Cut_Part_base::get_tree(momentum_configuration<T>& mc, const vector<int>& ind){return d_tree_ptr->get_value(mc,ind);}

template  complex<R> Cut_Part_base::get_tree(momentum_configuration<R>& mc, const vector<int>& ind);
template  complex<RHP> Cut_Part_base::get_tree(momentum_configuration<RHP>& mc, const vector<int>& ind);
template  complex<RVHP> Cut_Part_base::get_tree(momentum_configuration<RVHP>& mc, const vector<int>& ind);
#if BH_USE_GMP
template  complex<RGMP> Cut_Part_base::get_tree(momentum_configuration<RGMP>& mc, const vector<int>& ind);
#endif

template <class T> complex<T> Cut_Part_base::get_tree(const eval_param<T>& ep){
	// no mecanism yet to use caching, need to recompute...
	return d_tree_ptr->eval(ep);
	//return d_tree_ptr->get_value(ep);
}

template  complex<R> Cut_Part_base::get_tree(const eval_param<R>& ep);
template  complex<RHP> Cut_Part_base::get_tree(const eval_param<RHP>& ep);
template  complex<RVHP> Cut_Part_base::get_tree(const eval_param<RVHP>& ep);

#if BH_USE_GMP
template  complex<RGMP> Cut_Part_base::get_tree(const eval_param<RGMP>& ep);
#endif


}
