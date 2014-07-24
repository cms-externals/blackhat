/*
 * rec_tree_eval.cpp
 *
 *  Created on: Apr 7, 2009
 *      Author: dforde
 */

#ifndef BH_PUBLIC
#include "rec_tree.h"

#include "rational.h"
#include "BH_utilities.h"
#include "trees_eval/amplitudes_tree_eval.h"
#include "cut_part_factory.h"
#include "process_utils.h"
#include "process_utils_massive.h"
#include "tree_amp.h"
#include <stdlib.h>


using namespace std;

namespace BH {


/*
 *
 *
 * Code for known trees
 *
 *
 */

complex<R> Known_Rec_Tree_offset::eval(const eval_param<R>& ep){
	eval_param<R> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%_length));
	}
	return (*_eval_C_ep_ptr)(rotated_ep,*_masses);
}
complex<RHP> Known_Rec_Tree_offset::eval(const eval_param<RHP>& ep){
	eval_param<RHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%_length));
	}
	return (*_eval_CHP_ep_ptr)(rotated_ep,*_masses);
}
complex<RVHP> Known_Rec_Tree_offset::eval(const eval_param<RVHP>& ep){
	eval_param<RVHP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%_length));
	}
	return (*_eval_CVHP_ep_ptr)(rotated_ep,*_masses);
}

#if BH_USE_GMP
complex<RGMP> Known_Rec_Tree_offset::eval(const eval_param<RGMP>& ep){
	eval_param<RGMP> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(i,ep.p((_offset+(i-1))%_length));
	}
	return (*_eval_CGMP_ep_ptr)(rotated_ep,*_masses);
}
#endif

complex<R> Known_Rec_Tree_permutation::eval(const eval_param<R>& ep){
	eval_param<R> rotated_ep(ep.size());
	for(int i=0;i<ep.size();i++){
		rotated_ep.set(d_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_C_ep_ptr)(rotated_ep,*_masses);
}
complex<RHP> Known_Rec_Tree_permutation::eval(const eval_param<RHP>& ep){
	eval_param<RHP> rotated_ep(ep.size());
	for (int i=0;i<ep.size();i++){
		rotated_ep.set(d_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CHP_ep_ptr)(rotated_ep,*_masses);

}
complex<RVHP> Known_Rec_Tree_permutation::eval(const eval_param<RVHP>& ep){
	eval_param<RVHP> rotated_ep(ep.size());
	for (int i=0;i<ep.size();i++){
		rotated_ep.set(d_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CVHP_ep_ptr)(rotated_ep,*_masses);
}

#if BH_USE_GMP
complex<RGMP> Known_Rec_Tree_permutation::eval(const eval_param<RGMP>& ep){
	eval_param<RGMP> rotated_ep(ep.size());
	for (int i=0;i<ep.size();i++){
		rotated_ep.set(d_perm_ind[i]-1,ep.p(i));
	}
	return (*_eval_CGMP_ep_ptr)(rotated_ep,*_masses);
}
#endif


/*
 *
 *
 * Code for Unknown_Rec_Tree
 *
 *
 */


C Unknown_Rec_Tree::eval(const eval_param<R>& ep){
	static int depth=0;
	depth++;
	C res(0,0);
	int showme=depth; // needed to make it visible in gdb
	for (int i=0;i<daughters.size();i++){
		res+=daughters[i]->eval(ep);
	}
	depth--;
	return res;
}
CHP Unknown_Rec_Tree::eval(const eval_param<RHP>& ep){
	CHP res(0,0);
	for (int i=0;i<daughters.size();i++){
		res+=daughters[i]->eval(ep);
	}
	return res;
}
CVHP Unknown_Rec_Tree::eval(const eval_param<RVHP>& ep){
	CVHP res(0,0);
	for (int i=0;i<daughters.size();i++){
		res+=daughters[i]->eval(ep);
	}
	return res;
}

#if BH_USE_GMP
CGMP Unknown_Rec_Tree::eval(const eval_param<RGMP>& ep){
	CGMP res(0,0);
	for (int i=0;i<daughters.size();i++){
		res+=daughters[i]->eval(ep);
	}
	return res;
}
#endif

}
#endif
