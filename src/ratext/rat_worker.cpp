/*
 * rat_worker.cpp
 *
 *  Created on: 16-May-2009
 *      Author: daniel
 */

#include "ratext/rat_worker.h"
#include <iostream>
#include <cassert>
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "worker_utils.h"   // for read_process_from_stream

#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif

using namespace std;

namespace BH {

namespace ratext {


using BH::worker::read_process_from_stream;

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};

rat_worker::rat_worker(std::istream& is){
    TREE_FACTORY_TYPE TF;
    int nbr_corners;
    is >> nbr_corners;
    assert(nbr_corners > 1);
    assert(nbr_corners < 6);
    d_corners.reserve(nbr_corners);
    for (int i=0;i<nbr_corners;i++){
            d_corners.push_back(vector<int>());
            int nbr_entries;
            is >> nbr_entries;
            for (int j=0;j<nbr_entries;j++){
                    int entry;
                    is >> entry ;
                    assert(entry);
                    d_corners[i].push_back(entry);
            }
    }



	int nbr_decendants;

	is >> nbr_decendants;
	for (int i=0;i<nbr_decendants;i++){
		d_trees.push_back(vector<TREE_TYPE*>());
		for (int j=1;j<=nbr_corners;j++){
			process PRO;
			read_process_from_stream(PRO,is);
			d_trees[i].push_back(TF.new_tree(PRO));
		}
	}
	string title;
	is >> title;
	assert(title=="SF");
	int den,num;
	is >> num; is >> den;
	d_symmetry_factor=symmetry_factor(num,den);
};

std::ostream&  operator<<(std::ostream& os,const rat_worker& rw){
	switch(rw.corner_nbr()){
	case 2: 	return os<< "bubble rat_worker";
	case 3: 	return os<< "triangle rat_worker";
	case 4: 	return os<< "box rat_worker";
	case 5: 	return os<< "pentagon rat_worker";

	}
};


rat_worker::~rat_worker(){
	for (int i=0;i<d_trees.size();i++){
		for_each(d_trees[i].begin(),d_trees[i].end(),do_delete<TREE_TYPE>());
	}
}

std::complex<R> rat_worker::eval_tree(int n,int i,momentum_configuration<R>& mc, const std::vector<int>& ind){
	return (d_trees[n])[i]->eval(mc,ind);
}
std::complex<RHP> rat_worker::eval_tree(int n,int i,momentum_configuration<RHP>& mc, const std::vector<int>& ind){
	return (d_trees[n])[i]->eval(mc,ind);
}
std::complex<RVHP> rat_worker::eval_tree(int n,int i,momentum_configuration<RVHP>& mc, const std::vector<int>& ind){
	return (d_trees[n])[i]->eval(mc,ind);
}

std::complex<R> rat_worker::eval_tree(int descendant_n,int tree_n,const eval_param<R>& ep){
	return (d_trees[descendant_n])[tree_n]->eval(ep);
}
std::complex<RHP> rat_worker::eval_tree(int descendant_n,int tree_n,const eval_param<RHP>& ep){
	return (d_trees[descendant_n])[tree_n]->eval(ep);
}
std::complex<RVHP> rat_worker::eval_tree(int descendant_n,int tree_n,const eval_param<RVHP>& ep){
	return (d_trees[descendant_n])[tree_n]->eval(ep);
}

#if BH_USE_GMP
std::complex<RGMP> rat_worker::eval_tree(int n,int i,momentum_configuration<RGMP>& mc, const std::vector<int>& ind){
	return (d_trees[n])[i]->eval(mc,ind);
}
std::complex<RGMP> rat_worker::eval_tree(int descendant_n,int tree_n,const eval_param<RGMP>& ep){
	return (d_trees[descendant_n])[tree_n]->eval(ep);
}
#endif
}
}
