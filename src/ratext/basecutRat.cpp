/*
 * basecutRat.cpp
 *
 *  Created on: 7 Jul 2009
 *      Author: daniel
 */

#include "partitions.h"
#include "BH_utilities.h"
#include "rec_tree.h"
#include "rat_ext.h"
#include "ratext/basecutRat.h"
#include "settings.h"
#include "eval_param.h"
#include "mode_dependent_typedefs.h"
#include "BH_debug.h"
#include <cassert>
#include <set>
#ifdef BH_PUBLIC
#include "worker_tree.h"
#endif

using namespace std;

namespace BH {

#ifndef BH_PUBLIC

/*
 *
 *
 * base bubble, triangle, box and pentagon code
 *
 *
 */

bool baseCutD::is_match(cutD* cmpr){
	//If the number of corners is not the same then return false
	if(cmpr->nc()!=nc()) return false;

	for(size_t cor=1;cor<=nc();cor++){
		// If the corners are not of the same size return false
		if(cmpr->c(cor).size()!=c(cor).size()) return false;

		//If all the processes don't match then return false
		set<int> s1,s2;
		transform(c(cor).begin(),c(cor).end(),inserter(s1,s1.begin()),mem_fun_ref(&plabel::ind));
		transform(cmpr->c(cor).begin(),cmpr->c(cor).end(),inserter(s2,s2.begin()),mem_fun_ref(&plabel::ind));
		if(s1!=s2) return false;

		//If all the indexes don't match then return false
		vector<int> inds1,inds2;
		transform(c(cor).begin(),c(cor).end(),back_inserter(inds1),mem_fun_ref(&plabel::ind));
		transform(cmpr->c(cor).begin(),cmpr->c(cor).end(),back_inserter(inds2),mem_fun_ref(&plabel::ind));
		if(inds1!=inds2) return false;
	}

	// They must be the same so return true
	return true;
}
	
bool baseCutD::is_match(cutD* cmpr,size_t iclass){
	//If the class info is not the same then this is not the same
	if(iclass!=_iclass) return false;
	
	return is_match(cmpr);
}

baseCutD::~baseCutD(){
	for(int i=0;i<d_trees.size();i++){
		delete d_trees[i];
	}
}

std::complex<R> baseCutD::eval_tree(int descendant_n,int tree_n,momentum_configuration<R>& mc, const std::vector<int>& ind){
	return d_trees[(descendant_n*_width)+tree_n]->eval(mc,ind);
}
std::complex<RHP> baseCutD::eval_tree(int descendant_n,int tree_n,momentum_configuration<RHP>& mc, const std::vector<int>& ind){
	return d_trees[(descendant_n*_width)+tree_n]->eval(mc,ind);
}
std::complex<RVHP> baseCutD::eval_tree(int descendant_n,int tree_n,momentum_configuration<RVHP>& mc, const std::vector<int>& ind){
	return d_trees[(descendant_n*_width)+tree_n]->eval(mc,ind);
}

std::complex<R> baseCutD::eval_tree(int descendant_n,int tree_n,const eval_param<R>& ep){
#ifndef BH_PUBLIC
	BH_DEBUG_MESSAGE3(d_trees[(descendant_n*_width)+tree_n]->get_process()," dbl: ",d_trees[(descendant_n*_width)+tree_n]->eval(ep));
	BH_DEBUG_MESSAGE(d_trees[(descendant_n*_width)+tree_n]->get_process().p(1).mass_label());

	BH_DEBUG_MESSAGE(ep);
#endif
	return d_trees[(descendant_n*_width)+tree_n]->eval(ep);
}
std::complex<RHP> baseCutD::eval_tree(int descendant_n,int tree_n,const eval_param<RHP>& ep){
#ifndef BH_PUBLIC
	BH_DEBUG_MESSAGE3(d_trees[(descendant_n*_width)+tree_n]->get_process(),"  HP: ",d_trees[(descendant_n*_width)+tree_n]->eval(ep));
#endif
	return d_trees[(descendant_n*_width)+tree_n]->eval(ep);
}
std::complex<RVHP> baseCutD::eval_tree(int descendant_n,int tree_n,const eval_param<RVHP>& ep){
#ifndef BH_PUBLIC
	BH_DEBUG_MESSAGE3(d_trees[(descendant_n*_width)+tree_n]->get_process()," VHP: ",d_trees[(descendant_n*_width)+tree_n]->eval(ep));
#endif
	return d_trees[(descendant_n*_width)+tree_n]->eval(ep);
}

#if BH_USE_GMP
std::complex<RGMP> baseCutD::eval_tree(int descendant_n,int tree_n,momentum_configuration<RGMP>& mc, const std::vector<int>& ind){
	return d_trees[(descendant_n*_width)+tree_n]->eval(mc,ind);
}
std::complex<RGMP> baseCutD::eval_tree(int descendant_n,int tree_n,const eval_param<RGMP>& ep){
	return d_trees[(descendant_n*_width)+tree_n]->eval(ep);
}

#endif

/*
 *
 *
 * The base rational bubble classes
 *
 *
 */

basebubbleRat::basebubbleRat(const bubbleRat_comp* pcd) : baseCutD(pcd)
{
	// Store only the unique mass_param's for the cut legs
	_decendantbubbleRat.push_back(pcd);

	_width=2; // Set this up for a bubble
	// Add the new trees of the descendant bubble
#if __USE_NEW_REC_TREE_REC_RAT==1
#ifndef BH_PUBLIC
    BH_DEBUG_MESSAGE2("Constructing Bub (new) ",*this);
#endif
    Rec_tree_eval_factory TF;
#else
#ifndef BH_PUBLIC
    BH_DEBUG_MESSAGE2("Constructing Bub (old) ",*this);
#endif
    __RAT_TREE_TYPE TF;
#endif
    
    for(int i=0;i<2;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE3(i," : ",*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE5("Adding bubble at ",_decendantbubbleRat.size()," with ",pcd->get_process(1),pcd->get_process(2));
}

void basebubbleRat::add(const bubbleRat_comp* pcd)
{
	_decendantbubbleRat.push_back(pcd);

	// Add the new trees of the descendant bubble
//	_MESSAGE2("Constructing ",*this);
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<2;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE5("Adding bubble at ",_decendantbubbleRat.size()," with ",pcd->get_process(1),pcd->get_process(2));
}


basetriangleRat::basetriangleRat(const triangleRat_comp* pcd) : baseCutD(pcd)
{
	_decendanttriangleRat.push_back(pcd);

	_width=3; // Set this up for a triangle
	// Add the new trees of the descendant triangle
#ifndef BH_PUBLIC
	BH_DEBUG_MESSAGE2("Constructing Tri ",*this);
#endif
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<3;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE6("Adding triangle at ",_decendanttriangleRat.size()," with ",pcd->get_process(1),pcd->get_process(2),pcd->get_process(3));
}

void basetriangleRat::add(const triangleRat_comp* pcd)
{
	_decendanttriangleRat.push_back(pcd);

	// Add the new trees of the descendant triangle
//	_MESSAGE2("Constructing ",*this);
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<3;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE6("Adding triangle at ",_decendanttriangleRat.size()," with ",pcd->get_process(1),pcd->get_process(2),pcd->get_process(3));
}


baseboxRat::baseboxRat(const boxRat_comp* pcd) : baseCutD(pcd)
{
	_decendantboxRat.push_back(pcd);

	_width=4; // Set this up for a box
	// Add the new trees of the descendant box
#ifndef BH_PUBLIC
	BH_DEBUG_MESSAGE2("Constructing Box ",*this);
#endif
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<4;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE7("Adding box at ",_decendantboxRat.size()," with ",pcd->get_process(1),pcd->get_process(2),pcd->get_process(3),pcd->get_process(4));
}

void baseboxRat::add(const boxRat_comp* pcd)
{
	_decendantboxRat.push_back(pcd);

	// Add the new trees of the descendant box
//	_MESSAGE2("Adding ",*this);
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<4;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
    BH_DEBUG_MESSAGE7("Adding box at ",_decendantboxRat.size()," with ",pcd->get_process(1),pcd->get_process(2),pcd->get_process(3),pcd->get_process(4));
}


basepentagonRat::basepentagonRat(const pentagonRat_comp* pcd) : baseCutD(pcd)
{
	_decendantpentagonRat.push_back(pcd);

	_width=5; // Set this up for a pentagon
	// Add the new trees of the descendant pentagon
//	_MESSAGE2("Constructing ",*this);
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

    for(int i=0;i<5;i++){
//    	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//    	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//    		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//    	}
//    	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
}


void basepentagonRat::add(const pentagonRat_comp* pcd)
{
	_decendantpentagonRat.push_back(pcd);

	// Add the new trees of the descendant pentagon
#if __USE_NEW_REC_TREE_REC_RAT==1
    Rec_tree_eval_factory TF;
#else
    __RAT_TREE_TYPE TF;
#endif

//    _MESSAGE2("Adding to ",*this);
    for(int i=0;i<5;i++){
//       	Rec_Tree* nt=TF.new_tree(pcd->get_process(i+1));
//       	if(dynamic_cast<Known_Rec_Tree_base*>(nt)){
//       		_MESSAGE(*(static_cast<Known_Rec_Tree_base*>(nt)));
//       	}
//       	d_trees.push_back(nt);
		d_trees.push_back(TF.new_tree(pcd->get_process(i+1)));
    }
}

#endif

}
