/*
 * cached_EP.cpp
 *
 *  Created on: FEB 9, 2011
 *      Author: harald
 */

#include "cached_EP.h"
#include <map>
#include <vector>
#include <algorithm>
#include "settings.h"
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "mode_dependent_typedefs.h"
#include "index_vector.h"
#include "eval_param.h"
#include "BH_debug.h"
#ifndef BH_PUBLIC
#include "BH_grass.h"
#endif


#define _VERBOSE 0

using namespace std;

namespace BH {
namespace CachedTHA {


Cached_EP::Cached_EP(){}
Cached_EP::~Cached_EP(){
	for(size_t i=0;i<d_values.size();i++){
		if(d_values[i]!=0) { delete d_values[i]; d_values[i]=0; };
		if(d_values_HP[i]!=0) { delete d_values_HP[i]; d_values_HP[i]=0; };
		if(d_values_VHP[i]!=0) { delete d_values_VHP[i]; d_values_VHP[i]=0; };
	}
}

size_t Cached_EP::add(const std::vector<int>& indices){
	std::vector<std::vector<int> >::iterator it=find(d_index_vectors.begin(),d_index_vectors.end(),indices);
	if ( it == d_index_vectors.end() ){
		d_index_vectors.push_back(indices);
		size_t nbr(indices.size());
		d_load.push_back(1);
		/*
		d_values.push_back(0);
		d_values_HP.push_back(0);
		d_values_VHP.push_back(0);
		*/
		d_values.push_back(new eval_param<R>(nbr));
		d_values_HP.push_back(new eval_param<RHP>(nbr));
		d_values_VHP.push_back(new eval_param<RVHP>(nbr));
	
		d_mcIDs.push_back(0);
		d_mcIDs_HP.push_back(0);
		d_mcIDs_VHP.push_back(0);
		return d_index_vectors.size()-1;
	} else {
		size_t pos=it-d_index_vectors.begin();
		++d_load[pos];
		return pos;
	}
}

#ifndef BH_PUBLIC
size_t Cached_EP::add_gr(const std::vector<int>& indices){
	std::vector<std::vector<int> >::iterator it=find(d_index_vectors_gr.begin(),d_index_vectors_gr.end(),indices);
	if ( it == d_index_vectors_gr.end() ){
		d_index_vectors_gr.push_back(indices);
		d_load_gr.push_back(1);
		d_values_gr.push_back(0);
		d_values_gr_HP.push_back(0);
		d_values_gr_VHP.push_back(0);
		d_mcIDs_gr.push_back(0);
		d_mcIDs_gr_HP.push_back(0);
		d_mcIDs_gr_VHP.push_back(0);
		return d_index_vectors_gr.size()-1;
	} else {
		size_t pos=it-d_index_vectors_gr.begin();
		++d_load_gr[pos];
		return pos;
	}
}
#endif


/*
template <class T> eval_param<T>* Cached_EP::eval_fn(size_t n, momentum_configuration<T>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<T>(n) ){
		eval_param<T>* ep=new eval_param<T>(mc,d_index_vectors[n]);
		set_value(n, ep);
		set_mcID<T>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<T>(n);
}
*/


eval_param<R>* Cached_EP::eval(size_t n, momentum_configuration<R>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	BH_DEBUG_MESSAGE("\n=====================");
	BH_DEBUG_MESSAGE4("EP: double: ID",mc.get_ID()," ",n);
	if (mc.get_ID() != get_mcID<R>(n) ){

//		eval_param<R>* ep=new eval_param<R>(mc,d_index_vectors[n]);
//		set_value(n, ep);
		d_values[n]->update(mc,d_index_vectors[n]);

		set_mcID<R>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return d_values[n];
}


eval_param<RHP>* Cached_EP::eval(size_t n, momentum_configuration<RHP>& mc){
#if BH_USE_OMP
	get_lock();
#endif


	BH_DEBUG_MESSAGE("\n=====================");
	BH_DEBUG_MESSAGE4("EP: HP: ID",mc.get_ID()," ",n);
	if (mc.get_ID() != get_mcID<RHP>(n) ){

		d_values_HP[n]->update(mc,d_index_vectors[n]);
		set_mcID<RHP>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif

	return d_values_HP[n];
}


eval_param<RVHP>* Cached_EP::eval(size_t n, momentum_configuration<RVHP>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<RVHP>(n) ){
		//eval_param<RVHP>* ep=new eval_param<RVHP>(mc,d_index_vectors[n]);
		//set_value(n, ep);
		d_values_VHP[n]->update(mc,d_index_vectors[n]);
		set_mcID<RVHP>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return d_values_VHP[n];
}

#ifndef BH_PUBLIC
#if 1
BH_grass<R>* Cached_EP::eval_gr(size_t n, momentum_configuration<R>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID_gr<R>(n) ){
        size_t m=d_index_vectors_gr[n].size();
		BH_grass<R>* ep=new BH_grass<R>(m,mc,d_index_vectors_gr[n]);
		set_value_gr(n, ep);
		set_mcID_gr<R>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value_gr<R>(n);
}


BH_grass<RHP>* Cached_EP::eval_gr(size_t n, momentum_configuration<RHP>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID_gr<RHP>(n) ){
        size_t m=d_index_vectors_gr[n].size();
		BH_grass<RHP>* ep=new BH_grass<RHP>(m,mc,d_index_vectors_gr[n]);
		set_value_gr(n, ep);
		set_mcID_gr<RHP>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value_gr<RHP>(n);
}

BH_grass<RVHP>* Cached_EP::eval_gr(size_t n, momentum_configuration<RVHP>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID_gr<RVHP>(n) ){
        size_t m=d_index_vectors_gr[n].size();
		BH_grass<RVHP>* ep=new BH_grass<RVHP>(m,mc,d_index_vectors_gr[n]);
		set_value_gr(n, ep);
		set_mcID_gr<RVHP>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value_gr<RVHP>(n);
}
#endif


#endif

void Cached_EP::refresh(momentum_configuration<R>& mc) {
	for (int i=0;i<d_load.size();i++){
		eval(i,mc);
	}
}


/*
void Cached_EP::dry_run(){
	for (int n=0;n<d_index_vectors.size();n++){
		d_THA_p->dry_run(d_index_vectors[n]);
	}
}
*/

//Cached_EP_user_normal::~Cached_EP_user_normal(){};



template <class key,class T> struct do_delete_second : public std::unary_function<pair<key,T>&,void> {
	void operator()(pair<key,T>& p) { delete p.second;}
};


Cached_EP_factory::~Cached_EP_factory(){
	for_each(d_Cached_EP.begin(),d_Cached_EP.end(),do_delete_second<const size_t,Cached_EP*>());
}

Cached_EP* Cached_EP_factory::new_CEP(size_t n_mom_conf ){
    

    map<size_t,Cached_EP*>::iterator it = d_Cached_EP.find(n_mom_conf);
    if (it == d_Cached_EP.end()){
    	Cached_EP* new_CEP = new Cached_EP();
    	d_Cached_EP.insert(pair<size_t,Cached_EP*>(n_mom_conf,new_CEP));
        return new_CEP;
    } else {
        return (*it).second;
    }
}

void Cached_EP_factory::refresh(momentum_configuration<R>& mc){
	for ( map<size_t,Cached_EP*>::iterator it=d_Cached_EP.begin();it!=d_Cached_EP.end();it++){
		(*it).second->refresh(mc);
	};
}


Cached_EP_factory global_CEP;
Cached_EP_factory* Cached_EP_factory::default_CEP= &global_CEP;
}

}
