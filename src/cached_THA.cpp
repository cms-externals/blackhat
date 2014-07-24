/*
 * cached_THA.cpp
 *
 *  Created on: Sep 30, 2009
 *      Author: harald
 */

#include "cached_THA.h"
#include <map>
#include <vector>
#include <algorithm>
#include "settings.h"
#include "BH_typedefs.h"
#include "BH_utilities.h"
#include "mode_dependent_typedefs.h"
#ifndef BH_PUBLIC
#include "rec_tree.h"
#endif
#include "amplitudes.h"
#include "index_vector.h"
#include "eval_param.h"
#ifndef BH_PUBLIC
#include "BH_grass.h"
#include "grass_tree.h"
#endif
#include "cached_EP.h"

#include "timing.h"

#define _VERBOSE 0

using namespace std;
//using namespace CachedEP;

namespace BH {
namespace CachedTHA {


//Cached_THA::Cached_THA(const process& pro,color_structure cs){
//	d_THA_p=new OneLoopHelAmpl(pro,cs);
//}

//Cached_THA::Cached_THA(TreeHelAmpl* THAb) : d_THA_p(THAb){
//}

//Cached_THA::~Cached_THA(){
//	delete d_THA_p;
//}

void Cached_THA_ep::print_stat() {
	double sum=0;
	for (int i=0;i<d_load.size();i++){
		sum+=double(d_load[i]);
	}
	double av=sum/double(d_load.size());
	_MESSAGE6(d_THA_p->get_process(),": ",d_index_vectors.size()," index vectors used in average ", av ," times. Last values: ");
	for (int i=0;i<d_load.size();i++){
		_MESSAGE5(i,": ",d_index_vectors[i],": ",get_value<R>(i));
	}
}



size_t Cached_THA::add(const std::vector<int>& indices){
	std::vector<std::vector<int> >::iterator it=find(d_index_vectors.begin(),d_index_vectors.end(),indices);
	if ( it == d_index_vectors.end() ){
		d_index_vectors.push_back(indices);
		d_load.push_back(1);
		d_values.push_back(C(0,0));
		d_values_HP.push_back(CHP(0,0));
		d_values_VHP.push_back(CVHP(0,0));
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

/*
template <class T> std::complex<T> Cached_THA::eval(int n, momentum_configuration<T>& mc){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<T>(n) ){
// only safe (probably) if used with the thread safe(r) container used for the OMP case
#if BH_USE_OMP
		eval_param<T> ep(mc,d_index_vectors[n]);
		std::complex<T> treeVal=d_THA_p->eval(ep);
#else 
		std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);
#endif
		set_value(n, treeVal);
		set_mcID<T>(n, mc.get_ID()) ;
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<T>(n);
}
*/


template <class T> std::complex<T> Cached_THA_ep::eval_fn(int n, momentum_configuration<T>& mc,Cached_EP* CEP, size_t ep_ind){
#if BH_USE_OMP
	get_lock();
#endif

    
    //BH_START_TIMER(tree_eval_incl_lookup);
	if (mc.get_ID() != get_mcID<T>(n) ){
		eval_param<T>* ep=CEP->eval(ep_ind,mc);

        //BH_START_TIMER(ep_eval_single);
        //_MESSAGE(d_THA_p->get_process());
        std::complex<T> treeVal=d_THA_p->eval(*ep);
        //BH_STOP_TIMER(ep_eval_single);
		//std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);

		set_value(n, treeVal);
		set_mcID<T>(n, mc.get_ID());
	}
    //BH_STOP_TIMER(tree_eval_incl_lookup);
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<T>(n);
}

/*
template <class T> std::complex<T> Cached_THA_grass::eval_fn(int n, momentum_configuration<T>& mc,Cached_EP* CEP, size_t ep_ind){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<T>(n) ){
		BH_grass<T>* ep=CEP->eval_gr(ep_ind,mc);
	    std::complex<T> treeVal=d_THA_p->eval(*ep);
		//std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);

		set_value(n, treeVal);
		set_mcID<T>(n, mc.get_ID());
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<T>(n);
}
*/
#ifndef BH_PUBLIC
std::complex<R> Cached_THA_grass::eval(int n, momentum_configuration<R>& mc,Cached_EP* CEP, size_t ep_ind){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<R>(n) ){
		BH_grass<R>* ep=CEP->eval_gr(ep_ind,mc);
	    std::complex<R> treeVal=(*d_tree_eval)(*ep);
		//std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);

		set_value(n, treeVal);
		set_mcID<R>(n, mc.get_ID());
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<R>(n);
}


std::complex<RHP> Cached_THA_grass::eval(int n, momentum_configuration<RHP>& mc,Cached_EP* CEP, size_t ep_ind){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<RHP>(n) ){
		BH_grass<RHP>* ep=CEP->eval_gr(ep_ind,mc);
	    std::complex<RHP> treeVal=(*d_tree_eval_HP)(*ep);
		//std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);

		set_value(n, treeVal);
		set_mcID<RHP>(n, mc.get_ID());
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<RHP>(n);
}

std::complex<RVHP> Cached_THA_grass::eval(int n, momentum_configuration<RVHP>& mc,Cached_EP* CEP, size_t ep_ind){
#if BH_USE_OMP
	get_lock();
#endif
	if (mc.get_ID() != get_mcID<RVHP>(n) ){
		BH_grass<RVHP>* ep=CEP->eval_gr(ep_ind,mc);
	    std::complex<RVHP> treeVal=(*d_tree_eval_VHP)(*ep);
		//std::complex<T> treeVal=d_THA_p->eval(mc,d_index_vectors[n]);

		set_value(n, treeVal);
		set_mcID<RVHP>(n, mc.get_ID());
	}
#if BH_USE_OMP
	release_lock();
#endif
	return get_value<RVHP>(n);
}



#endif

/*
void Cached_THA::refresh(momentum_configuration<R>& mc) {
	for (int i=0;i<d_load.size();i++){
        eval(i,mc);
	}
}
*/


/*
void Cached_THA::dry_run(){
	for (int n=0;n<d_index_vectors.size();n++){
		d_THA_p->dry_run(d_index_vectors[n]);
	}
}
*/
const process& Cached_THA_ep::get_process() const {
	return d_THA_p->get_process();
}



Cached_THA_ep::~Cached_THA_ep(){
	delete d_THA_p;
	};
#ifndef BH_PUBLIC
Cached_THA_grass::~Cached_THA_grass(){};
#endif








Cached_THA_user_normal::~Cached_THA_user_normal(){};
Cached_THA_user_conjugate::~Cached_THA_user_conjugate(){};


template <class key,class T> struct do_delete_second : public std::unary_function<pair<key,T>&,void> {
	void operator()(std::pair<key,T>& p) { delete p.second;}
};


Cached_THA_factory::~Cached_THA_factory(){
	for_each(d_amplitudes.begin(),d_amplitudes.end(),do_delete_second<const process,Cached_THA*>());
}

void Cached_THA_factory::print_state(){
	_MESSAGE("=-=-=-=-=-=-=-=-=-=-= Cached_THA_factory =-=-=-=-=-=-=-=-=-=-= " );
	for ( map<const process,Cached_THA*>::iterator it=d_amplitudes.begin();it!=d_amplitudes.end();it++){
		(*it).second->print_stat();
	};
	_MESSAGE("=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= " );
}
/*
void Cached_THA_factory::fill_process_list(std::vector<process >& li){
	for ( map<const process,Cached_THA*>::iterator it=d_amplitudes.begin();it!=d_amplitudes.end();it++){
		li.push_back((*it).second->get_process());
	};
}
*/


void flip_helicity(process & pro){
       vector<particle_ID> v_pro;
       for(size_t i=0;i<pro.n();i++) v_pro.push_back(pro[i].helicity_conjugate());
       pro=process(v_pro);
       return;
}

Cached_THA_user* Cached_THA_factory::new_THA(process pro,const std::vector<int>& ind, short conj){

    /* conjugation infra structure
     *
     * switch: - conj==0:  no conjugation use normal user
     *         - conj==-1: conjugation with sign, use conjugate user and inset sign
     *         - conj==1: conjugation without sign, use conjugate user no sign
     */
    /*to switch off complex conj:*/ 
    //_MESSAGE("turned off conjugation");
    //conj=0;

     bool conjQ(false);
     int sign;
     switch(conj){
         case 0: {conjQ=false; sign=0;} break;
         case 1: {conjQ=true; sign=1;} break;
         case -1: {conjQ=true; sign=-1;} break;
     }
     //helicity conjugate particle IDs
     if(conjQ) flip_helicity(pro);
     /* end conjugation */

    
     /* begin eval param caching */

        size_t ep_n(0); //default for now
        Cached_EP* CEP=Cached_EP_factory::default_CEP->new_CEP(ep_n);
        size_t ep_ind(0); 
#ifndef BH_PUBLIC
        if(BH::settings::BH_interface_settings::s_use_grassman_trees){
            ep_ind=CEP->add_gr(ind);
        }
        else{
            ep_ind=CEP->add(ind);
        }
#else
        ep_ind=CEP->add(ind);
#endif

     /* end eval param caching */

    map<process,Cached_THA*>::iterator it = d_amplitudes.find(pro);
    if (it == d_amplitudes.end()){
        
       // new part
    	Cached_THA* new_CTHA(0);
        if(BH::settings::BH_interface_settings::s_use_grassman_trees){
#ifndef BH_PUBLIC
        	//TreeHelAmpl* d_THA_p;
            Tree_Fn_Ptr_eval_grass tree_eval=A_Tree_Ptr_Grass<R>(pro);
            Tree_Fn_Ptr_eval_grass_HP tree_eval_HP=A_Tree_Ptr_Grass<RHP>(pro);
            Tree_Fn_Ptr_eval_grass_VHP tree_eval_VHP=A_Tree_Ptr_Grass<RVHP>(pro);
            new_CTHA = new Cached_THA_grass(tree_eval,tree_eval_HP,tree_eval_VHP);
#endif
        }
        else{
            TreeHelAmpl* new_THAb = new TreeHelAmpl(pro);
    	    new_CTHA = new Cached_THA_ep(new_THAb);
        }
       //
        
        
        d_amplitudes.insert(pair<process,Cached_THA*>(pro,new_CTHA));
        if(conj!=0)   return new Cached_THA_user_conjugate(new_CTHA,new_CTHA->add(ind),CEP,ep_ind,sign);
        else        return new Cached_THA_user_normal(new_CTHA,new_CTHA->add(ind),CEP,ep_ind);
    } else {
    	if(conj!=0) return new Cached_THA_user_conjugate((*it).second,(*it).second->add(ind),CEP,ep_ind,sign);
        else      return new Cached_THA_user_normal((*it).second,(*it).second->add(ind),CEP,ep_ind);
    }
}

/*
void Cached_THA_factory::refresh(momentum_configuration<R>& mc){
	for ( map<const process,Cached_THA*>::iterator it=d_amplitudes.begin();it!=d_amplitudes.end();it++){
		(*it).second->refresh(mc);
	};
}
*/


Cached_THA_factory global_CTHA;
Cached_THA_factory* Cached_THA_factory::default_CTHA= &global_CTHA;


Cached_TA_factory::Cached_TA_factory(){};
Cached_TA_factory::~Cached_TA_factory(){
    
    for(size_t i=0;i<d_Cached_TA.size();i++){
        delete d_Cached_TA[i]; 
    }
    for(size_t i=0;i<d_eval_param.size();i++){
	delete d_eval_param[i];
    }
    for(size_t i=0;i<d_ind.size();i++){
	delete d_ind[i];
    }
 
};


class Cached_TA;
Cached_TA::Cached_TA(): d_complete(false) {};
Cached_TA::~Cached_TA(){};
 
size_t Cached_TA_factory::add_Cached_TA(){
    Cached_TA* CTA=new Cached_TA();
    if(d_Cached_TA.size()>0) d_Cached_TA.back()->complete_construction(); 
    d_Cached_TA.push_back(CTA);
    return d_Cached_TA.size()-1; 
};


size_t Cached_TA::add(size_t n_ep){
    std::map<size_t, size_t >::const_iterator it=d_known_eps.find(n_ep);
    if(it==d_known_eps.end()){
        eval_param<R>* ep=(default_CTA->d_eval_param)[n_ep];
        std::vector<int> * ind=(default_CTA->d_ind)[n_ep];
        d_eval_param.push_back(ep);
        d_ind.push_back(ind);
        d_known_eps.insert(make_pair(d_ind.size()-1 ,n_ep));
        return d_ind.size()-1;
    }
    return it->second;
}


void Cached_TA::add(TreeHelAmpl* THA, size_t n_ep, short conj){

    /* tree results, c.c. results and c.c. & minus sign dresses tree results are stored in p_val_c, p_val_mc
     * we use d_poll to store the final pointers to the tree-results:
     *      1) the pairs d_poll(n_tha,n_ep) are coordinates of trees and eval param pointers to read out results of the evaluation
     *      2) the pairs (d_vla_cc.size()-1,-1) are the position in d_val_c. similarly (d_vla_cc.size()-1,-2) are positions in d_val_cm
    */


    std::pair<size_t,size_t> p_val;
    map<std::pair<TreeHelAmpl*,size_t >, std::pair<size_t,size_t>  >::iterator it=d_map_vals.find(make_pair(THA,n_ep));
    if(it==d_map_vals.end()){
        size_t n_tha(0);
        std::map<TreeHelAmpl*, size_t >::const_iterator it_tha=d_known_tha.find(THA);
        if(d_known_tha.end()==it_tha){
            d_p_tree_hel_ampls.push_back(THA);
            n_tha=d_p_tree_hel_ampls.size()-1;
            d_known_tha.insert(make_pair(THA,n_tha));
        
            size_t n_loc_ep=this->add(n_ep);
            eval_param<R>* ep=(default_CTA->d_eval_param)[n_ep];
      
            std::vector<size_t > ep_vec;
            d_p_ep.push_back(ep_vec);
            d_p_ep[n_tha].push_back(n_loc_ep);
           
            p_val=make_pair(n_tha,0);
        }
        else{
            n_tha=it_tha->second;
        
            size_t n_loc_ep=this->add(n_ep);
            //eval_param<R>* ep=(default_CTA->d_eval_param)[n_ep];
            d_p_ep[n_tha].push_back(n_loc_ep);

            p_val=make_pair(n_tha,d_p_ep[n_tha].size()-1);
        }

        d_map_vals.insert(make_pair(make_pair(THA,n_ep),p_val));   
        
        // build conjugation
        switch(conj){
            case 0: {
                d_poll.push_back(p_val);
            } break;
            case 1: {
                d_poll_c.push_back(p_val);
                d_poll.push_back(make_pair(d_poll_c.size()-1,-1)); 
            } break;
            case -1: {
                d_poll_mc.push_back(p_val);
                d_poll.push_back(make_pair(d_poll_mc.size()-1,-2)); 
            } break;
        }
        return;
    }
    else{
        p_val=it->second;
        // build conjugation
        switch(conj){
            case 0: {
                d_poll.push_back(p_val);
            } break;
            case 1: {
                d_poll_c.push_back(p_val);
                d_poll.push_back(make_pair(d_poll_c.size()-1,-1)); 
            } break;
            case -1: {
                d_poll_mc.push_back(p_val);
                d_poll.push_back(make_pair(d_poll_mc.size()-1,-2)); 
            } break;
        }
        return;
    }
}

void Cached_TA::complete_construction(){

    if(!d_complete){
        p_vals=new complex<R> * [d_poll.size()]; 
        a_vals=new complex<R> * [d_p_ep.size()];
        
        if(d_poll_c.size()>0){
            d_p_vals_c=new complex<R> * [d_poll_c.size()];
            a_vals_c=new complex<R> [d_poll_c.size()];
        } 
        if(d_poll_mc.size()>0){
            d_p_vals_mc=new complex<R> * [d_poll_mc.size()];
            a_vals_mc=new complex<R> [d_poll_mc.size()];
        }


        for(size_t i=0;i<d_p_ep.size();i++){
            a_vals[i]= new complex<R>[d_p_ep[i].size()];
        }
  
//    _PRINT(d_poll.size());
//    _PRINT(d_poll_c.size());
//    _PRINT(d_poll_mc.size());

        size_t j(0),k(0);
        for(size_t i=0;i<d_poll.size();i++){
            if(d_poll[i].second>-1){
                //no sign or complex conjugation
                p_vals[i]=&a_vals[d_poll[i].first][d_poll[i].second];
            }
            else if(d_poll[i].second==-1){
                //point to cc result
                d_p_vals_c[j]=&a_vals[d_poll_c[j].first][d_poll_c[j].second];
                p_vals[i]=&a_vals_c[d_poll[i].first];
                ++j;
            }
            else{
                //_MESSAGE3(d_poll[i].first," ",d_poll[i].second);
                //point to minus cc result
                d_p_vals_mc[k]=&a_vals[d_poll_mc[k].first][d_poll_mc[k].second];
                p_vals[i]=&a_vals_mc[d_poll[i].first];
                ++k; 
            }
        }
        d_nbr_c=d_poll_c.size();
        d_nbr_mc=d_poll_mc.size();
        d_poll.clear();
        d_poll_c.clear();
        d_poll_mc.clear();
      
        d_map_vals.clear();
        d_known_tha.clear();
        d_known_eps.clear();
        
        d_complete=true;
    }

        //d_map_vals.clear();
        //d_known_eps.clear();
};

//template <class T> void Cached_TA::eval(momentum_configuration<T>& mc){
void Cached_TA::eval(momentum_configuration<R>& mc){

    it_ind=d_ind.begin();
    for(it_ep=d_eval_param.begin();it_ep!=d_eval_param.end();it_ep++){
        (*it_ep)->update(mc,**(it_ind++));
    }

    it_tha=d_p_tree_hel_ampls.begin();
    for(size_t i=0;i<d_p_ep.size();i++){
        it_n_ep=d_p_ep[i].begin();
        for(size_t j=0;j<d_p_ep[i].size();j++){
        //BH_START_TIMER(ep_eval_single);
           //a_vals[i][j]=(*it_tha)->eval(**it_ep++); 
           a_vals[i][j]=(*it_tha)->eval(*d_eval_param[(*it_n_ep++)]); 
        //BH_STOP_TIMER(ep_eval_single);
        }
        it_tha++;
    }
    for(size_t i=0;i<d_nbr_c;i++){
        a_vals_c[i]=conj(*(d_p_vals_c[i]));
    }
    for(size_t i=0;i<d_nbr_mc;i++){
        a_vals_mc[i]=-conj(*(d_p_vals_mc[i]));
    }
}

void Cached_TA::eval(momentum_configuration<RHP>& mc){

    for(it_ind=d_ind.begin();it_ind!=d_ind.end();it_ind++){
        d_eval_param_HP.push_back(new eval_param<RHP>(mc,**(it_ind)));
    }

    it_tha=d_p_tree_hel_ampls.begin();
    for(size_t i=0;i<d_p_ep.size();i++){
        it_n_ep=d_p_ep[i].begin();
        for(size_t j=0;j<d_p_ep[i].size();j++){
           a_vals[i][j]=to_double((*it_tha)->eval(*d_eval_param_HP[(*it_n_ep++)])); 
        }
        it_tha++;
    }
    for(size_t i=0;i<d_nbr_c;i++){
        a_vals_c[i]=conj(*(d_p_vals_c[i]));
    }
    for(size_t i=0;i<d_nbr_mc;i++){
        a_vals_mc[i]=-conj(*(d_p_vals_mc[i]));
    }
    
    for(size_t i=0;i<d_ind.size();i++){
        delete d_eval_param_HP[i];
    }
    d_eval_param_HP.clear();
}


void Cached_TA::eval(momentum_configuration<RVHP>& mc){

    for(it_ind=d_ind.begin();it_ind!=d_ind.end();it_ind++){
        d_eval_param_VHP.push_back(new eval_param<RVHP>(mc,**(it_ind)));
    }

    it_tha=d_p_tree_hel_ampls.begin();
    for(size_t i=0;i<d_p_ep.size();i++){
        it_n_ep=d_p_ep[i].begin();
        for(size_t j=0;j<d_p_ep[i].size();j++){
           a_vals[i][j]=to_double((*it_tha)->eval(*d_eval_param_VHP[(*it_n_ep++)])); 
        }
        it_tha++;
    }
    for(size_t i=0;i<d_nbr_c;i++){
        a_vals_c[i]=conj(*(d_p_vals_c[i]));
    }
    for(size_t i=0;i<d_nbr_mc;i++){
        a_vals_mc[i]=-conj(*(d_p_vals_mc[i]));
    }
    
    for(size_t i=0;i<d_ind.size();i++){
        delete d_eval_param_VHP[i];
    }
    d_eval_param_VHP.clear();
}



//template void Cached_TA::eval<R>(momentum_configuration<R>& mc);
//template void Cached_TA::eval<RHP>(momentum_configuration<RHP>& mc){};
//template void Cached_TA::eval<RVHP>(momentum_configuration<RVHP>& mc){};



void Cached_TA_factory::new_TreeAmpl(process pro, const std::vector<int>&  ind, short conj){

     //conjugation
     if(conj==1||conj==-1) flip_helicity(pro);
     //conjugation



    TreeHelAmpl* tha(0);
    size_t n_ep;

    std::map<const process, TreeHelAmpl* >::const_iterator it_ampl=d_amplitudes.find(pro);
    if(it_ampl==d_amplitudes.end()){
        tha=new TreeHelAmpl(pro); 
        d_amplitudes.insert(std::pair<const process,TreeHelAmpl*>(pro,tha));
        d_TreeHelAmpl.push_back(tha); 
    }
    else{
        tha=it_ampl->second;
    }

    std::map<const std::vector<int>, eval_param<R>* >::const_iterator it_ep=d_eval_params.find(ind);
    if(it_ep==d_eval_params.end()){
        eval_param<R>* ep=new eval_param<R>(ind.size());
        d_eval_param.push_back(ep);
        n_ep=d_eval_param.size()-1;
        
        d_ind.push_back(new std::vector<int>(ind));
        d_ep_n.insert(std::pair<eval_param<R>*,size_t >(ep,n_ep));

        d_eval_params.insert(std::pair<const std::vector<int>, eval_param<R>* >(ind,ep));
    }
    else{
        n_ep=d_ep_n[it_ep->second];
    }

    //complex<R>* res=d_Cached_TA.back()->add(tha,n_ep);
    d_Cached_TA.back()->add(tha,n_ep,conj);

    //return res; 
    return;
}



Cached_TA_factory global_CTA;
Cached_TA_factory* Cached_TA_factory::default_CTA= &global_CTA;
Cached_TA_factory* Cached_TA::default_CTA= &global_CTA;

}


}



