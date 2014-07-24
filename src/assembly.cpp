/*
 * assembly.cpp
 *
 *  Created on: Aug 25, 2008
 *      Author: daniel
 */

#include "assembly.h"
#include "cached_THA.h"
//#include "cached_OLHA.h"
#include <iostream>
#include "settings.h"
#include "index_vector.h"
#include "BH_debug.h"

#include "timing.h"

#if BH_USE_OMP
#include <omp.h>
#endif
#define _VERBOSE 0

using namespace std;


namespace BH {

//class Cached_OLHA_factory;
//class Cached_OLHA_user;

prop_hel_fn::prop_hel_fn(bool up_down_quark,int photonZW, int leading_vect_ax, int i, int j, const vector<ph_type> ql_ph_type): m_up_down_quark(up_down_quark) ,m_photonZW(photonZW), m_leading_vect_ax(leading_vect_ax), m_i(i), m_j(j), m_ql_ph_type(ql_ph_type) {
		R res=0.;
switch(leading_vect_ax){
	case 0:{		//leading case
	switch(m_up_down_quark){
		case 0:{m_Q=2;
			switch(m_photonZW){
				case 0: res=0.;break;  	//photon only
				case 1: {		//photon+Z
							//m_ql_ph_type [0,1] correspone to quark and electron ph_type
					if (m_ql_ph_type[0].helicity()==1 && !m_ql_ph_type[0].is_anti_particle())
						{if(m_ql_ph_type[1]==lp){res=constants::vupr*constants::ver;}
						else {res=constants::vupr*constants::vel;};}
		             		else
						{if(m_ql_ph_type[1]==lp){res=constants::vupl*constants::ver;}
						else {res=constants::vupl*constants::vel;};};
					}; break;
				case 2: {m_Q=0;		//m_ql_ph_type [0,1] correspone to quark and electron ph_type
					if (m_ql_ph_type[0].helicity()==1 && !m_ql_ph_type[0].is_anti_particle())
						{if(m_ql_ph_type[1]==lp){res=constants::vupr*constants::vnuer;}
						else {res=constants::vupr*constants::vnuel;};}
		             		else
						{if(m_ql_ph_type[1]==lp){res=constants::vupl*constants::vnuer;}
						else {res=constants::vupl*constants::vnuel;};};
					}; break;
				case 3: {		// W case
							//m_ql_ph_type [0,1] correspone to quark and electron ph_type
					if (!(m_ql_ph_type[0].helicity()==1 && !m_ql_ph_type[0].is_anti_particle()))
						{	if(!(m_ql_ph_type[1]==lp))	{res=1./constants::sin_th_2/2.;}
						}
					else {res=0.;}
					};break;
					};
		       };break;
		case 1:{m_Q=-1;
			switch(m_photonZW){
				case 0: res=0.;break;	//photon only
				case 1: {	 	//photon+Z
							//m_ql_ph_type [0,1] correspone to quark and electron ph_type
					if (m_ql_ph_type[0].helicity()==1 && !m_ql_ph_type[0].is_anti_particle())
						{if(m_ql_ph_type[1]==lp){res=constants::vdownr*constants::ver;}
						else {res=constants::vdownr*constants::vel;};}
		             		else
						{if(m_ql_ph_type[1]==lp){res=constants::vdownl*constants::ver;}
						else {res=constants::vdownl*constants::vel;};};
					} break;
				case 2: {m_Q=0;		//m_ql_ph_type [0,1] correspone to quark and electron ph_type
					if (m_ql_ph_type[0].helicity()==1 && !m_ql_ph_type[0].is_anti_particle())
						{if(m_ql_ph_type[1]==lp){res=constants::vdownr*constants::vnuer;}
						else {res=constants::vdownr*constants::vnuel;};}
		             		else
						{if(m_ql_ph_type[1]==lp){res=constants::vdownl*constants::vnuer;}
						else {res=constants::vdownl*constants::vnuel;};};
					}; break;
				case 3: {
					//_PRINT("No W-boson case"); 	//W-bosons
					res=0.;};break;
				};break;};
	};
	       };break;
	case 1:{		//vectorial case
		       		//sum over flavors included
		m_Q=2*2-3*1;
			switch(m_photonZW){
				case 0: res=0.;break;  	//photon only
				case 1: {		//photon+Z
							//m_ql_ph_type[1] corresponds to electron ph_type
					if(m_ql_ph_type[1]==lp){res=1./2.*(2.*(constants::vupr+constants::vupl)+3.*(constants::vdownr+constants::vdownl))*constants::ver;}
						else {res=1./2.*(2.*(constants::vupr+constants::vupl)+3.*(constants::vdownr+constants::vdownl))*constants::vel;};
					}; break;
				case 2: {m_Q=0;		//m_ql_ph_type[1] corresponds to electron ph_type
							// neutrinos, thus Z only
						if(m_ql_ph_type[1]==lp){res=1./2.*(2.*(constants::vupr+constants::vupl)+3.*(constants::vdownr+constants::vdownl))*constants::vnuer;}
						else {res=1./2.*(2.*(constants::vupr+constants::vupl)+3.*(constants::vdownr+constants::vdownl))*constants::vnuel;};
					}; break;
				case 3: {res=0.;  // no vectorial coupling for W-bosons
					};break;
					};
	};break;
	case 2:{		// axial case
		       		//sum over flavors included
		m_Q=0;
			switch(m_photonZW){
				case 0: res=0.;break;  	//photon only
				case 1: {		//photon+Z
							//m_ql_ph_type[1] corresponds to electron ph_type
					if(m_ql_ph_type[1]==lp){res=constants::ver/constants::sin_2th;}
						else {res=constants::vel/constants::sin_2th;};
					}; break;
				case 2: {m_Q=0;
							//m_ql_ph_type[1] corresponds to electron ph_type
							// neutrinos, thus Z only
						if(m_ql_ph_type[1]==lp){res=constants::vnuer/constants::sin_2th;}
						else {res=constants::vnuel/constants::sin_2th;};
					}; break;
				case 3: {res=0.;  // no axial coupling for W-bosons
					};break;
					};
	       };break;
};
	m_hel_coupling=res;
	}


prop_hel_fn_diphoton::prop_hel_fn_diphoton(int quark_pdb_label ){
	//only need charge squared
	switch(abs(quark_pdb_label)%2){
		case 0: m_Q2=4; break;//u-type
		case 1: m_Q2=1; break;//d-type
	}
}	

/*
subtraction& subtraction::operator=(const subtraction& sub){
	delete _tree;
	_tree=new TreeHelAmpl(sub._tree->get_process());
	_factor= sub._factor;
	_order= sub._order;
	return *this;
}
*/

template <class T> SeriesC<T> subtraction::eval(momentum_configuration<T>& mc){
	vector<complex<T> > val;
	val.push_back(_tree->eval(mc)*complex<T>(_factor,0));
for (int i=1;i<=-_order;i++){val.push_back(complex<T>(0,0));}
	SeriesC<T> res(_order,0,val);
	return res;
}



CTree_with_prefactor::CTree_with_prefactor(const process& pro, const std::vector<int>& ind):
	m_ind(ind), m_Cached_THA_user((*CachedTHA::Cached_THA_factory::default_CTHA).new_THA(pro,ind)) {
	multi_precision_constant one(1);
	m_prefactor= new constant_kinematic_function(one);};

Squared_ME::~Squared_ME(){
	for (int j=0;j<_one_loops.size();j++){
		delete _one_loops[j];
	}
	for (int j=0;j<_trees.size();j++){
		delete _trees[j];
	}

};



size_t Squared_ME::add(CTree_with_prefactor* CTwp){
	_ctrees_with_prefactor.push_back(CTwp);
	_tree_ind.push_back((*_ctrees_with_prefactor.back()).get_index_vector());
//	_tree_ind.push_back((*CTwp).get_index_vector());
	d_nbr_external=_tree_ind[0].size();
	return _ctrees_with_prefactor.size()-1;}


template <class T> SeriesC<T> Squared_ME::eval_fn(momentum_configuration<T>& mc,const vector<int>& ind, int mu_index){
	SeriesC<T> res(-2,0);
	for (int j=0;j<_cross_terms.size();j++){
		vector<int> new_ind;
		for (int k=0;k<_tree_ind[_cross_terms[j]._ind_2].size();k++){
			new_ind.push_back(ind[_tree_ind[_cross_terms[j]._ind_2][k]-1]);
		}
		BH_DEBUG_MESSAGE2("one-loop:",_one_loops[_cross_terms[j]._ind_1]->eval(mc,ind,mu_index));
		BH_DEBUG_MESSAGE2("(tree)*:",conj(_trees[_cross_terms[j]._ind_2]->eval(mc,new_ind)));
		res+= _one_loops[_cross_terms[j]._ind_1]->eval(mc,ind,mu_index)*
			conj(_trees[_cross_terms[j]._ind_2]->eval(mc,new_ind))*complex<T>(_cross_terms[j]._factor,0);
	}
	for (int j=0;j<_cross_terms_md.size();j++){
		vector<int> new_ind;
		for (int k=0;k<_tree_ind[_cross_terms_md[j]._ind_2].size();k++){
			new_ind.push_back(ind[_tree_ind[_cross_terms_md[j]._ind_2][k]-1]);
		}
		BH_DEBUG_MESSAGE2("one-loop:",_one_loops[_cross_terms[j]._ind_1]->eval(mc,ind,mu_index));
		BH_DEBUG_MESSAGE2("(tree)*:",conj(_trees[_cross_terms[j]._ind_2]->eval(mc,new_ind)));
		res+= _one_loops[_cross_terms_md[j]._ind_1]->eval(mc,ind,mu_index)*
			conj(_trees[_cross_terms_md[j]._ind_2]->eval(mc,new_ind))*
				complex<T>(_cross_terms_md[j]._factor,0)*
					std::complex<T>(T(_cross_terms_md[j].d_prefactor),0);
	}

// notice that indices are partially active here and thus should not be used.
	for (int j=0;j<_cached_cross_terms_md.size();j++){
		//ind={1,2,3...} expected

		BH_DEBUG_MESSAGE2("one-loop:",_one_loops[_cached_cross_terms_md[j]._ind_1]->eval(mc,ind,mu_index));
		BH_DEBUG_MESSAGE2("(tree)*:",conj(_ctrees_with_prefactor[_cached_cross_terms_md[j]._ind_2]->eval(mc)));
		BH_DEBUG_MESSAGE2("factors: ", std::complex<T>(T(_cached_cross_terms_md[j].d_prefactor),0));
		res+= _one_loops[_cached_cross_terms_md[j]._ind_1]->eval(mc,ind,mu_index)*
			conj(_ctrees_with_prefactor[_cached_cross_terms_md[j]._ind_2]->eval(mc))*
				complex<T>(_cached_cross_terms_md[j]._factor,0)*
					std::complex<T>(T(_cached_cross_terms_md[j].d_prefactor),0);
	}
	return res;
}

template <class T> complex<T> Squared_ME::eval_tree_fn(momentum_configuration<T>& mc,const vector<int>& ind){
	complex<T> res(0,0);
		for (int j=0;j<_cross_terms_tree.size();j++){
		vector<int> new_ind1,new_ind2;
		for (int k=0;k<_tree_ind[_cross_terms_tree[j]._ind_1].size();k++){
			new_ind1.push_back(ind[_tree_ind[_cross_terms_tree[j]._ind_1][k]-1]);
		}
		for (int k=0;k<_tree_ind[_cross_terms_tree[j]._ind_2].size();k++){
			new_ind2.push_back(ind[_tree_ind[_cross_terms_tree[j]._ind_2][k]-1]);
		}
		res+= _trees[_cross_terms_tree[j]._ind_1]->eval(mc,new_ind1)*
			conj(_trees[_cross_terms_tree[j]._ind_2]->eval(mc,new_ind2))*complex<T>(_cross_terms_tree[j]._factor,0);
	}
	for (int j=0;j<_cross_terms_tree_md.size();j++){
		vector<int> new_ind1,new_ind2;
		for (int k=0;k<_tree_ind[_cross_terms_tree_md[j]._ind_1].size();k++){
			new_ind1.push_back(ind[_tree_ind[_cross_terms_tree_md[j]._ind_1][k]-1]);
		}
		for (int k=0;k<_tree_ind[_cross_terms_tree_md[j]._ind_2].size();k++){
			new_ind2.push_back(ind[_tree_ind[_cross_terms_tree_md[j]._ind_2][k]-1]);
		}
		res+= _trees[_cross_terms_tree_md[j]._ind_1]->eval(mc,new_ind1)*
			conj(_trees[_cross_terms_tree_md[j]._ind_2]->eval(mc,new_ind2))*complex<T>(_cross_terms_tree_md[j]._factor,0)*T(_cross_terms_tree_md[j].d_prefactor);
	}

	for (int j=0;j<_cached_cross_terms_tree_md.size();j++){
		
		BH_DEBUG_MESSAGE2("tree:",_ctrees_with_prefactor[_cached_cross_terms_tree_md[j]._ind_1]->eval(mc));
		BH_DEBUG_MESSAGE2("tree:",conj(_ctrees_with_prefactor[_cached_cross_terms_tree_md[j]._ind_2]->eval(mc)));
		BH_DEBUG_MESSAGE2("factors: ", T(_cached_cross_terms_tree_md[j].d_prefactor));
		res+= _ctrees_with_prefactor[_cached_cross_terms_tree_md[j]._ind_1]->eval(mc)*
			conj(_ctrees_with_prefactor[_cached_cross_terms_tree_md[j]._ind_2]->eval(mc))*
				complex<T>(_cached_cross_terms_tree_md[j]._factor,0)*
					T(_cached_cross_terms_tree_md[j].d_prefactor);
	}
	return res;
}

void Squared_ME::dry_run(const vector<int>& ind){
	for (int j=0;j<_one_loops.size();j++){
		_one_loops[j]->dry_run(ind);
	}
}



Ampl_Info::Ampl_Info(const process& pro, const vector<int> & ind,double * re, double * im){
    real=re; 
    imag=im;
    int n=ind.size();
    int shift(0);
    for(int i=0;i<n;i++){
        m_hels.push_back(pro.p(i+1).helicity());
        shift=100*(pro.p(i+1).flavor()-1);
        m_perm.push_back(shift+ind[i]);
    }
};

/* 
void Squared_ME_Struct::print(){
 
    // for matrix element squared
    double weight_tree_right(0);
    size_t loop_left,tree_right;
    std::pair<int,int > kin_fn_pair;


    if(m_kin_fn.size()>0){
        //loop over prefactor
        for(it_ct=m_loop_cross_terms.begin();it_ct!=m_loop_cross_terms.end();it_ct++){
            //kin_fn_pair
            kin_fn_pair=it_ct->first;
            
            //loop over left-loop-amplitude
            for(it_ll=(++(it_ct->second).begin());it_ll!=(it_ct->second).end();it_ll++){
               
                loop_left=it_ll->first;
                process pro=m_loop_users[loop_left]->get_process(); 
                vector<int> ind(m_loop_users[loop_left]->get_index_vector()); 
                color_structure cs=m_loop_users[loop_left]->color_struct(); 
                _MESSAGE8("nsy: loopleft: ",pro," ",ind," ",cs, " entr ",loop_left);

                //tree_right
                for(it_rvtp=(it_ll->second).begin();it_rvtp!=(it_ll->second).end();it_rvtp++){
                    tree_right=it_rvtp->first;
                    weight_tree_right=it_rvtp->second;
                    vector<int> indp=m_tree_users[tree_right]->get_index_vector();
                    process prop(m_tree_users[tree_right]->get_process());
                    _MESSAGE8("nsy: righttree: \t\t",prop," ",indp," entr ",tree_right,":",weight_tree_right);
}}}}} 
  */ 



size_t Squared_ME_Struct::add(const process& pro, const std::vector<int>& ind,short conjQ, std::vector<kinematic_function*> kf){

    int p_tree_user,p_kin_fn;
    std::map<std::pair<process, std::vector<int> >, int >::const_iterator it0 = map_tree_users.find(make_pair(pro,ind));
    if(it0==map_tree_users.end()){
        m_tree_users.push_back((*CachedTHA::Cached_THA_factory::default_CTHA).new_THA(pro,ind,conjQ));
        m_tree_vals.push_back(complex<R>(0,0));
        p_tree_user=m_tree_users.size()-1;
        map_tree_users[make_pair(pro,ind)]=p_tree_user;
    }
    else p_tree_user=(*it0).second;

	std::vector<int> prod_kin_fn;
	if(kf.size()!=0){
		sort(kf.begin(),kf.end());
		for(size_t i=0;i<kf.size();i++){
			if(kf[i]!=0){
   				std::map<kinematic_function*,int >::const_iterator it1=map_kin_fn.find(kf[i]);
   		    	if(it1==map_kin_fn.end()){
     				m_kin_fn_vals.push_back(complex<R>(0,0));
           			m_kin_fn.push_back(kf[i]);
	            	p_kin_fn=m_kin_fn_vals.size()-1;
            		map_kin_fn[kf[i]]=p_kin_fn;
        		}
        		else{
            		p_kin_fn=(*it1).second;
        		}
    		}
    		else p_kin_fn=-1;
			prod_kin_fn.push_back(p_kin_fn);
		}
	}

    //m_trees.push_back(make_pair(p_kin_fn,p_tree_user));
    m_trees.push_back(make_pair(prod_kin_fn,p_tree_user));
    return (m_trees.size()-1);
};


/* begin symmetrize gluons*/

// returns false if gluon labels other than ( 1 or 2 )  are not ascendently ordered
bool is_ordered(const std::vector<int>& ind, const process & pro){

	int ref=0;
	for(size_t i=0; i<ind.size(); i++){
		if(pro[i].is_a(gluon) && ind[i]>2){
               if(ref > ind[i]) return false;
		       ref=ind[i];
		} 
    }
	    return true;
}

int symmetry_reweight(const std::vector<int> & ind, const process & pro){
	
	int sym_factor(1),incr(1);
	for(size_t i=0; i<ind.size(); i++){
		if(ind[i]>2 && pro[i].is_a(gluon) ) sym_factor*=(incr++);
	}
	return sym_factor;
}


/* end symmetrize gluons*/

size_t Squared_ME_Struct::add(const process& pro, const std::vector<int>& ind, color_structure cs, short conjQ, std::vector<kinematic_function*> kf){
    int p_loop_user,p_kin_fn;

    BH_DEBUG_MESSAGE4("pro: ",pro," cs: ",cs);

    //set first primitive to zero-premitive place holder    
    if(m_loop_users.size()==0){
        m_symmetry_factor=1;
	    if(!settings::BH_interface_settings::s_final_state_gluon_symmetrization){
            m_symmetry_factor=symmetry_reweight(ind,pro);
        }
        
        m_loop_users.push_back(0);
    	m_loop_finite.push_back(complex<R>(0,0));
        m_loop_single.push_back(complex<R>(0,0));
        m_loop_double.push_back(complex<R>(0,0));
       	p_loop_user=0;
       	map_loop_users[make_pair(CachedOLHA::pro_cs(),vector<int>())]=p_loop_user;
    }	

    /* begin final state symmetrization */
	
	    //return trivial primitive for non-ascending gluon momentum labels	
    	//will set all coefficients of this amplitude in partials to zero


    if(settings::BH_interface_settings::s_final_state_gluon_symmetrization
            || (!settings::BH_interface_settings::s_final_state_gluon_symmetrization && is_ordered(ind, pro) )){
        
        std::map<std::pair< CachedOLHA::pro_cs,std::vector<int> >, int >::const_iterator it0 = map_loop_users.find(make_pair(CachedOLHA::pro_cs(pro,cs),ind));
        if(it0==map_loop_users.end()){
            m_loop_users.push_back(CachedOLHA::Cached_OLHA_factory::default_COLHA->new_OLHA(pro,cs,ind,conjQ));
           
            m_loop_finite.push_back(complex<R>(0,0));
            m_loop_single.push_back(complex<R>(0,0));
            m_loop_double.push_back(complex<R>(0,0));
            p_loop_user=m_loop_users.size()-1;
            map_loop_users[make_pair(CachedOLHA::pro_cs(pro,cs),ind)]=p_loop_user;
        }
        else p_loop_user=(*it0).second;
    }
    else p_loop_user=0;
    
    /* end  final state symmetrization */

	vector<int> prod_kin_fn;
	if(kf.size()!=0){
		sort(kf.begin(),kf.end());
		for(size_t i=0;i<kf.size();i++){
		    if(kf[i]!=0){
	        	std::map<kinematic_function*,int >::const_iterator it1=map_kin_fn.find(kf[i]);
	        	if(it1==map_kin_fn.end()){
		            m_kin_fn_vals.push_back(complex<R>(0,0));
		            m_kin_fn.push_back(kf[i]);
		            p_kin_fn=m_kin_fn_vals.size()-1;
		            map_kin_fn[kf[i]]=p_kin_fn;
		        }
		        else{
		            p_kin_fn=(*it1).second;
		        }
	    	}
		    else p_kin_fn=-1;
			prod_kin_fn.push_back(p_kin_fn);
		}
	}
    m_loops.push_back(make_pair(prod_kin_fn,p_loop_user));
    return (m_loops.size()-1);
};



// (prefactor_left*prefactor_right)*(
//                                  left_tree1 ( right_tree1*prefactor+ ... )^*
//                                  left_tree2 ( right_tree1*prefactor+ ... )^*
//                                  left_tree3 ( right_tree1*prefactor+ ... )^*
//                                  )
//
void Squared_ME_Struct::add_tree(int tree_left, int tree_right, double pre_factor){

    if(pre_factor==0.) return;
    std::pair<vector<int>,vector<int> > kin_fn_pair(m_trees[tree_left].first,m_trees[tree_right].first);
    std::map<std::pair<std::vector<int>,std::vector<int> > , std::map<int,std::map<int,double> > >::iterator it=m_cross_terms.find(kin_fn_pair);
    if(it==m_cross_terms.end()){
        std::map<int,double> right_tree_and_prefactor;
        right_tree_and_prefactor[tree_right]=pre_factor;
        std::map<int,std::map<int,double> > m_left_tree_terms;
        m_left_tree_terms[tree_left]=right_tree_and_prefactor;
        m_cross_terms[kin_fn_pair]=m_left_tree_terms; 
    }
    else{
        std::map<int,std::map<int,double> > & m_left_tree_terms=(*it).second;
        std::map<int,std::map<int,double> >::iterator it_1=m_left_tree_terms.find(tree_left);
        if(it_1==m_left_tree_terms.end()){
            std::map<int,double> m_right_tree_and_prefactor;
            m_right_tree_and_prefactor[tree_right]=pre_factor;
            m_left_tree_terms[tree_left]=m_right_tree_and_prefactor;
        }
        else{
            std::map<int,double> & m_right_tree_and_prefactor=(*it_1).second;
            std::map<int,double>::iterator it_2=m_right_tree_and_prefactor.find(tree_right);
            if(it_2==m_right_tree_and_prefactor.end()){
                m_right_tree_and_prefactor[tree_right]=pre_factor;
            }
            else{
                m_right_tree_and_prefactor[tree_right]+=pre_factor;
            }
        }
    }
    return;
}

size_t Squared_ME_Struct::add_partial(std::vector< std::pair<int,double> >& partial_ampl){
    size_t p_partial_ampl;
    std::map<std::vector< std::pair<int,double> >, int >::iterator it=map_partial_ampls.find(partial_ampl);
    if(it==map_partial_ampls.end()){
        m_partial_ampls.push_back(partial_ampl);
        p_partial_ampl=m_partial_ampls.size()-1;
        map_partial_ampls[partial_ampl]=p_partial_ampl;
    }
    else p_partial_ampl=it->second;

    return p_partial_ampl;
}


void Squared_ME_Struct::add_loop(int partial_left, int tree_right, double pre_factor){

    std::vector< std::pair<int,double> > loop_left_vec(m_partial_ampls[partial_left]);
    double weight_loop_left(0);
    size_t loop_left;
    pre_factor*=double(m_symmetry_factor);
    
    for(size_t i=0;i<loop_left_vec.size();i++){//loop over primitives in partial
        loop_left=loop_left_vec[i].first;
        weight_loop_left=loop_left_vec[i].second;//weight of primitive in partial
	    if(weight_loop_left==0.) continue; //drop zero weight primitives
        if(m_loops[loop_left].second==0) continue; //drop zero primitive

    //-----for each primitive we sort as for tree amplitudes
    std::pair<std::vector<int>,std::vector<int> > kin_fn_pair(m_loops[loop_left].first,m_trees[tree_right].first);
    std::map<std::pair<std::vector<int>,std::vector<int> > , std::map<int,std::map<int,double> > >::iterator it=m_loop_cross_terms.find(kin_fn_pair);
    if(it==m_loop_cross_terms.end()){
        std::map<int,double> right_tree_and_prefactor;
        right_tree_and_prefactor[tree_right]=pre_factor*weight_loop_left;
        std::map<int,std::map<int,double> > m_left_loop_terms;
        m_left_loop_terms[m_loops[loop_left].second]=right_tree_and_prefactor;
        m_loop_cross_terms[kin_fn_pair]=m_left_loop_terms; 
    }
    else{
        std::map<int,std::map<int,double> > & m_left_loop_terms=(*it).second;
        std::map<int,std::map<int,double> >::iterator it_1=m_left_loop_terms.find(m_loops[loop_left].second);
        if(it_1==m_left_loop_terms.end()){
            std::map<int,double> m_right_tree_and_prefactor;
            m_right_tree_and_prefactor[tree_right]=pre_factor*weight_loop_left;
            m_left_loop_terms[m_loops[loop_left].second]=m_right_tree_and_prefactor;
        }
        else{
            std::map<int,double> & m_right_tree_and_prefactor=(*it_1).second;
            std::map<int,double>::iterator it_2=m_right_tree_and_prefactor.find(tree_right);
            if(it_2==m_right_tree_and_prefactor.end()){
                m_right_tree_and_prefactor[tree_right]=pre_factor*weight_loop_left;
            }
            else{
                m_right_tree_and_prefactor[tree_right]+=pre_factor*weight_loop_left;
            }
        }
    }
    //--- move on to next primitive in partial
    }
    return;
}


void Squared_ME_Struct::complete_construction(){
 
    if(m_cross_terms.size()==0) return;

    int left_label, right_label;
    if(m_kin_fn.size()>0){
        for(it_ct=m_cross_terms.begin();it_ct!=m_cross_terms.end();it_ct++){

			//v_kin_fn_vals_l.push_back(& m_kin_fn_vals[(it_ct->first).first]);
            //v_kin_fn_vals_r.push_back(& m_kin_fn_vals[(it_ct->first).second]); 
            v_kin_fn_vals_l.push_back(std::vector<complex<R>* >());
            v_kin_fn_vals_r.push_back(std::vector<complex<R>* >());
			for(size_t i=0; i<(it_ct->first).first.size();i++){
				BH_DEBUG_MESSAGE2("tree couplings_l: ",i);
				v_kin_fn_vals_l.back().push_back(& m_kin_fn_vals[((it_ct->first).first)[i]]);
			}
			for(size_t i=0; i<(it_ct->first).second.size();i++){
				BH_DEBUG_MESSAGE2("tree couplings_r: ",i);
				v_kin_fn_vals_r.back().push_back(& m_kin_fn_vals[((it_ct->first).second)[i]]);
			}


            v_tree_vals_l.push_back(std::vector<complex<R>* >());
            v_tree_vals_r.push_back(std::vector<std::vector<complex<R>* > >());
            v_color_factor.push_back(std::vector<std::vector<double > >());

            for(it_lt=(it_ct->second).begin();it_lt!=(it_ct->second).end();it_lt++){
                v_tree_vals_l.back().push_back(& m_tree_vals[it_lt->first]);
                
                std::vector<complex<R>* > loc_tree_vals_r;
                std::vector<double > loc_color_factor;
                for(it_rtp=(it_lt->second).begin();it_rtp!=(it_lt->second).end();it_rtp++){
                    loc_tree_vals_r.push_back(& m_tree_vals[it_rtp->first]);
                    loc_color_factor.push_back(it_rtp->second);
                }
                v_tree_vals_r.back().push_back(loc_tree_vals_r);
                v_color_factor.back().push_back(loc_color_factor);
            }
        }
    }
    else{    
        it_ct=m_cross_terms.begin();
            v_tree_vals_l.push_back(std::vector<complex<R>* >());


            v_tree_vals_r.push_back(std::vector<std::vector<complex<R>* > >());
            v_color_factor.push_back(std::vector<std::vector<double > >());

            

            for(it_lt=(it_ct->second).begin();it_lt!=(it_ct->second).end();it_lt++){
                v_tree_vals_l.back().push_back(& m_tree_vals[it_lt->first]);
                
                std::vector<complex<R>* > loc_tree_vals_r;
                std::vector<double > loc_color_factor;
                for(it_rtp=(it_lt->second).begin();it_rtp!=(it_lt->second).end();it_rtp++){
                    loc_tree_vals_r.push_back(& m_tree_vals[it_rtp->first]);
                    loc_color_factor.push_back(it_rtp->second);


                }
                v_tree_vals_r.back().push_back(loc_tree_vals_r);
                v_color_factor.back().push_back(loc_color_factor);
            }

    }

    //vectorize virtual if constructed
    this->complete_virt_construction();
    m_cross_terms.clear();
    return;
}


//transcribe into vectors
void Squared_ME_Struct::complete_virt_construction(){

	BH_DEBUG_MESSAGE("complete_virt_construction");
    if(m_loop_cross_terms.size()==0) return;

    if(m_kin_fn.size()>0){
        //loop over prefactor
        for(it_ct=m_loop_cross_terms.begin();it_ct!=m_loop_cross_terms.end();it_ct++){
          
			//vv_kin_fn_vals_l.push_back(& m_kin_fn_vals[(it_ct->first).first]);
            //vv_kin_fn_vals_r.push_back(& m_kin_fn_vals[(it_ct->first).second]); 
            vv_kin_fn_vals_l.push_back(std::vector<complex<R>* >());
            vv_kin_fn_vals_r.push_back(std::vector<complex<R>* >());
			for(size_t i=0; i<(it_ct->first).first.size();i++){
				BH_DEBUG_MESSAGE2("loop couplings_l: ",i);
				vv_kin_fn_vals_l.back().push_back(& m_kin_fn_vals[((it_ct->first).first)[i]]);
			}
			for(size_t i=0; i<(it_ct->first).second.size();i++){
				BH_DEBUG_MESSAGE2("loop couplings_r: ",i);
				vv_kin_fn_vals_r.back().push_back(& m_kin_fn_vals[((it_ct->first).second)[i]]);
			}


            vv_loop_finite_vals_l.push_back(std::vector<complex<R>* >());
            vv_loop_single_vals_l.push_back(std::vector<complex<R>* >());
            vv_loop_double_vals_l.push_back(std::vector<complex<R>* >());
            vv_tree_vals_r.push_back(std::vector<std::vector<complex<R>* > >());
            vv_color_factor.push_back(std::vector<std::vector<double > >());
            
            
            //loop over left-loop-amplitude
            for(it_ll=(it_ct->second).begin();it_ll!=(it_ct->second).end();it_ll++){
             
                vv_loop_finite_vals_l.back().push_back(& m_loop_finite[it_ll->first]);
                vv_loop_single_vals_l.back().push_back(& m_loop_single[it_ll->first]);
                vv_loop_double_vals_l.back().push_back(& m_loop_double[it_ll->first]);
                
                std::vector<complex<R>* > loc_tree_vals_r;
                std::vector<double > loc_color_factor;
                
                for(it_rvtp=(it_ll->second).begin();it_rvtp!=(it_ll->second).end();it_rvtp++){
                    loc_tree_vals_r.push_back(& m_tree_vals[it_rvtp->first]);
                    loc_color_factor.push_back(it_rvtp->second);
		    BH_DEBUG_MESSAGE2("write vv_color_factor: ",it_rvtp->second);
                }
                
                vv_tree_vals_r.back().push_back(loc_tree_vals_r);
                vv_color_factor.back().push_back(loc_color_factor);
        }
        }
    }
    else{    
        //without prefactor (e.g. pure QCD)
        it_ct=m_loop_cross_terms.begin();
        
            vv_loop_finite_vals_l.push_back(std::vector<complex<R>* >());
            vv_loop_single_vals_l.push_back(std::vector<complex<R>* >());
            vv_loop_double_vals_l.push_back(std::vector<complex<R>* >());


            vv_tree_vals_r.push_back(std::vector<std::vector<complex<R>* > >());
            vv_color_factor.push_back(std::vector<std::vector<double > >());

            //size_t n_rows=m_tree_vals.size();
            //double * color_matrix=new double[n_rows*n_rows];
        
            //loop over left-tree
            for(it_ll=(it_ct->second).begin();it_ll!=(it_ct->second).end();it_ll++){
            
                vv_loop_finite_vals_l.back().push_back(& m_loop_finite[it_ll->first]);
                vv_loop_single_vals_l.back().push_back(& m_loop_single[it_ll->first]);
                vv_loop_double_vals_l.back().push_back(& m_loop_double[it_ll->first]);
                
                std::vector<complex<R>* > loc_tree_vals_r;
                std::vector<double > loc_color_factor;
                
                for(it_rvtp=(it_ll->second).begin();it_rvtp!=(it_ll->second).end();it_rvtp++){
                    loc_tree_vals_r.push_back(& m_tree_vals[it_rvtp->first]);
                    loc_color_factor.push_back(it_rvtp->second);
                }
 
                vv_tree_vals_r.back().push_back(loc_tree_vals_r);
                vv_color_factor.back().push_back(loc_color_factor);
            }
    }

    m_loop_cross_terms.clear();
    return;
}


template <class T> complex<T> Squared_ME_Struct::eval_tree_fn(momentum_configuration<T>& mc, const std::vector<int>& ind){

    complex<R> res(0,0),prefactor,tree_product,tree_left,tree_sum_right;

   for(int i=0;i<m_tree_users.size();i++){
        m_tree_vals[i]=to_double(m_tree_users[i]->eval(mc));
//        BH_DEBUG_MESSAGE5(i/*," ",m_tree_users[i].get_index_vector()*/,":",m_tree_users[i]->get_process(),":",m_tree_vals[i]);
    }
    if(m_kin_fn.size()>0){
       for(int i=0;i<m_kin_fn.size();i++){
           m_kin_fn_vals[i]=to_double(m_kin_fn[i]->eval(mc));
       }
    }

    int left_label, right_label;
    BH_DEBUG_PRINT(v_kin_fn_vals_l.size());
    if(v_kin_fn_vals_l.size()>0){

        //it_kin_fn_vals_l=v_kin_fn_vals_l.begin();
        //it_kin_fn_vals_r=v_kin_fn_vals_r.begin();

        for(size_t i_kin=0;i_kin<v_kin_fn_vals_l.size();i_kin++){
           
			//prefactor=
            //    (**(it_kin_fn_vals_l++))
            //    *conj(**(it_kin_fn_vals_r++));
	    prefactor=complex<R>(1,0);
            it_kin_fn_vals_l=v_kin_fn_vals_l[i_kin].begin();
            it_kin_fn_vals_r=v_kin_fn_vals_r[i_kin].begin();
	    for(size_t h=0;h<v_kin_fn_vals_l[i_kin].size();h++){	
		prefactor*=(**(it_kin_fn_vals_l++));
	    }
	    for(size_t h=0;h<v_kin_fn_vals_r[i_kin].size();h++){	
                prefactor*=conj(**(it_kin_fn_vals_r++));
	    }

            tree_product=complex<R>(0,0);
            //loop over left-tree
            

            it_tree_vals_l=v_tree_vals_l[i_kin].begin();
	    for(size_t i_tree_l=0;i_tree_l<v_tree_vals_l[i_kin].size();i_tree_l++){
                tree_left=**(it_tree_vals_l++);
                tree_sum_right=complex<R>(0,0);
            	
		it_tree_vals_r=v_tree_vals_r[i_kin][i_tree_l].begin();
		it_color_factor=v_color_factor[i_kin][i_tree_l].begin();
                for(size_t i_tree_r=0;i_tree_r< v_color_factor[i_kin][i_tree_l].size();i_tree_r++){
                    //tree_sum_right+=m_tree_vals[it_rtp->first]*(it_rtp->second);  
                    tree_sum_right+=(**(it_tree_vals_r++))*(*(it_color_factor++));  
                }
                tree_product+=tree_left*conj(tree_sum_right);
            }
        res+=prefactor*tree_product;
        }
    }
    else{

	    size_t i_kin(0); 
        tree_product=complex<R>(0,0);
        //loop over left-tree
            
        it_tree_vals_l=v_tree_vals_l[i_kin].begin();
	    for(size_t i_tree_l=0;i_tree_l<v_tree_vals_l[i_kin].size();i_tree_l++){
            tree_left=**(it_tree_vals_l++);
            tree_sum_right=complex<R>(0,0);
            	
		    it_tree_vals_r=v_tree_vals_r[i_kin][i_tree_l].begin();
    		it_color_factor=v_color_factor[i_kin][i_tree_l].begin();
            for(size_t i_tree_r=0;i_tree_r< v_color_factor[i_kin][i_tree_l].size();i_tree_r++){
                tree_sum_right+=(**(it_tree_vals_r++))*(*(it_color_factor++));  
            }
            tree_product+=tree_left*conj(tree_sum_right);
        }
        res+=tree_product;
    }
    return complex<T>(res);
}
#if 0

//should be functional up to update for multi V-bosons not done here.

template <class T> complex<T> Squared_ME_Struct::eval_tree_vect_fn(momentum_configuration<T>& mc, const std::vector<int>& ind){
   
    BH_START_TIMER(tree_eval_fn);
    complex<R> res(0,0),prefactor,tree_product,tree_left,tree_sum_right;
    
    BH_START_TIMER(eval_tree_fn);
    for(int i=0;i<m_tree_users.size();i++){
        m_tree_vals[i]=to_double(m_tree_users[i]->eval(mc));
//        BH_DEBUG_MESSAGE5(i/*," ",m_tree_users[i].get_index_vector()*/,":",m_tree_users[i]->get_process(),":",m_tree_vals[i]);
    }
    if(m_kin_fn.size()>0){
       for(int i=0;i<m_kin_fn.size();i++){
           m_kin_fn_vals[i]=to_double(m_kin_fn[i]->eval(mc));
       }
    }
    BH_STOP_TIMER(eval_tree_fn);

    // for partial born we only need the trees
    BH_START_TIMER(eval_tree_vals_fn);
    if(m_tree_vals_re.size()>0){
        if(m_kin_fn.size()==0){
            for(int i=0;i<m_tree_users.size();i++){
                m_tree_vals_re[i]=m_couplings*m_tree_vals[i].real();
                m_tree_vals_im[i]=m_couplings*m_tree_vals[i].imag();
            }
            //compute full matrix elem squared as needed fro check from sherpa
            //return complex<R>(0,0);
        }
        else{
             complex<R> prefactor(1,0);
             for(int i=0;i<m_tree_users.size();i++){
                prefactor=m_kin_fn_vals[m_kin_ps_tree_val[i]]*m_couplings;
                
                m_tree_vals_re[i]=(prefactor*m_tree_vals[i]).real();
                m_tree_vals_im[i]=(prefactor*m_tree_vals[i]).imag();
            }
            //compute full matrix elem squared as needed fro check from sherpa
            //return complex<R>(0,0);
        }

    }
    BH_STOP_TIMER(eval_tree_vals_fn);
    
    // for matrix element squared
    int left_label, right_label;
    if(m_kin_fn.size()>0){
        //loop over prefactor
        for(it_ct=m_cross_terms.begin();it_ct!=m_cross_terms.end();it_ct++){
            prefactor=m_kin_fn_vals[(it_ct->first).first]*
                conj(m_kin_fn_vals[(it_ct->first).second]);

            tree_product=complex<R>(0,0);
            //loop over left-tree
            for(it_lt=(it_ct->second).begin();it_lt!=(it_ct->second).end();it_lt++){
                tree_left=m_tree_vals[it_lt->first];
                tree_sum_right=complex<R>(0,0);
                for(it_rtp=(it_lt->second).begin();it_rtp!=(it_lt->second).end();it_rtp++){
                    tree_sum_right+=m_tree_vals[it_rtp->first]*(it_rtp->second);  
                }
                tree_product+=tree_left*conj(tree_sum_right);
            }
        res+=prefactor*tree_product;
        }
    }
    else{    
        //no kin functions 
        it_ct=m_cross_terms.begin();
            tree_product=complex<R>(0,0);
            //loop over left-tree
            for(it_lt=(it_ct->second).begin();it_lt!=(it_ct->second).end();it_lt++){
                tree_left=m_tree_vals.at(it_lt->first);
                tree_sum_right=complex<R>(0,0);
    BH_START_TIMER(tree_sum_right);
                for(it_rtp=(it_lt->second).begin();it_rtp!=(it_lt->second).end();it_rtp++){
                   tree_sum_right+=m_tree_vals.at(it_rtp->first)*(it_rtp->second);  
                }
    BH_STOP_TIMER(tree_sum_right);
                tree_product+=tree_left*conj(tree_sum_right);
            }
        res+=tree_product;
    }

    BH_STOP_TIMER(tree_eval_fn);
    
    //BH_START_TIMER(tree_eval_fn_new);
    //_PRINT(this->eval_tree_vect_fn(mc,ind));
    //BH_STOP_TIMER(tree_eval_fn_new);

    return complex<T>(res);
}
#endif

/* old
template <class T> SeriesC<T> Squared_ME_Struct::eval_fn(momentum_configuration<T>& mc, const std::vector<int>& ind,int mu_index){
  
    //BH_DEBUG(this->print());

    complex<R> prefactor,tree_sum_right;
    complex<R> tree_product[3], loop_left[3], res[3];
    SeriesC<T> res_loc;
  
    for(int i=0;i<m_tree_users.size();i++){
        m_tree_vals[i]=to_double(m_tree_users[i]->eval(mc));
    }
    for(int i=1;i<m_loop_users.size();i++){
        res_loc=m_loop_users[i]->eval(mc,mu_index);
        m_loop_finite[i]=to_double(res_loc[0]);
        m_loop_single[i]=to_double(res_loc[-1]);
        m_loop_double[i]=to_double(res_loc[-2]);
    }

    if(m_kin_fn.size()>0){
       for(int i=0;i<m_kin_fn.size();i++){
           m_kin_fn_vals[i]=to_double(m_kin_fn[i]->eval(mc));
       }
    }
   
    // for matrix element squared
    int left_label, right_label;
    if(m_kin_fn.size()>0){
        //loop over prefactor
        for(it_ct=m_loop_cross_terms.begin();it_ct!=m_loop_cross_terms.end();it_ct++){
            prefactor=m_kin_fn_vals[(it_ct->first).first]*
                conj(m_kin_fn_vals[(it_ct->first).second]);
            tree_product[0]=complex<R>(0,0);
            tree_product[1]=complex<R>(0,0);
            tree_product[2]=complex<R>(0,0);
            
            //loop over left-loop-amplitude
            for(it_ll=(it_ct->second).begin();it_ll!=(it_ct->second).end();it_ll++){
                loop_left[0]=m_loop_finite[it_ll->first];
                loop_left[1]=m_loop_single[it_ll->first];
                loop_left[2]=m_loop_double[it_ll->first];
               tree_sum_right=complex<R>(0,0);
                for(it_rvtp=(it_ll->second).begin();it_rvtp!=(it_ll->second).end();it_rvtp++){
                    tree_sum_right+=m_tree_vals[it_rvtp->first]*(it_rvtp->second);  
                }
                tree_product[0]+=loop_left[0]*conj(tree_sum_right);
                tree_product[1]+=loop_left[1]*conj(tree_sum_right);
                tree_product[2]+=loop_left[2]*conj(tree_sum_right);
        }
        res[0]+=prefactor*tree_product[0];
        res[1]+=prefactor*tree_product[1];
        res[2]+=prefactor*tree_product[2];
        }
    }
    else{    
        //without prefactor (e.g. pure QCD)
        it_ct=m_loop_cross_terms.begin();
            tree_product[0]=complex<R>(0,0);
            tree_product[1]=complex<R>(0,0);
            tree_product[2]=complex<R>(0,0);
            //loop over left-tree
            for(it_ll=(it_ct->second).begin();it_ll!=(it_ct->second).end();it_ll++){
                loop_left[0]=m_loop_finite[it_ll->first];
                loop_left[1]=m_loop_single[it_ll->first];
                loop_left[2]=m_loop_double[it_ll->first];
                tree_sum_right=complex<R>(0,0);
                for(it_rvtp=(it_ll->second).begin();it_rvtp!=(it_ll->second).end();it_rvtp++){
                   tree_sum_right+=m_tree_vals[it_rvtp->first]*(it_rvtp->second);  
                }
                tree_product[0]+=loop_left[0]*conj(tree_sum_right);
                tree_product[1]+=loop_left[1]*conj(tree_sum_right);
                tree_product[2]+=loop_left[2]*conj(tree_sum_right);
            }
        res[0]+=tree_product[0];
        res[1]+=tree_product[1];
        res[2]+=tree_product[2];
    }

    vector<complex<T> >res_v;
    res_v.push_back(complex<T>(res[2]));
    res_v.push_back(complex<T>(res[1]));
    res_v.push_back(complex<T>(res[0]));

    return SeriesC<T>(-2,0,res_v);
    }
*/


template <class T> SeriesC<T> Squared_ME_Struct::eval_fn(momentum_configuration<T>& mc, const std::vector<int>& ind,int mu_index){
  
    //BH_DEBUG(this->print());

    complex<R> prefactor,tree_sum_right;
    complex<R> tree_product[3], loop_left[3], res[3];
    SeriesC<T> res_loc;
  
    for(int i=0;i<m_tree_users.size();i++){
        m_tree_vals[i]=to_double(m_tree_users[i]->eval(mc));
	BH_DEBUG_MESSAGE2("tree_vals: ",m_tree_vals[i]);
    }
    for(int i=1;i<m_loop_users.size();i++){
        res_loc=m_loop_users[i]->eval(mc,mu_index);
        m_loop_finite[i]=to_double(res_loc[0]);
        m_loop_single[i]=to_double(res_loc[-1]);
        m_loop_double[i]=to_double(res_loc[-2]);
	BH_DEBUG_MESSAGE4("loop_finite: ",i,":",m_loop_finite[i]);
    }
/*
    for(int i=1;i<m_loop_users.size();i++){
	_MESSAGE7(i,":",m_loop_finite[i]," 1/e:",m_loop_single[i]," 1/e^2",m_loop_double[i]);
	}
*/

    if(m_kin_fn.size()>0){
       for(int i=0;i<m_kin_fn.size();i++){
           m_kin_fn_vals[i]=to_double(m_kin_fn[i]->eval(mc));
	   BH_DEBUG_MESSAGE2("kin_fn: ",m_kin_fn_vals[i]);
       }
    }
   
    // for matrix element squared
    int left_label, right_label;
    if(vv_kin_fn_vals_l.size()>0){

        //itv_kin_fn_vals_l=vv_kin_fn_vals_l.begin();
        //itv_kin_fn_vals_r=vv_kin_fn_vals_r.begin();

        for(size_t i_kin=0;i_kin<vv_kin_fn_vals_l.size();i_kin++){
            //prefactor=
            //    (**(itv_kin_fn_vals_l++))
            //    *conj(**(itv_kin_fn_vals_r++));

			prefactor=complex<R>(1,0);
	        	//itv_kin_fn_vals_l=vv_kin_fn_vals_l[i_kin].begin();
        		//itv_kin_fn_vals_r=vv_kin_fn_vals_r[i_kin].begin();
			for(size_t h=0;h<vv_kin_fn_vals_l[i_kin].size();h++){	
				//prefactor*=(**(itv_kin_fn_vals_l++));
				prefactor*= *vv_kin_fn_vals_l[i_kin][h];	
			}
			for(size_t h=0;h<vv_kin_fn_vals_r[i_kin].size();h++){	
                		//prefactor*=conj(**(itv_kin_fn_vals_r++));
				prefactor*= conj(*vv_kin_fn_vals_r[i_kin][h]);	
			}


            tree_product[0]=complex<R>(0,0);
            tree_product[1]=complex<R>(0,0);
            tree_product[2]=complex<R>(0,0);
            //loop over left-tree

            //itv_loop_finite_vals_l=vv_loop_finite_vals_l[i_kin].begin();
            //itv_loop_single_vals_l=vv_loop_single_vals_l[i_kin].begin();
            //itv_loop_double_vals_l=vv_loop_double_vals_l[i_kin].begin();

    	    for(size_t i_loop_l=0;i_loop_l<vv_loop_finite_vals_l[i_kin].size();i_loop_l++){
		loop_left[0]=*vv_loop_finite_vals_l[i_kin][i_loop_l];
		loop_left[1]=*vv_loop_single_vals_l[i_kin][i_loop_l];
		loop_left[2]=*vv_loop_double_vals_l[i_kin][i_loop_l];
                //loop_left[0]=**(itv_loop_finite_vals_l++);
                //loop_left[1]=**(itv_loop_single_vals_l++);
                //loop_left[2]=**(itv_loop_double_vals_l++);
                tree_sum_right=complex<R>(0,0);
            	
		        //itv_tree_vals_r=vv_tree_vals_r[i_kin][i_loop_l].begin();
        		//itv_color_factor=vv_color_factor[i_kin][i_loop_l].begin();
                for(size_t i_tree_r=0;i_tree_r< vv_color_factor[i_kin][i_loop_l].size();i_tree_r++){
                    //tree_sum_right+=(**(itv_tree_vals_r++))*(*(itv_color_factor++));  
		    //BH_DEBUG_MESSAGE2("tree_color_factor: ",(*(itv_color_factor++)));
                    tree_sum_right+=(*vv_tree_vals_r[i_kin][i_loop_l][i_tree_r])*(vv_color_factor[i_kin][i_loop_l][i_tree_r]);  
			BH_DEBUG_MESSAGE2("tree_color_factor: ",vv_color_factor[i_kin][i_loop_l][i_tree_r]);

                }
                tree_product[0]+=loop_left[0]*conj(tree_sum_right);
                tree_product[1]+=loop_left[1]*conj(tree_sum_right);
                tree_product[2]+=loop_left[2]*conj(tree_sum_right);
            }
            res[0]+=prefactor*tree_product[0];
            res[1]+=prefactor*tree_product[1];
            res[2]+=prefactor*tree_product[2];
        }
    }
    else{    
 
	    size_t i_kin(0); 

            tree_product[0]=complex<R>(0,0);
            tree_product[1]=complex<R>(0,0);
            tree_product[2]=complex<R>(0,0);
            //loop over left-tree

            itv_loop_finite_vals_l=vv_loop_finite_vals_l[i_kin].begin();
            itv_loop_single_vals_l=vv_loop_single_vals_l[i_kin].begin();
            itv_loop_double_vals_l=vv_loop_double_vals_l[i_kin].begin();

    	    for(size_t i_loop_l=0;i_loop_l<vv_loop_finite_vals_l[i_kin].size();i_loop_l++){
                loop_left[0]=**(itv_loop_finite_vals_l++);
                loop_left[1]=**(itv_loop_single_vals_l++);
                loop_left[2]=**(itv_loop_double_vals_l++);
                tree_sum_right=complex<R>(0,0);
            	
		        itv_tree_vals_r=vv_tree_vals_r[i_kin][i_loop_l].begin();
        		itv_color_factor=vv_color_factor[i_kin][i_loop_l].begin();
                for(size_t i_tree_r=0;i_tree_r< vv_color_factor[i_kin][i_loop_l].size();i_tree_r++){
                    tree_sum_right+=(**(itv_tree_vals_r++))*(*(itv_color_factor++));  
                }
                tree_product[0]+=loop_left[0]*conj(tree_sum_right);
                tree_product[1]+=loop_left[1]*conj(tree_sum_right);
                tree_product[2]+=loop_left[2]*conj(tree_sum_right);
            }
            res[0]+=tree_product[0];
            res[1]+=tree_product[1];
            res[2]+=tree_product[2];
    }

    vector<complex<T> >res_v;
    res_v.push_back(complex<T>(res[2]));
    res_v.push_back(complex<T>(res[1]));
    res_v.push_back(complex<T>(res[0]));

    return SeriesC<T>(-2,0,res_v);
}
    
    
#if 1
void Squared_ME_Struct::set_partial_born(){

    for(int i=0;i<m_tree_vals.size();i++){
        m_tree_vals_re.push_back(0);
        m_tree_vals_im.push_back(0);
        m_kin_ps_tree_val.push_back(std::vector<int>());
    }

   for( std::map<std::pair<process, std::vector<int> >, int >::iterator iter= map_tree_users.begin(); iter != map_tree_users.end(); ++iter){
        this->add(new Ampl_Info((iter->first).first,
                                (iter->first).second,
                                 &(m_tree_vals_re[iter->second]),
                                 &(m_tree_vals_im[iter->second]))); 
   };

   //redundant but easy 
   for(int i=0;i<m_trees.size();i++) m_kin_ps_tree_val[m_trees[i].second]=m_trees[i].first;

    return;
}
#endif
   
Squared_ME_Struct::~Squared_ME_Struct(){

    for (int j=0;j<partial_tree_ampls.size();j++){
	delete partial_tree_ampls[j];
	}
    //cleared within cache
    //for (int j=0;j<m_loop_users.size();j++){
    //	if(m_loop_users[j]!=0) delete m_loop_users[j];
    //	}
    for (int j=0;j<m_tree_users.size();j++){
        delete m_tree_users[j];
    }
    for (int j=0;j<m_kin_fn.size();j++){
	if(m_kin_fn[j]!=0) delete m_kin_fn[j];
    }
};


//Virtual_SME::Virtual_SME(const Virtual_SME& vm){
//	for (int j=0;j<vm._MEs.size();j++){
//		_MEs.push_back(new Squared_ME(*vm._MEs[j]));
//	}
//}
//
//Virtual_SME& Virtual_SME::operator=(const Virtual_SME& vm){
//	for (int j=0;j<_MEs.size();j++){
//		delete _MEs[j];
//	}
//	_MEs.clear();
//
//	for (int j=0;j<vm._MEs.size();j++){
//		_MEs.push_back(new Squared_ME(*vm._MEs[j]));
//	}
//	return *this;
//
//}


void Virtual_SME::set_partial_born(){

		for(int j=0;j<_MEs.size();j++){
			_MEs[j]->set_partial_born();
            std::vector<Ampl_Info* > ampl_info(_MEs[j]->get_partial_born());
            for(int i=0;i<ampl_info.size();i++){
    			partial_tree_ampls.push_back(ampl_info[i]);
            }
        }
return;
}

void Virtual_SME::get_partial_born_map(vector<vector<int> >& permutation,
                vector<vector<int> >& helicity){
    
    permutation.clear();
    helicity.clear();

    for(int i=0;i<partial_tree_ampls.size();i++){
        permutation.push_back(partial_tree_ampls[i]->m_perm); 
        helicity.push_back(partial_tree_ampls[i]->m_hels); 
    }

return;
}

void Virtual_SME::get_vals_partial_born(vector<double* >& re_born,
    vector<double* >& im_born){
    
    re_born.clear();
    im_born.clear();

    for(int i=0;i<partial_tree_ampls.size();i++){
        re_born.push_back(partial_tree_ampls[i]->real); 
        im_born.push_back(partial_tree_ampls[i]->imag); 
    }

    return;

}



Virtual_SME::~Virtual_SME(){
	for (int j=0;j<_MEs.size();j++){
		delete _MEs[j];
	}
	_MEs.clear();

}


template <class T> SeriesC<T> Virtual_SME::eval_fn(momentum_configuration<T>& mc, T mu){

	SeriesC<T> res(-2,0);
	int mu_index = DefineMu<T>(mc,mu);
	for (int j=0;j<_MEs.size();j++){
		res+= _MEs[j]->eval(mc,d_index_vector,mu_index);
		BH_DEBUG_MESSAGE2("MEs:",_MEs[j]->eval(mc,d_index_vector,mu_index));
	}
	return res;
}

void Virtual_SME::set_couplings(R couplings){
    for (int j=0;j<_MEs.size();j++){
        _MEs[j]->set_couplings(couplings);
    }
}


template <class T> complex<T> Virtual_SME::eval_tree_fn(momentum_configuration<T>& mc){
	complex<T> res(0,0);
#if BH_USE_OMP
	std::vector<std::complex<T> > tempres(_MEs.size());
	std::vector<Squared_ME*> shuffle(_MEs);
	std::random_shuffle(shuffle.begin(),shuffle.end());
	BH_DEBUG_MESSAGE2("Number of trees in Virtual_SME::eval_tree: ",_MEs.size());
#pragma omp parallel
	{
	#pragma omp for 
	for (int j=0;j<_MEs.size();j++){
		std::complex<T> temp =shuffle[j]->eval_tree(mc,d_index_vector);
		tempres[j]=temp;
	}
	}

	for (int j=0;j<_MEs.size();j++){
		res+= tempres[j];
	}

#else
		for (int j=0;j<_MEs.size();j++){
			res+= _MEs[j]->eval_tree(mc,d_index_vector);
		}
#endif
return res;
}

void Virtual_SME_with_prefactor::add_amplitude_prefactor(kinematic_function* kin){
	d_prefactor_p=kin;
}



template <class T> SeriesC<T> Virtual_SME_with_prefactor::eval_fn(momentum_configuration<T>& mc, T mu){

	SeriesC<T> res(-2,0);
	res=Virtual_SME::eval_fn<T>(mc,mu);
	complex<T> prefactor=d_prefactor_p->eval(mc);
	return conj(prefactor)*prefactor*res;
}

template <class T> complex<T> Virtual_SME_with_prefactor::eval_tree_fn(momentum_configuration<T>& mc){
	complex<T> res = Virtual_SME::eval_tree_fn<T>(mc);
	complex<T> prefactor = d_prefactor_p->eval(mc);
	return conj(prefactor)*prefactor*res;
return res;
}

Virtual_SME_with_prefactor::~Virtual_SME_with_prefactor(){
}



void Virtual_SME::dry_run(){

	for (int j=0;j<_MEs.size();j++){
		_MEs[j]->dry_run(d_index_vector);
	}
}


void Virtual_SME::add(Squared_ME* sm){

	_MEs.push_back(sm);
	if (d_index_vector.size() == 0){
		for (int i=1;i<=sm->get_nbr_external();i++){d_index_vector.push_back(i);};
	}
}

template <class T> complex<T> prop_hel_fn_diphoton::eval_fn(momentum_configuration<T>& mc){
		return (complex<T>(m_Q2,0)/complex<T>(9,0));
}

template <class T> complex<T> prop_hel_fn::eval_fn(momentum_configuration<T>& mc){
	switch(m_photonZW){
		case 0:	{return(-complex<T>(m_Q,0)/complex<T>(3,0));} break;
		case 1: {complex<T> s=mc.s(m_i,m_j);
			return(-complex<T>(m_Q)/complex<T>(3,0)+complex<T>(m_hel_coupling,0)*s/(s
						-T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(constants::MZ*constants::MZ,0)
						+T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(0,constants::GZ*constants::MZ)));
			}	break;
		case 2: {complex<T> s=mc.s(m_i,m_j);
			return(complex<T>(m_hel_coupling,0)*s/(s
						-T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(constants::MZ*constants::MZ,0)
						+T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(0,constants::GZ*constants::MZ)));
			}	break;
		case 3: {complex<T> s=mc.s(m_i,m_j);
			return(complex<T>(m_hel_coupling,0)*s/(s
						-T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(constants::MW*constants::MW,0)
						+T(constants::s_GeV)*T(constants::s_GeV)*complex<T>(0,constants::GW*constants::MW)));
			}	break;
	}
}


template complex<R> prop_hel_fn::eval_fn(momentum_configuration<R>& mc);
template complex<RHP> prop_hel_fn::eval_fn(momentum_configuration<RHP>& mc);
template complex<RVHP> prop_hel_fn::eval_fn(momentum_configuration<RVHP>& mc);


template complex<R> prop_hel_fn_diphoton::eval_fn(momentum_configuration<R>& mc);
template complex<RHP> prop_hel_fn_diphoton::eval_fn(momentum_configuration<RHP>& mc);
template complex<RVHP> prop_hel_fn_diphoton::eval_fn(momentum_configuration<RVHP>& mc);

template SeriesC<R> subtraction::eval(momentum_configuration<R>& mc);
template SeriesC<RHP> subtraction::eval(momentum_configuration<RHP>& mc);
template SeriesC<RVHP> subtraction::eval(momentum_configuration<RVHP>& mc);


}
