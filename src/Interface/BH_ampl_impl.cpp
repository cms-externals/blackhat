/*
 * BH_Ampl_impl.cpp
 *
 *  Created on: 28 Jul 2009
 *      Author: daniel
 */

#include "Interface/BH_ampl_impl.h"
#include "assembly.h"
#include "settings.h"
#include "Interface/BH_interface_impl.h"
#include <algorithm>

//if HP result for tree mathches double-result up to 1/BORNPRECISION no VHP is called even if trees are large.
#define BORNPRECISION 10000 
#define _MEASURE_TIME 0

using namespace std;

namespace BH {


void BH_Ampl_impl::print_virt_info(){
	
	mom_conf& mc=*d_parent_p->get_mc();
	R couplings=
			std::pow(12.56637061435917*constants::alpha_S,d_PowersOfAlphaS)*
			std::pow(12.56637061435917*constants::alpha_QED,d_PowersOfAlphaQED)*
		std::pow(constants::s_GeV,d_GeVdim);
	std::cout<<"*****************************************"<<std::endl;
	std::cout<<"PS point with large ME2 (multi-precision):"<<std::endl;
	for(int i=0;i<d_NbrExtParticles;i++){
	std::cout<<i+1<<": ("<<1./constants::s_GeV*mc.p(i+1).E().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).Z().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).X().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).Y().real() <<");"<<std::endl;
	}
	std::cout<<"process particle IDs: ";print_process();
	std::cout<<"normalized virt. ME2: "<<std::endl;
	std::cout<<"   1/eps^2: "<< d_factor*d_last_result_2<<std::endl;
	std::cout<<"   1/eps^1: "<< d_factor*d_last_result_1<<std::endl;
	std::cout<<"   finite:  "<< d_factor*d_last_result_0<<std::endl;
	std::cout<<"born ME2 (couplings included): "<< d_factor*couplings*d_last_result_tree<<std::endl;
	std::cout<<"*****************************************"<<std::endl;
}

void BH_Ampl_impl::print_born_info(){
	
	mom_conf& mc=*d_parent_p->get_mc();
	R couplings=
			std::pow(12.56637061435917*constants::alpha_S,d_PowersOfAlphaS)*
			std::pow(12.56637061435917*constants::alpha_QED,d_PowersOfAlphaQED)*
		std::pow(constants::s_GeV,d_GeVdim);
	std::cout<<"*****************************************"<<std::endl;
	std::cout<<"PS point with large ME2 (VHP):"<<std::endl;
	for(int i=0;i<d_NbrExtParticles;i++){
	std::cout<<i+1<<": ("<<1./constants::s_GeV*mc.p(i+1).E().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).Z().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).X().real() <<","
			   <<1./constants::s_GeV*mc.p(i+1).Y().real() <<");"<<std::endl;
	}
	std::cout<<"process particle IDs: ";print_process();
	std::cout<<"born ME2 (couplings included): "<<std::endl;
	std::cout<<"   born:            "<< d_factor*couplings*d_last_result_tree<<std::endl;
	std::cout<<"   born(HP):            "<< d_factor*couplings*to_double(d_last_result_HP_tree)<<std::endl;
	std::cout<<"   born(VHP):           "<< d_factor*couplings*to_double(d_last_result_VHP_tree)<<std::endl;
	std::cout<<"   born(no couplings):  "<< d_last_result_tree<<std::endl;
	std::cout<<"   catch_bound:         "<< d_catch_large_born_ME2<<std::endl;
	std::cout<<"*****************************************"<<std::endl;
}
	
double BH_Ampl_impl::get_finite(){

	mom_conf& mc=*d_parent_p->get_mc();

	if(d_mc_ID_virt==(mc).get_ID()){return (d_factor*d_last_result_0).real();}
	else {d_mc_ID_virt=(mc).get_ID();d_mc_ID_born=d_mc_ID_virt;};

#if _MEASURE_TIME
clock_t before, after;
before=clock();
#endif

	d_last_result_tree =d_VSM_p->eval_tree(mc);
	if(((std::abs(d_last_result_tree)>d_catch_large_born_ME2)||
                (d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points))&& settings::general::s_use_higher_precision){
		vector<int> indices(d_NbrExtParticles);
		for (int ii=1;ii<=d_NbrExtParticles;ii++){
			indices[ii-1]=ii;
		}
		mom_conf_HP mcHP=mc.extend<RHP>(indices);

		d_last_result_HP_tree=d_VSM_p->eval_tree(mcHP);
		C last_result_tree=d_last_result_tree;
		d_last_result_tree=to_double(d_last_result_HP_tree);
			
			if(std::abs(last_result_tree/d_last_result_tree-1.)<1./double(BORNPRECISION)){
				if(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points){
					d_eval_id=d_eval_id+1;
					d_catch_large_born_ME2=d_catch_large_born_ME2+std::abs(d_last_result_tree)/double(BH::settings::BH_interface_settings::s_number_of_warmup_points)*d_catch_large_ME2;
				}
			}
			else if((std::abs(d_last_result_tree)>d_catch_large_born_ME2)||(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points)){
			vector<int> indices(d_NbrExtParticles);
			for (int ii=1;ii<=d_NbrExtParticles;ii++){
				indices[ii-1]=ii;
			}
			mom_conf_VHP mcVHP=mc.extend<RVHP>(indices);

		d_last_result_VHP_tree=d_VSM_p->eval_tree(mcVHP);
		C new_result_tree=to_double(d_last_result_VHP_tree);
			
			if(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points){
				d_eval_id=d_eval_id+1;
				d_catch_large_born_ME2=d_catch_large_born_ME2+std::abs(new_result_tree)/double(BH::settings::BH_interface_settings::s_number_of_warmup_points)*d_catch_large_ME2;
			}
			if((std::abs(new_result_tree)>d_catch_large_born_ME2)&&(d_eval_id>=BH::settings::BH_interface_settings::s_number_of_warmup_points)){
				print_born_info();
			}
		d_last_result_tree=new_result_tree;
		}
	}	

	// just regular computation
	SeriesC<R> res(-2,0);
	res=d_VSM_p->eval(mc,d_parent_p->get_mu());

	d_last_result_2 = res[-2]/d_last_result_tree/2.;
	d_last_result_1 = res[-1]/d_last_result_tree/2.;
	d_last_result_0 = res[0]/d_last_result_tree/2.;

	if(std::abs(d_last_result_2)>d_catch_large_ME2||
		std::abs(d_last_result_1)>d_catch_large_ME2||
		std::abs(d_last_result_0)>d_catch_large_ME2){
		print_virt_info();
	}

#if _MEASURE_TIME
after=clock();
cout<< "total evaluation time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;
#endif
 
 counter_manager::s_cm.print();
	return (d_factor*d_last_result_0).real()-d_scheme_shift;
}

double BH_Ampl_impl::get_born(){

#if _MEASURE_TIME
clock_t before, after;
before=clock();
#endif

	mom_conf& mc=*d_parent_p->get_mc();


	R couplings=
			std::pow(12.56637061435917*constants::alpha_S,d_PowersOfAlphaS)*
			std::pow(12.56637061435917*constants::alpha_QED,d_PowersOfAlphaQED)*
		std::pow(constants::s_GeV,d_GeVdim);

	if(d_mc_ID_born==(mc).get_ID()){return (d_factor*d_last_result_tree).real()*couplings;}
	else {d_mc_ID_born=(mc).get_ID();};

    //improtant for rescaling primitive amplitudes for subtraction terms
    d_VSM_p->set_couplings(sqrt(couplings));
	
	d_last_result_tree =d_VSM_p->eval_tree(mc);
	
	//if tree is large call higher precision trees
	if(((std::abs(d_last_result_tree)>d_catch_large_born_ME2)||(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points))&& settings::general::s_use_higher_precision){
		vector<int> indices(d_NbrExtParticles);
		for (int ii=1;ii<=d_NbrExtParticles;ii++){
			indices[ii-1]=ii;
		}
		mom_conf_HP mcHP=mc.extend<RHP>(indices);
	
	
		d_last_result_HP_tree =d_VSM_p->eval_tree(mcHP);
		C last_result_tree=d_last_result_tree;
		d_last_result_tree =to_double(d_last_result_HP_tree);

		if(std::abs(last_result_tree/d_last_result_tree-1.)<1./double(BORNPRECISION)){
			if(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points){
				d_eval_id=d_eval_id+1;
				d_catch_large_born_ME2=d_catch_large_born_ME2+std::abs(d_last_result_tree)/double(BH::settings::BH_interface_settings::s_number_of_warmup_points)*d_catch_large_ME2;
			}
		}
		else if((std::abs(d_last_result_tree)>d_catch_large_born_ME2)||(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points)){
			vector<int> indices(d_NbrExtParticles);
			for (int ii=1;ii<=d_NbrExtParticles;ii++){
				indices[ii-1]=ii;
			}
			mom_conf_VHP mcVHP=mc.extend<RVHP>(indices);
	
		d_last_result_VHP_tree =d_VSM_p->eval_tree(mcVHP);
		C new_result_tree=to_double(d_last_result_VHP_tree);
			
			if(d_eval_id<BH::settings::BH_interface_settings::s_number_of_warmup_points){
				d_eval_id=d_eval_id+1;
				d_catch_large_born_ME2=d_catch_large_born_ME2+std::abs(new_result_tree)/double(BH::settings::BH_interface_settings::s_number_of_warmup_points)*d_catch_large_ME2;}

			if((std::abs(new_result_tree)>d_catch_large_born_ME2)&&(d_eval_id>BH::settings::BH_interface_settings::s_number_of_warmup_points)){
				print_born_info();
			}
		d_last_result_tree=new_result_tree;
		}
	}


#if _MEASURE_TIME
print_born_info();
after=clock();
cout<< "total evaluation time:"<< double(after-before)/double(CLOCKS_PER_SEC) <<endl;
time_t seconds;
seconds = time (NULL);
cout<< "absolute time since January 1, 1970: "<< seconds <<endl;
#endif


	return (d_factor*d_last_result_tree).real()*couplings;
}


double BH_Ampl_impl::get_single_pole(){
	return (d_factor*d_last_result_1).real()-d_renormalization_shift;
}

double BH_Ampl_impl::get_double_pole(){

	return (d_factor*d_last_result_2).real();
}

double BH_Ampl_impl::get_finite_HP(){

	mom_conf& mc=*d_parent_p->get_mc();

	vector<int> indices(d_NbrExtParticles);
	for (int ii=1;ii<=d_NbrExtParticles;ii++){
		indices[ii-1]=ii;
	}
	mom_conf_HP mcHP=mc.extend<RHP>(indices);


// just regular computation with HP

	CHP treeHP =d_VSM_p->eval_tree(mcHP);

	SeriesC<RHP> resHP(-2,0);
	resHP=d_VSM_p->eval(mcHP,RHP(d_parent_p->get_mu()));


	d_last_result_HP_2 = to_double(resHP[-2]/treeHP)/2.;
	d_last_result_HP_1 = to_double(resHP[-1]/treeHP)/2.;
	d_last_result_HP_0 = to_double(resHP[0]/treeHP)/2.;

	return (d_factor*d_last_result_HP_0).real()-d_scheme_shift;
}

double BH_Ampl_impl::get_born_HP(){

	mom_conf& mc=*d_parent_p->get_mc();

	vector<int> indices(d_NbrExtParticles);
	for (int ii=1;ii<=d_NbrExtParticles;ii++){
		indices[ii-1]=ii;
	}
	mom_conf_HP mcHP=mc.extend<RHP>(indices);
	
	R couplings=
			std::pow(12.56637061435917*constants::alpha_S,d_PowersOfAlphaS)*
			std::pow(12.56637061435917*constants::alpha_QED,d_PowersOfAlphaQED)*
		std::pow(constants::s_GeV,d_GeVdim);

    //improtant for rescaling primitive amplitudes for subtraction terms
    d_VSM_p->set_couplings(sqrt(couplings));
	
	CHP treeHP =d_VSM_p->eval_tree(mcHP);
	d_last_result_tree =to_double(treeHP);

	return (d_factor*d_last_result_tree).real()*couplings;
}

double BH_Ampl_impl::get_single_pole_HP(){
	return (d_factor*d_last_result_HP_1).real()-d_renormalization_shift;
}

double BH_Ampl_impl::get_double_pole_HP(){

	return (d_factor*d_last_result_HP_2).real();
}

double BH_Ampl_impl::get_finite_VHP(){

	mom_conf& mc=*d_parent_p->get_mc();
	vector<int> indices(d_NbrExtParticles);
	for (int ii=1;ii<=d_NbrExtParticles;ii++){
		indices[ii-1]=ii;
	}
	mom_conf_VHP mcVHP=mc.extend<RVHP>(indices);

// just regular computation with VHP
	CVHP treeVHP =d_VSM_p->eval_tree(mcVHP);
	
	SeriesC<RVHP> resVHP(-2,0);
	resVHP=d_VSM_p->eval(mcVHP,RVHP(d_parent_p->get_mu()));


	d_last_result_VHP_2 = to_double(resVHP[-2]/treeVHP)/2.;
	d_last_result_VHP_1 = to_double(resVHP[-1]/treeVHP)/2.;
	d_last_result_VHP_0 = to_double(resVHP[0]/treeVHP)/2.;

	return (d_factor*d_last_result_VHP_0).real()-d_scheme_shift;
}


double BH_Ampl_impl::get_born_VHP(){

	mom_conf& mc=*d_parent_p->get_mc();

	vector<int> indices(d_NbrExtParticles);
	for (int ii=1;ii<=d_NbrExtParticles;ii++){
		indices[ii-1]=ii;
	}
	mom_conf_VHP mcVHP=mc.extend<RVHP>(indices);
	
	R couplings=
			std::pow(12.56637061435917*constants::alpha_S,d_PowersOfAlphaS)*
			std::pow(12.56637061435917*constants::alpha_QED,d_PowersOfAlphaQED)*
		std::pow(constants::s_GeV,d_GeVdim);
    //improtant for rescaling primitive amplitudes for subtraction terms
    d_VSM_p->set_couplings(sqrt(couplings));
	
	CVHP treeVHP =d_VSM_p->eval_tree(mcVHP);
	d_last_result_tree =to_double(treeVHP);

	return (d_factor*d_last_result_tree).real()*couplings;
}



double BH_Ampl_impl::get_single_pole_VHP(){
	return (d_factor*d_last_result_VHP_1).real()-d_renormalization_shift;
}

double BH_Ampl_impl::get_double_pole_VHP(){

	return (d_factor*d_last_result_VHP_2).real();
}


void BH_Ampl_impl::set_partial_born(){
    d_VSM_p->set_partial_born(); 
    return;
};


void BH_Ampl_impl::get_map(
        vector<vector<int> >& permutation,
        vector<vector<int> >& helicity){d_VSM_p->get_partial_born_map(permutation,helicity);return;};


void BH_Ampl_impl::get_vals(
        vector<double* >& re_ampl,
        vector<double* >& im_ampl){d_VSM_p->get_vals_partial_born(re_ampl,im_ampl);return;};



double BH_Ampl_impl::getScaleVariationCoefficient(int logMuOrder){

	switch (logMuOrder) {
	case 0: return get_finite();
	case 1:	return get_single_pole()+(d_PowersOfAlphaS /*NbrPowersOfAlphaS*/* 23/6. /*beta_0*/);
	case 2:	return get_double_pole();
	}


};




}
