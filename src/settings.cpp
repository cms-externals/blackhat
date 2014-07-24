/*
 * settings.cpp
 *
 *  Created on: 13-Nov-2008
 *      Author: daniel
 */

#include "settings.h"
#include "BH_utilities.h"
#include <algorithm>

using namespace std;

namespace BH {



template <class T> bool settings_table::apply_setting(const string& name,T& param){
	vector<setting_base*>::iterator it=find_if(d_settings_list.begin(),d_settings_list.end(),setting_name_is(name));
	if ( it != d_settings_list.end() ) {
		do_apply_setting( (*it), param);
		return true;
	} else {
		_WARNING3("Setting \"",name,"\" not found." );
		return false;
	}
}

template <class T> void settings_table::add(const std::string& name,T value){
	vector<setting_base*>::iterator it=find_if(d_settings_list.begin(),d_settings_list.end(),setting_name_is(name));
	if ( it != d_settings_list.end() ) {
		delete (*it);*it=new setting<T>(name,value);
	} else {
		d_settings_list.push_back(new setting<T>(name,value));

	}
}

template <class T> bool settings_table::do_set(const std::string& name,T value){
	vector<setting_base*>::iterator it=find_if(d_settings_list.begin(),d_settings_list.end(),setting_name_is(name));
	if ( it != d_settings_list.end() ) {
		delete (*it);*it=new setting<T>(name,value);
		return true;
	} else {
		_WARNING3("Setting ",name," is not in the list of settings.");
		return false;
	}
}

settings_table::~settings_table(){
	for ( std::vector<setting_base*>::iterator it= d_settings_list.begin();it!=d_settings_list.end();++it){
		delete *it;
	}
}


struct print_setting {
	ostream& d_ostream;
	std::string d_head;
	print_setting(ostream& os,const string& head): d_ostream(os), d_head(head){};
	bool operator()(setting_base* se){ d_ostream << d_head << se->get_name()<<  ": " ; se->print_value(d_ostream); d_ostream << endl;};
};

void settings_table::display(ostream& os, const string& head){
		for_each(d_settings_list.begin(),d_settings_list.end(),print_setting(os,head));
}

std::string settings_table::s_default_head;

template void settings_table::add(const std::string& name,double value);
template void settings_table::add(const std::string& name,int value);
template void settings_table::add(const std::string& name,bool value);
template void settings_table::add(const std::string& name,std::string value);


template bool settings_table::apply_setting(const string& name,int& param);
template bool settings_table::apply_setting(const string& name,bool& param);
template bool settings_table::apply_setting(const string& name,double& param);
template bool settings_table::apply_setting(const string& name,string& param);


namespace settings {

bool skeleton_settings::s_use_skeleton_in_rat_ext=false;
spurious_poles_settings::cut_part_type spurious_poles_settings::s_cut_part_type=spurious_poles_settings::Darren;
bool rational_settings::s_skip_sub_leading_color=false;
bool rational_settings::s_use_known_formulae_in_ratext=true;
bool rational_settings::s_use_IR_in_ratext=true;
bool rational_settings::s_use_data_files=true;
bool rational_settings::s_set_all_zero=false;
double rational_settings::s_rat_ext_precision=4.;

BH_interface_settings::BH_interface_mode BH_interface_settings::s_BH_interface_mode=BH_interface_settings::normal;
BH_interface_settings::BH_interface_eval_mode BH_interface_settings::s_BH_interface_eval_mode=BH_interface_settings::normal_eval;
std::string BH_interface_settings::s_PS_collection_filename="PSpoints.dat";
std::string BH_interface_settings::s_echo_input_filename="MEs.";
int BH_interface_settings::s_file_base=200;
#ifndef BH_PUBLIC
general::cut_type general::s_cut_type=general::Darren;
general::rat_type general::s_rat_type=general::ratext;
#else
general::cut_type general::s_cut_type=general::worker;
general::rat_type general::s_rat_type=general::ratext_worker;
#endif
string general::s_data_path="not set";
string general::s_parent_data_path="not set";
string general::s_assembly_data_path="not set";

double general::s_bub_cut_precision=4.;
bool general::s_use_check_in_cut_part=true;

bool general::s_use_known_formulae=true;
bool general::s_use_cached_integrals=true;
bool general::s_use_parent_files=false;
bool general::s_show_parent_diagrams=false;
bool general::s_generate_worker_tree=false;
bool general::s_use_higher_precision=true;

bool general::s_use_new_rec_tree=false;    
bool general::s_use_ep_only=true;
bool general::s_use_g3_coupling=false;
bool general::s_use_only_BG_trees=false;
bool general::s_normalise_debug_coeff_output=true;

BH_interface_settings::BH_color_mode BH_interface_settings::s_BH_color_mode=full_color;
BH_interface_settings::BH_color_mode BH_interface_settings::s_BH_tree_color_mode=full_color;
int BH_interface_settings::s_nc=3;
int BH_interface_settings::s_nf=5;
int BH_interface_settings::s_photon_only=1;   // 0 photon only; 1 for all other
double BH::settings::BH_interface_settings::s_catch_large_ME2=100000.;//to catch about 10% effect on 10^6 PS points
int BH::settings::BH_interface_settings::s_number_of_warmup_points=4;//to set number of warmup points to adapt catch large ME2
int BH::settings::BH_interface_settings::s_w5jet_number_of_threads=3;
bool BH_interface_settings::s_born_final_state_gluon_symmetrization=true;
bool BH_interface_settings::s_final_state_gluon_symmetrization=true;
bool BH::settings::BH_interface_settings::s_use_automated_assembly=false;
bool BH::settings::BH_interface_settings::s_generate_assembly_files=false;
bool BH::settings::BH_interface_settings::s_use_symmetrized_assembly_files=false;
bool BH::settings::BH_interface_settings::s_use_W_polarization_A567_assembly=false;
int BH::settings::BH_interface_settings::s_W_polarization_A=-1;
bool BH::settings::BH_interface_settings::s_use_grassman_trees=false;
int BH::settings::BH_interface_settings::s_same_helicity_projection=-1; //[-1,n] = [normal mode, keep only helicity configurations in virtual with sum_i h_i = n]

};




}

