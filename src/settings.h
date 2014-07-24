/*
 * settings.h
 *
 *  Created on: 13-Nov-2008
 *      Author: daniel
 */

#ifndef SETTINGS_H_
#define SETTINGS_H_

#include <string>
#include <iostream>
#include <vector>
#include <string>
#include "setable.h"

namespace BH {

namespace settings {

struct skeleton_settings {
	static bool s_use_skeleton_in_rat_ext;
	static void use_use_skeleton_in_rat_ext(){s_use_skeleton_in_rat_ext=true;};
	static void dont_use_use_skeleton_in_rat_ext(){s_use_skeleton_in_rat_ext=false;};
};

struct general {
	enum cut_type { Darren, FHZ, Darren_wS, FHZ_wS, worker };
	enum rat_type { ratext, recursive, ratext_worker/*, ratmix*/ };
	static cut_type s_cut_type;
	static rat_type s_rat_type;
	static bool s_use_known_formulae;
	static bool s_use_cached_integrals;
	static bool s_use_parent_files;
	static bool s_show_parent_diagrams;
	static bool s_generate_worker_tree;
	static std::string s_data_path;
	static std::string s_parent_data_path;
	static std::string s_assembly_data_path;
	static double s_bub_cut_precision;
	static bool s_use_check_in_cut_part;
	static bool s_use_ep_only;
    static bool s_use_new_rec_tree;
	static bool s_use_higher_precision;
	static bool s_use_g3_coupling;
	static bool s_normalise_debug_coeff_output;
  static bool s_use_only_BG_trees;
};

struct spurious_poles_settings {
	enum cut_part_type { Darren, FHZ };
	static cut_part_type s_cut_part_type;
};

struct rational_settings {
	static bool s_skip_sub_leading_color;
	static double s_rat_ext_precision;
	static bool s_use_known_formulae_in_ratext;
	static bool s_use_IR_in_ratext;
	static bool s_use_data_files;
	static bool s_set_all_zero;

};

struct BH_interface_settings {
	enum BH_color_mode { full_color, leading_color, full_minus_leading_color };
	static BH_color_mode s_BH_color_mode;
	static BH_color_mode s_BH_tree_color_mode;
	static int s_nc;
	static int s_nf;
	enum BH_interface_mode { normal, gridWarmup, collectPS, echo, cached };
	enum BH_interface_eval_mode { normal_eval, polarization_eval };
	static BH_interface_mode s_BH_interface_mode;
	static BH_interface_eval_mode s_BH_interface_eval_mode;
	static std::string s_PS_collection_filename;
	static std::string s_echo_input_filename;
	static int s_file_base;
	static int s_photon_only;
	static double s_catch_large_ME2;
	static int s_number_of_warmup_points;
	static int s_w5jet_number_of_threads;
	static bool s_final_state_gluon_symmetrization;
	static bool s_born_final_state_gluon_symmetrization;
	static bool s_use_automated_assembly;
	static bool s_use_symmetrized_assembly_files;
	static bool s_generate_assembly_files;
	static bool s_use_W_polarization_A567_assembly;
    static int s_W_polarization_A;
    static bool s_use_grassman_trees;
    static int s_same_helicity_projection;
};


}



class setting_base {
	void* d_value_p;
	std::string d_name;
public:
	const std::string& get_name() const { return d_name;}
	void* get_address(){return d_value_p;};
	virtual void  print_value(std::ostream& )=0;
	template <class T> void set_address(T& ref){d_value_p=&ref;}
	setting_base(const std::string& name): d_name(name) {};
	virtual ~setting_base(){}
};

template <class T> class setting : public setting_base {
	T d_value;
public:
	setting(std::string name, T& ref) : setting_base(name), d_value(ref)  {set_address(d_value); } ;
	void print_value(std::ostream& os ){ os << d_value; };
};

template <class T> void do_apply_setting(setting_base* se, T& param) {
	param=(*static_cast<T*>(se->get_address()));
}



class settings_table : public setable {
	std::vector<setting_base*> d_settings_list;
public:
	template <class T> bool apply_setting(const std::string& name, T& param);
	//! adds a setting object (by making a copy)
	/** if the entry exists, it wulll be replaced */
	template <class T> void add(const std::string&,T value);
	virtual bool set(const std::string& name,int value){ return do_set(name,value);}
	virtual bool set(const std::string& name,bool value){ return do_set(name,value);}
	virtual bool set(const std::string& name,double value){ return do_set(name,value);}
	virtual bool set(const std::string& name,std::string value){ return do_set(name,value);}
	void display(std::ostream& os = std::cout , const std::string& head = s_default_head );
	~settings_table();
private:
	//! sets the value of an existing setting
	/** if the entry doesn't exist, a warning is issued, but nothing is done */
	template <class T> bool do_set(const std::string&,T value);
	static std::string s_default_head;
};

struct setting_name_is : public std::unary_function<setting_base*,bool> {
	std::string d_name;
	setting_name_is(const std::string& name): d_name(name){};
	bool operator()(setting_base* se){return se->get_name() == d_name;};
};




}
#endif /* SETTINGS_H_ */
