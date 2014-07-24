#include "settings_reader.h"
#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include "BH_utilities.h"
#include "settings.h"
//#include "skeleton/skeleton.h"
#include "cached_OLHA.h"
#include "cut_part_factory.h"
#include <limits>

using namespace std;

namespace BH{

namespace settings {

enum setting_name {
	use_skeleton_in_ratext ,
	comment,
	use_known_formulae,
	use_cached_integrals,
	unknown_setting,
	spur_cut_type,
	co_cut_type,
	co_rat_type,
	all_rat_zero,
	skip_subleading,
	BH_color_mode,
	BH_tree_color_mode,
	Nc,
	Nf,
	photon_only,
	BH_interface_mode,
	BH_interface_eval_mode,
	ps_collection_filename,
	me_echo_filename,
	psp_nbr,
	rat_ext_precision,
	bub_cut_precision,
	use_known_formulae_in_ratext,
	use_IR_in_ratext,
	use_check_in_cut_part,
	ge_cut_type,
	ge_rat_type,
	data_path,
	gen_worker_tree,
	use_ep_only,
	parent_data_path,
	assembly_data_path,
	use_parent_files,
	show_parent_diagrams,
	use_automated_assembly,
	catch_large_ME2,
	number_of_warmup_points,
	w5jet_number_of_threads,
	final_state_gluon_symmetrization,
	born_final_state_gluon_symmetrization,
	use_higher_precision,
	use_symmetrized_assembly_files,
	generate_assembly_files,
    use_g3_coupling,
    use_only_BG_trees,
    use_W_polarization_A567_assembly,
    normalise_debug_coeff_output,
	W_polarization_A,
	use_grassman_trees,
    same_helicity_projection
};

map<string,setting_name> Init_settings_list(	){
	std::map<string,setting_name> set_map;
	set_map.insert(std::pair<string,setting_name>("USE_SKELETON_IN_RATEXT",use_skeleton_in_ratext));
	set_map.insert(std::pair<string,setting_name>("USE_KNOWN_FORMULAE",use_known_formulae));
	set_map.insert(std::pair<string,setting_name>("USE_CACHED_INTEGRALS",use_cached_integrals));
	set_map.insert(std::pair<string,setting_name>("COLHA_CUT_TYPE",co_cut_type));
	set_map.insert(std::pair<string,setting_name>("COLHA_RAT_TYPE",co_rat_type));
	set_map.insert(std::pair<string,setting_name>("SET_ALL_RAT_TO_ZERO",all_rat_zero));
	set_map.insert(std::pair<string,setting_name>("RAT_SKIP_SUBLEADING",skip_subleading));
	set_map.insert(std::pair<string,setting_name>("SPURIOUS_CUT_TYPE",spur_cut_type));
	set_map.insert(std::pair<string,setting_name>("COLOR_MODE",BH_color_mode));
	set_map.insert(std::pair<string,setting_name>("TREE_COLOR_MODE",BH_tree_color_mode));
	set_map.insert(std::pair<string,setting_name>("Nc",Nc));
	set_map.insert(std::pair<string,setting_name>("Nf",Nf));
	set_map.insert(std::pair<string,setting_name>("PHOTON_ONLY",photon_only));
	set_map.insert(std::pair<string,setting_name>("INTERFACE_MODE",BH_interface_mode));
	set_map.insert(std::pair<string,setting_name>("INTERFACE_EVAL_MODE",BH_interface_eval_mode));
	set_map.insert(std::pair<string,setting_name>("PS_COLLECTION_FILENAME",ps_collection_filename));
	set_map.insert(std::pair<string,setting_name>("ME_ECHO_FILENAME",me_echo_filename));
	set_map.insert(std::pair<string,setting_name>("NBR_PS_POINTS_PER_FILE",psp_nbr));
	set_map.insert(std::pair<string,setting_name>("SET_RATIONAL_PRECISION",rat_ext_precision));
	set_map.insert(std::pair<string,setting_name>("SET_CUT_PART_PRECISION",bub_cut_precision));
	set_map.insert(std::pair<string,setting_name>("USE_KNOWN_FORMULAE_IN_RATEXT",use_known_formulae_in_ratext));
	set_map.insert(std::pair<string,setting_name>("USE_IR_IN_RATEXT",use_IR_in_ratext));
	set_map.insert(std::pair<string,setting_name>("USE_CHECK_IN_CUT_PART",use_check_in_cut_part));
	set_map.insert(std::pair<string,setting_name>("CUT_TYPE",ge_cut_type));
	set_map.insert(std::pair<string,setting_name>("RAT_TYPE",ge_rat_type));
	set_map.insert(std::pair<string,setting_name>("DATA_PATH",data_path));
	set_map.insert(std::pair<string,setting_name>("GENERATE_WORKER_TREE",gen_worker_tree));
	set_map.insert(std::pair<string,setting_name>("USE_EVAL_PARAM_ONLY",use_ep_only));
	set_map.insert(std::pair<string,setting_name>("PARENT_DATA_PATH",parent_data_path));
	set_map.insert(std::pair<string,setting_name>("ASSEMBLY_DATA_PATH",assembly_data_path));
	set_map.insert(std::pair<string,setting_name>("USE_PARENT_FILES",use_parent_files));
	set_map.insert(std::pair<string,setting_name>("SHOW_PARENT_DIAGRAMS",show_parent_diagrams));
	set_map.insert(std::pair<string,setting_name>("CATCH_LARGE_ME2",catch_large_ME2));
	set_map.insert(std::pair<string,setting_name>("NUMBER_OF_WARMUP_POINTS",number_of_warmup_points));
	set_map.insert(std::pair<string,setting_name>("USE_AUTOMATED_ASSEMBLY",use_automated_assembly));
	set_map.insert(std::pair<string,setting_name>("USE_HIGHER_PRECISION",use_higher_precision));
	set_map.insert(std::pair<string,setting_name>("W5JET_NUMBER_OF_THREADS",w5jet_number_of_threads));
	set_map.insert(std::pair<string,setting_name>("FINAL_STATE_GLUON_SYMMETRIZATION",final_state_gluon_symmetrization));
	set_map.insert(std::pair<string,setting_name>("BORN_FINAL_STATE_GLUON_SYMMETRIZATION",born_final_state_gluon_symmetrization));
	set_map.insert(std::pair<string,setting_name>("USE_SYMMETRIZED_ASSEMBLY_FILES",use_symmetrized_assembly_files));
	set_map.insert(std::pair<string,setting_name>("GENERATE_ASSEMBLY_FILES",generate_assembly_files));
	set_map.insert(std::pair<string,setting_name>("USE_G3_COUPLING",use_g3_coupling));
	set_map.insert(std::pair<string,setting_name>("USE_ONLY_BG_TREES",use_only_BG_trees));
	set_map.insert(std::pair<string,setting_name>("USE_W_POLARIZATION_A567_ASSEMBLY",use_W_polarization_A567_assembly));
	set_map.insert(std::pair<string,setting_name>("W_POLARIZATION_A",W_polarization_A));
	set_map.insert(std::pair<string,setting_name>("USE_GRASSMAN_TREES",use_grassman_trees));
	set_map.insert(std::pair<string,setting_name>("SAME_HELICITY_PROJECTION",same_helicity_projection));

	return set_map;
}


template <class charT, class traits>
inline
std::basic_istream<charT,traits>&
ignoreLine (std::basic_istream<charT,traits>& strm)
{
    // skip until end-of-line
    strm.ignore(std::numeric_limits<int>::max(),strm.widen('\n'));

    // return stream for concatenation
    return strm;
}



setting_name get_command(istream& in){
	static std::map<string,setting_name> set_list=Init_settings_list();

	string the_command;
//	cout << "> ";
	in >> the_command;

	if (the_command.c_str()[0] == '#' ){
		in >> ignoreLine;
		return comment;

	}

	map<string,setting_name>::iterator pos;
	pos= set_list.find(the_command);

	if (pos!=set_list.end()){
		return pos->second;
	}
	else return unknown_setting;

}


bool read_answer(istream& in){
	string the_command;

	in >> the_command;
	static string yes_answers[]={"yes","YES","Yes","ON","on"};
	static string no_answers[]={"no","NO","No","OFF","off"};
	string* yes_answers_end=yes_answers+5;
	string* no_answers_end=no_answers+5;
	string* the_answer=find(yes_answers,yes_answers_end,the_command);
	if ( the_answer != yes_answers_end ){
		return true;
	}
	the_answer=find(no_answers,no_answers_end,the_command);
	if ( the_answer != no_answers_end ){
		return false;
	}
	_MESSAGE3("Sorry, could not understand your answer: ",the_command," assuming no.");
	return false;
}
// could provide more safety here
void read_string_answer(istream& in,string& the_answer){
	in >> the_answer;
}
// could provide more safety here
void read_int_answer(istream& in,int& the_answer){
	in >> the_answer;
}
// could provide more safety here
void read_double_answer(istream& in,double& the_answer){
	in >> the_answer;
}

//return 0 based index in the list of answers provided.
int read_answer(istream& in,vector<string>& possible_answers){
	string the_answer;

	in >> the_answer;
	vector<string>::iterator it=find(possible_answers.begin(),possible_answers.end(),the_answer);
	if ( it != possible_answers.end()){
		return it-possible_answers.begin();
	} else {
		cerr << "Choice " << the_answer << " is not available. Possible choices are:";
		copy(possible_answers.begin(),possible_answers.end(),ostream_iterator<string>(cerr,"\n"));
		return -1;
	}
};


bool read_from_stream(istream& in) {


while (!in.fail() && !in.eof()){
	setting_name the_command=get_command(in);
	if (!in.fail() && !in.eof()){
		switch (the_command) {
		case  use_skeleton_in_ratext : {
//			if ( read_answer(in) ){
//				BH::settings::skeleton_settings::s_use_skeleton_in_rat_ext=true;
//				_MESSAGE("USE_SKELETON_IN_RATEXT set to YES");
//			} else {
//				BH::settings::skeleton_settings::s_use_skeleton_in_rat_ext=false;
//				_MESSAGE("USE_SKELETON_IN_RATEXT set to NO");
//			};
			_MESSAGE("Deprecated option!");
			break;
		}
		case  use_known_formulae : {
			if ( read_answer(in) ){
				settings::general::s_use_known_formulae=true;
				_MESSAGE("USE_KNOWN_FORMULAE set to YES");
			} else {
				settings::general::s_use_known_formulae=false;
				_MESSAGE("USE_KNOWN_FORMULAE set to NO");
			};
			break;
		}
		case  use_cached_integrals : {
			if ( read_answer(in) ){
				settings::general::s_use_cached_integrals=true;
				_MESSAGE("USE_CACHED_INTEGRALS set to YES");
			} else {
				settings::general::s_use_cached_integrals=false;
				_MESSAGE("USE_CACHED_INTEGRALS set to NO");
			};
			break;
		}
		case  use_parent_files : {
			if ( read_answer(in) ){
				settings::general::s_use_parent_files=true;
				_MESSAGE("USE_PARENT_FILES set to YES");
			} else {
				settings::general::s_use_parent_files=false;
				_MESSAGE("USE_PARENT_FILES set to NO");
			};
			break;
		}
		case  show_parent_diagrams : {
			if ( read_answer(in) ){
				settings::general::s_show_parent_diagrams=true;
				_MESSAGE("SHOW_PARENT_DIAGRAMS set to YES");
			} else {
				settings::general::s_show_parent_diagrams=false;
				_MESSAGE("SHOW_PARENT_DIAGRAMS set to NO");
			};
			break;
		}
		case  all_rat_zero : {
			if ( read_answer(in) ){
				BH::settings::rational_settings::s_set_all_zero=true;
				_MESSAGE("SET_ALL_RAT_TO_ZERO set to YES");
			} else {
				BH::settings::rational_settings::s_set_all_zero=false;
				_MESSAGE("SET_ALL_RAT_TO_ZERO set to NO");
			};
			break;
		}
		case  skip_subleading : {
			if ( read_answer(in) ){
				BH::settings::rational_settings::s_skip_sub_leading_color=true;
				_MESSAGE("RAT_SKIP_SUBLEADING set to YES");
			} else {
				BH::settings::rational_settings::s_skip_sub_leading_color=false;
				_MESSAGE("RAT_SKIP_SUBLEADING set to NO");
			};
			break;
		}
		case  co_cut_type : {
			_WARNING("COLHA_CUT_TYPE is deprecated, please use CUT_TYPE instead! ");

			break;
		}
		case  co_rat_type : {
			_WARNING("COLHA_RAT_TYPE is deprecated, please use RAT_TYPE instead! ");
			//			static string possible_answers[]={"ratext","recursive"};
//			vector<string> ci_cut_type_choices(possible_answers,possible_answers+2);
//			switch ( read_answer(in,ci_cut_type_choices) ){
//			case 0:
//				BH::CachedOLHA::Cached_OLHA_factory_impl<OneLoopHelAmpl>::s_use_rational_type=BH::CachedOLHA::Cached_OLHA_factory::ratext;
//				BH::CachedOLHA::Cached_OLHA_factory_impl<IR_checked_OLHA>::s_use_rational_type=BH::CachedOLHA::Cached_OLHA_factory::ratext;
//				_MESSAGE("COLHA_RAT_TYPE set to: ratext"); break;
//			case 1:
//				BH::CachedOLHA::Cached_OLHA_factory_impl<OneLoopHelAmpl>::s_use_rational_type=BH::CachedOLHA::Cached_OLHA_factory::recursive;
//				BH::CachedOLHA::Cached_OLHA_factory_impl<IR_checked_OLHA>::s_use_rational_type=BH::CachedOLHA::Cached_OLHA_factory::recursive;
//				_MESSAGE("COLHA_RAT_TYPE set to: recursive"); break;
//			default:
//				_MESSAGE("nothing set for COLHA_RAT_TYPE."); break;
//			};
			break;
		}
		case  ge_cut_type : {
			static string possible_answers[]={"Darren","FHZ","Darren_wS","FHZ_wS","worker"};
			vector<string> ci_cut_type_choices(possible_answers,possible_answers+5);
			switch ( read_answer(in,ci_cut_type_choices) ){
			case 0:
				BH::settings::general::s_cut_type=BH::settings::general::Darren;
				_MESSAGE("CUT_TYPE set to: Darren"); break;
			case 1:
				BH::settings::general::s_cut_type=BH::settings::general::FHZ;
				_MESSAGE("CUT_TYPE set to: FHZ"); break;
			case 2:
//				BH::settings::general::s_cut_type=BH::settings::general::Darren_wS;
//				_MESSAGE("CUT_TYPE set to: Darren (with skeletons)"); break;
				_MESSAGE("Deprecated option!");
			case 3:
//				BH::settings::general::s_cut_type=BH::settings::general::FHZ_wS;
//				_MESSAGE("CUT_TYPE set to: FHZ (with skeletons)"); break;
				_MESSAGE("Deprecated option!");
			case 4:
				BH::settings::general::s_cut_type=BH::settings::general::worker;
				_MESSAGE("CUT_TYPE set to: worker"); break;
			default:
				_MESSAGE("nothing set for CUT_TYPE."); break;
			};
			break;
		}
		case  ge_rat_type : {
			static string possible_answers[]={"ratext","recursive","ratext_worker"};
			vector<string> ci_cut_type_choices(possible_answers,possible_answers+3);
			switch ( read_answer(in,ci_cut_type_choices) ){
			case 0:
				BH::settings::general::s_rat_type=BH::settings::general::ratext;
				_MESSAGE("RAT_TYPE set to: ratext"); break;
			case 1:
				BH::settings::general::s_rat_type=BH::settings::general::recursive;
				_MESSAGE("RAT_TYPE set to: recursive"); break;
			case 2:
				BH::settings::general::s_rat_type=BH::settings::general::ratext_worker;
				_MESSAGE("RAT_TYPE set to: ratext_worker"); break;
			default:
				_MESSAGE("nothing set for RAT_TYPE."); break;
			};
			break;
		}

		case  spur_cut_type : {
			static string possible_answers[]={"Darren","FHZ"};
			vector<string> ci_cut_type_choices(possible_answers,possible_answers+2);
			switch ( read_answer(in,ci_cut_type_choices) ){
			case 0:
				settings::spurious_poles_settings::s_cut_part_type=settings::spurious_poles_settings::Darren;
				_MESSAGE("SPURIOUS_CUT_TYPE set to: Darren"); break;
			case 1:
				settings::spurious_poles_settings::s_cut_part_type=settings::spurious_poles_settings::FHZ;
				_MESSAGE("SPURIOUS_CUT_TYPE set to: FHZ"); break;
			default:
				_MESSAGE("nothing set for SPURIOUS_CUT_TYPE."); break;
			};
			break;
		}
		case  BH_interface_mode : {
			static string possible_answers[]={"normal","gridWarmup","collectPS","echo","cached"};
			vector<string> choices(possible_answers,possible_answers+5);
			switch ( read_answer(in,choices) ){
			case 0:
				settings::BH_interface_settings::s_BH_interface_mode=settings::BH_interface_settings::normal;
				_MESSAGE("INTERFACE_MODE set to: normal"); break;
			case 1:
				settings::BH_interface_settings::s_BH_interface_mode=settings::BH_interface_settings::gridWarmup;
				_MESSAGE("INTERFACE_MODE set to: gridWarmup"); break;
			case 2:
				settings::BH_interface_settings::s_BH_interface_mode=settings::BH_interface_settings::collectPS;
				_MESSAGE("INTERFACE_MODE set to: collectPS"); break;
			case 3:
				settings::BH_interface_settings::s_BH_interface_mode=settings::BH_interface_settings::echo;
				_MESSAGE("INTERFACE_MODE set to: echo"); break;
			case 4:
				settings::BH_interface_settings::s_BH_interface_mode=settings::BH_interface_settings::cached;
				_MESSAGE("INTERFACE_MODE set to: cached"); break;
			default:
				_MESSAGE("nothing set for INTERFACE_MODE."); break;
			};
			break;
		}

		case  BH_interface_eval_mode : {
			static string possible_answers[]={"normal_eval","polarization_eval"};
			vector<string> choices(possible_answers,possible_answers+2);
			switch ( read_answer(in,choices) ){
			case 0:
				settings::BH_interface_settings::s_BH_interface_eval_mode=settings::BH_interface_settings::normal_eval;
				_MESSAGE("INTERFACE_EVAL_MODE set to: normal"); break;
			case 1:
				settings::BH_interface_settings::s_BH_interface_eval_mode=settings::BH_interface_settings::polarization_eval;
				_MESSAGE("INTERFACE_EVAL_MODE set to: gridWarmup"); break;
			default:
				_MESSAGE("nothing set for INTERFACE_EVAL_MODE."); break;
			};
			break;
		}
		case  BH_tree_color_mode : {
			static string possible_answers[]={"full_color","leading_color"};
			vector<string> choices(possible_answers,possible_answers+2);
			switch ( read_answer(in,choices) ){
			case 0:
				settings::BH_interface_settings::s_BH_tree_color_mode=settings::BH_interface_settings::full_color;
				_MESSAGE("TREE_COLOR_MODE set to: full_color"); break;
			case 1:
				settings::BH_interface_settings::s_BH_tree_color_mode=settings::BH_interface_settings::leading_color;
				_MESSAGE("TREE_COLOR_MODE set to: leading_color"); break;
			default:
				_MESSAGE("nothing set for TREE_COLOR_MODE."); break;
			};
			break;
		}

		case  BH_color_mode : {
			static string possible_answers[]={"full_color","leading_color","full_minus_leading_color"};
			vector<string> choices(possible_answers,possible_answers+3);
			switch ( read_answer(in,choices) ){
			case 0:
				settings::BH_interface_settings::s_BH_color_mode=settings::BH_interface_settings::full_color;
				_MESSAGE("COLOR_MODE set to: full_color"); break;
			case 1:
				settings::BH_interface_settings::s_BH_color_mode=settings::BH_interface_settings::leading_color;
				_MESSAGE("COLOR_MODE set to: leading_color"); break;
			case 2:
				settings::BH_interface_settings::s_BH_color_mode=settings::BH_interface_settings::full_minus_leading_color;
				_MESSAGE("COLOR_MODE set to: full_minus_leading_color"); break;
			default:
				_MESSAGE("nothing set for COLOR_MODE."); return false; break;
			};
			break;
		}

		case  Nc : {
			int nbr;
			read_int_answer(in,nbr);
			BH::settings::BH_interface_settings::s_nc=nbr;
			_MESSAGE2("Nc set to: ",nbr); break;
		};	break;

		case  Nf : {
			int nbr;
			read_int_answer(in,nbr);
			BH::settings::BH_interface_settings::s_nf=nbr;
			_MESSAGE2("Nf set to: ",nbr); break;
		};	break;

		case  photon_only : {
			if ( read_answer(in) ){
				BH::settings::BH_interface_settings::s_photon_only=0;
				_MESSAGE("PHOTON_ONLY set to YES");
			} else {
				BH::settings::BH_interface_settings::s_photon_only=1;
				_MESSAGE("PHOTON_ONLY set to NO");
			};
		};	break;

		case  ps_collection_filename : {
			string filename;
			read_string_answer(in,filename);
			settings::BH_interface_settings::s_PS_collection_filename=filename;
			_MESSAGE2("PS_COLLECTION_FILENAME set to: ",filename); break;
		};	break;
		case  me_echo_filename : {
			string filename;
			read_string_answer(in,filename);
			settings::BH_interface_settings::s_echo_input_filename=filename;
			_MESSAGE2("ME_ECHO_FILENAME set to: ",filename); break;
		};	break;
		case  data_path : {
			string filename;
			read_string_answer(in,filename);
			settings::general::s_data_path=filename;
			_MESSAGE2("DATA_PATH set to: ",filename); break;
		};	break;
		case  parent_data_path : {
			string filename;
			read_string_answer(in,filename);
			settings::general::s_parent_data_path=filename;
			_MESSAGE2("PARENT_DATA_PATH set to: ",filename); break;
		};	break;
		case  assembly_data_path : {
			string filename;
			read_string_answer(in,filename);
			settings::general::s_assembly_data_path=filename;
			_MESSAGE2("ASSEMBLY_DATA_PATH set to: ",filename); break;
		};	break;
		case  psp_nbr : {
			int nbr;
			read_int_answer(in,nbr);
			BH::settings::BH_interface_settings::s_file_base=nbr;
			_MESSAGE2("NBR_PS_POINTS_NBR set to: ",nbr); break;
		};	break;

		case  rat_ext_precision : {
			double nbr;
			read_double_answer(in,nbr);
			BH::settings::rational_settings::s_rat_ext_precision=nbr;
			_MESSAGE2("SET_RATIONAL_PRECISION set to: ",nbr); break;
		};	break;

		case  bub_cut_precision : {
			double nbr;
			read_double_answer(in,nbr);
			BH::settings::general::s_bub_cut_precision=nbr;
			_MESSAGE2("SET_CUT_PART_PRECISION set to: ",nbr); break;
		};	break;


		case  use_known_formulae_in_ratext : {
			if ( read_answer(in) ){
				BH::settings::rational_settings::s_use_known_formulae_in_ratext=true;
				_MESSAGE("USE_KNOWN_FORMULAE_IN_RATEXT set to YES");
			} else {
				BH::settings::rational_settings::s_use_known_formulae_in_ratext=false;
				_MESSAGE("USE_KNOWN_FORMULAE_IN_RATEXT set to NO");
			};
			break;
		}
		case  gen_worker_tree : {
			if ( read_answer(in) ){
				settings::general::s_generate_worker_tree=true;
				_MESSAGE("GENERATE_WORKER_TREE set to YES");
			} else {
				settings::general::s_generate_worker_tree=false;
				_MESSAGE("GENERATE_WORKER_TREE set to NO");
			};
			break;
		}
		case  use_IR_in_ratext : {
			if ( read_answer(in) ){
				BH::settings::rational_settings::s_use_IR_in_ratext=true;
				_MESSAGE("USE_IR_IN_RATEXT set to YES");
			} else {
				BH::settings::rational_settings::s_use_IR_in_ratext=false;
				_MESSAGE("USE_IR_IN_RATEXT set to NO");
			};
			break;
		}

		case  use_check_in_cut_part : {
			if ( read_answer(in) ){
				BH::settings::general::s_use_check_in_cut_part=true;
				_MESSAGE("USE_CHECK_IN_CUT_PART set to YES");
			} else {
				BH::settings::general::s_use_check_in_cut_part=false;
				_MESSAGE("USE_CHECK_IN_CUT_PART set to NO");
			};
			break;
		}
		case  use_higher_precision : {
			if ( read_answer(in) ){
				BH::settings::general::s_use_higher_precision=true;
				_MESSAGE("USE_HIGHER_PRECISION set to YES");
			} else {
				BH::settings::general::s_use_higher_precision=false;
				_MESSAGE("USE_HIGHER_PRECISION set to NO");
			};
			break;
		}

		case  use_ep_only : {
			if ( read_answer(in) ){
				settings::general::s_use_ep_only=true;
				_MESSAGE("USE_EVAL_PARAM_ONLY set to YES");
			} else {
				settings::general::s_use_ep_only=false;
				_MESSAGE("USE_EVAL_PARAM_ONLY set to NO");
			};
			break;
		}
		case  use_g3_coupling : {
			if ( read_answer(in) ){
				settings::general::s_use_g3_coupling=true;
				_MESSAGE("USE_G3_COUPLING set to YES");
			} else {
				settings::general::s_use_g3_coupling=false;
				_MESSAGE("USE_G3_COUPLING set to NO");
			};
			break;
		}

		case  use_only_BG_trees : {
			if ( read_answer(in) ){
				settings::general::s_use_only_BG_trees=true;
				_MESSAGE("USE_ONLY_BG_TREES set to YES");
			} else {
				settings::general::s_use_only_BG_trees=false;
				_MESSAGE("USE_ONLY_BG_TREES set to NO");
			};
			break;
		}

		case  catch_large_ME2 : {
			double nbr;
			read_double_answer(in,nbr);
			BH::settings::BH_interface_settings::s_catch_large_ME2=nbr;
			_MESSAGE2("CATCH_LARGE_ME2 set to: ",nbr); break;
		};	break;

        case  number_of_warmup_points : {
			int nbr;
			read_int_answer(in,nbr);
			BH::settings::BH_interface_settings::s_number_of_warmup_points=nbr;
			_MESSAGE2("NUMBER_OF_WARMUP_POINTS set to: ",nbr); break;
		};	break;

		case  w5jet_number_of_threads : {
			int nbr;
			read_int_answer(in,nbr);
			BH::settings::BH_interface_settings::s_w5jet_number_of_threads=nbr;
			_MESSAGE2("W5JET_NUMBER_OF_THREADS set to: ",nbr); break;
		}	break;
		case  use_automated_assembly: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_use_automated_assembly=true;
				_MESSAGE("USE_AUTOMATED_ASSEMBLY set to YES");
			} else {
			settings::BH_interface_settings::s_use_automated_assembly=false;
				_MESSAGE("USE_AUTOMATED_ASSEMBLY set to NO");
			};
		}	break;
        case  final_state_gluon_symmetrization: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_final_state_gluon_symmetrization=true;
				_MESSAGE("FINAL_STATE_GLUON_SYMMETRIZATION set to YES");
			} else {
			settings::BH_interface_settings::s_final_state_gluon_symmetrization=false;
				_MESSAGE("FINAL_STATE_GLUON_SYMMETRIZATION set to NO");
			};
		}	break; 	
        case  born_final_state_gluon_symmetrization: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_born_final_state_gluon_symmetrization=true;
				_MESSAGE("BORN_FINAL_STATE_GLUON_SYMMETRIZATION set to YES");
			} else {
			settings::BH_interface_settings::s_born_final_state_gluon_symmetrization=false;
				_MESSAGE("BORN_FINAL_STATE_GLUON_SYMMETRIZATION set to NO");
			};
		}	break; 	
        case  use_grassman_trees: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_use_grassman_trees=true;
				_MESSAGE("USE_GRASSMAN_TREES set to YES");
			} else {
			settings::BH_interface_settings::s_use_grassman_trees=false;
				_MESSAGE("USE_GRASSMAN_TREES set to NO");
			};
		}	break; 	
        case  use_symmetrized_assembly_files: {
			if ( read_answer(in) ){
				settings::BH_interface_settings::s_use_symmetrized_assembly_files=true;
				_MESSAGE("USE_SYMMETRIZED_ASSEMBLY_FILES set to YES");
			} else {
				settings::BH_interface_settings::s_use_symmetrized_assembly_files=false;
				_MESSAGE("USE_SYMMETRIZED_ASSEMBLY_FILES set to NO");
			};
		}	break;
		case  generate_assembly_files: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_generate_assembly_files=true;
				_MESSAGE("GENERATE_ASSEMBLY_FILES set to YES");
			} else {
			settings::BH_interface_settings::s_generate_assembly_files=false;
				_MESSAGE("GENERATE_ASSEMBLY_FILES set to NO");
			};
		} break;
		case  use_W_polarization_A567_assembly: {
			if ( read_answer(in) ){
			settings::BH_interface_settings::s_use_W_polarization_A567_assembly=true;
				_MESSAGE("USE_W_POLARIZATION_A567_ASSEMBLY set to YES (presently ONLY implemented for W-bosons from 2q1g2l assembly.)");
			} else {
			settings::BH_interface_settings::s_use_W_polarization_A567_assembly=false;
				_MESSAGE("USE_W_POLARIZATION_A567_ASSEMBLY set to NO");
			};
		} break;
		case  W_polarization_A : {
			int nbr;
			read_int_answer(in,nbr);
			settings::BH_interface_settings::s_W_polarization_A=nbr;
			_MESSAGE2("W_POLARIZATION_A set to: ",nbr); break;
		};	break;
		case  same_helicity_projection : {
			int nbr;
			read_int_answer(in,nbr);
			settings::BH_interface_settings::s_same_helicity_projection=nbr;
			_MESSAGE2("SAME_HELICITY_PROJECTION set to: ",nbr); break;
		};	break;
		case  unknown_setting : {
			_MESSAGE("Unknown setting");
			return false;
		}	break;
	}
}
}
//if I have not returned false so far, then I should return true
return true;
}

void read_settings_from_file(const string& path,bool printWarning) {
	std::ifstream file;
	file.open(path.c_str());
	if (file.fail()){
		if (printWarning){
			_WARNING3("Could not open ",path,": done nothing. ");
		}
	} else {
		_MESSAGE3("#-#-#-#-# Reading settings from file ",path," #-#-#-#-#");
		read_from_stream(file);
		_MESSAGE("#-#-#-#-# Done #-#-#-#-#");
	}
}

void use_setting(const string& option) {

	std::istringstream opt(option);
	if (opt.fail()){
		_WARNING3("Could not understand ",option,": done nothing. ");
	} else {
		//_MESSAGE("#-#-#-#-# applying setting  #-#-#-#-#");
		read_from_stream(opt);
		//_MESSAGE("#-#-#-#-# Done #-#-#-#-#");
	}
}


}
}

