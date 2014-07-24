/*
 * calculator.cpp
 *
 *  Created on: Aug 1, 2008
 *      Author: daniel
 */

#include "config.h"

#include "scheme.h"
#include "info.h"
#include <typeinfo>
#include "print_cutD.h"
#include "settings.h"
#include "settings_reader.h"
#include "Interface/BH_interface.h"
#include "Interface/BH_Ampl.h"
#include "OneLoopHelAmpl.h"
#include "IR_checked.h"
#include "cut_part_normal.h"
#include "ratext/ratext_part_normal.h"
#include "polylog.h" //for pi
#include "BH_debug.h"

#include <stdio.h>
#include <iostream>
#if USE_READLINE
#include <dlfcn.h>
#include <readline/readline.h>
#include <readline/history.h>
#endif
#include <string.h>
#include <cassert>

#define GREEN "\e[32m"
#define PURPLE "\e[35m"
#define NOCOLOR "\e[0m"


#define _BHC_MESSAGE(X) std::cout << PURPLE <<  (X) << NOCOLOR << std::endl
#define _BHC_MESSAGE2(X,Y) std::cout << PURPLE << (X)  << (Y) << NOCOLOR << std::endl
#define _BHC_MESSAGE3(X,Y,Z) std::cout << PURPLE << (X) << (Y) <<  (Z) << NOCOLOR << std::endl
#define _BHC_MESSAGE4(W,X,Y,Z) std::cout << PURPLE << (W) << (X) << (Y) <<  (Z) << NOCOLOR << std::endl
#define _BHC_MESSAGE5(V,W,X,Y,Z) std::cout << PURPLE << (V) << (W) << (X) << (Y) <<  (Z) << NOCOLOR << std::endl




using namespace std;

#if USE_READLINE
/*
    dupstr: duplicate a string
*/
char *dupstr(char *s)
{
  char *Result = NULL;
  size_t slen = 0;

  /* sanity check */
  assert(NULL != s);

  /* get string length */
  slen = strlen(s);

  /* allocate enough storage */
  Result = (char *)(malloc(slen + 1));

  /* populate string */
  if(NULL != Result)
  {
    memcpy(Result, s, slen);
    *(Result + slen) = '\0';
  }

  return Result;
}




char* command_names[]={
	"ampl",
	"quit",
	"set",
    "mom_conf",
    "tree",
    "help",
    "cut",
    "print_graph",
    "print_rat_graph",
    "vector",
    "mass",
    "vector_sum",
    "helcode",
    "bhi",
    "repeat",
    "endrepeat",
    "read",
    "use",
    "ordering_constraint",
    (char *)NULL
};

char* option_names[]={
    "use_analytic",
    "info",
    "scheme",
    "precision",
    "mu",
    "output_precision",
    "IR_checked",
    "mom_conf",
    "index",
//
    "leading_color",
    "sub_leading_color",
	"nf",
	"nf_top",
    "LT",
    "RT",
    "LLT",
    "RLT",
    "LRT",
    "RRT",
    "LLLT",
    "RLLT",
    "LRLT",
    "RRLT",
    "LLRT",
    "RLRT",
    "LRRT",
    "RRRT",
    "nfLT",
    "nfRT",
    "nfLLT",
    "nfRLT",
    "nfLRT",
    "nfRRT",
    "nfLLLT",
    "nfRLLT",
    "nfLRLT",
    "nfRRLT",
    "nfLLRT",
    "nfRLRT",
    "nfLRRT",
    "nfRRRT",
    "glue",
    "AX",
    "VECT",
    "AXSL",
    "slc_q",
    "slc_G",
    (char *)NULL
};



/* Generator function for command completion.  STATE lets us know whether
   to start from scratch; without any state (i.e. STATE == 0), then we
   start at the top of the list. */
char * command_generator (const char *text, int state)
{
  static int list_index, len;
  char *name;

  /* If this is a new word to complete, initialize now.  This includes
     saving the length of TEXT for efficiency, and initializing the index
     variable to 0. */
  if (!state)
    {
      list_index = 0;
      len = strlen (text);
    }

  /* Return the next name which partially matches from the command list. */
  while (name = command_names[list_index])
    {
      list_index++;
      if (strncmp (name, text, len) == 0)
        return dupstr(name);
    }

  /* If no names matched, then return NULL. */
  return ((char *)NULL);
}

char * option_generator (const char *text, int state)
{
  static int list_index, len;
  char *name;

  /* If this is a new word to complete, initialize now.  This includes
     saving the length of TEXT for efficiency, and initializing the index
     variable to 0. */
  if (!state)
    {
      list_index = 0;
      len = strlen (text);
    }

  /* Return the next name which partially matches from the command list. */
  while (name = option_names[list_index])
    {
      list_index++;
      if (strncmp (name, text, len) == 0)
        return dupstr(name);
    }

  /* If no names matched, then return NULL. */
  return ((char *)NULL);
}

//rl_completion_matches;



char ** fileman_completion (const char *text, int start, int end)
{
  char **matches;

  matches = (char **)NULL;

  /* If this word is at the start of the line, then it is a command
     to complete.  Otherwise it is the name of a file in the current
     directory. */
  if (start == 0){
    matches = rl_completion_matches (text, command_generator);
	} else {
    matches = rl_completion_matches (text, option_generator);
        }
  return (matches);
}


typedef struct {
  char *name;			/* User printable name of the function. */
  rl_icpfunc_t *func;		/* Function to call to do the job. */
  char *doc;			/* Documentation for this function.  */
} COMMAND;


void initialize_readline ()
{
  /* Allow conditional parsing of the ~/.inputrc file. */
	//  rl_readline_name = "FileMan";


}


#endif

#if __SUNPRO_CC
	#include "typeinfo.h"
#endif

namespace BH {
long helcode(const process& pro);
}
using namespace std;
using namespace BH;

enum bhc_status { BHC_CONTINUE, BHC_QUIT };

template <class T> struct do_delete : public std::unary_function<T*,void> {
	void operator()(T* ptr) { delete ptr;}
};



void read_process(vector<particle_ID>& pps,color_structure& cs,std::istream& is){
	string new_p;
	bool is_particle=true;
	bool has_parenthesis=false;

	while ( is_particle && 	is >> new_p ) {
	particle_ID np;
	if (new_p == "(")  {has_parenthesis=true;}
	else if (new_p == ")")  {is_particle=false;}
	else if (new_p == "qbp")  {np=qbp;	pps.push_back(np);}
	else if (new_p == "qbm")  {np=qbm;	pps.push_back(np);}
	else if (new_p == "qp")  {np=qp;	pps.push_back(np);}
	else if (new_p == "qm")  {np=qm;	pps.push_back(np);}
	else if (new_p == "qb2p")  {np=qb2p;	pps.push_back(np);}
	else if (new_p == "qb2m")  {np=qb2m;	pps.push_back(np);}
	else if (new_p == "q2p")  {np=q2p;	pps.push_back(np);}
	else if (new_p == "q2m")  {np=q2m;	pps.push_back(np);}
	else if (new_p == "qb3p")  {np=qb3p;	pps.push_back(np);}
	else if (new_p == "qb3m")  {np=qb3m;	pps.push_back(np);}
	else if (new_p == "q3p")  {np=q3p;	pps.push_back(np);}
	else if (new_p == "q3m")  {np=q3m;	pps.push_back(np);}
	else if (new_p == "qb21p")  {np=particle_ID(quark,+1,21,true);	pps.push_back(np);}
	else if (new_p == "qb21m")  {np=particle_ID(quark,-1,21,true);	pps.push_back(np);}
	else if (new_p == "q21p")  {np=particle_ID(quark,+1,21,false);	pps.push_back(np);}
	else if (new_p == "q21m")  {np=particle_ID(quark,-1,21,false);	pps.push_back(np);}
	else if (new_p == "qb11p")  {np=particle_ID(quark,+1,11,true);	pps.push_back(np);}
	else if (new_p == "qb11m")  {np=particle_ID(quark,-1,11,true);	pps.push_back(np);}
	else if (new_p == "q11p")  {np=particle_ID(quark,+1,11,false);	pps.push_back(np);}
	else if (new_p == "q11m")  {np=particle_ID(quark,-1,11,false);	pps.push_back(np);}
	else if (new_p == "qb12p")  {np=particle_ID(quark,+1,12,true);	pps.push_back(np);}
	else if (new_p == "qb12m")  {np=particle_ID(quark,-1,12,true);	pps.push_back(np);}
	else if (new_p == "q12p")  {np=particle_ID(quark,+1,12,false);	pps.push_back(np);}
	else if (new_p == "q12m")  {np=particle_ID(quark,-1,12,false);	pps.push_back(np);}
	else if (new_p == "qb22p")  {np=particle_ID(quark,+1,22,true);	pps.push_back(np);}
	else if (new_p == "qb22m")  {np=particle_ID(quark,-1,22,true);	pps.push_back(np);}
	else if (new_p == "q22p")  {np=particle_ID(quark,+1,22,false);	pps.push_back(np);}
	else if (new_p == "q22m")  {np=particle_ID(quark,-1,22,false);	pps.push_back(np);}
	else if (new_p == "p")  {np=p;	pps.push_back(np);}
	else if (new_p == "m")  {np=m;	pps.push_back(np);}
	else if (new_p == "lp")  {np=lp;	pps.push_back(np);}
	else if (new_p == "lm")  {np=lm;	pps.push_back(np);}
	else if (new_p == "lbp")  {np=lbp;	pps.push_back(np);}
	else if (new_p == "lbm")  {np=lbm;	pps.push_back(np);}
	else if (new_p == "yp")  {np=yp;	pps.push_back(np);}
	else if (new_p == "ym")  {np=ym;	pps.push_back(np);}
	else if (new_p == "y2p")  {np=y2p;	pps.push_back(np);}
	else if (new_p == "y2m")  {np=y2m;	pps.push_back(np);}
	else if (new_p == "Gbp")  {np=Gbp;	pps.push_back(np);}
	else if (new_p == "Gbm")  {np=Gbm;	pps.push_back(np);}
	else if (new_p == "Gp")  {np=Gp;	pps.push_back(np);}
	else if (new_p == "Gm")  {np=Gm;	pps.push_back(np);}
	else if (new_p == "Gb2p")  {np=particle_ID(gluino,+1,2,true);	pps.push_back(np);}
	else if (new_p == "Gb2m")  {np=particle_ID(gluino,-1,2,true);	pps.push_back(np);}
	else if (new_p == "G2p")  {np=particle_ID(gluino,+1,2,false);	pps.push_back(np);}
	else if (new_p == "G2m")  {np=particle_ID(gluino,-1,2,false);	pps.push_back(np);}
	else if (new_p == "Gb3p")  {np=particle_ID(gluino,+1,3,true);	pps.push_back(np);}
	else if (new_p == "Gb3m")  {np=particle_ID(gluino,-1,3,true);	pps.push_back(np);}
	else if (new_p == "G3p")  {np=particle_ID(gluino,+1,3,false);	pps.push_back(np);}
	else if (new_p == "G3m")  {np=particle_ID(gluino,-1,3,false);	pps.push_back(np);}
    else if (new_p == "nf")  {cs=nf;		is_particle=false;	}
	else if (new_p == "nf_top")  {cs=nf_top;		is_particle=false;	}
    else if (new_p == "LT")  {cs=LT;		is_particle=false;	}
	else if (new_p == "RT")  {cs=RT;		is_particle=false;	}
	else if (new_p == "LLT")  {cs=LLT;		is_particle=false;	}
	else if (new_p == "LRT")  {cs=LRT;		is_particle=false;	}
	else if (new_p == "RLT")  {cs=RLT;		is_particle=false;	}
	else if (new_p == "RRT")  {cs=RRT;		is_particle=false;	}
	else if (new_p == "LLLT")  {cs=LLLT;		is_particle=false;	}
	else if (new_p == "LRLT")  {cs=LRLT;		is_particle=false;	}
	else if (new_p == "RLLT")  {cs=RLLT;		is_particle=false;	}
	else if (new_p == "RRLT")  {cs=RRLT;		is_particle=false;	}
	else if (new_p == "LLRT")  {cs=LLRT;		is_particle=false;	}
	else if (new_p == "LRRT")  {cs=LRRT;		is_particle=false;	}
	else if (new_p == "RLRT")  {cs=RLRT;		is_particle=false;	}
	else if (new_p == "RRRT")  {cs=RRRT;		is_particle=false;	}
    else if (new_p == "nfLT")  {cs=nfLT;		is_particle=false;	}
	else if (new_p == "nfRT")  {cs=nfRT;		is_particle=false;	}
	else if (new_p == "nfLLT")  {cs=nfLLT;		is_particle=false;	}
	else if (new_p == "nfLRT")  {cs=nfLRT;		is_particle=false;	}
	else if (new_p == "nfRLT")  {cs=nfRLT;		is_particle=false;	}
	else if (new_p == "nfRRT")  {cs=nfRRT;		is_particle=false;	}
	else if (new_p == "nfLLLT")  {cs=nfLLLT;		is_particle=false;	}
	else if (new_p == "nfLRLT")  {cs=nfLRLT;		is_particle=false;	}
	else if (new_p == "nfRLLT")  {cs=nfRLLT;		is_particle=false;	}
	else if (new_p == "nfRRLT")  {cs=nfRRLT;		is_particle=false;	}
	else if (new_p == "nfLLRT")  {cs=nfLLRT;		is_particle=false;	}
	else if (new_p == "nfLRRT")  {cs=nfLRRT;		is_particle=false;	}
	else if (new_p == "nfRLRT")  {cs=nfRLRT;		is_particle=false;	}
	else if (new_p == "nfRRRT")  {cs=nfRRRT;		is_particle=false;	}
    else if (new_p == "glue")  {cs=glue;	is_particle=false;	}
	else if (new_p == "leading_color")  {cs=leading_color;	is_particle=false;	}
	else if (new_p == "sub_leading_color")  {cs=sub_leading_color;	is_particle=false;	}
	else if (new_p == "AX")  {cs=AX;	is_particle=false;	}
	else if (new_p == "AXSL")  {cs=AXSL;	is_particle=false;	}
	else if (new_p == "VECT")  {cs=VECT;	is_particle=false;	}
	else if (new_p == "slc_q")  {cs=slc_q;  is_particle=false;      }
	else if (new_p == "slc_G")  {cs=slc_G;  is_particle=false;      }

	else  {_MESSAGE2("]> !! Unknown particle or color structure ",new_p);}
}
}



void read_index_vector(vector<int>* ind,std::istream& is){
	string new_p;
	bool is_ok=true;
	bool missed_comma=false;
	bool has_parenthesis=false;
	is >> new_p;
	if (new_p != "{")  { throw BHerror("index vector has to start with a {");}
int next_index;
	while ( is_ok && is >> next_index && is >> new_p ) {
	if ( new_p != "," && new_p != "}"){
		missed_comma=true;
	}
	ind->push_back(next_index);
	if (new_p == "}")  {is_ok=false;}
	}
	if ( missed_comma){
		cout << " Warning: set index expects input as { N1 , N2 , ... , NS }" << endl;
	}
}


void read_bhi_process(vector<int>& pdg,istream& is){
	// gluon 21
	// photon 22
	// d 1
	// u 2
	//s c b t
	//e- 11
	//ne 12
	//mu- 13
	//nmu 14
	//tau- 15
	//ntau- 16
	string new_p;
	bool is_particle=true;
	bool has_parenthesis=false;

	while ( is_particle && 	is >> new_p ) {
	particle_ID np;
	if (new_p == "(")  {has_parenthesis=true;}
	else if (new_p == ")")  {is_particle=false;}
	else if (new_p == "g")  {pdg.push_back(21);}
	else if (new_p == "y")  {pdg.push_back(22);}
	else if (new_p == "d")  {pdg.push_back(1);}
	else if (new_p == "u")  {pdg.push_back(2);}
	else if (new_p == "db")  {pdg.push_back(-1);}
	else if (new_p == "ub")  {pdg.push_back(-2);}
	else if (new_p == "e-")  {pdg.push_back(11);}
	else if (new_p == "e+")  {pdg.push_back(-11);}
	else if (new_p == "ne")  {pdg.push_back(12);}
	else if (new_p == "mu-")  {pdg.push_back(13);}
	else if (new_p == "nmu")  {pdg.push_back(14);}
	else if (new_p == "tau-")  {pdg.push_back(15);}
	else if (new_p == "ntau")  {pdg.push_back(16);}

	else  {_MESSAGE2("]> !! Unknown particle: ",new_p);}
}
}


enum prec_setting { NP, HP, VHP,multi,not_known};
enum command {
	ampl_cmd,
	set_cmd,
	quit_cmd,
	mom_conf_cmd,
	not_known_cmd,
	tree_cmd,
	help_cmd,
	cut_cmd,
	print_graph_cmd,
	print_rat_graph_cmd,
	mass_cmd,
	vector_cmd,
	vector_sum_cmd,
	helcode_cmd,
	bhi_cmd,
	repeat_cmd,
	read_cmd,
	use_cmd,
	ord_constr_cmd
};

map<string,command> Init_command_list(	){
	std::map<string,command> cmd_map;
	cmd_map.insert(std::pair<string,command>("ampl",ampl_cmd));
	cmd_map.insert(std::pair<string,command>("quit",quit_cmd));
	cmd_map.insert(std::pair<string,command>("set",set_cmd));
	cmd_map.insert(std::pair<string,command>("mom_conf",mom_conf_cmd));
	cmd_map.insert(std::pair<string,command>("tree",tree_cmd));
	cmd_map.insert(std::pair<string,command>("help",help_cmd));
	cmd_map.insert(std::pair<string,command>("cut",cut_cmd));
	cmd_map.insert(std::pair<string,command>("print_graph",print_graph_cmd));
	cmd_map.insert(std::pair<string,command>("print_rat_graph",print_rat_graph_cmd));
	cmd_map.insert(std::pair<string,command>("vector",vector_cmd));
	cmd_map.insert(std::pair<string,command>("mass",mass_cmd));
	cmd_map.insert(std::pair<string,command>("vector_sum",vector_sum_cmd));
	cmd_map.insert(std::pair<string,command>("helcode",helcode_cmd));
	cmd_map.insert(std::pair<string,command>("bhi",bhi_cmd));
	cmd_map.insert(std::pair<string,command>("repeat",repeat_cmd));
	cmd_map.insert(std::pair<string,command>("endrepeat",repeat_cmd));
	cmd_map.insert(std::pair<string,command>("read",read_cmd));
	cmd_map.insert(std::pair<string,command>("use",use_cmd));
	cmd_map.insert(std::pair<string,command>("ordering_constraint",ord_constr_cmd));
return cmd_map;
}


command get_command(istream& is, bool prompt=true){
	static std::map<string,command> cmd_list=Init_command_list();

	string the_command;
	if (prompt) { cout << "]> ";};
	is >> the_command;

	map<string,command>::iterator pos;
	pos= cmd_list.find(the_command);

	if (pos!=cmd_list.end()){
		return pos->second;
	}
	else return not_known_cmd;

}

struct session_info {
	size_t old_n;
	vector<vector<int>*> all_new_index_vectors;
	vector<int>* ind;
	vector<int>* ind_ptr4;
	vector<int>* ind_ptr5;
	vector<int>* ind_ptr6;
	vector<int>* ind_ptr7;
	vector<int>* ind_ptr8;
	vector<int>* ind_ptr9;
	vector<int>* ind_ptr10;
	mom_conf* mc;
	mom_conf_HP* mc_HP;
	mom_conf_VHP* mc_VHP;
	bool do_delete_mc;
	bool print_info;
	bool ir_checked;
	multi_precision_constant mu;
	prec_setting precision;
	scheme scheme_setting;
	BH_interface BHI;

	session_info() :old_n(0), mc(0), mc_HP(0), mc_VHP(0), do_delete_mc(false), print_info(false), ir_checked(false), mu(1), precision(NP),scheme_setting(FDH) {
			ind_ptr4=0;
			ind_ptr5=0;
			ind_ptr6=0;
			ind_ptr7=0;
			ind_ptr8=0;
			ind_ptr9=0;
			ind_ptr10=0;
	};
};

std::vector <vector<int>* > local_ind(11);


std::vector <momentum_configuration<R>* > local_mc;
std::vector <momentum_configuration<RHP>* > local_mc_HP;
std::vector <momentum_configuration<RVHP>* > local_mc_VHP;


void set_mc(std::istream& is,session_info& si){
	string value;
	is >> value;
	si.mc=0;
	if ( value == "egz"){ si.mc=local_mc[0];si.mc_HP=local_mc_HP[0];si.mc_VHP=local_mc_VHP[0];si.old_n=6;};
	if ( value == "gkm6"){ si.mc=local_mc[1];si.mc_HP=local_mc_HP[1];si.mc_VHP=local_mc_VHP[1];si.old_n=6;};
	if ( value == "bbdfk7pt"){ si.mc=local_mc[2];si.mc_HP=local_mc_HP[2];si.mc_VHP=local_mc_VHP[2];si.old_n=7;};
	if ( value == "mc4"){ si.mc=local_mc[3];si.old_n=4;};
	if ( value == "mc5"){ si.mc=local_mc[4];si.old_n=5;};
	if ( value == "mc6"){ si.mc=local_mc[5];si.old_n=6;};
	if ( value == "mc7"){ si.mc=local_mc[6];si.old_n=7;};
	if ( value == "mc8"){ si.mc=local_mc[7];si.old_n=8;};
	if ( value == "mc9"){ si.mc=local_mc[8];si.old_n=9;};
	if ( value == "gz7pt"){si.mc=local_mc[9];si.old_n=7;};
	if ( value == "gz8pt"){ si.mc=local_mc[10];si.old_n=8;};
if ( value == "bbdfk8pt"){ si.mc=local_mc[11];si.mc_HP=local_mc_HP[3];si.mc_VHP=local_mc_VHP[3];si.old_n=8;}

	if (si.mc == 0){ // mc not set yet
		int number;
		is >> number;
		switch (si.precision){
		case multi: {
			multi_precision_reader* mpr=new multi_precision_reader(value.c_str(),number);
			mpr->next();
			si.mc=mpr;	break;
		}
		case NP: {
			mc_reader* mpr=new mc_reader(value.c_str(),number);
			mpr->next();
			si.mc=mpr;	break;
		}
		case HP: {
			mc_reader_HP* mpr=new mc_reader_HP(value.c_str(),number);
			mpr->next();
			si.mc_HP=mpr;	break;
		}
		case VHP: {
			mc_reader_VHP* mpr=new mc_reader_VHP(value.c_str(),number);
			mpr->next();
			si.mc_VHP=mpr;	break;
		}

		}
		si.old_n=number;
	}

}

void set_index(int n, session_info& si){
	switch (n){
	case 4: if (si.ind_ptr4 == 0) {si.ind=local_ind[4];} else {si.ind=si.ind_ptr4;};   break;
	case 5: if (si.ind_ptr5 == 0) {si.ind=local_ind[5];} else {si.ind=si.ind_ptr5;};   break;
	case 6: if (si.ind_ptr6 == 0) {si.ind=local_ind[6];} else {si.ind=si.ind_ptr6;};   break;
	case 7: if (si.ind_ptr7 == 0) {si.ind=local_ind[7];} else {si.ind=si.ind_ptr7;};   break;
	case 8: if (si.ind_ptr8 == 0) {si.ind=local_ind[8];} else {si.ind=si.ind_ptr8;};   break;
	case 9: if (si.ind_ptr9 == 0) {si.ind=local_ind[9];} else {si.ind=si.ind_ptr9;};   break;
	case 10: if (si.ind_ptr10 == 0) {si.ind=local_ind[10];} else {si.ind=si.ind_ptr10;};   break;
	}
}

bhc_status execute(istream& is,session_info& si){
	while (!is.fail() && !is.eof()){
		command the_command=get_command(is,false);
if (!is.fail() ){
	switch (the_command) {
	case set_cmd : {
		string param;
		is >> param;
		if (param == "mu"){
			double value;
			is >> value;
			si.mu.set(value);
		}
		if (param == "info"){
			string value;
			is >> value;
			if (value == "on") {si.print_info=true; _BHC_MESSAGE("Info set on"); };
			if (value == "off") {si.print_info=false;_BHC_MESSAGE("Info set off");};
		}
		if (param == "use_analytic"){
			string value;
			is >> value;
			if (value == "on") {settings::general::s_use_known_formulae=true;_BHC_MESSAGE("use_analytic set on");};
			if (value == "off") {settings::general::s_use_known_formulae=false;_BHC_MESSAGE("use_analytic set off");};
		}
		if (param == "use_cached_integrals"){
			string value;
			is >> value;
			if (value == "on") {settings::general::s_use_cached_integrals=true;};
			if (value == "off") {settings::general::s_use_cached_integrals=false;};
		}
		if (param == "scheme"){
			string value;
			is >> value;
			if (value == "FHD") {si.scheme_setting=FDH;}
			else if (value == "HV") {si.scheme_setting=HV;}
			else {_MESSAGE3("Unknown scheme: ", value, " using FDH."); si.scheme_setting=FDH; };
		}
		if (param == "precision"){
			string setting;
			is >> setting;
			si.precision= not_known;
			if (setting == "NP") {si.precision = NP;}
			if (setting == "HP") {si.precision = HP;}
			if (setting == "VHP") {si.precision = VHP;}
			if (setting == "multi") {si.precision = multi;}
			if (si.precision == not_known) { throw BHerror("Unknown precision setting."); }
		}
		if (param == "output_precision"){
			int value;
			is >> value;
			cout << setprecision(value);
		}
		if (param == "IR_checked"){
			string value;
			is >> value;
			if (value == "on") {si.ir_checked=true;};
			if (value == "off") {si.ir_checked=false;};
		}
		if (param == "mom_conf"){
			set_mc(is,si);
		}
		if (param == "index"){
			vector<int>* indices= new vector<int> ;
			read_index_vector(indices,is);
			si.all_new_index_vectors.push_back(indices);
			switch (indices->size()){
			case 4: si.ind_ptr4=indices;   break;
			case 5: si.ind_ptr5=indices;   break;
			case 6: si.ind_ptr6=indices;   break;
			case 7: si.ind_ptr7=indices;   break;
			case 8: si.ind_ptr8=indices;   break;
			case 9: si.ind_ptr9=indices;   break;
			case 10: si.ind_ptr10=indices;   break;
			}
		}

		break;
	}

	case ampl_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		process PRO(pps);
		if (si.mc == 0 || PRO.n() != si.old_n) {
			switch (PRO.n()){
			case 4: si.mc=local_mc[3];   break;
			case 5: si.mc=local_mc[4];   break;
			case 6: si.mc=local_mc[5];   break;
			case 7: si.mc=local_mc[6];   break;
			case 8: si.mc=local_mc[7];   break;
			case 9: si.mc=local_mc[8];   break;
			}
		}
		set_index(PRO.n(),si);
		_MESSAGE3("Creating amplitude ",PRO," ...");

		OneLoopAmplitude_base* A;
		switch (si.ir_checked){
			case true :	A= new IR_checked_OLHA(PRO,cs); break;
			case false: A= new OneLoopHelAmpl(PRO,cs);break;
		}

		A->set_mu(si.mu);
		A->set_scheme(si.scheme_setting);

		_MESSAGE("Computing, please be patient ...");

		switch (si.precision){
		case NP: case multi: A->eval(*si.mc,*si.ind); break;
		case HP:  A->eval(*si.mc_HP,*si.ind);break;
		case VHP:  A->eval(*si.mc_VHP,*si.ind);break;
		}
		_MESSAGE("-----------------------------------------------");
		_MESSAGE4("Process: ",PRO," ",cs);
		_MESSAGE("-----------------------------------------------");
		_MESSAGE("Result not normalized");
		_MESSAGE("-----------------------------------------------");
		switch (si.precision){
		case NP: case multi: {
			_MESSAGE2("cut: ",A->get_cut(*si.mc,*si.ind));
			_MESSAGE2("rational: ",A->get_rational(*si.mc,*si.ind));
			_MESSAGE2("amplitude: ",A->get_value(*si.mc,*si.ind));
			break;
		}
		case HP: {
			_MESSAGE2("cut: ",A->get_cut(*si.mc_HP,*si.ind));
			_MESSAGE2("rational: ",A->get_rational(*si.mc_HP,*si.ind));
			_MESSAGE2("amplitude: ",A->get_value(*si.mc_HP,*si.ind));
			break;
		}
		case VHP: {
			_MESSAGE2("cut: ",A->get_cut(*si.mc_VHP,*si.ind));
			_MESSAGE2("rational: ",A->get_rational(*si.mc_VHP,*si.ind));
			_MESSAGE2("amplitude: ",A->get_value(*si.mc_VHP,*si.ind));
			break;
		};
		}

		TreeHelAmpl TREE(PRO);

		if (! TREE.is_zero()){
			_MESSAGE("-----------------------------------------------");
			_MESSAGE("Result normalized");
			_MESSAGE("-----------------------------------------------");
			switch (si.precision){
			case NP: case multi: {
				C tree=TREE.eval(*si.mc,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",A->get_cut(*si.mc,*si.ind)/tree);
				_MESSAGE2("rational: ",A->get_rational(*si.mc,*si.ind)/tree);
				_MESSAGE2("amplitude: ",A->get_value(*si.mc,*si.ind)/tree);
				break;
			}
			case HP:  {
				CHP tree=TREE.eval(*si.mc_HP,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",A->get_cut(*si.mc_HP,*si.ind)/tree);
				_MESSAGE2("rational: ",A->get_rational(*si.mc_HP,*si.ind)/tree);
				_MESSAGE2("amplitude: ",A->get_value(*si.mc_HP,*si.ind)/tree);
				break;
			}
			case VHP:  {
				CVHP tree=TREE.eval(*si.mc_VHP,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",A->get_cut(*si.mc_VHP,*si.ind)/tree);
				_MESSAGE2("rational: ",A->get_rational(*si.mc_VHP,*si.ind)/tree);
				_MESSAGE2("amplitude: ",A->get_value(*si.mc_VHP,*si.ind)/tree);
				break;
			}
			}
			_MESSAGE("-----------------------------------------------");

		}
		if (si.print_info){
			switch (si.ir_checked){
			case true : {
				IR_checked_OLHA* IRCO = dynamic_cast<IR_checked_OLHA*>(A);
				switch (si.precision){
				case NP: case multi: {
					info(IRCO->cut_part(),*si.mc,*si.ind);
					info(IRCO->rational_part(),*si.mc,*si.ind);
					break;
				}
				case HP: {
					info(IRCO->rational_part(),*si.mc_HP,*si.ind);
					break;
				}
				case VHP: {
					info(IRCO->rational_part(),*si.mc_VHP,*si.ind);
					break;
				};
				}
			}
			case false:{
				OneLoopHelAmpl* IRCO = dynamic_cast<OneLoopHelAmpl*>(A);
				switch (si.precision){
				case NP: case multi: {
					info(IRCO->cut_part(),*si.mc,*si.ind);
					info(IRCO->rational_part(),*si.mc,*si.ind);
					break;
				}
				case HP: {
					info(IRCO->cut_part(),*si.mc,*si.ind);
					info(IRCO->rational_part(),*si.mc_HP,*si.ind);
					break;
				}
				case VHP: {
					info(IRCO->rational_part(),*si.mc_VHP,*si.ind);
					break;
				};
				}
			}
			}
		}
		delete A;

		break;
	}
	case cut_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		process PRO(pps);
		if (si.mc == 0 || PRO.n() != si.old_n) {
			switch (PRO.n()){
			case 4: si.mc=local_mc[3];   break;
			case 5: si.mc=local_mc[4];   break;
			case 6: si.mc=local_mc[5];   break;
			case 7: si.mc=local_mc[6];   break;
			case 8: si.mc=local_mc[7];   break;
			case 9: si.mc=local_mc[8];   break;
			}
		}
		set_index(PRO.n(),si);
		_MESSAGE3("Creating amplitude ",PRO," ...");

		Cut_Part_base* CP;
		CP = cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(PRO)->new_cut_part(PRO,cs);


		CP->set_mu(si.mu);
//		A->set_scheme(scheme_setting);

		_MESSAGE("Computing, please be patient ...");

		switch (si.precision){
		case NP: case multi: CP->eval(*si.mc,*si.ind); break;
		case HP:  CP->eval(*si.mc_HP,*si.ind);break;
		case VHP:  CP->eval(*si.mc_VHP,*si.ind);break;
		}
		_MESSAGE("-----------------------------------------------");
		_MESSAGE4("Process: ",PRO," ",cs);
		_MESSAGE("-----------------------------------------------");
		_MESSAGE("Result not normalized");
		_MESSAGE("-----------------------------------------------");
		switch (si.precision){
		case NP: case multi: {
			_MESSAGE2("cut: ",CP->get_value(*si.mc,*si.ind));
			break;
		}
		case HP: {
			_MESSAGE2("cut: ",CP->get_value(*si.mc_HP,*si.ind));
			break;
		}
		case VHP: {
			_MESSAGE2("cut: ",CP->get_value(*si.mc_VHP,*si.ind));
			break;
		};
		}

		TreeHelAmpl TREE(PRO);

		if (! TREE.is_zero()){
			_MESSAGE("-----------------------------------------------");
			_MESSAGE("Result normalized");
			_MESSAGE("-----------------------------------------------");
			switch (si.precision){
			case NP: case multi: {
				C tree=TREE.eval(*si.mc,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",CP->get_value(*si.mc,*si.ind)/tree);
				break;
			}
			case HP:  {
				CHP tree=TREE.eval(*si.mc_HP,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",CP->get_value(*si.mc_HP,*si.ind)/tree);
				break;
			}
			case VHP:  {
				CVHP tree=TREE.eval(*si.mc_VHP,*si.ind);
				_MESSAGE2("tree: ",tree);
				_MESSAGE2("cut: ",CP->get_value(*si.mc_VHP,*si.ind)/tree);
				break;
			}
			}
			_MESSAGE("-----------------------------------------------");

		}
		if (si.print_info){
				switch (si.precision){
				case NP: case multi: {
					info(CP,*si.mc,*si.ind);
					break;
				}
				case HP: {
					info(CP,*si.mc_HP,*si.ind);
					break;
				}
				case VHP: {
					info(CP,*si.mc_VHP,*si.ind);
					break;
				};
				}


			}

		delete CP;

		break;
	}

	case tree_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		process PRO(pps);
		if (si.mc == 0 || PRO.n() != si.old_n) {
			switch (PRO.n()){
			case 4: si.mc=local_mc[3];   break;
			case 5: si.mc=local_mc[4];   break;
			case 6: si.mc=local_mc[5];   break;
			case 7: si.mc=local_mc[6];   break;
			case 8: si.mc=local_mc[7];   break;
			case 9: si.mc=local_mc[8];   break;
			}
		}
		set_index(PRO.n(),si);
		_MESSAGE3("Creating amplitude ",PRO," ...");

		TreeHelAmpl A(PRO);

		_MESSAGE("-----------------------------------------------");
		switch (si.precision){
		case NP: case multi: 	_MESSAGE4("Process: ",PRO," tree: ",A.eval(*si.mc,*si.ind));break;
		case HP:  				_MESSAGE4("Process: ",PRO," tree: ",A.eval(*si.mc_HP,*si.ind));break;
		case VHP:  				_MESSAGE4("Process: ",PRO," tree: ",A.eval(*si.mc_VHP,*si.ind));break;
		}
		_MESSAGE("-----------------------------------------------");
		if (si.print_info){
			info(A.pointee(),*si.mc,*si.ind);
		}
		break;
	}
	case use_cmd:{
		string param;
		is >> param;
		if (param != "setting"){
			_WARNING("always use \"use setting\" for the \"use\" command!"); break;
		}
		getline(is,param) ;
		settings::use_setting(param);
		break;
	}
	case helcode_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		process PRO(pps);
		_MESSAGE("-----------------------------------------------");
		_MESSAGE4("helcode for process ",PRO," : ",helcode(PRO) );
		_MESSAGE("-----------------------------------------------");


		break;
	}
	case ord_constr_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		process PRO(pps);
		ordering_constraint oc = get_ordering_constraint(PRO);
		_MESSAGE("-----------------------------------------------");
		_MESSAGE3("ordering_constraint for process ",PRO," : \n" );
		_MESSAGE("strong:");
		for (int j=0;j<oc.strong.size();j++){
			copy(oc.strong[j].begin(),oc.strong[j].end(),ostream_iterator<int>(cout," ")); cout << "\n";
		}
		_MESSAGE("weak:");
		for (int j=0;j<oc.weak.size();j++){
			copy(oc.weak[j].begin(),oc.weak[j].end(),ostream_iterator<int>(cout," ")); cout << "\n";
		}
		_MESSAGE("-----------------------------------------------");


		break;
	}

	case vector_cmd: {
		int index;
		is >> index;
		switch (si.precision){
		case NP: case multi:  if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("vector ",index,": ",si.mc->p(index));}; break;
		case HP: if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("vector ",index,": ",si.mc_HP->p(index));}; break;
		case VHP: if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("vector ",index,": ",si.mc_VHP->p(index));}; break;
		}
		break;
	}

	case mass_cmd: {
		int index;
		is >> index;
		switch (si.precision){
		case NP: case multi:  if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("mass^2 of vector ",index,": ",si.mc->m2(index));}; break;
		case HP: if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("mass^2 of vector ",index,": ",si.mc_HP->m2(index));}; break;
		case VHP: if (index > si.mc->n()) { _MESSAGE("Index too large for the current mom_conf");}  else { _MESSAGE4("mass^2 of vector ",index,": ",si.mc_VHP->m2(index));}; break;
		}
		break;
	}
	case vector_sum_cmd: {

		vector<int>* indices= new vector<int> ;
		read_index_vector(indices,is);
		si.all_new_index_vectors.push_back(indices);
		switch (si.precision){
		case NP: case multi:  {

			size_t P=si.mc->Sum(*indices);
			 _MESSAGE2("sum :",si.mc->p(P));
		}; break;
		case HP: {

			size_t P=si.mc_HP->Sum(*indices);
			 _MESSAGE2("sum :",si.mc_HP->p(P));
		}; break;
		case VHP: {

			size_t P=si.mc_VHP->Sum(*indices);
			 _MESSAGE2("sum :",si.mc_VHP->p(P));
		}; break;
	}; break;
	}
	case print_graph_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		string path;
		is >> path;

		process PRO(pps);

		_MESSAGE3("Creating amplitude ",PRO," ...");
		_MESSAGE( "FeynDiagram by Bill Dimm, bdimm@hotneuron.com\n");
		_MESSAGE( "Version $Revision: 1.27 $   $Date: 2003/08/13 06:23:38");
		Cut_Part_base* CP;
		CP= cut_part_factory<Cut_Part_base>::s_default_cut_part_factory(PRO)->new_cut_part(PRO,cs);

		_MESSAGE("Generating graph, please be patient ...");
		typedef cut::normal_cut_part<cut::Darren_CutD_Factory> CutType;

		CutType* CT = dynamic_cast<CutType*>(CP);

		if ( CT ){
			print_tree_graph(CT,path.c_str());
		} else {
			_MESSAGE4("Wrong type for cut type: ",typeid(*CP).name()," expected: ",typeid(CutType *).name());
		}



		delete CP;

		break;
	}
	case print_rat_graph_cmd: {
		vector<particle_ID> pps;
		color_structure cs;
		read_process(pps,cs,is);
		string path;
		is >> path;

		process PRO(pps);

		_MESSAGE3("Creating rational part ",PRO," ...");
		_MESSAGE( "FeynDiagram by Bill Dimm, bdimm@hotneuron.com\n");
		_MESSAGE( "Version $Revision: 1.27 $   $Date: 2003/08/13 06:23:38");

		ratext::normal_ratext_factory explicit_NRF;

		Rational_base* RB = explicit_NRF.new_rational(PRO,cs);
		ratext::normal_ratext* NR=dynamic_cast<ratext::normal_ratext*>(RB);


		if ( NR ){
			Cut_Part_D_Dims* CPDD=NR->get_rat_structure();
			print_tree_graph(CPDD,path.c_str());
		} else {
			_MESSAGE4("Wrong type for rat type: ",typeid(*RB).name()," expected: ",typeid(ratext::normal_ratext).name());
		}



		delete RB;

		break;
	}
	case mom_conf_cmd: {
		string param;
		is >> param;
		mom_conf_reader_base* mcrb;
		switch(si.precision){
		case multi: case NP : mcrb=dynamic_cast<mom_conf_reader_base*>(si.mc); break;
		case HP : mcrb=dynamic_cast<mom_conf_reader_base*>(si.mc_HP);break;
		case VHP : mcrb=dynamic_cast<mom_conf_reader_base*>(si.mc_VHP);break;
		}
		if ( mcrb != 0 ){
			if (param == "next"){
				if (! mcrb->next()) {_MESSAGE("Could not perform operation NEXT.");};
			}
			if (param == "goto"){
				size_t n;
				is >> n;
				mcrb->go_to(n);
			}
		}
		else {
			_PRINT(typeid(*si.mc).name());
			throw BHerror("Command next or goto applied for a non mc_reader type.");
		}
		break;
	}
	case repeat_cmd:{
		string all_commands,next_line;
		int nbr_of_iteration;
		is >> nbr_of_iteration;
		//ignore anything after that
		getline(is,next_line) ;
		bool found_endrepeat=false;
		while( (!found_endrepeat) && getline(is,next_line) ) {
		  if( next_line.size() >= 3 && next_line[0] == 'e' && next_line[1] == 'n' && next_line[2] == 'd' ) {
			  found_endrepeat=true;
		  } else {
			  all_commands += "\n";
			  all_commands+=next_line;
		  }
		}
		_MESSAGE("repeating:");
		cout << all_commands << endl;
		cout << nbr_of_iteration << " times" << endl;

		for (int i=1;i<=nbr_of_iteration;i++){
			istringstream ss(all_commands);
			execute(ss,si);
		}


	} break;

	case read_cmd : {
		string filename;
		is >> filename;
		ifstream ifile;
		ifile.open(filename.c_str());
		if (ifile){
			_MESSAGE2("Reading from file ",filename);
			execute(ifile,si);
		} else {
			_WARNING2("Could not open file ",filename);
		}

	} break;

	case  bhi_cmd: {
		vector<int> pdgs;
		read_bhi_process(pdgs,is);
		BH_Ampl* ampl=si.BHI.new_ampl(pdgs);
		if (si.mc == 0 || pdgs.size() != si.old_n) {
			switch (pdgs.size()){
			case 4: si.mc=local_mc[3];   break;
			case 5: si.mc=local_mc[4];   break;
			case 6: si.mc=local_mc[5];   break;
			case 7: si.mc=local_mc[6];   break;
			case 8: si.mc=local_mc[7];   break;
			case 9: si.mc=local_mc[8];   break;
			}
		}

		vector<vector <double> > momenta;
		for (int i=1;i<=pdgs.size();i++){
			vector<double> p;
			if (i<=2){
				p.push_back(-si.mc->mom(i).E().real());
				p.push_back(-si.mc->mom(i).X().real());
				p.push_back(-si.mc->mom(i).Y().real());
				p.push_back(-si.mc->mom(i).Z().real());
			} else {
				p.push_back(si.mc->mom(i).E().real());
				p.push_back(si.mc->mom(i).X().real());
				p.push_back(si.mc->mom(i).Y().real());
				p.push_back(si.mc->mom(i).Z().real());
			}
			momenta.push_back(p);
		}

		BHinput input(momenta,si.mu);

		si.BHI(input);
		ampl->get_finite();
		cout << ampl->get_double_pole() << " e^-2  + " <<  ampl->get_single_pole() << " e^-1 + " << ampl->get_finite() << endl;

	} break;

	case help_cmd: {
		_MESSAGE(" =============================================================================== \n");
		_MESSAGE(" This is the BlackHat calculator.  \n");
		_MESSAGE(" list of commands:  \n");
		_MESSAGE(" set [parameter] [value] \n");
		_MESSAGE("     [parameter] || [value] can be \n"
				 "    mu  NUM    : \n"
				 "    mom_conf FILE_NAME NUM \t : opens FILE_NAME with a mom_conf reader for NUM momenta \n"
				 "    precision [NP|HP|VHP|multi]               \n"
				 "    scheme  [HV|FDH] \t\t : sets the scheme \n"
				 "    output_precision NUM \t : sets the precision \n"
				 "    IR_checked [on|off] \t : sets IR checks on or off \n"
				 "    info [on|off] \t\t : sets information on or off \n"
				 "    use_analytic [on|off] \t : sets whether to use analytic formulae for cut_parts \n"
				 "\n"
				 "  mom_conf next \t\t : makes the mom_conf reader go to the next point \n"
				 "\n"
				 "  ampl PROCESS COLOR_STRUCTURE \t : computes the one-loop amplitude for PROCESS \n"
				 "\n"
				 "  tree ( PROCESS ) \t\t : computes the tree amplitude for PROCESS \n"
				 "\n"
				 "  vector INDEX \t\t\t : prints vector INDEX of the current mom_conf \n"
				 "  mass INDEX \t\t\t :  prints the mass of vector INDEX of the current mom_conf  \n"
				 "  tree INDEX_VECTOR \t\t : computes the sum of the vectors INDEX_VECTOR \n"
				 "\n"
				 "  help \t\t\t : displays this help\n"
				 "\n"
				 "  quit \t\t\t : terminates the program\n"
				);
		_MESSAGE(" =============================================================================== \n");
		break;
	}
	case quit_cmd: {
		if (si.do_delete_mc) {delete si.mc;}
		return BHC_QUIT;
	}
	case not_known_cmd: {
		_MESSAGE("Command unknown");
		break;
	}
	}
	}
}
	return BHC_CONTINUE;

}


int main() {

	   unsigned int old_cw;

	fpu_fix_start(&old_cw);

	#include "std_momenta.hpp"

	local_mc_HP.push_back(&egzHP);
	local_mc_HP.push_back(&gkm6_HP);
	local_mc_HP.push_back(&bbdfk7pt_HP);
	local_mc_VHP.push_back(&egzVHP);
	local_mc_VHP.push_back(&gkm6_VHP);
	local_mc_VHP.push_back(&bbdfk7pt_VHP);


	local_mc.push_back(&egz);
	local_mc.push_back(&gkm6);
	local_mc.push_back(&bbdfk7pt);
	local_mc.push_back(&mc4);
	local_mc.push_back(&mc5);
	local_mc.push_back(&mc6);
	local_mc.push_back(&mc7);
	local_mc.push_back(&mc8);
	local_mc.push_back(&mc9);

	local_mc.push_back(&gz7pt);
	local_mc.push_back(&gz8pt);

	local_ind[4]=&ind4;
	local_ind[5]=&ind5;
	local_ind[6]=&ind6;
	local_ind[7]=&ind7;
	local_ind[8]=&ind8;
	local_ind[9]=&ind9;
	local_ind[10]=&ind10;

	size_t old_n=0;
	vector<vector<int>*> all_new_index_vectors;
	vector<int>* ind=0;
	vector<int>* ind_ptr4=&ind4;
	vector<int>* ind_ptr5=&ind5;
	vector<int>* ind_ptr6=&ind6;
	vector<int>* ind_ptr7=&ind7;
	vector<int>* ind_ptr8=&ind8;
	vector<int>* ind_ptr9=&ind9;
	vector<int>* ind_ptr10=&ind10;

#if 0
	     // open the library
//	    cout << "Opening readline.so...\n";
	    void* handle = dlopen("libreadline.so.5", RTLD_LAZY);

	    if (!handle) {
	        cerr << "Cannot open library: " << dlerror() << '\n';
	        return 1;
	    }

	    // load the symbol
//	    cout << "Loading symbol readline...\n";
	    typedef char* (*ff_t)(char*);
	    ff_t readline = (ff_t) dlsym(handle, "readline");
	    if (!readline) {
	        cerr << "Cannot load symbol 'readline': " << dlerror() <<
	            '\n';
	        dlclose(handle);
	        return 1;
	    }

	    typedef void (*fff_t)(char*);
	    fff_t add_history = (fff_t) dlsym(handle, "add_history");
	    if (!add_history) {
	        cerr << "Cannot load symbol 'add_history': " << dlerror() <<
	            '\n';
	        dlclose(handle);
	        return 1;
	    }



	    fileman_completion_type rl_attempted_completion_function_from_readline = reinterpret_cast<fileman_completion_type>(dlsym(handle, "rl_attempted_completion_function"));
	    if (!rl_attempted_completion_function_from_readline) {
	        cerr << "Cannot load symbol 'rl_attempted_completion_function': " << dlerror() <<
	            '\n';
	        dlclose(handle);
	        return 1;
	    } else {
	    	  /* Tell the completer that we want a crack first. */
	    	  rl_attempted_completion_function_from_readline = fileman_completion;
	    	  /* Tell the completer that we want a crack first. */
//	    	  rl_attempted_completion_function = fileman_completion;
	    }

	    rl_completion_matches_from_readline = reinterpret_cast<rl_completion_matches_type>(dlsym(handle, "rl_completion_matches"));
	    if (!rl_completion_matches_from_readline) {
	        cerr << "Cannot load symbol 'rl_completion_matches': " << dlerror() <<
	            '\n';
	        dlclose(handle);
	        return 1;
	    }


#endif


#if USE_READLINE
   	  rl_attempted_completion_function = fileman_completion;
#endif
   	bool history_file_exist = false;
   	fstream fin;
   	fin.open(".BH_history",ios::in);
   	if( fin.is_open() )
   	{
   		history_file_exist=true;
   	}
   	fin.close();
#if USE_HISTORY
   	if (history_file_exist){
   	   	read_history_range (".BH_history", 0,100);
   	}
#endif

	session_info si;

	string new_line;
	bhc_status status = BHC_CONTINUE;
	try {
	while ( status == BHC_CONTINUE ){
#if USE_READLINE
		char *rl_line;
		cout << GREEN;
		rl_line=readline( "]> " );
		string rl_line_str(rl_line);
		stringstream ss(rl_line);
#else
		cout << GREEN "]> ";
		getline(cin,new_line);
		stringstream ss(new_line);
#endif
			cout << NOCOLOR ;
			status = execute(ss,si);
#if USE_READLINE
		    if (rl_line){
			    add_history (rl_line);
		    }
			free(rl_line);
#endif
			}
#if 0
//    cout << "Closing readline library ...\n";
    dlclose(handle);
#endif
	}
	catch (BHerror& be){
		be.print();
		throw;
	}

	catch (...){
#if USE_HISTORY
	ofstream File;
	File.open(".BH_history",ios::out|ios::trunc);
	append_history (100, ".BH_history");
#endif
	}

	for_each(all_new_index_vectors.begin(),all_new_index_vectors.end(),do_delete<vector<int> >());

	fpu_fix_end(&old_cw);

#if USE_HISTORY
	ofstream File;
	File.open(".BH_history",ios::out|ios::trunc);
	append_history (100, ".BH_history");
#endif

}
