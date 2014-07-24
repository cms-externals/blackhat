/*
 * known_rec_tree.cpp
 *
 *  Created on: 24 Jul 2009
 *      Author: daniel
 */

#include "process.h"
#include "settings.h"
#include "BH_error.h"
#include "BH_debug.h"
#include "tree_amp.h"
#include "helcode.h"
#include "process_utils.h"
#include "trees_eval/amplitudes_tree_eval.h"
#if BH_USE_GMP
#include "gmp_r.h"
#endif

#define _DISTINGUISH_QUARKS_GLUINOS_REC 0 // Set to 1 if you want to treat gluinos differently to quarks in the recursion relations, otherwise set to 0 (so far this makes no difference)

using namespace std;

namespace BH {


string string_name(const particle_ID& p,size_t flavor){
	if ( p.type() == gluon ) {
		switch (p.helicity() ){
		case 1: return string("gp");
		case -1: return string("gm");
		}
	}
	if ( p.type() == quark
#if _DISTINGUISH_QUARKS_GLUINOS_REC==0
			|| p.type() == gluino
#endif
		) {
		string res;
		switch (flavor) {
		case 1: res="q"; break;
		case 2: res="u"; break;
		case 3: res="d"; break;
		case 4: res="s"; break;
		case 5: res="c"; break;
		case 6: res="b"; break;
		}
		switch (p.helicity() ){

		case 1: res+="p";break;
		case -1: res+="m";break;
		}
		return res;
	}
#if _DISTINGUISH_QUARKS_GLUINOS_REC==1
	if (p.type() == gluino) {
		string res;
		switch (flavor) {
		case 1: res="z"; break;
		case 2: res="x"; break;
		case 3: res="v"; break;
		case 4: res="n"; break;
		case 5: res="a"; break;
		case 6: res="f"; break;
		}
		switch (p.helicity() ){

		case 1: res+="p";break;
		case -1: res+="m";break;
		}
		return res;
	}
#endif
	if ( p.type() == quark_massive
#if _DISTINGUISH_QUARKS_GLUINOS_REC==0
			|| p.type() == gluino_massive
#endif
	) {
		string res;
		switch (flavor) {
		case 1: res="Q"; break;
		case 2: res="U"; break;
		case 3: res="D"; break;
		case 4: res="S"; break;
		case 5: res="C"; break;
		case 6: res="B"; break;
		}
		switch (p.helicity() ){

		case 1: res+="p";break;
		case -1: res+="m";break;
		}
		return res;
	}
#if _DISTINGUISH_QUARKS_GLUINOS_REC==1
	if (p.type() == gluino_massive ) {
		string res;
		switch (flavor) {
		case 1: res="Z"; break;
		case 2: res="X"; break;
		case 3: res="V"; break;
		case 4: res="N"; break;
		case 5: res="A"; break;
		case 6: res="F"; break;
		}
		switch (p.helicity() ){

		case 1: res+="p";break;
		case -1: res+="m";break;
		}
		return res;
	}
#endif
	if ( p.type() == lepton ) {
		switch (p.helicity() ){
		case 1: return string("lp");
		case -1: return string("lm");
		}
	}
	if ( p.type() == scalar) {
		switch (p.helicity() ){
		case 0: return string("s0");
		case 1: return string("s0");
		case -1: return string("s0");
		}
	}
	if ( p.type() == gluon_massive_scalar ||  p.type() == gluon_massive ||  p.type() == scalar_massive) {
		switch (p.helicity() ){
            case 0: return string("S0");
            case 1: return string("Sp");
            case -1: return string("Sm");
		}
	}
	if ( p.type() == photon) {
        string res;
		switch (flavor%2) { //maximal 2 photons
		case 0: res="y"; break;
		//case 1: res="y"; break;
		case 1: res="e"; break;
		}
		switch (p.helicity() ){
		case 1: res+="p";break;
		case -1: res+="m";break;
		}
	    return res;
    }
	if ( p.type() == higgs) {
		switch (p.helicity() ){
		case 0: return string("h0");
		}
	}
	_WARNING2("Failed to find string for particle of type ",p);
	return string("");
}


string string_name(const process& pro){
	string res;
	map<int,int> flavor_index, flavor_index_massive, flavor_index_gluino, flavor_index_massive_gluino;

	for(int i=1;i<=pro.n();i++){
		if(pro.p(i).is_a(quark)){
			map<int,int>::iterator pos;
			pos=flavor_index.find(pro.p(i).flavor());
			if ( pos == flavor_index.end() ) {
				int new_flavour=flavor_index.size()+flavor_index_gluino.size()+1;
				flavor_index.insert(pair<int,int>(pro.p(i).flavor(),new_flavour));
				res+=string_name(pro.p(i),new_flavour);
			} else {
				res+=string_name(pro.p(i),pos->second);
			}
		}
		else if(pro.p(i).is_a(gluino)){
			map<int,int>::iterator pos;
			pos=flavor_index_gluino.find(pro.p(i).flavor());
			if ( pos == flavor_index_gluino.end() ) {
				int new_flavour=flavor_index.size()+flavor_index_gluino.size()+1;
				flavor_index_gluino.insert(pair<int,int>(pro.p(i).flavor(),new_flavour));
				res+=string_name(pro.p(i),new_flavour);
			} else {
				res+=string_name(pro.p(i),pos->second);
			}
		}
		else if(pro.p(i).is_a(quark_massive)){
			map<int,int>::iterator pos;
			pos=flavor_index_massive.find(pro.p(i).flavor());
			if ( pos == flavor_index_massive.end() ) {
				int new_flavour=flavor_index_massive.size()+flavor_index_massive_gluino.size()+1;
				flavor_index_massive.insert(pair<int,int>(pro.p(i).flavor(),new_flavour));
				res+=string_name(pro.p(i),new_flavour);
			} else {
				res+=string_name(pro.p(i),pos->second);
			}
		}
		else if(pro.p(i).is_a(gluino_massive)){
			map<int,int>::iterator pos;
			pos=flavor_index_massive_gluino.find(pro.p(i).flavor());
			if ( pos == flavor_index_massive_gluino.end() ) {
				int new_flavour=flavor_index_massive.size()+flavor_index_massive_gluino.size()+1;
				flavor_index_massive_gluino.insert(pair<int,int>(pro.p(i).flavor(),new_flavour));
				res+=string_name(pro.p(i),new_flavour);
			} else {
				res+=string_name(pro.p(i),pos->second);
			}
		}
        	else if(pro.p(i).is_a(photon)){
			map<int,int>::iterator pos;
			pos=flavor_index.find(pro.p(i).flavor());
			if ( pos == flavor_index.end() ) {
				int new_flavour=flavor_index.size()+flavor_index_gluino.size()+1;
				flavor_index.insert(pair<int,int>(pro.p(i).flavor(),new_flavour));
				res+=string_name(pro.p(i),new_flavour);
			} else {
				res+=string_name(pro.p(i),pos->second);
			}
		}
		else res+=string_name(pro.p(i),0);

	}
			return res;
}


#ifdef USE_MC_ALSO
template <class T> complex<T> (*A_Tree_Ptr(const process& pro))(momentum_configuration<T>&,const  vector<int>&){

	if (Tree_is_zero(pro)) {return &ZeroF;}
	switch(pro.pcode()){

		case 3: return A3g_Tree_Ptr<T>(helcode_g(pro));
        case 4:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A4g_G3_Tree_Ptr<T>(helcode_g(pro)); //add G3 insteraction for pure gluons
                    case false: return A4g_Tree_Ptr<T>(helcode_g(pro));}}
        case 5:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A5g_G3_Tree_Ptr<T>(helcode_g(pro)); //add G3 insteraction for pure gluons
                    case false: return A5g_Tree_Ptr<T>(helcode_g(pro));}}
		case 6: return A6g_Tree_Ptr<T>(helcode_g(pro));
		case 7: return A7g_Tree_Ptr<T>(helcode_g(pro)); //REMOVED TO MATCH EVAL_PARAM CODE
		case 8: return A8g_Tree_Ptr<T>(helcode_g(pro)); //REMOVED TO MATCH EVAL_PARAM CODE
		case 21:return A2q1g_Tree_Ptr<T>(helcode_2q(fix_flavors(pro)));
		case 22:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A2q2g_G3_Tree_Ptr<T>(helcode_2q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q2g_Tree_Ptr<T>(helcode_2q(fix_flavors(pro)));}}
        case 23:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A2q3g_G3_Tree_Ptr<T>(helcode_2q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q3g_Tree_Ptr<T>(helcode_2q(fix_flavors(pro)));}}
        case 24:return A2q4g_Tree_Ptr<T>(helcode_2q(fix_flavors(pro))); //REMOVED TO MATCH EVAL_PARAM CODE
		case 40:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return A4q_Tree_Ptr<T>(helcode_4q(fix_flavors(pro)));
			case 2:	{switch(settings::general::s_use_g3_coupling){
                        //distinct flavors only
                    case true: return A2q2Q_G3_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q2Q_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));}}
			default: throw BHerror("wrong number of quark flavors");
			}
		case 41:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return A4q1g_Tree_Ptr<T>(helcode_4q(fix_flavors(pro)));
            case 2:{switch(settings::general::s_use_g3_coupling){
                       //distinct flavors only
                    case true: return A2q1g2Q_G3_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q1g2Q_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));}}
			default:throw BHerror("wrong number of quark flavors");
			}
//		case 42:
//			switch (nbr_of_flavors(fix_flavors(pro),quark)){
//			case 1: return A4q2g_Tree_Ptr<T>(helcode_4q(fix_flavors(pro)));
//			case 2: return A2q2g2Q_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
//			default: throw BHerror("wrong number of quark flavors");
//			}
//		case 43:
//			switch (nbr_of_flavors(fix_flavors(pro),quark)){
//			case 1: return A4q3g_Tree_Ptr<T>(helcode_4q(fix_flavors(pro)));
//			case 2: return A2q3g2Q_Tree_Ptr<T>(helcode_2q2Q(fix_flavors(pro)));
//			default: throw BHerror("wrong number of quark flavors");
//			}
		case 220: return A2q2l_Tree_Ptr<T>(helcode_2q2l(pro));
		case 221: return A2q1g2l_Tree_Ptr<T>(helcode_2q2l(pro));
		case 222: return A2q2g2l_Tree_Ptr<T>(helcode_2q2l(pro)); //posiible problem in future: points to full color trees unlike eval-param
		case 223: return A2q3g2l_Tree_Ptr<T>(helcode_2q2l(pro)); //posiible problem in future: points to full color trees unlike eval-param
		case 224: return A2q4g2l_Tree_Ptr<T>(helcode_2q2l(pro));
		case 225: return A2q5g2l_Tree_Ptr<T>(helcode_2q2l(pro));
		// 2 scalars and n-2 gluon amplitudes
		case 2001: return A2s1g_Tree_Ptr<T>(helcode_2s(pro));
		case 2002: return A2s2g_Tree_Ptr<T>(helcode_2s(pro));
		case 2003: return A2s3g_Tree_Ptr<T>(helcode_2s(pro));
		case 2004: return A2s4g_Tree_Ptr<T>(helcode_2s(pro));
		// 2 scalars and n-2 quarks
		case 2010: return &ZeroF;
		case 2020: return A2s2q_Tree_Ptr<T>(helcode_2Q2qs_massive(fix_flavors(pro))); // We do not care about the flavours here, at least for now
		case 2200: return A2s2l_Tree_Ptr<T>(helcode_2Q2qs_lepton_massive(fix_flavors(pro)));
//		case 2220: return A2s2q2l_Tree_Ptr<T>(helcode_2Q2qs_lepton_massive(fix_flavors_massive(pro)));// DOES NOT WORK

		// Set these to zero for now in the future we may need them
		case 4000: case 4001: case 4002: case 4003: case 4004: case 4005: return &ZeroF;

		// 2 Massive quarks and n gluons
		case 20001: return A2QM1g_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(pro)));
		case 20002: return A2QM2g_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(pro)));
		// 2 Massive quarks and n massless quarks
		case 20010: return &ZeroF;
		case 20020: return A2QM2q_Tree_Ptr<T>(helcode_2Q2qs_massive(fix_flavors(pro)));
		case 20200: return A2QM2l_Tree_Ptr<T>(helcode_2Q2qs_lepton_massive(fix_flavors(pro)));
		// 1 Massive quark and 1 massless quark and a scalar
		case 11010: return A1QM1qs_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(pro)));
		// 1 Massive quark and 1 massless quark and a scalar and a gluon
		case 11011:  return A1QM1q1gs_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(pro)));

		// 2 Massive quarks and two massless quarks and a photon (which couples to the massless quark)
//		case 120020: return A2QM2q1y_Tree_Ptr<T>(helcode_2Q2qs1y_massive(fix_flavors(pro))); // DOES NOT WORK

		case 100020: return A2q1g_Tree_Ptr<T>(helcode_2q(replace_photon_with_gluon(pro)));
		case 100200: return A2q1g_Tree_Ptr<T>(helcode_2q(replace_lepton_with_quark(replace_photon_with_gluon(pro))));

		// 1 photon and two massive quarks
		case 120000: return A2QM1g_Tree_Ptr<T>(helcode_2qs_massive(replace_photon_with_gluon(pro)));
		// 1 photon, 1 Massive quark, 1 Massive scalar and a quark
		case 111010: return A1QM1q1gs_Tree_Ptr<T>(helcode_2qs_massive(replace_photon_with_gluon(pro)));
		// 1 photon, 2 Massive quarks and a gluon
		case 120001: return A2QM1g1y_Tree_Ptr<T>(helcode_2Q1g1y_massive(pro));
		// 2 leptons, 2 Massive quarks and a gluon
		case 20201: return A2QM1g2l_Tree_Ptr<T>(helcode_2Q1g1y_massive(pro)); //INCOMPLETE

		case 240: switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q2l_Tree_Ptr<T>(helcode_2q2l2Q(pro));
			default: throw BHerror("wrong number of quark flavors");
		}
		case 241:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q1g2l_Tree_Ptr<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		case 242:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q2g2l_Tree_Ptr<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		case 243:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q3g2l_Tree_Ptr<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		// 2q1g1y
		case 100021: return A2q1g1y_Tree_Ptr<T>(helcode_2q1y(pro));
		// 2q1g1y
		case 100022: return A2q2g1y_Tree_Ptr<T>(helcode_2q1y(pro));
		// 2q2y	
		case 200020: return A2q1g1y_Tree_Ptr<T>(helcode_2q1y(replace_one_photon_with_gluon(pro)));
		// 4q1y
		case 100040:{
			BH_DEBUG_PRINT(fix_flavors(pro).pcode());
			switch(fix_flavors(pro).pcode()){
				case 100040: return A2q2Q1y_Tree_Ptr<T>(helcode(fix_flavors(pro)));
				case 2100020: return A_Tree_Ptr<T>(fix_flavors(pro));
				default: throw BHerror("bad flavor arrangement");
			}
		}
		// 2 gluinos
		case 2000001: case 2000002: case 2000003: case 2000004: case 2000005: case 2000006: return A_Tree_Ptr<T>(replace_gluino_with_quark(pro));
		case 4000000: case 4000001: case 4000002: case 4000003: case 4000004: case 4000005: case 4000006: return A_Tree_Ptr<T>(replace_gluino_with_quark(pro));

		case 2000020: case 2000021: case 2000022: case 2000023:case 2000024:{
			return A_Tree_Ptr<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1));
		}
		case 2100020: case 2100021: case 2100022: case 2100023: case 2100024: {
			return A2q2G1y_Tree_Ptr<T>(helcode_2q2G1y(pro));
		}
		case 2000220:
			 return A2q2l2G_Tree_Ptr<T>(helcode_2q2l2G(pro));
		case 2000221: case 2000222: case 2000223: {
			return A_Tree_Ptr<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1));
		}

		//2 gluinos and 2 massive scalars
		case 2002000: return A2s2q_Tree_Ptr<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(fix_flavors(pro)))); // We do not care about the flavours here, at least for now

		//2 gluinos and 2 massive quarks
		case 2020000:  return A2QM2q_Tree_Ptr<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark_massive)[0]+1)));

		//1 Massive gluino, 1 Massless gluino, 1 Massive quark and 1 massless quark
		case 11010010: return A2QM2q_Tree_Ptr<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark)[0]+1)));

		// 2 Massive gluinos and two massless quarks
		case 20000020: return A2QM2q_Tree_Ptr<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark)[0]+1)));
//			return A2LM2q_Tree_Ptr<T>(helcode_2L2qs_massive(fix_flavors(pro))); // Contains a bug
		// 2 Massive gluinos and two massless quarks and a photon
//		case 20100020: return A2LM2q1y_Tree_Ptr<T>(helcode_2L2qs1y_massive(fix_flavors(pro)));//DOES NOT WORK

		// 1 photon and two massive gluinos
		case 20100000: return &ZeroF;
		// 1 photon, 1 Massive gluino, 1 Massive scalar and a gluino
		case 11101000: return &ZeroF;
		// 1 photon, 2 Massive gluinos and a gluon
		case 20100001: return &ZeroF;

		// Massive gluino
		// 2 Massive gluinos and n gluons
		case 20000001: return A2QM1g_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		case 20000002: return A2QM2g_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 2 Massive gluinos and n massless gluinos
		case 21000000: return &ZeroF;
		case 22000000: return A2QM2q_Tree_Ptr<T>(helcode_2Q2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 1 Massive gluino and 1 massless gluino and a scalar
		case 11001000: return A1QM1qs_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 1 Massive gluino and 1 massless gluino and a scalar and a gluon
		case 11001001: return A1QM1q1gs_Tree_Ptr<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
						//A1LM1G1gs_Tree_Ptr<T>(helcode_2Ls_massive(fix_flavors(pro)));

		// Trees with the complex higgs scalar phi or phid
		case 100000002: return A1ph2g_Tree_Ptr<T>(helcode_phi_1q(pro));
//		case 100000003: return A1ph3g_Tree_Ptr<T>(helcode_phi_1q(pro));
		case 100000020: return &ZeroF;
//		case 100000021: return A1ph1g2q_Tree_Ptr<T>(helcode_phi_1q(fix_flavors(pro)));
//		case 102000001: return A1ph1g2q_Tree_Ptr<T>(helcode_phi_1q(replace_gluino_with_quark(fix_flavors(pro))));
//		case 100000022: return A1ph2g2q_Tree_Ptr<T>(helcode_phi_1q(pro)); //NOT IMPEMENTED
//		case 100000040: return A1ph2q2Q_Tree_Ptr<T>(helcode_phi_2q(fix_flavors(pro)));
//		case 104000000: return A1ph2q2Q_Tree_Ptr<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
//		case 100020020: return A1ph2q2QM_Tree_Ptr<T>(helcode_phi_2q(fix_flavors(pro)));
//		case 122000000: return A1ph2q2QM_Tree_Ptr<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		// Trees with the complex higgs scalar phi or phid and scalars/massive quarks
		case 100002000: return A1ph2sc_Tree_Ptr<T>(helcode_phi_1q(pro));
		case 100010010: return &ZeroF;
		case 100020000: return &ZeroF;

		default:
//			_WARNING3("Tree ",pro," is not known in Known_Tree_Rec constructor.");
			return 0;
		}
}
#endif
/*
 *
 *
 * Code for Known_Rec_Tree
 *
 *
 */

template <class T> complex<T> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<T>&, const mass_param_coll&){
	if (Tree_is_zero(pro)) {return &ZeroF_eval;}

	BH_DEBUG_MESSAGE4("A_Tree_Ptr_eval: ",pro.pcode()," ",pro);
    switch(pro.pcode()){

		case 3: return A3g_Tree_Ptr_eval<T>(helcode_g(pro));
        case 4:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A4g_G3_Tree_Ptr_eval<T>(helcode_g(pro)); //add G3 insteraction for pure gluons
                    case false: return A4g_Tree_Ptr_eval<T>(helcode_g(pro));}}
        case 5:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A5g_G3_Tree_Ptr_eval<T>(helcode_g(pro)); //add G3 insteraction for pure gluons
                    case false: return A5g_Tree_Ptr_eval<T>(helcode_g(pro));}}
		case 6: return A6g_Tree_Ptr_eval<T>(helcode_g(pro));
//		case 7: return A7g_Tree_Ptr_eval<T>(helcode_g(pro));
//		case 8: return A8g_Tree_Ptr_eval<T>(helcode_g(pro));
		case 21:return A2q1g_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro)));
        case 22:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A2q2g_G3_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q2g_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro)));}}
        case 23:	{switch(settings::general::s_use_g3_coupling){
                    case true: return A2q3g_G3_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q3g_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro)));}}
//		case 24:return A2q4g_Tree_Ptr_eval<T>(helcode_2q(fix_flavors(pro)));
		case 40:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return A4q_Tree_Ptr_eval<T>(helcode_4q(fix_flavors(pro)));
			case 2:	{switch(settings::general::s_use_g3_coupling){
                        //distinct flavors only
                    case true: return A2q2Q_G3_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q2Q_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));}}
			default: throw BHerror("wrong number of quark flavors");
			}
		case 41:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return A4q1g_Tree_Ptr_eval<T>(helcode_4q(fix_flavors(pro)));
            case 2:{switch(settings::general::s_use_g3_coupling){
                       //distinct flavors only
                    case true: return A2q1g2Q_G3_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro))); //add G3 insteraction for pure gluons
                    case false: return A2q1g2Q_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));}}
            default:throw BHerror("wrong number of quark flavors");
			}
//		case 42:
//			switch (nbr_of_flavors(fix_flavors(pro),quark)){
//			case 1: return A4q2g_Tree_Ptr_eval<T>(helcode_4q(fix_flavors(pro)));
//			case 2: return A2q2g2Q_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
//			default: throw BHerror("wrong number of quark flavors");
//			}
//		case 43:
//			switch (nbr_of_flavors(fix_flavors(pro),quark)){
////			case 1: return A4q3g_Tree_Ptr_eval<T>(helcode_4q(fix_flavors(pro)));
////			case 2: return A2q3g2Q_Tree_Ptr_eval<T>(helcode_2q2Q(fix_flavors(pro)));
//			default: throw BHerror("wrong number of quark flavors");
//			}
		case 220: return A2q2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		case 221: return A2q1g2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		case 222: return A2q2g2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		case 223: return A2q3g2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		case 224: return A2q4g2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		case 225: return A2q5g2l_Tree_Ptr_eval<T>(helcode_2q2l(pro));
		// 2 scalars and n-2 gluon amplitudes
		case 2001: return A2s1g_Tree_Ptr_eval<T>(helcode_2s(pro));
		case 2002: return A2s2g_Tree_Ptr_eval<T>(helcode_2s(pro));
		case 2003: return A2s3g_Tree_Ptr_eval<T>(helcode_2s(pro));
		case 2004: return A2s4g_Tree_Ptr_eval<T>(helcode_2s(pro));
		// 2 scalars and n-2 quarks
		case 2010: return &ZeroF_eval;
		case 2020: return A2s2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(fix_flavors(pro))); // We do not care about the flavours here, at least for now
		case 2200: return A2s2l_Tree_Ptr_eval<T>(helcode_2Q2qs_lepton_massive(fix_flavors(pro)));

		// Set these to zero for now in the future we may need them
		case 4000: case 4001: case 4002: case 4003: case 4004: case 4005: return &ZeroF_eval;

		// 2 Massive quarks and n gluons
		case 20001: return A2QM1g_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(pro)));
		case 20002: return A2QM2g_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(pro)));
		// 2 Massive quarks and n massless quarks
		case 20010: return &ZeroF_eval;
		case 20020: return A2QM2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(fix_flavors(pro)));
		case 20200: return A2QM2l_Tree_Ptr_eval<T>(helcode_2Q2qs_lepton_massive(fix_flavors(pro)));
		// 1 Massive quark and 1 massless quark and a scalar
		case 11010: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(pro)));
		// 1 Massive quark and 1 massless quark and a scalar and a gluon
		case 11011:  return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(pro)));

		case 100020: return A2q1g_Tree_Ptr_eval<T>(helcode_2q(replace_photon_with_gluon(pro)));
		case 100200: return A2q1g_Tree_Ptr_eval<T>(helcode_2q(replace_lepton_with_quark(replace_photon_with_gluon(pro))));

		// 1 photon and two massive quarks
		case 120000: return A2QM1g_Tree_Ptr_eval<T>(helcode_2qs_massive(replace_photon_with_gluon(pro)));
		// 1 photon, 1 Massive quark, 1 Massive scalar and a quark
		case 111010: return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(replace_photon_with_gluon(pro)));
		// 1 photon, 2 Massive quarks and a gluon
		case 120001: return A2QM1g1y_Tree_Ptr_eval<T>(helcode_2Q1g1y_massive(pro));
		// 2 photon, 2 Massive quarks (same a one photon and a gluon)
		case 220000: return A2QM1g1y_Tree_Ptr_eval<T>(helcode_2Q1g1y_massive(replace_one_photon_with_gluon(pro)));
		// 2 leptons, 2 Massive quarks and a gluon
		case 20201: return A2QM1g2l_Tree_Ptr_eval<T>(helcode_2Q1g1y_massive(pro)); //INCOMPLETE

		case 240: switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q2l_Tree_Ptr_eval<T>(helcode_2q2l2Q(pro));
			default: throw BHerror("wrong number of quark flavors");
		}
		case 241:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q1g2l_Tree_Ptr_eval<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		case 242:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q2g2l_Tree_Ptr_eval<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		case 243:
			switch (nbr_of_flavors(fix_flavors(pro),quark)){
			case 1: return 0;
			case 2: return A2q2Q3g2l_Tree_Ptr_eval<T>(helcode_2q2l2Q(fix_flavors(pro)));
			default: throw BHerror("wrong number of quark flavors");
			}
		// 2q1g1y
		case 100021: return A2q1g1y_Tree_Ptr_eval<T>(helcode_2q1y(pro));
		// 2q2g1y
		case 100022: return A2q2g1y_Tree_Ptr_eval<T>(helcode_2q1y(pro));
		//2q2y
		case 200020: return A2q1g1y_Tree_Ptr_eval<T>(helcode_2q1y(replace_one_photon_with_gluon(pro)));
		// 4q1y
		case 100040:{
			switch(fix_flavors(pro).pcode()){
				case 100040: return A2q2Q1y_Tree_Ptr_eval<T>(helcode_2q2G1y(fix_flavors(pro)));
				case 2100020: return A_Tree_Ptr_eval<T>(fix_flavors(pro));
				default: throw BHerror("bad flavor arrangement");
			}
		}
		// 2 gluinos
		case 2000001: case 2000002: case 2000003: case 2000004: case 2000005: case 2000006: return A_Tree_Ptr_eval<T>(replace_gluino_with_quark(pro));
		case 4000000: case 4000001: case 4000002: case 4000003: case 4000004: case 4000005: case 4000006: return A_Tree_Ptr_eval<T>(fix_flavors(replace_gluino_with_quark(pro)));
		case 2000020: case 2000021: case 2000022: case 2000023:case 2000024:{
			return A_Tree_Ptr_eval<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1));
		}
		case 2100020: case 2100021: case 2100022: case 2100023: case 2100024: {
			return A2q2G1y_Tree_Ptr_eval<T>(helcode_2q2G1y(pro));
		}
		case 2000220:
			 return A2q2l2G_Tree_Ptr_eval<T>(helcode_2q2l2G(pro));
		case 2000221: case 2000222: case 2000223: {
			return A_Tree_Ptr_eval<T>(replace_gluino_with_quark(pro,find_all_flavors(pro,quark)[0]+1));
		}

		//2 gluinos and 2 massive scalars
		case 2002000: return A2s2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(fix_flavors(pro)))); // We do not care about the flavours here, at least for now

		//2 gluinos and 2 massive quarks
		case 2020000:  return A2QM2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark_massive)[0]+1)));

		//1 Massive gluino, 1 Massless gluino, 1 Massive quark and 1 massless quark
		case 11010010: return A2QM2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		// 2 Massive gluinos and two massless quarks
		case 20000020:  return A2QM2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(pro,find_all_flavors(fix_flavors(pro),quark)[0]+1)));

		// Massive gluino
		// 2 Massive gluinos and n gluons
		case 20000001: return A2QM1g_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		case 20000002: return A2QM2g_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 2 Massive gluinos and n massless gluinos
		case 21000000: return &ZeroF_eval;
		case 22000000: return A2QM2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 1 Massive gluino and 1 massless gluino and a scalar
		case 11001000: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));
		// 1 Massive gluino and 1 massless gluino and a scalar and a gluon
		case 11001001: return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(pro))));

		// 1 photon and two massive gluinos
		case 20100000: return &ZeroF_eval;
		// 1 photon, 1 Massive gluino, 1 Massive scalar and a gluino
		case 11101000: return &ZeroF_eval;
		// 1 photon, 2 Massive gluinos and a gluon
		case 20100001: return &ZeroF_eval;


		// Trees with the complex higgs scalar phi or phid
		case 100000002: return A1ph2g_Tree_Ptr_eval<T>(helcode_phi_1q(pro));
//		case 100000003: return A1ph3g_Tree_Ptr_eval<T>(helcode_phi_1q(pro));

		case 100000040: return A1ph2q2Q_Tree_Ptr_eval<T>(helcode_phi_2q(fix_flavors(pro))); // Broken
		case 102000020: return A1ph2q2Q_Tree_Ptr_eval<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		case 104000000: return A1ph2q2Q_Tree_Ptr_eval<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		case 100020020: return A1ph2q2QM_Tree_Ptr_eval<T>(helcode_phi_2q(fix_flavors(pro)));
		case 111010010: return A1ph2q2QM_Tree_Ptr_eval<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		case 120000020: return A1ph2q2QM_Tree_Ptr_eval<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		case 122000000: return A1ph2q2QM_Tree_Ptr_eval<T>(helcode_phi_2q(replace_gluino_with_quark(fix_flavors(pro),find_all_flavors(fix_flavors(pro),quark)[0]+1)));
		case 100000020: return &ZeroF_eval;

		case 102000000: return &ZeroF_eval;
		case 111000000: return &ZeroF_eval;
		case 120000000: return &ZeroF_eval;
            
		// Trees with the complex higgs scalar phi or phid and massive gluon scalars/massive quarks
		case 20100000000: return A1ph2sc_Tree_Ptr_eval<T>(helcode_phi_1q(pro));
		case 20100000001: return A1ph2sc1g_Tree_Ptr_eval<T>(helcode_phi_1q(pro));
		case 100020001: return A1ph2QM1g_Tree_Ptr_eval<T>(helcode_phi_1q(pro));
		case 120000001: return A1ph2QM1g_Tree_Ptr_eval<T>(helcode_phi_1q(replace_gluino_with_quark(pro)));
		case 100010010: return &ZeroF_eval;
		case 100020000: return &ZeroF_eval;
			
			
		// Trees with the a pair of massless scalars and gluons
		case 2000000001: return A2sc1g_Tree_Ptr_eval<T>(helcode_2qs_massless(pro));
		//case 2000000002: return A2sc2g_Tree_Ptr_eval<T>(helcode_2qs_massless(pro)); // NEED TO CHECK THESE
			
		// Trees with the a pair of massless scalars, gluons and a higgs
		case 2100000000: return &ZeroF_eval;
		case 2100000001: return A1ph2scm1g_Tree_Ptr_eval<T>(helcode_2qs_massless(pro));
			
//		// Trees with the a pair of massive gluon scalars and gluons
//		case 20000000001: return A2s1g_Tree_Ptr_eval<T>(helcode_2s(replace_scalar_massive_with_gluon(pro)));

			
		// Trees with massive scalars, gluons and a higgs
//		case 100002000: return A1ph2SM_Tree_Ptr_eval<T>(helcode_2qs_massless(pro));
		case 100002000: return &ZeroF_eval;
		case 100002001: return A1ph2SM1g_Tree_Ptr_eval<T>(helcode_2qs_massless(pro)); 
        case 100011010: return A1ph1QM1SM1q_Tree_Ptr_eval<T>(helcode_phi_SM_1q(pro));
		case 100002020: return A1ph2SM2q_Tree_Ptr_eval<T>(helcode_phi_SM_1q(pro));
		case 120002000: return A1ph2SM2q_Tree_Ptr_eval<T>(helcode_phi_SM_1q(replace_gluino_with_quark(pro)));

            
        // Massive gluons and massless gluons
        case 20000000001: return A2Gsc1g_Tree_Ptr_eval<T>(helcode_2Gsc(pro));
               
		// Massive gluons and quarks
		case 20000000010: return &ZeroF_eval;
		case 20000000020: return A2s2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(fix_flavors(replace_gluon_massive_with_scalar(pro)))); // We do not care about the flavours here, at least for now
		case 20000000200: return A2s2l_Tree_Ptr_eval<T>(helcode_2Q2qs_lepton_massive(fix_flavors(replace_gluon_massive_with_scalar(pro))));
			
		// 1 Massive quark, a massless quark and a massive gluon scalar
		case 10000010010: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluon_massive_with_scalar(pro))));
		case 10000010011: return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluon_massive_with_scalar(pro))));
		case 20100000020: return A1ph2sc2q_Tree_Ptr_eval<T>(helcode_phi_SM_1q(replace_gluon_massive_with_scalar(pro)));
		case 20120000000: return A1ph2sc2q_Tree_Ptr_eval<T>(helcode_phi_SM_1q(replace_gluino_with_quark(replace_gluon_massive_with_scalar(pro))));
		case 10100010010: return A1ph1QM1sc1q_Tree_Ptr_eval<T>(helcode_phi_1q(replace_gluon_massive_with_scalar(pro)));

		
		// 1 Massive quark, a massless quark and a massive gluon scalar and a photon
		//case 10000110010: return A_Tree_Ptr_eval<T>(replace_photon_with_gluon(pro));
		// careful cannot do with colored scalar around
		//case 10000110010: case 11110010: return A_Tree_Ptr_eval<T>(replace_photon_with_gluon(pro));

		case 100000010010: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluon_massive_with_scalar(pro))));
		case 100000010011: return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluon_massive_with_scalar(pro))));
			 
		case 100011000000: // These are all the same
		case 110000000001: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(replace_gluon_massive_with_scalar(pro)))));
			// 1 Massive gluino and 1 massless gluino and a scalar and a gluon
		case 100011000001: return A1QM1q1gs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(replace_gluon_massive_with_scalar(pro)))));
			
		case 200000000001: return A2s1g_Tree_Ptr_eval<T>(helcode_2s(replace_gluon_massive_with_scalar(pro)));
		case 200000000002: return A2s2g_Tree_Ptr_eval<T>(helcode_2s(replace_gluon_massive_with_scalar(pro)));
		case 200000000003: return A2s3g_Tree_Ptr_eval<T>(helcode_2s(replace_gluon_massive_with_scalar(pro)));
		case 200000000004: return A2s4g_Tree_Ptr_eval<T>(helcode_2s(replace_gluon_massive_with_scalar(pro)));
		case 200000002000: return A2s2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(fix_flavors(replace_gluon_massive_with_scalar(pro))))); // We do not care about the flavours here, at least for now
			
			
		case 10011000000: return A1QM1qs_Tree_Ptr_eval<T>(helcode_2qs_massive(fix_flavors(replace_gluino_with_quark(replace_gluon_massive_with_scalar(pro)))));
		case 20002000000: return A2s2q_Tree_Ptr_eval<T>(helcode_2Q2qs_massive(replace_gluino_with_quark(fix_flavors(replace_gluon_massive_with_scalar(pro))))); // We do not care about the flavours here, at least for now

		default:
//			_WARNING3("Tree ",pro," is not known in Known_Tree_Rec constructor.");
			return 0;
		}

	return 0;
}

// Explicit Instantiation
#ifdef USE_MC_ALSO

template complex<R> (*A_Tree_Ptr(const process& pro))(momentum_configuration<R>&,const  vector<int>&);
template complex<RHP> (*A_Tree_Ptr(const process& pro))(momentum_configuration<RHP>&,const  vector<int>&);
template complex<RVHP> (*A_Tree_Ptr(const process& pro))(momentum_configuration<RVHP>&,const  vector<int>&);

#if BH_USE_GMP
template complex<RGMP> (*A_Tree_Ptr(const process& pro))(momentum_configuration<RGMP>&,const  vector<int>&);
#endif
#endif
template complex<R> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<R>&,const mass_param_coll&);
template complex<RHP> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<RHP>&,const mass_param_coll&);
template complex<RVHP> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<RVHP>&,const mass_param_coll&);

#if BH_USE_GMP
template complex<RGMP> (*A_Tree_Ptr_eval(const process& pro))(const eval_param<RGMP>&,const mass_param_coll&);
#endif


}
