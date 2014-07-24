/*
 * parent_diagrams_from_file.cpp
 *
 *  Created on: Oct 4, 2009
 *      Author: daniel
 */

#define _VERBOSE 1

#ifndef BH_PUBLIC
#include "parent_diagrams_from_file.h"
#include "parent_diagrams.h"
#endif
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <map>
#include "particles.h"
#include "BH_utilities.h"
#include "partitions.h"
#include "iterators.h"
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <cassert>
#include "BH_debug.h"
#include "settings.h"
#include "eval_param.h"
#include "BHpath.h"
#include <sys/stat.h>
#include <openssl/md5.h>
#include <stdlib.h>

using namespace std;

namespace BH {

string NoSpaces(string s)
{
  while(s.find(" ") != string::npos)
  {
    s.replace(s.find(" "), 1, "");
  }
  return s;
}





std::map<std::string,particle*> init_particle_map(){
	std::map<std::string,particle*> m;
	m.insert(std::make_pair("g",&gluon));
	m.insert(std::make_pair("q",&quark));
	m.insert(std::make_pair("qb",&quark));
	m.insert(std::make_pair("l",&lepton));
	m.insert(std::make_pair("y",&photon));
	m.insert(std::make_pair("Gl",&gluino));
	m.insert(std::make_pair("Glb",&gluino));
	m.insert(std::make_pair("sc",&scalar_massive));
	return m;
}

particle_ID map_massless_to_massive(particle_ID pID){
	size_t hel(pID.helicity());
	particle* type(pID.type());
        int flav(pID.flavor());
	if(type==&gluon) {type=&scalar_massive; hel=0;}
	else if(type==&scalar) {type=&scalar_massive; hel=0;}
	// use the quark flavor 105 for nf-terms with external gluons only
	// we then compute the rational term with a scalari_massive 
	// in order to reduce the helicity sum. 
	// trick not yet possible for external quarks, as we do not have a 
	// conserved massive_scalar, but rather it mimics a D-dim gluon mode
	else if(type==&quark && flav==105) {type=&scalar_massive; hel=0; flav=0; }
	else if(type==&quark) type=&quark_massive;
	else if(type==&gluino) type=&gluino_massive;
//	else if(type==&gluino_massive) type=&gluino_massive;
//	else if(type==&scalar_massive){type=&scalar_massive;hel=0;}
//	else if(type==&quark_massive) type=&quark_massive;
//	else if(type==&gluino_massive) type=&gluino_massive;
	else{_MESSAGE("Check consistency in map_massless_to_massive.");};
	
	return particle_ID(type,hel,flav,pID.is_anti_particle());
	
}

std::map<color_structure,std::string> init_cs_map(){
	std::map<color_structure,std::string> m;
	m.insert(std::make_pair(nf,"nf"));
	m.insert(std::make_pair(leading_color,"leading_color"));
	m.insert(std::make_pair(glue,"glue"));
	m.insert(std::make_pair(LT,"L"));
	m.insert(std::make_pair(RT,"R"));
	m.insert(std::make_pair(LLT,"LL"));
	m.insert(std::make_pair(LRT,"LR"));
	m.insert(std::make_pair(RLT,"RL"));
	m.insert(std::make_pair(RRT,"RR"));
	m.insert(std::make_pair(LLLT,"LLL"));
	m.insert(std::make_pair(RLLT,"RLL"));
	m.insert(std::make_pair(LRLT,"LRL"));
	m.insert(std::make_pair(LLRT,"LLR"));
	m.insert(std::make_pair(RRLT,"RRL"));
	m.insert(std::make_pair(RLRT,"RLR"));
	m.insert(std::make_pair(LRRT,"LRR"));
	m.insert(std::make_pair(RRRT,"RRR"));
	m.insert(std::make_pair(nfLT,"nfL"));
	m.insert(std::make_pair(nfRT,"nfR"));
	m.insert(std::make_pair(nfLLT,"nfLL"));
	m.insert(std::make_pair(nfLRT,"nfLR"));
	m.insert(std::make_pair(nfRLT,"nfRL"));
	m.insert(std::make_pair(nfRRT,"nfRR"));
//
	m.insert(std::make_pair(nfLLLT,"nfLLL"));
	m.insert(std::make_pair(nfRLLT,"nfRLL"));
	m.insert(std::make_pair(nfLRLT,"nfLRL"));
	m.insert(std::make_pair(nfLLRT,"nfLLR"));
	m.insert(std::make_pair(nfRRLT,"nfRRL"));
	m.insert(std::make_pair(nfRLRT,"nfRLR"));
	m.insert(std::make_pair(nfLRRT,"nfLRR"));
	m.insert(std::make_pair(nfRRRT,"nfRRR"));
	return m;
}




particle* particle_type_from_string(const std::string& s){
	static std::map<std::string,particle*> s_map=init_particle_map();
	// to get rid of leading and trailing spaces
	string partic = NoSpaces(s);


	std::map<std::string,particle*>::iterator it=s_map.find(partic);

	if (it != s_map.end()){
		return (*it).second;
	}
	else {
		_MESSAGE2("problems reading particle type ",partic);
		throw BHerror("Syntax error");
	}
}


particle_ID get_particle_ID_from_string(const std::string& str){
	bool ap=false;
	int flavor=1;
	int helicity=1;
	BH_DEBUG_MESSAGE3("get_particle_ID_from_string for string: '",str,"'");
	int end_of_type=-1;

	int pos_b=str.find('b',0);
	if ( pos_b !=string::npos ){
		ap=true;
		end_of_type=pos_b;
		BH_DEBUG_PRINT(pos_b);
		BH_DEBUG_MESSAGE("anti-particle");
	} else {
		BH_DEBUG_MESSAGE("particle");
	}

	int pos_f=1000;
	for (int j=0;j<10;j++){
		int pos_f_new=str.find(char(int('0')+j),0);
//		_PRINT(char(int('0')+j));
//		_PRINT(pos_f_new);
		if ( pos_f_new != string::npos ){
			BH_DEBUG_MESSAGE2("found ",j+1);
			if ( pos_f_new < pos_f ){
				pos_f=pos_f_new;
			}
		}
	}
	// treats the case of negative flavor number
	int pos_f_new=str.find('-',0);
	if ( pos_f_new != string::npos ){
		//			_MESSAGE2("found ",j+1);
		if ( pos_f_new < pos_f ){
			pos_f=pos_f_new;
		}
	}

	if (pos_f != 1000){
		if (end_of_type == -1){
			end_of_type = pos_f;
		}
	}

	int end_of_flavor=str.size();

	int pos_m=str.find('m',0);
	if ( pos_m !=string::npos ){
		helicity = -1 ;
		end_of_flavor--;
		if (end_of_type == -1){
			end_of_type=pos_m;
		}
	}
	int pos_p=str.find('p',0);
	if ( pos_p !=string::npos ){
		helicity = +1 ;
		end_of_flavor--;
		if (end_of_type == -1){
			end_of_type=pos_p;
		}
	}

	if (pos_f != 1000){
		string flavor_str=str.substr(pos_f,end_of_flavor);
		BH_DEBUG_MESSAGE2("flavor: ",flavor_str);
		flavor=atoi(flavor_str.c_str());
	}
	if (end_of_type == -1){
		end_of_type=str.size();
	}

	BH_DEBUG_MESSAGE2("type: ",str.substr(0,end_of_type));

	particle* p_type = particle_type_from_string(str.substr(0,end_of_type));

	if ( *p_type == scalar || *p_type == gluon_massive_scalar || *p_type == higgs ){
		helicity=0;
	}
	return particle_ID(p_type,helicity,flavor,ap);
}

plabel plabel_from_string(const std::string& s){
	int pos=s.find('(',0);
	int pos2=s.find(')',0);

	std::string part=s.substr(0,pos);
	std::string index=s.substr(pos+1,pos2-pos-1);
	int ind; std::stringstream(index) >> ind;

	return plabel(get_particle_ID_from_string(part),ind);

}



void read_processes(const std::string& input,	vector<vector<plabel> >& labels){
	int bracket_open=0, bracket_end=0;


	vector<particle_ID> props;

	int expr_start = input.find('{',0);
	int expr_end = input.find('}',expr_start);
	if (expr_start == string::npos){
		_WARNING("Missing \'{\'");
		throw BHerror("Syntax error");
	}
	if (expr_end == string::npos){
		_WARNING("Missing \'}\'");
		throw BHerror("Syntax error");
	}

	bracket_open = expr_start;

	while ( bracket_open  != expr_end ){

		bracket_end = input.find('|',bracket_open+1);
		if (bracket_end == string::npos){
			bracket_end = expr_end;
		}
		labels.push_back(std::vector<plabel>());

		std::string corner_str=input.substr(bracket_open+1,bracket_end-bracket_open-1);

//		_PRINT(bracket_end);
//		_PRINT(bracket_open);
//		_MESSAGE3("corner_str: \'",corner_str,"\'");
		std::stringstream ss_corner(corner_str);

		while (ss_corner.good()){
			std::string partic;
			ss_corner >> partic;
			if ( !partic.empty() ){
//				_MESSAGE3("partic: \'",partic,"\'");
				labels.back().push_back(plabel_from_string(partic));
			}
		}


		bracket_open = bracket_end;

	}

//	_PRINT(labels.size());

//	for (int i=0;i< labels.size();i++){
//		copy(labels[i].begin(),labels[i].end(),ostream_iterator<plabel>(cout," ")); cout << endl;
//	}
}

//! Returns the directory where the parent data files are stored.
/**
The order of priority is the following:
	 1) a directory specified through a BHsettings file
	 2) the local source directory ${abs_srcdir}/datafiles/parents/
	 3) the installation directory ${abs_top_src_dir}/datafiles/parents/
*/
std::string GetParentDataDirectory(){
	if (settings::general::s_parent_data_path != string("not set")){
		BH_DEBUG_MESSAGE3("Using ",settings::general::s_parent_data_path," as the partent data path.");
		return settings::general::s_parent_data_path;
	} else {
		struct stat st;
		std::string datafiles="/datafiles/parents/";
		if( stat(BH_SOURCE_PATH,&st) == 0){
			BH_DEBUG_MESSAGE3("Using ",BH_SOURCE_PATH+datafiles," as the parent data path.");
			return BH_SOURCE_PATH+datafiles;
		} else {
			BH_DEBUG_MESSAGE3("Using ",BH_INSTALL_PATH+datafiles," as the parent data path.");
			return BH_INSTALL_PATH+datafiles;
		}
	}
}


bool getPathFromEnv(string& path) {

	char * pPath;
    pPath = getenv ("WORKER_DATA_PATH");
    char path_points[50];
    char path_points_HP[50];
    if (pPath!=NULL)
    {
  	  path= string(pPath);
  	   return true;
    } else {
    	return false;
    }

}

/* First use the envirenoment, then BHsettings and then default to $prefix/share/datafiles
 * the environment is useful for
 * */

std::string get_worker_dir(const std::string& subdir){
		struct stat st;
		static std::string pathFromEnv;
		static bool has_env=getPathFromEnv(pathFromEnv);
	    if (has_env)
	    {
	    	return pathFromEnv +std::string("/") + subdir  +std::string("/") ;
	    }


		if (settings::general::s_data_path != "not set"){
			return settings::general::s_data_path  +std::string("/") + subdir  +std::string("/") ;
			//_WARNING("Data path not set. Please set it using the DATA_PATH setting in your BHsettings file. ");
		} else {
			std::string datafiles="/share/blackhat/datafiles/";
			if ( stat(BH_INSTALL_PATH,&st) == 0){
				BH_DEBUG_MESSAGE3("Using ",BH_INSTALL_PATH+datafiles," as the tree data path.");
				return std::string(BH_INSTALL_PATH)+datafiles+subdir;
//			} else if ( stat(BH_SOURCE_PATH,&st) == 0){
//				BH_DEBUG_MESSAGE3("Using ",BH_SOURCE_PATH+datafiles," as the tree data path.");
//				return std::string(BH_SOURCE_PATH) +std::string("/datafiles/")+subdir;
			} else {
				throw BHerror("No path for data!");
			}
		}
}

char letters[]="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";

string myHash(const string& input,int maxLength){
//	locale loc;                 // the "C" locale
//	const collate<char>& coll = use_facet<collate<char> >(loc);
//	unsigned long myhash = coll.hash(input.data(),input.data()+input.length());
//	unsigned long rest=myhash;
//	char buffer[10];
//	int index=0;
//	while (rest != 0 and index<maxLength){
//		short mod=rest%52;
//		rest/=52;
//		buffer[index]=letters[mod];
//		index++;
//	}
//	buffer[index]='\0';

	unsigned char md[33];
	MD5((unsigned char *)input.data(),input.size(),&md[0]);

	char md5string[33];
	for(int i = 0; i < maxLength ; ++i){
		md5string[i]=letters[(unsigned int)md[i]%62];
	}
	md5string[maxLength]='\0';
//	    sprintf(&md5string[i*2], "%02x", (unsigned int)md[i]);

//	cout << md5string << std::endl;

	return string(md5string);
}


}
